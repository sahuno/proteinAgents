#!/usr/bin/env python
# coding: utf-8

### single agent for A1 block only
from langgraph.graph import MessagesState
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.memory import MemorySaver

from langgraph.graph import START, StateGraph
from langgraph.prebuilt import tools_condition, ToolNode
from IPython.display import Image, display
from Bio import Entrez

import os, getpass
def _set_env(var: str):
    if not os.environ.get(var):
        os.environ[var] = getpass.getpass(f"Please enter your {var}: ")

_set_env("OPENAI_API_KEY")
_set_env("LANGSMITH_API_KEY")
_set_env("ANTHROPIC_API_KEY")
_set_env("NCBI_API_KEY")


Entrez.email = "ekwame001@gmail.com"


def search_pubmed(query, retmax=200):
    """Search PubMed and return results."""
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax,
        usehistory="y"
    )
    record = Entrez.read(handle)
    handle.close()
    return {
        "count": int(record["Count"]),
        "ids": record["IdList"],
        "webenv": record["WebEnv"],
        "query_key": record["QueryKey"],
    }

### Tool 2: concatenate pubmed queries
def concatenate_pubmed_queries(A1_block):
    """
    Concatenate PubMed queries with a base query using a logical operator .
    Args:
        A1_block (str): Terms for Inclusion criteria (concepts) that would concatenated with O_block terms.
    Returns:
        str: The combined PubMed query.

    Examples:
    Example 1: User wants to create an advanced pubmed query to retrieve research papers on protein protein interactions. Below are examples of A1 blocks that would be concatenated with O_block to create the final query.
    -----------
        # A1: INCLUSION CRITERIA - CONCEPTS (protein complexes, interactions)
A1_block = (
    '("nucleoproteins"[MeSH Terms] '
    'OR "protein interaction mapping"[MeSH Terms] '
    'OR (("nucleoprotein"[All Fields] OR "nucleoproteins"[All Fields] '
    'OR "multiprotein"[All Fields] OR "multiproteins"[All Fields] '
    'OR "proteins"[MeSH Terms] OR "protein"[All Fields] '
    'OR "proteins"[All Fields] OR "enzyme"[All Fields]) '
    'AND ("interact"[All Fields] OR "interacted"[All Fields] '
    'OR "interacting"[All Fields] OR "interaction"[All Fields] '
    'OR "interactions"[All Fields] OR "interactivity"[All Fields] '
    'OR "interacts"[All Fields])) '
    'OR "protein interaction"[All Fields] '
    'OR "protein interactions"[All Fields] '
    'OR "interacting protein"[All Fields] '
    'OR "interacting proteins"[All Fields] '
    'OR "multiprotein complexes"[MeSH Terms] '
    'OR (("nucleoprotein"[All Fields] OR "nucleoproteins"[All Fields] '
    'OR "multiprotein"[All Fields] OR "multiproteins"[All Fields] '
    'OR "proteins"[MeSH Terms] OR "protein"[All Fields] '
    'OR "proteins"[All Fields] OR "enzyme"[All Fields]) '
    'AND ("complex"[All Fields] OR "complexes"[All Fields] '
    'OR "heteromer"[All Fields] OR "heteromers"[All Fields] '
    'OR "homomer"[All Fields] OR "homomers"[All Fields] '
    'OR "heteromeric"[All Fields] OR "homomeric"[All Fields] '
    'OR "subunit"[All Fields] OR "subunits"[All Fields])) '
    'OR "protein complex"[All Fields] '
    'OR "protein complexes"[All Fields] '
    'OR ("protein"[All Fields] '
    'AND ("RNA"[All Fields] OR "DNA"[All Fields] '
    'OR "ribonucleic"[All Fields] OR "deoxyribonucleic"[All Fields]) '
    'AND ("interaction"[All Fields] OR "interactions"[All Fields])))'
)
    """
    # O_block (always constant): Terms to exclude unwanted publication types and non-English articles
    O_block = (
    '("english"[Language] '
    'NOT "meta-analysis"[Publication Type] '
    'NOT "review"[Publication Type] '
    'NOT "retracted publication"[Publication Type] '
    'NOT "retraction of publication"[Publication Type] '
    'NOT "published erratum"[Publication Type] '
    'NOT "controlled clinical trial"[Publication Type] '
    'NOT "clinical study"[Publication Type] '
    'NOT "clinical trial"[Publication Type] '
    'NOT "clinical trial protocol"[Publication Type] '
    'NOT "clinical trial, phase i"[Publication Type] '
    'NOT "clinical trial, phase ii"[Publication Type] '
    'NOT "clinical trial, phase iii"[Publication Type] '
    'NOT "clinical trial, phase iv"[Publication Type] '
    'NOT "clinical trial, veterinary"[Publication Type])'
)
    return f"({O_block} AND ({A1_block}))" 


llm = ChatOpenAI(model="gpt-4o")
tools = [concatenate_pubmed_queries, search_pubmed]
llm_with_tools = llm.bind_tools(tools)

# System message
sys_msg = SystemMessage(content="""You are a Molecular Biologist and an expert in constructing Advanced Pubmed Queries. 
                        you perform your job through a dialogue with scientists. 
                        Users give you a few examples of CONCEPTS they are interested in and your task is to understand and suggest similar CONCEPTS that are similar to the provided CONCEPTS list.
                        You will then use this understanding to create a refined pubmed query that can be used to retrieve those papers.
                        This means your role is to both help the scientists refine what they want to search for and then construct the actual PubMed Advanced query that captures the papers of interest.
                        
                        Understand what the user is for studying then begin by creating pubmed query terms we call `blocks` that capture different aspects of the user's interest.
                        Specifically, you will always create the following blocks:
                        1. A1_block: A block of terms (for example, MeSH terms) that captures the methods/approaches of interest called `A1_block`
                        Note that, There is always a block of terms `O_block` that excludes non-English articles and unwanted publication types such as reviews, meta-analyses, clinical trials etc. that should always be included in the final query.
                        You will then concatenate all the blocks into a final PubMed query using the following logic: `(O_block AND (A1_block))` which you can achieve by using the tool `concatenate_pubmed_queries()`.
                        Always call `concatenate_pubmed_queries()` tool to generate the final query for user and explicity say "THIS IS SUGGESTED FINAL QUERY".
                        Pls use `search_pubmed()` tool for searching PubMed after generating the final query.
                        """, 
                        name="system")

# If there there is human feedback it should be options, add this, remove this or that etc...

# Node
def assistant(state: MessagesState):
   return {"messages": [llm_with_tools.invoke([sys_msg] + state["messages"])]}


# Graph
builder = StateGraph(MessagesState)
# Define nodes: these do the work, i don't see why we should define `ToolNode(tools)` if we have already bound the tools to the llm; i guess this explicity uses tools as a node that could we used or ignored by the llm
builder.add_node("assistant", assistant)
builder.add_node("tools", ToolNode(tools))
# Define edges: these determine how the control flow moves
builder.add_edge(START, "assistant")
builder.add_conditional_edges(
    "assistant",
    # If the latest message (result) from assistant is a tool call -> tools_condition routes to tools
    # If the latest message (result) from assistant is a not a tool call -> tools_condition routes to END
    tools_condition,
)
builder.add_edge("tools", "assistant")
react_graph = builder.compile()

# Show
display(Image(react_graph.get_graph(xray=True).draw_mermaid_png()))

#### WITH MEMORY CASE:
memory = MemorySaver()
react_graph_with_memory = builder.compile(checkpointer=memory)
config = {"configurable": {"thread_id":"1"}}

pmed_humanMsg= "nucleoprotein, multiproteins, proteins, interacting, multiprotein complexes, complexes, RNA"
messages = [HumanMessage(content=pmed_humanMsg)]
messages = react_graph_with_memory.invoke({"messages": messages},config)
for m in messages['messages']:
    m.pretty_print()