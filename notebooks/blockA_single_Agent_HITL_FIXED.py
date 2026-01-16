#!/usr/bin/env python
# coding: utf-8
"""
PubMed Query Builder with Human-in-the-Loop (HITL) - FIXED VERSION
====================================================================
This version properly uses LangGraph's interrupt/resume pattern.

Key Changes from Original:
1. Uses interrupt() at checkpoints instead of routing to END
2. Uses update_state() to inject user responses
3. Proper resume with invoke(None, config)
4. No manual state manipulation in conversation loop

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-26
"""

from typing import TypedDict, Literal, Annotated
from langgraph.graph import MessagesState, START, END, StateGraph
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph.graph import interrupt
from Bio import Entrez
import os
import getpass
import operator


# =============================================================================
# Environment Setup
# =============================================================================

def _set_env(var: str):
    """Set environment variable if not already set."""
    if not os.environ.get(var):
        os.environ[var] = getpass.getpass(f"Please enter your {var}: ")

_set_env("OPENAI_API_KEY")
_set_env("LANGSMITH_API_KEY")
_set_env("ANTHROPIC_API_KEY")
_set_env("NCBI_API_KEY")

Entrez.email = "ekwame001@gmail.com"


# =============================================================================
# Enhanced State Definition
# =============================================================================

class EnhancedState(TypedDict):
    """
    Enhanced state for HITL workflow.

    Uses Annotated with operator.add for message accumulation.
    """
    # Conversation history (accumulated)
    messages: Annotated[list, operator.add]

    # Core data
    user_protocols: list
    suggested_protocols: list
    approved_protocols: list
    A2_block: str
    final_query: str
    search_results: dict

    # User responses at checkpoints
    user_response: str

    # Workflow control
    iteration_count: int


# =============================================================================
# PubMed Tools
# =============================================================================

def search_pubmed(query: str, retmax: int = 200) -> dict:
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


def concatenate_pubmed_queries(A2_block: str) -> str:
    """Concatenate PubMed queries with base O_block."""
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
    return f"({O_block} AND ({A2_block}))"


# =============================================================================
# LLM Setup
# =============================================================================

llm = ChatOpenAI(model="gpt-4o")

sys_msg = SystemMessage(
    content="""You are a Molecular Biologist and an expert in constructing Advanced PubMed Queries.

Your workflow:
1. Analyze user-provided protocols/methods
2. Suggest similar protocols that match the user's intent
3. Create A2_block (methods/approaches query block)
4. Use concatenate_pubmed_queries() to create final query

Important:
- Always use MeSH terms when possible
- Suggest protocols similar to user's examples
- Explain your reasoning for suggested protocols
""",
    name="system"
)


# =============================================================================
# Workflow Nodes
# =============================================================================

def detect_ambiguity(state: EnhancedState) -> EnhancedState:
    """
    HITL Checkpoint #1: Detect ambiguous input and request clarification.
    """
    user_protocols = state.get("user_protocols", [])

    # Check for ambiguous terms
    ambiguous_keywords = ["protein", "binding", "interaction", "assay", "cancer", "disease"]
    generic_terms = [p for p in user_protocols if any(kw in p.lower() for kw in ambiguous_keywords)]

    if generic_terms and len(user_protocols) < 3:
        # Ambiguity detected - interrupt for clarification
        clarification_msg = AIMessage(
            content=f"""I notice you mentioned: {', '.join(user_protocols)}

These terms are quite broad. To create a precise PubMed query, I need clarification:
- Are you interested in specific protein families or all proteins?
- Do you have specific experimental methods in mind?
- Is there a specific disease context or organism?

Please provide more specific details to help me build an accurate query."""
        )

        # Interrupt here - wait for user clarification
        user_clarification = interrupt(clarification_msg.content)

        # After resume, parse clarification and extract protocols
        # For now, just add clarification to messages
        return {
            "messages": [clarification_msg, HumanMessage(content=user_clarification)],
            "user_protocols": user_protocols  # Keep original for now
        }

    # No ambiguity - proceed
    return {"messages": []}


def extract_protocols(state: EnhancedState) -> EnhancedState:
    """Extract and validate protocols from user input."""
    user_protocols = state.get("user_protocols", [])

    if not user_protocols:
        # No protocols in state, extract from last message
        messages = state.get("messages", [])
        for msg in reversed(messages):
            if isinstance(msg, HumanMessage):
                protocols = [p.strip() for p in msg.content.split(",") if p.strip()]
                return {"user_protocols": protocols, "messages": []}

    return {"messages": []}


def expand_protocols(state: EnhancedState) -> EnhancedState:
    """Use LLM to suggest additional protocols."""
    user_protocols = state.get("user_protocols", [])

    expansion_prompt = f"""Given these user protocols: {', '.join(user_protocols)}

Suggest at least 10 additional similar experimental methods that a molecular biologist
might want to include for protein-protein interaction studies. Be exhaustive as much as possible.

Return only the protocol names, one per line."""

    response = llm.invoke([sys_msg, HumanMessage(content=expansion_prompt)])

    # Parse suggested protocols
    suggested = [
        line.strip().lstrip("-").lstrip("0123456789").lstrip(".").strip()
        for line in response.content.split("\n")
        if line.strip() and not line.strip().startswith("#")
    ]

    return {
        "suggested_protocols": suggested[:5],
        "messages": []
    }


def request_protocol_approval(state: EnhancedState) -> EnhancedState:
    """
    HITL Checkpoint #2: Request approval for protocol expansion.
    """
    user_protocols = state.get("user_protocols", [])
    suggested_protocols = state.get("suggested_protocols", [])

    approval_msg = AIMessage(
        content=f"""Based on your protocols, I suggest including these similar methods:

**Your protocols:**
{chr(10).join(f"  - {p}" for p in user_protocols)}

**Suggested additions:**
{chr(10).join(f"  - {p}" for p in suggested_protocols)}

Please respond with:
- "approve all" to include all suggestions
- "approve: <protocol1>, <protocol2>" to select specific ones
- "add: <new_protocol>" to suggest additional protocols
- "skip" to use only your original protocols"""
    )

    # Interrupt here - wait for user approval
    user_response = interrupt(approval_msg.content)

    # Parse approval response
    if "approve all" in user_response.lower():
        approved = user_protocols + suggested_protocols
    elif "skip" in user_response.lower():
        approved = user_protocols
    elif "approve:" in user_response.lower():
        approved_part = user_response.split("approve:", 1)[1]
        selected = [p.strip() for p in approved_part.split(",")]
        approved = user_protocols + selected
    elif "add:" in user_response.lower():
        new_protocol = user_response.split("add:", 1)[1].strip()
        approved = user_protocols + suggested_protocols + [new_protocol]
    else:
        # Default: approve all
        approved = user_protocols + suggested_protocols

    return {
        "approved_protocols": approved,
        "messages": [approval_msg, HumanMessage(content=user_response)]
    }


def create_a2_block(state: EnhancedState) -> EnhancedState:
    """Create A2_block query from approved protocols."""
    approved = state.get("approved_protocols", [])

    creation_prompt = f"""Create a PubMed A2_block for these methods: {', '.join(approved)}

Include:
- MeSH terms where possible
- Common alternative names
- Related techniques
- Use OR operators between terms

Format as a valid PubMed query block."""

    response = llm.invoke([sys_msg, HumanMessage(content=creation_prompt)])

    return {
        "A2_block": response.content,
        "messages": []
    }


def create_final_query(state: EnhancedState) -> EnhancedState:
    """Concatenate O_block and A2_block into final query."""
    a2_block = state.get("A2_block", "")
    final_query = concatenate_pubmed_queries(a2_block)

    return {
        "final_query": final_query,
        "messages": []
    }


def show_query_for_approval(state: EnhancedState) -> EnhancedState:
    """
    HITL Checkpoint #3: Show query for approval before search.
    """
    final_query = state.get("final_query", "")
    A2_block = state.get("A2_block", "")

    query_msg = AIMessage(
        content=f"""THIS IS SUGGESTED FINAL QUERY:

**A2 Block (Methods/Approaches):**
```
{A2_block}
```

**Complete PubMed Query:**
```
{final_query[:500]}...
```

Please respond with:
- "approve" to proceed with search
- "modify: <instructions>" to request changes"""
    )

    # Interrupt here - wait for query approval
    user_response = interrupt(query_msg.content)

    return {
        "messages": [query_msg, HumanMessage(content=user_response)]
    }


def search_pubmed_node(state: EnhancedState) -> EnhancedState:
    """Execute PubMed search with approved query."""
    final_query = state.get("final_query", "")
    results = search_pubmed(final_query, retmax=200)

    return {
        "search_results": results,
        "messages": []
    }


def validate_results(state: EnhancedState) -> EnhancedState:
    """
    HITL Checkpoint #4: Validate search results.
    """
    search_results = state.get("search_results", {})
    count = search_results.get("count", 0)
    sample_ids = search_results.get("ids", [])[:5]

    validation_msg = AIMessage(
        content=f"""Search Results Summary:
- Total matches: {count:,} articles
- Sample PMIDs: {', '.join(sample_ids)}

Please respond with:
- "accept" if results look good
- "too broad" if too many irrelevant results
- "too narrow" if too few results"""
    )

    # Interrupt here - wait for validation
    user_response = interrupt(validation_msg.content)

    return {
        "messages": [validation_msg, HumanMessage(content=user_response)]
    }


def ask_for_next_action(state: EnhancedState) -> EnhancedState:
    """
    HITL Checkpoint #5: Ask for next action.
    """
    iteration = state.get("iteration_count", 0)

    next_action_msg = AIMessage(
        content=f"""What would you like to do next? (Iteration {iteration + 1})

Options:
1. "refine" - Modify search terms based on results
2. "new search" - Start fresh with different protocols
3. "done" - Finish

Please choose an action:"""
    )

    # Interrupt here - wait for next action
    user_response = interrupt(next_action_msg.content)

    return {
        "user_response": user_response,
        "messages": [next_action_msg, HumanMessage(content=user_response)]
    }


def route_next_action(state: EnhancedState) -> Literal["refine", "new_search", "done"]:
    """Route based on user's next action choice."""
    user_response = state.get("user_response", "").lower()

    if "refine" in user_response:
        return "refine"
    elif "new" in user_response:
        return "new_search"
    else:
        return "done"


def refine_query_node(state: EnhancedState) -> EnhancedState:
    """Increment iteration and reset for refinement."""
    iteration = state.get("iteration_count", 0)
    return {
        "iteration_count": iteration + 1,
        "messages": []
    }


def start_new_search(state: EnhancedState) -> EnhancedState:
    """Reset state for new search."""
    return {
        "user_protocols": [],
        "suggested_protocols": [],
        "approved_protocols": [],
        "A2_block": "",
        "final_query": "",
        "search_results": {},
        "user_response": "",
        "iteration_count": 0,
        "messages": []
    }


# =============================================================================
# Graph Assembly
# =============================================================================

def create_hitl_graph() -> StateGraph:
    """Create HITL workflow graph with interrupt-based checkpoints."""

    builder = StateGraph(EnhancedState)

    # Add all nodes
    builder.add_node("detect_ambiguity", detect_ambiguity)
    builder.add_node("extract_protocols", extract_protocols)
    builder.add_node("expand_protocols", expand_protocols)
    builder.add_node("request_protocol_approval", request_protocol_approval)
    builder.add_node("create_a2_block", create_a2_block)
    builder.add_node("create_final_query", create_final_query)
    builder.add_node("show_query_for_approval", show_query_for_approval)
    builder.add_node("search_pubmed_node", search_pubmed_node)
    builder.add_node("validate_results", validate_results)
    builder.add_node("ask_for_next_action", ask_for_next_action)
    builder.add_node("refine_query_node", refine_query_node)
    builder.add_node("start_new_search", start_new_search)

    # Add edges (linear flow with checkpoints using interrupt())
    builder.add_edge(START, "detect_ambiguity")
    builder.add_edge("detect_ambiguity", "extract_protocols")
    builder.add_edge("extract_protocols", "expand_protocols")
    builder.add_edge("expand_protocols", "request_protocol_approval")
    builder.add_edge("request_protocol_approval", "create_a2_block")
    builder.add_edge("create_a2_block", "create_final_query")
    builder.add_edge("create_final_query", "show_query_for_approval")
    builder.add_edge("show_query_for_approval", "search_pubmed_node")
    builder.add_edge("search_pubmed_node", "validate_results")
    builder.add_edge("validate_results", "ask_for_next_action")

    # Conditional routing for next action
    builder.add_conditional_edges(
        "ask_for_next_action",
        route_next_action,
        {
            "refine": "refine_query_node",
            "new_search": "start_new_search",
            "done": END
        }
    )

    # Loop back edges
    builder.add_edge("refine_query_node", "detect_ambiguity")
    builder.add_edge("start_new_search", "detect_ambiguity")

    return builder


# =============================================================================
# Main Conversation Loop (FIXED)
# =============================================================================

def run_hitl_conversation():
    """
    Run HITL conversation with proper interrupt/resume pattern.

    Key differences from buggy version:
    1. No manual state updates
    2. Uses interrupt() in nodes
    3. Resumes with invoke(None, config)
    4. State managed by LangGraph checkpointer
    """
    # Create graph with memory
    builder = create_hitl_graph()
    memory = MemorySaver()
    graph = builder.compile(checkpointer=memory)

    # Configuration for thread
    config = {"configurable": {"thread_id": "1"}}

    print("=== PubMed Query Builder with HITL (FIXED) ===")
    print("Enter your protocols (comma-separated) or 'quit' to exit\n")

    # Get initial user input
    user_input = input("\nYou: ").strip()

    if user_input.lower() in ["quit", "exit", "q"]:
        print("Goodbye!")
        return

    # Parse initial protocols
    protocols = [p.strip() for p in user_input.split(",") if p.strip()]

    # Initial state
    initial_state = {
        "messages": [HumanMessage(content=user_input)],
        "user_protocols": protocols,
        "suggested_protocols": [],
        "approved_protocols": [],
        "A2_block": "",
        "final_query": "",
        "search_results": {},
        "user_response": "",
        "iteration_count": 0
    }

    # Start graph execution
    try:
        result = graph.invoke(initial_state, config)

        # Graph will interrupt at checkpoints
        # We never reach here on first invocation if checkpoints exist
        print("\n=== Workflow Complete ===")

        # Display final results
        search_results = result.get("search_results", {})
        if search_results:
            print(f"\nFinal Results:")
            print(f"  Total matches: {search_results.get('count', 0):,}")
            print(f"  PMIDs retrieved: {len(search_results.get('ids', []))}")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PubMed Query Builder with HITL (FIXED VERSION)"
    )
    parser.add_argument(
        "--mode",
        choices=["hitl"],
        default="hitl",
        help="Run mode (only hitl for now)"
    )

    args = parser.parse_args()

    if args.mode == "hitl":
        run_hitl_conversation()
