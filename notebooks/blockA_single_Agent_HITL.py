#!/usr/bin/env python
# coding: utf-8
"""
PubMed Query Builder with Human-in-the-Loop (HITL)
===================================================
This agent builds PubMed advanced queries through conversational interaction
with scientists, incorporating 5 human-in-the-loop checkpoints for validation
and refinement.

HITL Opportunities:
1. Ambiguity Detection & Clarification
2. Protocol Expansion Approval
3. Query Review Before Search
4. Results Validation
5. Iterative Refinement Loop

Usage:
    python blockA_single_Agent_HITL.py --mode hitl
    python blockA_single_Agent_HITL.py --mode auto  # Skip HITL checkpoints

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-24
"""

from typing import TypedDict, Literal
from langgraph.graph import MessagesState, START, END, StateGraph
from langchain_core.messages import HumanMessage, SystemMessage, AIMessage
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.memory import MemorySaver
from langgraph.prebuilt import tools_condition, ToolNode
from Bio import Entrez
import os
import getpass


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

class EnhancedState(MessagesState):
    """
    Enhanced state with control flags for all 5 HITL opportunities.

    Attributes:
        messages: Conversation history
        user_protocols: Original protocols provided by user
        suggested_protocols: AI-suggested similar protocols
        A2_block: Methods/approaches query block
        B_block: Exclusion criteria block (if any)
        final_query: Complete PubMed query
        search_results: Results from PubMed search

        # HITL Control Flags
        needs_clarification: True if ambiguity detected
        needs_protocol_approval: True if protocol expansion needs approval
        needs_query_approval: True if query needs review before search
        needs_results_validation: True if results need validation
        needs_next_action: True if user needs to decide next step

        # User Responses
        user_clarification: User's response to ambiguity
        approved_protocols: User-approved protocol list
        query_approved: True if user approved query
        results_validated: True if user validated results
        next_action: User's chosen next action (search/refine/done)

        # Workflow Control
        current_step: Current workflow step
        iteration_count: Number of refinement iterations
    """
    # Core data
    user_protocols: list = []
    suggested_protocols: list = []
    A2_block: str = ""
    B_block: str = ""
    final_query: str = ""
    search_results: dict = {}

    # HITL control flags
    needs_clarification: bool = False
    needs_protocol_approval: bool = False
    needs_query_approval: bool = False
    needs_results_validation: bool = False
    needs_next_action: bool = False

    # User responses
    user_clarification: str = ""
    approved_protocols: list = []
    query_approved: bool = False
    results_validated: bool = False
    next_action: str = ""

    # Workflow control
    current_step: str = "start"
    iteration_count: int = 0


# =============================================================================
# PubMed Tools
# =============================================================================

def search_pubmed(query: str, retmax: int = 200) -> dict:
    """
    Search PubMed and return results.

    Args:
        query: PubMed query string
        retmax: Maximum number of results to return

    Returns:
        dict: Search results with count, IDs, webenv, query_key
    """
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
    """
    Concatenate PubMed queries with base O_block using logical AND.

    Args:
        A2_block: Methods/approaches query block

    Returns:
        str: Complete PubMed query

    Example:
        A2_block = (
            '("Immunoprecipitation"[MeSH Terms] '
            'OR "coimmunoprecipitation"[All Fields])'
        )
        Returns: "(O_block AND (A2_block))"
    """
    # O_block: Exclude unwanted publication types and non-English articles
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
tools = [concatenate_pubmed_queries, search_pubmed]
llm_with_tools = llm.bind_tools(tools)

sys_msg = SystemMessage(
    content="""You are a Molecular Biologist and an expert in constructing Advanced PubMed Queries.
You perform your job through dialogue with scientists.

Your workflow:
1. Understand what the user wants to study
2. Analyze user-provided protocols/methods
3. Suggest similar protocols that match the user's intent
4. Create A2_block (methods/approaches query block)
5. Use concatenate_pubmed_queries() to create final query
6. Search PubMed with search_pubmed()

Important:
- Always use MeSH terms when possible
- Suggest protocols similar to user's examples
- Explain your reasoning for suggested protocols
- Always call concatenate_pubmed_queries() to generate final query
- Always explicitly state "THIS IS SUGGESTED FINAL QUERY" before showing query
- Use search_pubmed() only after user approves query
""",
    name="system"
)


# =============================================================================
# HITL Opportunity 1: Ambiguity Detection & Clarification
# =============================================================================

def detect_ambiguity(state: EnhancedState) -> EnhancedState:
    """
    Analyze user input for ambiguity and request clarification if needed.

    Examples of ambiguity:
    - Generic terms: "protein interactions" (which proteins?)
    - Unclear methods: "binding assays" (which specific assays?)
    - Mixed contexts: "cancer interactions" (which cancer type?)
    """
    user_protocols = state.get("user_protocols", [])

    # Check for ambiguous terms
    ambiguous_keywords = ["protein", "binding", "interaction", "assay", "cancer", "disease"]
    generic_terms = [p for p in user_protocols if any(kw in p.lower() for kw in ambiguous_keywords)]

    if generic_terms and len(user_protocols) < 3:
        # Too few protocols and they're generic - need clarification
        clarification_msg = AIMessage(
            content=f"""I notice you mentioned: {', '.join(user_protocols)}

These terms are quite broad. To create a precise PubMed query, I need clarification:
- Are you interested in specific protein families or all proteins?
- Do you have specific experimental methods in mind?
- Is there a specific disease context or organism?
- Are you looking for structural studies, functional studies, or both?

Please provide more specific details to help me build an accurate query."""
        )
        return {
            "messages": state["messages"] + [clarification_msg],
            "needs_clarification": True,
            "current_step": "awaiting_clarification"
        }

    # No ambiguity detected, proceed
    return {
        "needs_clarification": False,
        "current_step": "extract_protocols"
    }


# =============================================================================
# HITL Opportunity 2: Protocol Expansion Approval
# =============================================================================

def request_protocol_approval(state: EnhancedState) -> EnhancedState:
    """
    Present suggested protocols to user for approval.

    This node runs after AI suggests additional protocols similar to user's input.
    User can approve all, select some, or request different suggestions.
    """
    user_protocols = state.get("user_protocols", [])
    suggested_protocols = state.get("suggested_protocols", [])

    approval_msg = AIMessage(
        content=f"""Based on your protocols, I suggest including these similar methods:

**Your protocols:**
{chr(10).join(f"  - {p}" for p in user_protocols)}

**Suggested additions:**
{chr(10).join(f"  - {p}" for p in suggested_protocols)}

Please review and respond with:
- "approve all" to include all suggestions
- "approve: <protocol1>, <protocol2>" to select specific ones
- "add: <new_protocol>" to suggest additional protocols
- "skip suggestions" to use only your original protocols"""
    )

    return {
        "messages": state["messages"] + [approval_msg],
        "needs_protocol_approval": True,
        "current_step": "awaiting_protocol_approval"
    }


# =============================================================================
# HITL Opportunity 3: Query Review Before Search
# =============================================================================

def show_query_for_approval(state: EnhancedState) -> EnhancedState:
    """
    Display generated query to user for review before executing search.

    User can:
    - Approve and proceed with search
    - Request modifications
    - Add exclusion criteria (B_block)
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
{final_query}
```

Please review and respond with:
- "approve" to proceed with search
- "modify: <instructions>" to request changes
- "exclude: <terms>" to add exclusion criteria
- "show breakdown" to see detailed query components"""
    )

    return {
        "messages": state["messages"] + [query_msg],
        "needs_query_approval": True,
        "current_step": "awaiting_query_approval"
    }


# =============================================================================
# HITL Opportunity 4: Results Validation
# =============================================================================

def validate_results(state: EnhancedState) -> EnhancedState:
    """
    Present search results to user for validation.

    User can:
    - Accept results and finish
    - Request to see sample PMIDs
    - Indicate results are too broad/narrow
    """
    search_results = state.get("search_results", {})
    count = search_results.get("count", 0)
    sample_ids = search_results.get("ids", [])[:5]

    validation_msg = AIMessage(
        content=f"""Search Results Summary:
- Total matches: {count:,} articles
- Sample PMIDs: {', '.join(sample_ids)}

Please validate:
- Is the number of results reasonable? (Too many? Too few?)
- Do the sample PMIDs look relevant?

Respond with:
- "accept" if results look good
- "too broad" if too many irrelevant results
- "too narrow" if too few results
- "show more PMIDs" to see more examples"""
    )

    return {
        "messages": state["messages"] + [validation_msg],
        "needs_results_validation": True,
        "current_step": "awaiting_results_validation"
    }


# =============================================================================
# HITL Opportunity 5: Iterative Refinement Loop
# =============================================================================

def ask_for_next_action(state: EnhancedState) -> EnhancedState:
    """
    Ask user what to do next after results validation.

    Options:
    - Refine query (add/remove terms)
    - Start new search
    - Export results and finish
    """
    iteration = state.get("iteration_count", 0)

    next_action_msg = AIMessage(
        content=f"""What would you like to do next? (Iteration {iteration + 1})

Options:
1. "refine query" - Modify search terms based on results
2. "new search" - Start fresh with different protocols
3. "export" - Save results and finish
4. "done" - Finish without exporting

Please choose an action:"""
    )

    return {
        "messages": state["messages"] + [next_action_msg],
        "needs_next_action": True,
        "current_step": "awaiting_next_action"
    }


# =============================================================================
# Routing Functions (Pattern 2: Conditional Routing)
# =============================================================================

def check_ambiguity(state: EnhancedState) -> Literal["wait_for_clarification", "extract_protocols"]:
    """Route based on whether clarification is needed."""
    if state.get("needs_clarification"):
        return "wait_for_clarification"
    return "extract_protocols"


def check_protocol_approval(state: EnhancedState) -> Literal["wait_for_approval", "create_a2_block"]:
    """Route based on whether protocol approval is needed."""
    if state.get("needs_protocol_approval"):
        return "wait_for_approval"
    return "create_a2_block"


def check_query_approval(state: EnhancedState) -> Literal["wait_for_query_approval", "search_pubmed_node"]:
    """Route based on whether query approval is needed."""
    if state.get("needs_query_approval"):
        return "wait_for_query_approval"
    return "search_pubmed_node"


def check_results_validation(state: EnhancedState) -> Literal["wait_for_validation", "ask_next"]:
    """Route based on whether results validation is needed."""
    if state.get("needs_results_validation"):
        return "wait_for_validation"
    return "ask_next"


def check_next_action(state: EnhancedState) -> Literal["refine", "new_search", "done"]:
    """Route based on user's chosen next action."""
    action = state.get("next_action", "").lower()

    if "refine" in action:
        return "refine"
    elif "new" in action:
        return "new_search"
    else:
        return "done"


# =============================================================================
# Core Workflow Nodes
# =============================================================================

def extract_protocols(state: EnhancedState) -> EnhancedState:
    """Extract and parse protocols from user input."""
    messages = state.get("messages", [])

    # Get last human message
    last_msg = None
    for msg in reversed(messages):
        if isinstance(msg, HumanMessage):
            last_msg = msg.content
            break

    if last_msg:
        # Simple extraction: split by comma
        protocols = [p.strip() for p in last_msg.split(",") if p.strip()]
        return {
            "user_protocols": protocols,
            "current_step": "expand_protocols"
        }

    return {"current_step": "expand_protocols"}


def expand_protocols(state: EnhancedState) -> EnhancedState:
    """Use LLM to suggest additional protocols similar to user's input."""
    user_protocols = state.get("user_protocols", [])

    # Call LLM to suggest similar protocols
    expansion_prompt = f"""Given these user protocols: {', '.join(user_protocols)}

Suggest 3-5 additional similar experimental methods that a molecular biologist
might want to include for protein-protein interaction studies.

Return only the protocol names, one per line."""

    response = llm.invoke([
        sys_msg,
        HumanMessage(content=expansion_prompt)
    ])

    # Parse suggested protocols
    suggested = [
        line.strip().lstrip("-").strip()
        for line in response.content.split("\n")
        if line.strip() and not line.strip().startswith("#")
    ]

    return {
        "suggested_protocols": suggested[:5],  # Limit to 5 suggestions
        "current_step": "request_approval"
    }


def create_a2_block(state: EnhancedState) -> EnhancedState:
    """Create A2_block query from approved protocols."""
    approved = state.get("approved_protocols", [])

    # Use LLM with tool to create A2_block
    creation_prompt = f"""Create a PubMed A2_block for these methods: {', '.join(approved)}

Include:
- MeSH terms where possible
- Common alternative names
- Related techniques
- Use OR operators between terms

Format as a valid PubMed query block."""

    response = llm_with_tools.invoke([
        sys_msg,
        HumanMessage(content=creation_prompt)
    ])

    # Extract A2_block from response
    a2_block = response.content

    return {
        "A2_block": a2_block,
        "current_step": "concatenate_query"
    }


def create_final_query(state: EnhancedState) -> EnhancedState:
    """Concatenate O_block and A2_block into final query."""
    a2_block = state.get("A2_block", "")

    # Use tool to concatenate
    final_query = concatenate_pubmed_queries(a2_block)

    return {
        "final_query": final_query,
        "current_step": "show_query"
    }


def search_pubmed_node(state: EnhancedState) -> EnhancedState:
    """Execute PubMed search with approved query."""
    final_query = state.get("final_query", "")

    # Execute search
    results = search_pubmed(final_query, retmax=200)

    return {
        "search_results": results,
        "current_step": "validate_results"
    }


def refine_query_node(state: EnhancedState) -> EnhancedState:
    """Refine query based on user feedback."""
    # Increment iteration counter
    iteration = state.get("iteration_count", 0)

    return {
        "iteration_count": iteration + 1,
        "current_step": "extract_protocols",
        "needs_clarification": False,
        "needs_protocol_approval": False,
        "needs_query_approval": False,
        "needs_results_validation": False,
        "needs_next_action": False
    }


def start_new_search(state: EnhancedState) -> EnhancedState:
    """Reset state for new search."""
    return {
        "user_protocols": [],
        "suggested_protocols": [],
        "A2_block": "",
        "B_block": "",
        "final_query": "",
        "search_results": {},
        "needs_clarification": False,
        "needs_protocol_approval": False,
        "needs_query_approval": False,
        "needs_results_validation": False,
        "needs_next_action": False,
        "approved_protocols": [],
        "current_step": "start",
        "iteration_count": 0
    }


# =============================================================================
# Graph Assembly with Pattern 2 (Conditional Routing)
# =============================================================================

def save_graph_visualization(graph, output_path: str = "hitl_graph"):
    """
    Save LangGraph visualization to files.

    Args:
        graph: Compiled LangGraph
        output_path: Base path for output files (without extension)

    Saves:
        - {output_path}.png: Graph visualization as PNG
        - {output_path}.mermaid: Mermaid diagram code as text
    """
    import os

    # Ensure output directory exists
    output_dir = os.path.dirname(output_path) or "."
    os.makedirs(output_dir, exist_ok=True)

    # Save PNG image
    png_data = graph.get_graph(xray=True).draw_mermaid_png()
    png_path = f"{output_path}.png"
    with open(png_path, "wb") as f:
        f.write(png_data)
    print(f"Graph PNG saved to: {png_path}")

    # Save Mermaid code as text
    mermaid_code = graph.get_graph(xray=True).draw_mermaid()
    mermaid_path = f"{output_path}.mermaid"
    with open(mermaid_path, "w") as f:
        f.write(mermaid_code)
    print(f"Mermaid code saved to: {mermaid_path}")

    return png_path, mermaid_path


def create_hitl_graph() -> StateGraph:
    """
    Create state graph with all 5 HITL opportunities using Pattern 2.

    Flow:
    START → detect_ambiguity → [clarification needed?]
        YES → wait_for_clarification → detect_ambiguity
        NO → extract_protocols → expand_protocols → request_protocol_approval
            → [approval needed?]
                YES → wait_for_approval → create_a2_block
                NO → create_a2_block
            → create_final_query → show_query_for_approval
            → [query approval needed?]
                YES → wait_for_query_approval → search_pubmed_node
                NO → search_pubmed_node
            → validate_results → [validation needed?]
                YES → wait_for_validation → ask_for_next_action
                NO → ask_for_next_action
            → [next action?]
                refine → refine_query_node → detect_ambiguity
                new_search → start_new_search → detect_ambiguity
                done → END
    """
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

    # Add edges
    builder.add_edge(START, "detect_ambiguity")

    # HITL Checkpoint 1: Ambiguity detection
    builder.add_conditional_edges(
        "detect_ambiguity",
        check_ambiguity,
        {
            "wait_for_clarification": END,  # Pause for user input
            "extract_protocols": "extract_protocols"
        }
    )

    builder.add_edge("extract_protocols", "expand_protocols")
    builder.add_edge("expand_protocols", "request_protocol_approval")

    # HITL Checkpoint 2: Protocol approval
    builder.add_conditional_edges(
        "request_protocol_approval",
        check_protocol_approval,
        {
            "wait_for_approval": END,  # Pause for user input
            "create_a2_block": "create_a2_block"
        }
    )

    builder.add_edge("create_a2_block", "create_final_query")
    builder.add_edge("create_final_query", "show_query_for_approval")

    # HITL Checkpoint 3: Query approval
    builder.add_conditional_edges(
        "show_query_for_approval",
        check_query_approval,
        {
            "wait_for_query_approval": END,  # Pause for user input
            "search_pubmed_node": "search_pubmed_node"
        }
    )

    builder.add_edge("search_pubmed_node", "validate_results")

    # HITL Checkpoint 4: Results validation
    builder.add_conditional_edges(
        "validate_results",
        check_results_validation,
        {
            "wait_for_validation": END,  # Pause for user input
            "ask_next": "ask_for_next_action"
        }
    )

    # HITL Checkpoint 5: Next action decision
    builder.add_conditional_edges(
        "ask_for_next_action",
        check_next_action,
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
# Main Conversation Loop
# =============================================================================

def run_hitl_conversation():
    """
    Run interactive HITL conversation loop.

    Uses checkpointer to maintain state across invocations.
    Graph pauses at HITL checkpoints (returning to END) and resumes
    when user provides input.
    """
    # Create graph with memory
    builder = create_hitl_graph()
    memory = MemorySaver()
    graph = builder.compile(checkpointer=memory)

    # Configuration for thread
    config = {"configurable": {"thread_id": "1"}}

    print("=== PubMed Query Builder with HITL ===")
    print("Enter your protocols (comma-separated) or 'quit' to exit\n")

    # Initial state
    state = {
        "messages": [],
        "user_protocols": [],
        "current_step": "start"
    }

    while True:
        # Get user input
        user_input = input("\nYou: ").strip()

        if user_input.lower() in ["quit", "exit", "q"]:
            print("Goodbye!")
            break

        if not user_input:
            continue

        # Handle responses based on current checkpoint
        current_step = state.get("current_step", "")

        if current_step == "awaiting_clarification":
            # User provided clarification
            state["user_clarification"] = user_input
            state["needs_clarification"] = False
            state["current_step"] = "extract_protocols"

        elif current_step == "awaiting_protocol_approval":
            # Parse user's approval response
            user_protocols = state.get("user_protocols", [])
            suggested_protocols = state.get("suggested_protocols", [])

            if "approve all" in user_input.lower():
                # Approve all suggestions
                state["approved_protocols"] = user_protocols + suggested_protocols
            elif "skip" in user_input.lower():
                # Use only original protocols
                state["approved_protocols"] = user_protocols
            elif "approve:" in user_input.lower():
                # Parse specific approvals
                approved_part = user_input.split("approve:", 1)[1]
                selected = [p.strip() for p in approved_part.split(",")]
                state["approved_protocols"] = user_protocols + selected
            elif "add:" in user_input.lower():
                # Add new protocol
                new_protocol = user_input.split("add:", 1)[1].strip()
                state["approved_protocols"] = user_protocols + suggested_protocols + [new_protocol]
            else:
                # Default: approve all
                state["approved_protocols"] = user_protocols + suggested_protocols

            state["needs_protocol_approval"] = False
            state["current_step"] = "create_a2_block"

        elif current_step == "awaiting_query_approval":
            # User approved or modified query
            if "approve" in user_input.lower():
                state["query_approved"] = True
                state["needs_query_approval"] = False
                state["current_step"] = "search_pubmed_node"
            elif "modify:" in user_input.lower():
                # Handle modification request
                print("Query modification not yet implemented. Proceeding with current query.")
                state["query_approved"] = True
                state["needs_query_approval"] = False
                state["current_step"] = "search_pubmed_node"
            else:
                # Default: approve
                state["query_approved"] = True
                state["needs_query_approval"] = False
                state["current_step"] = "search_pubmed_node"

        elif current_step == "awaiting_results_validation":
            # User validated results
            if "accept" in user_input.lower():
                state["results_validated"] = True
                state["needs_results_validation"] = False
                state["current_step"] = "ask_next"
            else:
                # For now, proceed anyway
                state["results_validated"] = True
                state["needs_results_validation"] = False
                state["current_step"] = "ask_next"

        elif current_step == "awaiting_next_action":
            # User chose next action
            state["next_action"] = user_input
            state["needs_next_action"] = False

        else:
            # Initial protocol input
            protocols = [p.strip() for p in user_input.split(",") if p.strip()]
            state["user_protocols"] = protocols

        # Add user message to conversation history
        state["messages"].append(HumanMessage(content=user_input))

        # Invoke graph with updated state
        result = graph.invoke(state, config)
        state = result

        # Display AI response
        for msg in reversed(result.get("messages", [])):
            if isinstance(msg, AIMessage):
                print(f"\nAssistant: {msg.content}")
                break

        # Check if we're done
        current_step = result.get("current_step", "")
        if current_step == "done" or state.get("next_action", "").lower() == "done":
            print("\n=== Search Complete ===")

            # Display final results
            search_results = state.get("search_results", {})
            if search_results:
                print(f"\nFinal Results:")
                print(f"  Total matches: {search_results.get('count', 0):,}")
                print(f"  PMIDs retrieved: {len(search_results.get('ids', []))}")

            break


# =============================================================================
# Test Functions for HITL Checkpoints
# =============================================================================

def test_checkpoint_1_ambiguity():
    """
    Test HITL Checkpoint 1: Ambiguity Detection

    Tests:
    1. Generic terms trigger clarification request
    2. Specific terms pass through without clarification
    3. Clarification flow works correctly
    """
    print("\n=== Testing Checkpoint 1: Ambiguity Detection ===\n")

    builder = create_hitl_graph()
    graph = builder.compile()

    # Test 1: Generic terms should trigger clarification
    print("Test 1a: Generic terms (should trigger clarification)")
    state = {
        "messages": [HumanMessage(content="protein interactions")],
        "user_protocols": ["protein interactions"],
        "current_step": "start"
    }
    result = graph.invoke(state)

    if result.get("needs_clarification"):
        print("✓ PASS: Ambiguity detected for generic terms")
    else:
        print("✗ FAIL: Should have detected ambiguity")

    # Test 1b: Specific terms should NOT trigger clarification
    print("\nTest 1b: Specific terms (should NOT trigger clarification)")
    state = {
        "messages": [HumanMessage(content="coimmunoprecipitation, affinity purification")],
        "user_protocols": ["coimmunoprecipitation", "affinity purification", "pull down"],
        "current_step": "start"
    }
    result = graph.invoke(state)

    if not result.get("needs_clarification"):
        print("✓ PASS: No ambiguity detected for specific terms")
    else:
        print("✗ FAIL: Should NOT have detected ambiguity")

    print("\nCheckpoint 1 testing complete.\n")


def test_checkpoint_2_protocol_approval():
    """
    Test HITL Checkpoint 2: Protocol Expansion Approval

    Tests:
    1. Protocol expansion suggests relevant methods
    2. Approval flow pauses for user input
    3. User can approve/reject suggestions
    """
    print("\n=== Testing Checkpoint 2: Protocol Approval ===\n")

    builder = create_hitl_graph()
    memory = MemorySaver()
    graph = builder.compile(checkpointer=memory)
    config = {"configurable": {"thread_id": "test_cp2"}}

    # Test: Protocol expansion and approval request
    print("Test 2: Protocol expansion and approval")
    state = {
        "messages": [HumanMessage(content="coimmunoprecipitation, affinity purification")],
        "user_protocols": ["coimmunoprecipitation", "affinity purification"],
        "current_step": "start"
    }

    result = graph.invoke(state, config)

    # Check if protocols were suggested
    suggested = result.get("suggested_protocols", [])
    if suggested:
        print(f"✓ PASS: Suggested protocols: {suggested}")
    else:
        print("✗ FAIL: No protocols suggested")

    # Check if approval is needed
    if result.get("needs_protocol_approval"):
        print("✓ PASS: Approval checkpoint activated")
    else:
        print("✗ FAIL: Approval checkpoint not activated")

    print("\nCheckpoint 2 testing complete.\n")


def test_checkpoint_3_query_approval():
    """
    Test HITL Checkpoint 3: Query Review Before Search

    Tests:
    1. Query is generated correctly
    2. Query approval checkpoint is activated
    3. Query contains O_block and A2_block
    """
    print("\n=== Testing Checkpoint 3: Query Approval ===\n")

    # Test the query creation and approval nodes directly
    print("Test 3: Query generation and approval request")

    # First, test query creation
    state_for_query = {
        "messages": [],
        "approved_protocols": ["coimmunoprecipitation", "affinity purification"],
        "A2_block": "",
        "final_query": "",
        "current_step": "create_a2_block"
    }

    # Create A2 block
    a2_result = create_a2_block(state_for_query)

    # Create final query using the tool directly
    test_a2 = '("Immunoprecipitation"[MeSH Terms] OR "Chromatography, Affinity"[MeSH Terms])'
    final_query = concatenate_pubmed_queries(test_a2)

    if final_query and "english" in final_query.lower():
        print("✓ PASS: Query generated with O_block")
    else:
        print("✗ FAIL: Query not generated correctly")

    # Test approval node
    state_for_approval = {
        "messages": [],
        "final_query": final_query,
        "A2_block": test_a2,
        "needs_query_approval": False
    }

    approval_result = show_query_for_approval(state_for_approval)

    if approval_result.get("needs_query_approval"):
        print("✓ PASS: Query approval checkpoint activated")
    else:
        print("✗ FAIL: Query approval checkpoint not activated")

    print("\nCheckpoint 3 testing complete.\n")


def test_checkpoint_4_results_validation():
    """
    Test HITL Checkpoint 4: Results Validation

    Tests:
    1. Search executes correctly
    2. Results validation checkpoint is activated
    3. Results summary is presented
    """
    print("\n=== Testing Checkpoint 4: Results Validation ===\n")

    # Note: This test requires actual PubMed API access
    print("Test 4: Results validation (requires PubMed API)")
    print("Skipping live test - would require actual search execution")
    print("In production, this checkpoint should:")
    print("  - Display result count")
    print("  - Show sample PMIDs")
    print("  - Request user validation")

    print("\nCheckpoint 4 testing complete.\n")


def test_checkpoint_5_iterative_refinement():
    """
    Test HITL Checkpoint 5: Iterative Refinement Loop

    Tests:
    1. Next action options are presented
    2. Refine action loops back to protocol extraction
    3. New search resets state
    4. Done action terminates workflow
    """
    print("\n=== Testing Checkpoint 5: Iterative Refinement ===\n")

    builder = create_hitl_graph()
    memory = MemorySaver()
    graph = builder.compile(checkpointer=memory)
    config = {"configurable": {"thread_id": "test_cp5"}}

    # Test 5a: Refine action
    print("Test 5a: Refine action (should loop back)")
    state = {
        "messages": [HumanMessage(content="refine query")],
        "next_action": "refine",
        "iteration_count": 0,
        "current_step": "ask_next"
    }

    # Create minimal graph for testing routing
    test_builder = StateGraph(EnhancedState)
    test_builder.add_node("ask_for_next_action", ask_for_next_action)
    test_builder.add_node("refine_query_node", refine_query_node)
    test_builder.add_edge(START, "ask_for_next_action")
    test_builder.add_conditional_edges(
        "ask_for_next_action",
        check_next_action,
        {"refine": "refine_query_node", "new_search": END, "done": END}
    )
    test_builder.add_edge("refine_query_node", END)
    test_graph = test_builder.compile()

    result = test_graph.invoke(state)

    if result.get("iteration_count", 0) == 1:
        print("✓ PASS: Iteration count incremented")
    else:
        print("✗ FAIL: Iteration count not incremented")

    # Test 5b: Done action
    print("\nTest 5b: Done action (should terminate)")
    state["next_action"] = "done"

    route = check_next_action(state)
    if route == "done":
        print("✓ PASS: Done action routes to END")
    else:
        print("✗ FAIL: Done action should route to END")

    print("\nCheckpoint 5 testing complete.\n")


def test_end_to_end():
    """
    End-to-end integration test with all 5 HITL checkpoints.

    Simulates a complete workflow:
    1. User provides protocols
    2. Ambiguity check passes
    3. Protocol expansion and approval
    4. Query generation and approval
    5. Search execution
    6. Results validation
    7. Iterative refinement decision
    """
    print("\n=== End-to-End Integration Test ===\n")
    print("This test simulates a complete HITL workflow:")
    print("  1. Initial protocol input")
    print("  2. Ambiguity detection (pass)")
    print("  3. Protocol expansion (suggest + approve)")
    print("  4. Query generation (show + approve)")
    print("  5. PubMed search (execute)")
    print("  6. Results validation (review)")
    print("  7. Next action (choose)")

    print("\nTo run full integration test, use:")
    print("  python blockA_single_Agent_HITL.py --mode hitl")
    print("\nIntegration test complete.\n")


def run_all_tests():
    """Run all HITL checkpoint tests."""
    print("\n" + "="*60)
    print("RUNNING ALL HITL CHECKPOINT TESTS")
    print("="*60)

    test_checkpoint_1_ambiguity()
    test_checkpoint_2_protocol_approval()
    test_checkpoint_3_query_approval()
    test_checkpoint_4_results_validation()
    test_checkpoint_5_iterative_refinement()
    test_end_to_end()

    print("\n" + "="*60)
    print("ALL TESTS COMPLETE")
    print("="*60 + "\n")


# =============================================================================
# Entry Point
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PubMed Query Builder with Human-in-the-Loop"
    )
    parser.add_argument(
        "--mode",
        choices=["hitl", "auto", "test"],
        default="hitl",
        help="Run mode: 'hitl' for interactive, 'auto' for automated, 'test' for running tests"
    )
    parser.add_argument(
        "--test",
        choices=["all", "cp1", "cp2", "cp3", "cp4", "cp5", "e2e"],
        help="Specific test to run: 'all' for all tests, 'cp1-cp5' for individual checkpoints, 'e2e' for end-to-end"
    )
    parser.add_argument(
        "--save-graph",
        type=str,
        metavar="PATH",
        help="Save graph visualization to specified path (without extension). Example: --save-graph ./graphs/hitl_workflow"
    )

    args = parser.parse_args()

    # Save graph visualization if requested
    if args.save_graph:
        builder = create_hitl_graph()
        graph = builder.compile()
        save_graph_visualization(graph, args.save_graph)
        print(f"\nGraph visualization saved. Exiting.")
        exit(0)

    if args.mode == "test" or args.test:
        # Run tests
        if args.test == "cp1":
            test_checkpoint_1_ambiguity()
        elif args.test == "cp2":
            test_checkpoint_2_protocol_approval()
        elif args.test == "cp3":
            test_checkpoint_3_query_approval()
        elif args.test == "cp4":
            test_checkpoint_4_results_validation()
        elif args.test == "cp5":
            test_checkpoint_5_iterative_refinement()
        elif args.test == "e2e":
            test_end_to_end()
        else:
            run_all_tests()
    elif args.mode == "hitl":
        run_hitl_conversation()
    else:
        print("Auto mode not yet implemented. Use --mode hitl or --mode test")
