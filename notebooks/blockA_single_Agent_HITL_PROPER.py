#!/usr/bin/env python3
"""
HITL Workflow - Proper Implementation with Separate Graph Invocations

This implementation uses separate graphs for each workflow stage, connected by
simple Python orchestration with CLI input() calls. This approach:
- Works perfectly with CLI
- No complex state management
- Easy to test each stage independently
- Production-ready

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-01-11
"""

import os
import json
import warnings
from typing import TypedDict, List, Optional, Dict, Any
from datetime import datetime

# LangGraph imports
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

# LangChain imports
from langchain_openai import ChatOpenAI
from langchain_core.messages import AIMessage, HumanMessage, SystemMessage

# PubMed tools
from Bio import Entrez

# Suppress warnings
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

# Set your email for NCBI Entrez
Entrez.email = "ekwame001@gmail.com"

# API Key setup
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
if not OPENAI_API_KEY:
    raise ValueError("OPENAI_API_KEY environment variable not set")

# Initialize LLM
llm = ChatOpenAI(
    model="gpt-4",
    temperature=0.3,
    openai_api_key=OPENAI_API_KEY
)

# O_block: Base filter to exclude unwanted publication types and non-English articles
O_BLOCK = (
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

# =============================================================================
# State Definitions (Minimal, focused states for each stage)
# =============================================================================

class ProtocolExpansionInput(TypedDict):
    """State for Stage 1: Protocol Expansion Graph"""
    user_protocols: List[str]
    suggested_protocols: Optional[List[str]]
    expansion_reasoning: Optional[str]


class QueryCreationInput(TypedDict):
    """State for Stage 2: Query Creation Graph"""
    approved_protocols: List[str]
    a2_block: Optional[str]
    final_query: Optional[str]
    query_components: Optional[Dict[str, str]]


class SearchInput(TypedDict):
    """State for Stage 3: PubMed Search Graph"""
    final_query: str
    search_results: Optional[Dict[str, Any]]
    result_count: Optional[int]
    articles: Optional[List[Dict[str, str]]]


class SessionState(TypedDict):
    """Session state for multi-iteration workflows"""
    iteration: int
    protocol_history: List[List[str]]
    query_history: List[str]
    search_history: List[Dict[str, Any]]
    timestamp: str


# =============================================================================
# Stage 1: Protocol Expansion Graph
# =============================================================================

def expand_protocols_node(state: ProtocolExpansionInput) -> ProtocolExpansionInput:
    """
    Node: Expand user protocols by suggesting similar methods.
    Uses LLM to suggest at least 10 related protein-protein interaction techniques.
    """
    user_protocols = state["user_protocols"]

    # Validate input
    if not user_protocols or len(user_protocols) == 0:
        return {
            "user_protocols": user_protocols,
            "suggested_protocols": [],
            "expansion_reasoning": "No protocols provided for expansion"
        }

    try:
        prompt = f"""You are a biomedical research expert specializing in protein-protein interaction (PPI) studies.

The user has provided these experimental protocols:
{', '.join(user_protocols)}

Please suggest 3-5 additional similar or complementary methods commonly used for PPI studies that the user might want to include in their literature search.

Focus on widely-used, well-established techniques. Provide brief reasoning for your suggestions.

Format your response as:
SUGGESTED PROTOCOLS:
- [protocol 1]
- [protocol 2]
- [protocol 3]

REASONING:
[Brief explanation of why these protocols are relevant]
"""

        messages = [
            SystemMessage(content="You are a biomedical research expert."),
            HumanMessage(content=prompt)
        ]

        response = llm.invoke(messages)
        response_text = response.content

        # Parse suggested protocols
        suggested = []
        reasoning = ""

        if "SUGGESTED PROTOCOLS:" in response_text:
            parts = response_text.split("REASONING:")
            protocol_section = parts[0].replace("SUGGESTED PROTOCOLS:", "").strip()
            reasoning = parts[1].strip() if len(parts) > 1 else ""

            # Extract protocols from bullet points
            for line in protocol_section.split("\n"):
                line = line.strip()
                if line.startswith("-"):
                    protocol = line[1:].strip()
                    if protocol:
                        suggested.append(protocol)

        return {
            "user_protocols": user_protocols,
            "suggested_protocols": suggested,
            "expansion_reasoning": reasoning
        }

    except Exception as e:
        print(f"⚠️ Error during protocol expansion: {e}")
        return {
            "user_protocols": user_protocols,
            "suggested_protocols": [],
            "expansion_reasoning": f"Error during expansion: {str(e)}"
        }


def format_suggestions_node(state: ProtocolExpansionInput) -> ProtocolExpansionInput:
    """
    Node: Format expansion results for display.
    This is a simple transformation node that could add metadata or logging.
    """
    # For now, just pass through (could add logging, validation, etc.)
    return state


def create_protocol_expansion_graph() -> StateGraph:
    """Create the graph for Stage 1: Protocol Expansion"""
    builder = StateGraph(ProtocolExpansionInput)

    # Add nodes
    builder.add_node("expand_protocols", expand_protocols_node)
    builder.add_node("format_suggestions", format_suggestions_node)

    # Add edges
    builder.add_edge(START, "expand_protocols")
    builder.add_edge("expand_protocols", "format_suggestions")
    builder.add_edge("format_suggestions", END)

    return builder


# =============================================================================
# Stage 2: Query Creation Graph
# =============================================================================

def create_a2_block_node(state: QueryCreationInput) -> QueryCreationInput:
    """
    Node: Create A2 block (methods section) from approved protocols.
    """
    approved_protocols = state["approved_protocols"]

    # Validate input
    if not approved_protocols or len(approved_protocols) == 0:
        return {
            "approved_protocols": approved_protocols,
            "a2_block": "",
            "final_query": state.get("final_query"),
            "query_components": state.get("query_components")
        }

    try:
        prompt = f"""You are creating a PubMed advanced search query.

Create the A2 block (methods/techniques section) for these protocols:
{', '.join(approved_protocols)}

Use proper PubMed search syntax with:
- [Title/Abstract] tags for each method
- OR operators between methods
- Parentheses for grouping

Example format:
("co-immunoprecipitation"[Title/Abstract] OR "affinity purification"[Title/Abstract] OR "pull-down assay"[Title/Abstract])

Provide ONLY the A2 block query text, no explanation.
"""

        messages = [
            SystemMessage(content="You are a PubMed query expert."),
            HumanMessage(content=prompt)
        ]

        response = llm.invoke(messages)
        a2_block = response.content.strip()

        return {
            "approved_protocols": approved_protocols,
            "a2_block": a2_block,
            "final_query": state.get("final_query"),
            "query_components": state.get("query_components")
        }

    except Exception as e:
        print(f"⚠️ Error during A2 block creation: {e}")
        return {
            "approved_protocols": approved_protocols,
            "a2_block": "",
            "final_query": state.get("final_query"),
            "query_components": state.get("query_components")
        }


def create_final_query_node(state: QueryCreationInput) -> QueryCreationInput:
    """
    Node: Combine A2 block with O_block to create final query.
    Format: (O_block AND (A2_block))
    """
    a2_block = state["a2_block"]

    # Validate input
    if not a2_block or len(a2_block.strip()) == 0:
        print("⚠️ Warning: A2 block is empty. Query creation may fail.")
        a2_block = ""

    # Combine O_block (base filters) AND A2_block (methods)
    # Format: (O_block AND (A2_block))
    if a2_block:
        final_query = f"({O_BLOCK} AND ({a2_block}))"
    else:
        final_query = O_BLOCK

    query_components = {
        "a2_block": a2_block,
        "o_block": O_BLOCK,  # Base publication type and language filters
        "b_block": None   # Could add exclusion criteria later
    }

    return {
        "approved_protocols": state["approved_protocols"],
        "a2_block": a2_block,
        "final_query": final_query,
        "query_components": query_components
    }


def create_query_creation_graph() -> StateGraph:
    """Create the graph for Stage 2: Query Creation"""
    builder = StateGraph(QueryCreationInput)

    # Add nodes
    builder.add_node("create_a2_block", create_a2_block_node)
    builder.add_node("create_final_query", create_final_query_node)

    # Add edges
    builder.add_edge(START, "create_a2_block")
    builder.add_edge("create_a2_block", "create_final_query")
    builder.add_edge("create_final_query", END)

    return builder


# =============================================================================
# Stage 3: PubMed Search Graph
# =============================================================================

def search_pubmed_node(state: SearchInput) -> SearchInput:
    """
    Node: Execute PubMed search with the final query.
    """
    final_query = state["final_query"]

    try:
        # Search PubMed
        handle = Entrez.esearch(
            db="pubmed",
            term=final_query,
            retmax=100,
            sort="relevance"
        )
        search_results = Entrez.read(handle)
        handle.close()

        result_count = int(search_results["Count"])
        id_list = search_results["IdList"]

        # Fetch article details for top 10 results
        articles = []
        if id_list:
            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=id_list[:10],
                rettype="medline",
                retmode="text"
            )
            fetch_results = fetch_handle.read()
            fetch_handle.close()

            # Parse results (simplified - just store raw for now)
            # In production, would parse into structured article objects
            articles = [{"pmid": pmid, "raw": fetch_results} for pmid in id_list[:10]]

        return {
            "final_query": final_query,
            "search_results": search_results,
            "result_count": result_count,
            "articles": articles
        }

    except Exception as e:
        print(f"Error during PubMed search: {e}")
        return {
            "final_query": final_query,
            "search_results": {"Error": str(e)},
            "result_count": 0,
            "articles": []
        }


def create_search_graph() -> StateGraph:
    """Create the graph for Stage 3: PubMed Search"""
    builder = StateGraph(SearchInput)

    # Add nodes
    builder.add_node("search_pubmed", search_pubmed_node)

    # Add edges
    builder.add_edge(START, "search_pubmed")
    builder.add_edge("search_pubmed", END)

    return builder


# =============================================================================
# Helper Functions for HITL Checkpoints
# =============================================================================

def detect_ambiguity(protocols: List[str]) -> bool:
    """
    HITL Checkpoint #1: Detect ambiguous protocol input.
    Returns True if input needs clarification.
    """
    # Check for generic keywords
    generic_keywords = [
        "protein interaction", "ppi", "interactions",
        "assay", "experiment", "method", "technique"
    ]

    # Check if user provided too few protocols
    if len(protocols) < 3:
        # Check if any protocol contains generic keywords
        for protocol in protocols:
            protocol_lower = protocol.lower()
            if any(keyword in protocol_lower for keyword in generic_keywords):
                return True

    return False


def display_expansion_results(user_protocols: List[str], suggested_protocols: List[str], reasoning: str):
    """Display protocol expansion results in a formatted way."""
    print("\n" + "="*80)
    print("PROTOCOL EXPANSION RESULTS")
    print("="*80)

    print("\nYour protocols:")
    for i, p in enumerate(user_protocols, 1):
        print(f"  {i}. {p}")

    print("\nSuggested additional protocols:")
    for i, p in enumerate(suggested_protocols, 1):
        print(f"  {i}. {p}")

    print(f"\nReasoning:\n{reasoning}")
    print("="*80)


def display_query(query: str, components: Dict[str, str]):
    """Display the created query in a formatted way."""
    print("\n" + "="*80)
    print("QUERY CREATED")
    print("="*80)

    print("\nO Block (Base Filters - excludes reviews, clinical trials, non-English):")
    print(f"  {components['o_block']}")

    print("\nA2 Block (Methods):")
    print(f"  {components['a2_block']}")

    print("\nFinal Query (O_block AND A2_block):")
    print(f"  {query}")
    print("="*80)


def display_search_results(result_count: int, articles: List[Dict[str, str]]):
    """Display search results in a formatted way."""
    print("\n" + "="*80)
    print("SEARCH RESULTS")
    print("="*80)

    print(f"\nTotal results found: {result_count}")

    if articles:
        print(f"\nTop {len(articles)} articles:")
        for i, article in enumerate(articles, 1):
            print(f"  {i}. PMID: {article['pmid']}")
    else:
        print("\nNo articles retrieved.")

    print("="*80)


# =============================================================================
# Main CLI Orchestrator
# =============================================================================

def run_hitl_workflow():
    """
    Main CLI orchestrator that connects all stages with HITL checkpoints.

    This function:
    1. Creates the three stage graphs
    2. Runs them in sequence
    3. Inserts HITL checkpoints between stages
    4. Supports iterative refinement
    """
    print("="*80)
    print("HITL PUBMED QUERY BUILDER - PROPER IMPLEMENTATION")
    print("="*80)
    print("\nThis workflow helps you build PubMed queries for protein-protein interaction studies.")
    print("You'll be asked to approve suggestions at each stage.\n")

    # Compile the three stage graphs
    expansion_graph = create_protocol_expansion_graph().compile()
    query_graph = create_query_creation_graph().compile()
    search_graph = create_search_graph().compile()

    # Initialize session state
    session = SessionState(
        iteration=0,
        protocol_history=[],
        query_history=[],
        search_history=[],
        timestamp=datetime.now().isoformat()
    )

    # Main workflow loop (supports iterative refinement)
    while True:
        session["iteration"] += 1
        print(f"\n{'='*80}")
        print(f"ITERATION {session['iteration']}")
        print(f"{'='*80}\n")

        # =================================================================
        # HITL Checkpoint #1: Protocol Input & Ambiguity Detection
        # =================================================================
        while True:
            print("Enter your experimental protocols (comma-separated):")
            print("Example: co-immunoprecipitation, affinity purification, pull-down assay")
            user_input = input("\nYou: ").strip()

            if not user_input:
                print("No input provided. Please enter at least one protocol.")
                continue

            # Parse protocols
            user_protocols = [p.strip() for p in user_input.split(",") if p.strip()]

            # Check for ambiguity
            if detect_ambiguity(user_protocols):
                print("\n⚠️  Your input seems generic. Please be more specific.")
                print("Tip: Provide specific method names like 'co-immunoprecipitation', 'yeast two-hybrid', etc.")
                retry = input("\nTry again? (yes/no): ").strip().lower()
                if retry != "yes":
                    break
            else:
                break

        # Save to session
        session["protocol_history"].append(user_protocols)

        # =================================================================
        # Stage 1: Protocol Expansion
        # =================================================================
        print("\n🔄 Expanding protocols with LLM suggestions...")
        expansion_input = ProtocolExpansionInput(
            user_protocols=user_protocols,
            suggested_protocols=None,
            expansion_reasoning=None
        )
        expansion_result = expansion_graph.invoke(expansion_input)

        suggested_protocols = expansion_result["suggested_protocols"]
        reasoning = expansion_result["expansion_reasoning"]

        # =================================================================
        # HITL Checkpoint #2: Protocol Approval
        # =================================================================
        display_expansion_results(user_protocols, suggested_protocols, reasoning)

        print("\nOptions:")
        print("  1. approve all - Include all suggested protocols")
        print("  2. skip - Use only your original protocols")
        print("  3. custom - Manually select which to include")

        approval = input("\nYour choice: ").strip().lower()

        if approval == "approve all" or approval == "1":
            approved_protocols = user_protocols + suggested_protocols
            print(f"\n✓ Approved all. Total protocols: {len(approved_protocols)}")
        elif approval == "skip" or approval == "2":
            approved_protocols = user_protocols
            print(f"\n✓ Using original protocols only. Total: {len(approved_protocols)}")
        elif approval == "custom" or approval == "3":
            print("\nEnter the numbers of suggested protocols to include (comma-separated):")
            for i, p in enumerate(suggested_protocols, 1):
                print(f"  {i}. {p}")
            selections = input("\nNumbers: ").strip()
            selected_indices = [int(s.strip())-1 for s in selections.split(",") if s.strip().isdigit()]
            selected = [suggested_protocols[i] for i in selected_indices if 0 <= i < len(suggested_protocols)]
            approved_protocols = user_protocols + selected
            print(f"\n✓ Approved {len(selected)} suggestions. Total protocols: {len(approved_protocols)}")
        else:
            # Default to user protocols only
            approved_protocols = user_protocols
            print(f"\n✓ Using original protocols only. Total: {len(approved_protocols)}")

        # =================================================================
        # Stage 2: Query Creation
        # =================================================================
        print("\n🔄 Creating PubMed query...")
        query_input = QueryCreationInput(
            approved_protocols=approved_protocols,
            a2_block=None,
            final_query=None,
            query_components=None
        )
        query_result = query_graph.invoke(query_input)

        final_query = query_result["final_query"]
        query_components = query_result["query_components"]

        # =================================================================
        # HITL Checkpoint #3: Query Review
        # =================================================================
        display_query(final_query, query_components)

        query_approval = input("\nApprove this query? (yes/no): ").strip().lower()

        if query_approval != "yes":
            print("\n❌ Query not approved. Let's refine the protocols.")
            refine = input("Refine protocols? (yes/no): ").strip().lower()
            if refine == "yes":
                continue  # Go back to protocol input
            else:
                print("\nWorkflow cancelled.")
                break

        # Save to session
        session["query_history"].append(final_query)

        # =================================================================
        # Stage 3: PubMed Search
        # =================================================================
        print("\n🔄 Searching PubMed...")
        search_input = SearchInput(
            final_query=final_query,
            search_results=None,
            result_count=None,
            articles=None
        )
        search_result = search_graph.invoke(search_input)

        result_count = search_result["result_count"]
        articles = search_result["articles"]

        # =================================================================
        # HITL Checkpoint #4: Results Validation
        # =================================================================
        display_search_results(result_count, articles)

        validation = input("\nAccept these results? (yes/no): ").strip().lower()

        if validation == "yes":
            # Save to session
            session["search_history"].append({
                "query": final_query,
                "result_count": result_count,
                "timestamp": datetime.now().isoformat()
            })

            # =================================================================
            # HITL Checkpoint #5: Next Action
            # =================================================================
            print("\nWhat would you like to do next?")
            print("  1. done - Save results and exit")
            print("  2. refine - Refine the query with different protocols")
            print("  3. new - Start a completely new search")

            next_action = input("\nYour choice: ").strip().lower()

            if next_action == "done" or next_action == "1":
                # Save session data
                session_file = f"hitl_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
                with open(session_file, 'w') as f:
                    json.dump(session, f, indent=2)

                print(f"\n✓ Session saved to {session_file}")
                print("\n🎉 Workflow complete! Thank you for using the HITL Query Builder.")
                break

            elif next_action == "refine" or next_action == "2":
                print("\n🔄 Starting refinement...")
                continue

            elif next_action == "new" or next_action == "3":
                # Reset session
                session = SessionState(
                    iteration=0,
                    protocol_history=[],
                    query_history=[],
                    search_history=[],
                    timestamp=datetime.now().isoformat()
                )
                print("\n🔄 Starting new search...")
                continue

            else:
                print("\nInvalid choice. Exiting.")
                break

        else:
            print("\n❌ Results not accepted.")
            refine = input("Refine the query? (yes/no): ").strip().lower()
            if refine == "yes":
                continue
            else:
                print("\nWorkflow cancelled.")
                break


# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    try:
        run_hitl_workflow()
    except KeyboardInterrupt:
        print("\n\n⚠️  Workflow interrupted by user. Exiting...")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
