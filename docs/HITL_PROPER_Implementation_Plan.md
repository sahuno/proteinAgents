# Implementation Plan: blockA_single_Agent_HITL_PROPER.py

**Date:** 2025-12-26
**Author:** Samuel Ahuno (ekwame001@gmail.com)
**Strategy:** Option 2 - Separate Graph Invocations per Stage
**Status:** Planning Phase

---

## Executive Summary

This plan outlines the implementation of a production-ready HITL workflow using **separate, focused graphs** for each workflow stage, connected by simple Python orchestration code. This approach:

- ✅ Works perfectly with CLI `input()` calls
- ✅ No complex state management bugs
- ✅ Easy to test each stage independently
- ✅ Maintainable and extensible
- ✅ **Actually works** (proven pattern)

---

## Architecture Overview

### High-Level Structure

```
┌─────────────────────────────────────────────────────────────┐
│                    CLI Orchestrator                          │
│  (Simple Python function with input() calls)                │
└───────┬─────────────────────────────────────────────────────┘
        │
        ├─> Stage 1: Protocol Expansion Graph
        │   Input: user_protocols
        │   Output: suggested_protocols
        │   [HITL Checkpoint #2: Approve protocols]
        │
        ├─> Stage 2: Query Creation Graph
        │   Input: approved_protocols
        │   Output: A2_block, final_query
        │   [HITL Checkpoint #3: Approve query]
        │
        ├─> Stage 3: PubMed Search Graph
        │   Input: final_query
        │   Output: search_results
        │   [HITL Checkpoint #4: Validate results]
        │
        └─> [HITL Checkpoint #5: Refine/New/Done]
```

### Key Principles

1. **One graph per workflow stage** (not one monolithic graph)
2. **HITL checkpoints are simple Python `input()` calls** (not graph routing)
3. **State passed explicitly between stages** (no shared global state)
4. **Each stage is independently testable**

---

## Detailed Component Design

### Component 1: State Classes

**Purpose:** Type-safe data structures for each stage

```python
# File: blockA_single_Agent_HITL_PROPER.py (lines 50-120)

from typing import TypedDict, List, Dict

class ProtocolExpansionInput(TypedDict):
    """Input for protocol expansion stage."""
    user_protocols: List[str]

class ProtocolExpansionOutput(TypedDict):
    """Output from protocol expansion stage."""
    user_protocols: List[str]
    suggested_protocols: List[str]
    reasoning: str  # Why these protocols were suggested

class QueryCreationInput(TypedDict):
    """Input for query creation stage."""
    approved_protocols: List[str]

class QueryCreationOutput(TypedDict):
    """Output from query creation stage."""
    A2_block: str
    final_query: str
    query_breakdown: Dict[str, str]  # For display

class SearchInput(TypedDict):
    """Input for PubMed search stage."""
    final_query: str
    retmax: int

class SearchOutput(TypedDict):
    """Output from PubMed search stage."""
    count: int
    ids: List[str]
    webenv: str
    query_key: str

class WorkflowSession(TypedDict):
    """Complete workflow state across iterations."""
    iteration: int
    user_protocols: List[str]
    suggested_protocols: List[str]
    approved_protocols: List[str]
    A2_block: str
    final_query: str
    search_results: SearchOutput
```

**Benefits:**
- Type safety
- Clear contracts between stages
- Easy to serialize for logging
- Self-documenting

---

### Component 2: Stage 1 - Protocol Expansion Graph

**Purpose:** Suggest similar protocols using LLM

**File:** `blockA_single_Agent_HITL_PROPER.py` (lines 150-250)

```python
def create_protocol_expansion_graph() -> StateGraph:
    """
    Create graph for protocol expansion stage.

    Nodes:
    - expand_protocols: Call LLM to suggest similar methods
    - format_suggestions: Format for display

    Flow:
    START → expand_protocols → format_suggestions → END
    """
    builder = StateGraph(ProtocolExpansionInput)

    builder.add_node("expand_protocols", expand_protocols_node)
    builder.add_node("format_suggestions", format_suggestions_node)

    builder.add_edge(START, "expand_protocols")
    builder.add_edge("expand_protocols", "format_suggestions")
    builder.add_edge("format_suggestions", END)

    return builder

def expand_protocols_node(state: ProtocolExpansionInput) -> dict:
    """Use LLM to suggest similar protocols."""
    user_protocols = state["user_protocols"]

    prompt = f"""Given these protocols: {', '.join(user_protocols)}

Suggest 3-5 similar experimental methods for protein-protein interaction studies.
Include brief reasoning for each suggestion.

Format:
1. Protocol Name - Reason
2. Protocol Name - Reason
..."""

    response = llm.invoke([sys_msg, HumanMessage(content=prompt)])

    # Parse response
    suggested = parse_llm_suggestions(response.content)

    return {
        "suggested_protocols": suggested["protocols"],
        "reasoning": suggested["reasoning"]
    }

def format_suggestions_node(state: dict) -> dict:
    """Format suggestions for display."""
    # Add formatting logic
    return state
```

**Testing:**
```python
def test_protocol_expansion():
    graph = create_protocol_expansion_graph().compile()

    result = graph.invoke({
        "user_protocols": ["co-IP", "affinity purification"]
    })

    assert len(result["suggested_protocols"]) >= 3
    assert len(result["suggested_protocols"]) <= 5
    assert "reasoning" in result
```

---

### Component 3: Stage 2 - Query Creation Graph

**Purpose:** Generate PubMed A2_block and final query

**File:** `blockA_single_Agent_HITL_PROPER.py` (lines 280-380)

```python
def create_query_creation_graph() -> StateGraph:
    """
    Create graph for query creation stage.

    Nodes:
    - create_a2_block: LLM generates A2_block with MeSH terms
    - create_o_block: Generate O_block (publication filters)
    - concatenate_blocks: Combine into final query
    - validate_query: Check query syntax

    Flow:
    START → create_a2_block → create_o_block → concatenate_blocks
          → validate_query → END
    """
    builder = StateGraph(QueryCreationInput)

    builder.add_node("create_a2_block", create_a2_block_node)
    builder.add_node("create_o_block", create_o_block_node)
    builder.add_node("concatenate_blocks", concatenate_blocks_node)
    builder.add_node("validate_query", validate_query_node)

    builder.add_edge(START, "create_a2_block")
    builder.add_edge("create_a2_block", "create_o_block")
    builder.add_edge("create_o_block", "concatenate_blocks")
    builder.add_edge("concatenate_blocks", "validate_query")
    builder.add_edge("validate_query", END)

    return builder

def create_a2_block_node(state: QueryCreationInput) -> dict:
    """Generate A2_block with MeSH terms."""
    protocols = state["approved_protocols"]

    prompt = f"""Create a PubMed A2_block for: {', '.join(protocols)}

Requirements:
- Use official MeSH terms where possible
- Include common alternative names
- Include abbreviations
- Use OR operators
- Format as valid PubMed query syntax

Example format:
("Term1"[MeSH Terms] OR "alternative1"[All Fields] OR "abbrev1"[All Fields]
OR "Term2"[MeSH Terms] OR "alternative2"[All Fields])"""

    response = llm.invoke([sys_msg, HumanMessage(content=prompt)])

    return {
        "A2_block": response.content,
        "query_breakdown": {"A2": response.content}
    }

def create_o_block_node(state: dict) -> dict:
    """Generate O_block (publication filters)."""
    O_block = """("english"[Language]
NOT "meta-analysis"[Publication Type]
NOT "review"[Publication Type]
NOT "retracted publication"[Publication Type]
NOT "retraction of publication"[Publication Type]
NOT "published erratum"[Publication Type]
NOT "controlled clinical trial"[Publication Type]
NOT "clinical study"[Publication Type]
NOT "clinical trial"[Publication Type])"""

    state["query_breakdown"]["O"] = O_block
    return state

def concatenate_blocks_node(state: dict) -> dict:
    """Combine O_block and A2_block."""
    O_block = state["query_breakdown"]["O"]
    A2_block = state["A2_block"]

    final_query = f"({O_block} AND ({A2_block}))"

    return {
        "final_query": final_query
    }

def validate_query_node(state: dict) -> dict:
    """Validate query syntax (basic checks)."""
    query = state["final_query"]

    # Basic validation
    assert query.count("(") == query.count(")"), "Unbalanced parentheses"
    assert "AND" in query or "OR" in query, "Missing boolean operators"

    return state
```

**Testing:**
```python
def test_query_creation():
    graph = create_query_creation_graph().compile()

    result = graph.invoke({
        "approved_protocols": ["co-IP", "yeast two-hybrid"]
    })

    assert "A2_block" in result
    assert "final_query" in result
    assert "english" in result["final_query"].lower()  # O_block present
    assert result["final_query"].count("(") == result["final_query"].count(")")
```

---

### Component 4: Stage 3 - PubMed Search Graph

**Purpose:** Execute PubMed search and format results

**File:** `blockA_single_Agent_HITL_PROPER.py` (lines 410-480)

```python
def create_search_graph() -> StateGraph:
    """
    Create graph for PubMed search stage.

    Nodes:
    - execute_search: Call Entrez.esearch
    - format_results: Format for display
    - validate_results: Check result count is reasonable

    Flow:
    START → execute_search → format_results → validate_results → END
    """
    builder = StateGraph(SearchInput)

    builder.add_node("execute_search", execute_search_node)
    builder.add_node("format_results", format_results_node)
    builder.add_node("validate_results", validate_results_node)

    builder.add_edge(START, "execute_search")
    builder.add_edge("execute_search", "format_results")
    builder.add_edge("format_results", "validate_results")
    builder.add_edge("validate_results", END)

    return builder

def execute_search_node(state: SearchInput) -> dict:
    """Execute PubMed search."""
    query = state["final_query"]
    retmax = state.get("retmax", 200)

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
        "query_key": record["QueryKey"]
    }

def format_results_node(state: dict) -> dict:
    """Format results for display."""
    # Could add PMID detail fetching here
    return state

def validate_results_node(state: dict) -> dict:
    """Validate result count is reasonable."""
    count = state["count"]

    if count == 0:
        state["warning"] = "No results found - query may be too restrictive"
    elif count > 100000:
        state["warning"] = "Very many results - query may be too broad"

    return state
```

**Testing:**
```python
def test_search():
    graph = create_search_graph().compile()

    # Mock query that should return results
    test_query = '("Immunoprecipitation"[MeSH Terms])'

    result = graph.invoke({
        "final_query": test_query,
        "retmax": 10
    })

    assert "count" in result
    assert "ids" in result
    assert result["count"] > 0  # Should find some results
```

---

### Component 5: CLI Orchestrator (The Main Function)

**Purpose:** Connect all stages with HITL checkpoints

**File:** `blockA_single_Agent_HITL_PROPER.py` (lines 520-720)

```python
def run_hitl_workflow():
    """
    Main HITL workflow orchestrator.

    Connects all stages with HITL checkpoints using simple input() calls.
    """
    print("="*80)
    print("PubMed Query Builder with HITL (PROPER VERSION)")
    print("="*80)

    # Create all graphs (reusable across iterations)
    expansion_graph = create_protocol_expansion_graph().compile()
    query_graph = create_query_creation_graph().compile()
    search_graph = create_search_graph().compile()

    # Initialize session state
    session = WorkflowSession(
        iteration=0,
        user_protocols=[],
        suggested_protocols=[],
        approved_protocols=[],
        A2_block="",
        final_query="",
        search_results={}
    )

    # Main workflow loop (supports refinement iterations)
    while True:
        session["iteration"] += 1
        print(f"\n{'='*80}")
        print(f"Iteration {session['iteration']}")
        print(f"{'='*80}\n")

        # ===================================================================
        # STAGE 0: Get User Protocols
        # ===================================================================
        print("Step 1: Enter Protocols")
        print("-" * 40)
        user_input = input("Enter protocols (comma-separated) or 'quit': ").strip()

        if user_input.lower() in ["quit", "exit", "q"]:
            print("Goodbye!")
            break

        # Parse protocols
        protocols = [p.strip() for p in user_input.split(",") if p.strip()]
        session["user_protocols"] = protocols

        # ===================================================================
        # HITL CHECKPOINT #1: Ambiguity Detection
        # ===================================================================
        if detect_ambiguity(protocols):
            print("\n⚠️  Your protocols seem ambiguous.")
            clarify = input("Would you like to provide more specific terms? (yes/no): ")
            if clarify.lower() == "yes":
                continue  # Loop back to get new protocols

        # ===================================================================
        # STAGE 1: Protocol Expansion
        # ===================================================================
        print(f"\n{'='*80}")
        print("Step 2: Expanding Protocols (AI Suggestions)")
        print("-" * 40)

        expansion_result = expansion_graph.invoke({
            "user_protocols": protocols
        })

        session["suggested_protocols"] = expansion_result["suggested_protocols"]

        # Display suggestions
        print("\nYour protocols:")
        for i, p in enumerate(protocols, 1):
            print(f"  {i}. {p}")

        print("\nSuggested additions:")
        for i, p in enumerate(expansion_result["suggested_protocols"], 1):
            print(f"  {i}. {p}")

        print(f"\nReasoning:\n{expansion_result['reasoning']}")

        # ===================================================================
        # HITL CHECKPOINT #2: Protocol Approval
        # ===================================================================
        print(f"\n{'='*80}")
        print("CHECKPOINT #2: Protocol Approval")
        print("-" * 40)

        approval = input("\nOptions:\n"
                        "  1. 'approve all' - Include all suggestions\n"
                        "  2. 'skip' - Use only your protocols\n"
                        "  3. 'approve: 1,3,5' - Select specific numbers\n"
                        "\nYour choice: ").strip().lower()

        if "approve all" in approval:
            approved = protocols + expansion_result["suggested_protocols"]
        elif "skip" in approval:
            approved = protocols
        elif "approve:" in approval:
            # Parse selected indices
            try:
                indices = [int(x.strip()) - 1 for x in approval.split(":")[1].split(",")]
                selected = [expansion_result["suggested_protocols"][i] for i in indices
                           if 0 <= i < len(expansion_result["suggested_protocols"])]
                approved = protocols + selected
            except:
                print("Invalid selection, using all protocols")
                approved = protocols + expansion_result["suggested_protocols"]
        else:
            # Default: approve all
            approved = protocols + expansion_result["suggested_protocols"]

        session["approved_protocols"] = approved

        print(f"\nApproved protocols ({len(approved)} total):")
        for i, p in enumerate(approved, 1):
            print(f"  {i}. {p}")

        # ===================================================================
        # STAGE 2: Query Creation
        # ===================================================================
        print(f"\n{'='*80}")
        print("Step 3: Creating PubMed Query")
        print("-" * 40)

        query_result = query_graph.invoke({
            "approved_protocols": approved
        })

        session["A2_block"] = query_result["A2_block"]
        session["final_query"] = query_result["final_query"]

        # Display query
        print("\nA2 Block (Methods/Approaches):")
        print("-" * 40)
        print(query_result["A2_block"][:500] + "...")

        print("\nFinal PubMed Query:")
        print("-" * 40)
        print(query_result["final_query"][:500] + "...")

        # ===================================================================
        # HITL CHECKPOINT #3: Query Approval
        # ===================================================================
        print(f"\n{'='*80}")
        print("CHECKPOINT #3: Query Approval")
        print("-" * 40)

        query_approval = input("\nOptions:\n"
                              "  1. 'approve' - Proceed with search\n"
                              "  2. 'modify' - Go back to protocol selection\n"
                              "  3. 'show full' - See complete query\n"
                              "\nYour choice: ").strip().lower()

        if "show full" in query_approval:
            print("\nComplete Query:")
            print("="*80)
            print(query_result["final_query"])
            print("="*80)
            query_approval = input("\nApprove? (yes/no): ").strip().lower()

        if "modify" in query_approval or query_approval == "no":
            print("\nReturning to protocol selection...")
            continue  # Loop back

        # ===================================================================
        # STAGE 3: PubMed Search
        # ===================================================================
        print(f"\n{'='*80}")
        print("Step 4: Executing PubMed Search")
        print("-" * 40)

        search_result = search_graph.invoke({
            "final_query": query_result["final_query"],
            "retmax": 200
        })

        session["search_results"] = search_result

        # Display results
        print(f"\nSearch Results:")
        print(f"  Total matches: {search_result['count']:,} articles")
        print(f"  Retrieved: {len(search_result['ids'])} PMIDs")

        if "warning" in search_result:
            print(f"\n⚠️  {search_result['warning']}")

        if search_result["ids"]:
            print(f"\n  Sample PMIDs:")
            for i, pmid in enumerate(search_result["ids"][:5], 1):
                print(f"    {i}. {pmid}")

        # ===================================================================
        # HITL CHECKPOINT #4: Results Validation
        # ===================================================================
        print(f"\n{'='*80}")
        print("CHECKPOINT #4: Results Validation")
        print("-" * 40)

        validation = input("\nOptions:\n"
                          "  1. 'accept' - Results look good\n"
                          "  2. 'too broad' - Too many results\n"
                          "  3. 'too narrow' - Too few results\n"
                          "\nYour assessment: ").strip().lower()

        # ===================================================================
        # HITL CHECKPOINT #5: Next Action
        # ===================================================================
        print(f"\n{'='*80}")
        print("CHECKPOINT #5: Next Action")
        print("-" * 40)

        next_action = input("\nOptions:\n"
                           "  1. 'refine' - Modify protocols and try again\n"
                           "  2. 'new search' - Start fresh\n"
                           "  3. 'export' - Save results and finish\n"
                           "  4. 'done' - Finish without saving\n"
                           "\nYour choice: ").strip().lower()

        if "refine" in next_action:
            print("\nStarting refinement iteration...")
            continue  # Loop back with session state preserved

        elif "new" in next_action:
            print("\nStarting new search...")
            # Reset session except iteration count
            session = WorkflowSession(
                iteration=session["iteration"],
                user_protocols=[],
                suggested_protocols=[],
                approved_protocols=[],
                A2_block="",
                final_query="",
                search_results={}
            )
            continue

        elif "export" in next_action:
            # Export results
            export_results(session)
            print("\n✅ Results exported successfully!")
            break

        else:  # done
            print("\n✅ Workflow complete!")
            break

    # Final summary
    print(f"\n{'='*80}")
    print("Session Summary")
    print(f"{'='*80}")
    print(f"  Total iterations: {session['iteration']}")
    print(f"  Final query: {session['final_query'][:100]}...")
    print(f"  Results: {session['search_results'].get('count', 0):,} articles")
    print(f"{'='*80}\n")


def detect_ambiguity(protocols: List[str]) -> bool:
    """Check if protocols are too generic."""
    ambiguous_keywords = ["protein", "binding", "interaction", "assay"]
    generic_count = sum(1 for p in protocols
                       if any(kw in p.lower() for kw in ambiguous_keywords))

    return generic_count > 0 and len(protocols) < 3


def export_results(session: WorkflowSession):
    """Export session results to file."""
    import json
    from datetime import datetime

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"pubmed_query_results_{timestamp}.json"

    with open(filename, "w") as f:
        json.dump(session, f, indent=2)

    print(f"\nResults saved to: {filename}")
```

---

## Testing Strategy

### Unit Tests (Per Stage)

```python
# File: test_hitl_proper.py

def test_protocol_expansion_stage():
    """Test Stage 1 in isolation."""
    graph = create_protocol_expansion_graph().compile()

    result = graph.invoke({
        "user_protocols": ["co-IP", "affinity purification"]
    })

    assert len(result["suggested_protocols"]) >= 3
    assert "reasoning" in result

def test_query_creation_stage():
    """Test Stage 2 in isolation."""
    graph = create_query_creation_graph().compile()

    result = graph.invoke({
        "approved_protocols": ["co-IP", "yeast two-hybrid"]
    })

    assert "A2_block" in result
    assert "final_query" in result
    assert "MeSH" in result["A2_block"]

def test_search_stage():
    """Test Stage 3 in isolation."""
    graph = create_search_graph().compile()

    result = graph.invoke({
        "final_query": '("Immunoprecipitation"[MeSH Terms])',
        "retmax": 10
    })

    assert result["count"] > 0
    assert len(result["ids"]) > 0
```

### Integration Test

```python
def test_end_to_end_automated():
    """Test complete workflow with mocked user input."""
    # Mock user inputs
    inputs = [
        "co-IP, affinity purification",  # Initial protocols
        "approve all",                    # Checkpoint #2
        "approve",                        # Checkpoint #3
        "accept",                         # Checkpoint #4
        "done"                            # Checkpoint #5
    ]

    with mock_input(inputs):
        run_hitl_workflow()

    # Verify final state
    assert session["search_results"]["count"] > 0
```

---

## Implementation Timeline

### Phase 1: Core Infrastructure (30 min)
- [ ] Create file structure
- [ ] Define state classes
- [ ] Set up imports and utilities

### Phase 2: Stage Graphs (1 hour)
- [ ] Implement protocol expansion graph
- [ ] Implement query creation graph
- [ ] Implement search graph
- [ ] Write unit tests for each

### Phase 3: CLI Orchestrator (45 min)
- [ ] Implement main workflow loop
- [ ] Add all 5 HITL checkpoints
- [ ] Add input validation

### Phase 4: Testing & Refinement (45 min)
- [ ] Run end-to-end test
- [ ] Fix bugs
- [ ] Add error handling
- [ ] Create documentation

**Total Estimated Time: 3 hours**

---

## Success Criteria

✅ **Functional Requirements:**
1. User can input protocols without looping bugs
2. All 5 HITL checkpoints work correctly
3. PubMed search executes successfully
4. Results are displayed clearly
5. Refinement iterations work

✅ **Technical Requirements:**
1. Each stage graph independently testable
2. No shared mutable state
3. Clean separation of concerns
4. Type-safe state management
5. Error handling at each stage

✅ **User Experience:**
1. Clear prompts at each checkpoint
2. Helpful error messages
3. Progress indicators
4. Session state preserved across iterations

---

## Advantages over Original Implementation

| Aspect | Original (Buggy) | Proper Implementation |
|--------|------------------|----------------------|
| **Graph Structure** | One monolithic graph | 3 focused graphs |
| **HITL Checkpoints** | Conditional routing to END | Simple `input()` calls |
| **State Management** | Manual updates + re-invoke | Explicit passing between stages |
| **Resume Logic** | Broken (re-starts from START) | No resume needed |
| **Testing** | Hard (monolithic) | Easy (per-stage) |
| **Debugging** | Complex (state tracking) | Simple (linear flow) |
| **Works in CLI** | ❌ No (loops infinitely) | ✅ Yes |
| **Code Clarity** | Low (complex routing) | High (obvious flow) |
| **Maintainability** | Low (tight coupling) | High (loose coupling) |

---

## Next Steps

1. **Approve this plan** - Review and confirm approach
2. **Implementation** - Build according to plan
3. **Testing** - Run all test scenarios
4. **Documentation** - Create user guide
5. **Deployment** - Ready for production use

**Ready to proceed with implementation?**
