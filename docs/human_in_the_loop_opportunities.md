# Human-in-the-Loop (HITL) Opportunities for PubMed Query Agent

**Author:** Samuel Ahuno (ekwame001@gmail.com)
**Date:** 2025-12-08
**Based on:** `/notebooks/blockA_single_Agent.py`

---

## Overview

This document outlines opportunities to add Human-in-the-Loop (HITL) interactions to the PubMed query building agent. HITL allows users to review, approve, and refine the agent's work at critical decision points, ensuring accuracy and user satisfaction.

---

## Current Agent Behavior

The existing agent (`blockA_single_Agent.py`) operates in a **fully automated** manner:

1. User provides initial protocols (e.g., "co-IP, mass spectrometry, cryo-EM")
2. Agent suggests similar protocols
3. Agent creates A2_block automatically
4. Agent concatenates with O_block to create final query
5. Agent searches PubMed
6. Returns results

**Issue:** User has no control over agent decisions. If the agent misunderstands input or makes wrong assumptions, the user only discovers this at the end.

---

## 🎯 Five Key HITL Opportunities

### 1. Protocol Expansion Approval
**When:** After agent suggests similar protocols
**Why:** Ensures relevance, catches misunderstandings

#### Current Behavior
```python
# Agent automatically suggests protocols and creates A2_block
# No user review
```

#### Enhanced Behavior with HITL
```python
# Add to GraphState
class EnhancedState(MessagesState):
    suggested_protocols: list = []  # Agent's suggestions
    approved_protocols: list = []   # User's final selection
    needs_protocol_approval: bool = False
    user_feedback: str = ""

# New node: Request protocol approval
def request_protocol_approval(state: EnhancedState):
    """Agent suggests protocols and asks for approval."""
    suggestion_msg = (
        f"Based on your input, I suggest these protocols:\n"
        f"{', '.join(state['suggested_protocols'])}\n\n"
        f"Would you like to:\n"
        f"1. ✅ Approve these protocols\n"
        f"2. ➕ Add more protocols\n"
        f"3. ➖ Remove some protocols\n"
        f"Please respond with your choice and any modifications."
    )
    return {
        "messages": [AIMessage(content=suggestion_msg)],
        "needs_protocol_approval": True
    }

# Conditional edge: Wait for human approval
def check_approval_needed(state: EnhancedState):
    if state.get("needs_protocol_approval"):
        return "wait_for_human"  # Interrupt here
    else:
        return "create_a2_block"

# Add to graph
builder.add_node("request_approval", request_protocol_approval)
builder.add_conditional_edges(
    "assistant",
    check_approval_needed,
    {
        "wait_for_human": END,  # Pause and wait
        "create_a2_block": "create_a2_block"
    }
)
```

#### Usage
```python
# First invocation - agent suggests
messages = react_graph.invoke(
    {"messages": [HumanMessage(content="co-IP, mass spec")]},
    config
)

# Graph pauses, waits for human feedback
# User reviews suggestions in messages[-1]

# Second invocation - continue with feedback
messages = react_graph.invoke(
    {
        "messages": messages["messages"] + [
            HumanMessage(content="Remove FRET, add cryo-EM")
        ]
    },
    config
)
```

**Benefits:**
- User has final say on protocols
- Catches agent misunderstandings early
- Educational - user sees what agent is thinking

---

### 2. Query Review Before Search
**When:** After final query generation
**Why:** Prevents wasted API calls, educational

#### Current Behavior
```python
# Agent creates query and immediately searches PubMed
# User sees query only after search completes
```

#### Enhanced Behavior with HITL
```python
# Add node: Show query for approval
def show_query_for_approval(state: EnhancedState):
    final_query = state.get("final_query")
    estimated_count = state.get("estimated_count", "unknown")

    return {
        "messages": [AIMessage(content=f"""
THIS IS THE SUGGESTED FINAL QUERY:
{final_query}

Expected results: ~{estimated_count} papers

Would you like to:
1. ✅ Run this query
2. ✏️ Modify the query
3. 🔄 Start over

Please respond with your choice.
""")],
        "needs_query_approval": True
    }

# Routing function
def route_after_query_generation(state: EnhancedState):
    last_message = state["messages"][-1]

    # Check if user approved
    if "run" in last_message.content.lower() or "✅" in last_message.content:
        return "search_pubmed"
    elif "modify" in last_message.content.lower():
        return "assistant"  # Back to agent for refinement
    else:
        return "wait_for_approval"

builder.add_node("show_query", show_query_for_approval)
builder.add_conditional_edges(
    "show_query",
    route_after_query_generation,
    {
        "search_pubmed": "search_pubmed",
        "assistant": "assistant",
        "wait_for_approval": END
    }
)
```

**Benefits:**
- Avoids unnecessary PubMed API calls
- User learns PubMed query syntax
- Opportunity to refine before executing expensive search

---

### 3. Results Validation
**When:** After PubMed search returns results
**Why:** Ensures satisfaction, allows refinement

#### Current Behavior
```python
# Returns results and ends
# No feedback mechanism
```

#### Enhanced Behavior with HITL
```python
def validate_results(state: EnhancedState):
    result_count = state.get("search_results", {}).get("count", 0)
    sample_pmids = state.get("search_results", {}).get("ids", [])[:5]

    return {
        "messages": [AIMessage(content=f"""
📊 SEARCH RESULTS:
Found {result_count:,} papers

Sample PMIDs: {', '.join(sample_pmids)}
Preview: https://pubmed.ncbi.nlm.nih.gov/{sample_pmids[0] if sample_pmids else ''}

Are these results what you expected?
- Type 'good' to save and finish
- Type 'too many' to narrow down the search
- Type 'too few' to broaden the search
- Type 'refine' to adjust protocols
""")]
    }

def route_based_on_results_feedback(state: EnhancedState):
    last_msg = state["messages"][-1].content.lower()

    if "good" in last_msg:
        return END
    elif "too many" in last_msg or "narrow" in last_msg:
        return "add_exclusions"  # New node to add B_block
    elif "too few" in last_msg or "broaden" in last_msg:
        return "remove_filters"
    else:
        return "assistant"  # Refine with agent

builder.add_node("validate", validate_results)
builder.add_conditional_edges("validate", route_based_on_results_feedback)
```

**Benefits:**
- Immediate feedback on result quality
- Iterative refinement until satisfied
- Prevents downloading irrelevant papers

---

### 4. Ambiguity Resolution
**When:** User input contains unclear terms
**Why:** Prevents errors from agent assumptions

#### Current Behavior
```python
# Agent makes assumptions about ambiguous terms
# May misinterpret user intent
```

#### Enhanced Behavior with HITL
```python
def detect_ambiguity(state: EnhancedState):
    """Agent checks if user input is ambiguous."""
    user_input = state["messages"][-1].content

    # Check for ambiguous terms
    ambiguous_terms = []

    # Common ambiguities in biomedical literature
    if "IP" in user_input and "co-IP" not in user_input.lower():
        ambiguous_terms.append("IP (did you mean co-IP/immunoprecipitation or intellectual property?)")

    if "MS" in user_input and "mass" not in user_input.lower():
        ambiguous_terms.append("MS (mass spectrometry or multiple sclerosis?)")

    if "NMR" in user_input and "nuclear" not in user_input.lower():
        ambiguous_terms.append("NMR (nuclear magnetic resonance or neuronal?)")

    # Check for vague time references
    if "recent" in user_input.lower() or "latest" in user_input.lower():
        ambiguous_terms.append("Time range (what year range do you consider recent?)")

    if ambiguous_terms:
        return {
            "messages": [AIMessage(content=f"""
I found some potentially ambiguous terms in your input:
{chr(10).join(f'- {term}' for term in ambiguous_terms)}

Please clarify what you meant.
""")],
            "needs_clarification": True
        }

    return {"needs_clarification": False}

def check_needs_clarification(state: EnhancedState):
    if state.get("needs_clarification"):
        return "wait_for_clarification"
    else:
        return "assistant"

builder.add_node("detect_ambiguity", detect_ambiguity)
builder.add_conditional_edges("detect_ambiguity", check_needs_clarification)
```

**Benefits:**
- Prevents misinterpretation
- Forces explicit communication
- Improves query accuracy

---

### 5. Iterative Refinement Loop
**When:** After each major step
**Why:** Allows continuous improvement

#### Current Behavior
```python
# One-shot interaction
# User must restart to make changes
```

#### Enhanced Behavior with HITL
```python
def ask_for_next_action(state: EnhancedState):
    return {
        "messages": [AIMessage(content="""
Query created successfully! What would you like to do next?

1. 🔍 Run the search
2. ➕ Add more protocols to include
3. ➖ Add protocols to exclude (create B_block)
4. 🔄 Start over with different criteria
5. 💾 Save query without searching
6. 📊 Show me the query again

Type the number or describe your action.
""")]
    }

def route_next_action(state: EnhancedState):
    last_msg = state["messages"][-1].content.lower()

    if "1" in last_msg or "run" in last_msg or "search" in last_msg:
        return "search_pubmed"
    elif "2" in last_msg or "add more" in last_msg:
        return "assistant"  # Add more to A2
    elif "3" in last_msg or "exclude" in last_msg:
        return "create_exclusions"  # Create B_block
    elif "4" in last_msg or "start over" in last_msg:
        return START
    elif "5" in last_msg or "save" in last_msg:
        return "save_query"
    elif "6" in last_msg or "show" in last_msg:
        return "show_query"
    else:
        return "ask_for_next_action"  # Ask again

builder.add_node("ask_next_action", ask_for_next_action)
builder.add_conditional_edges("ask_next_action", route_next_action)
```

**Benefits:**
- Multi-turn refinement
- User maintains control throughout
- Natural conversation flow

---

## 🏗️ Complete Enhanced Graph Architecture

```python
from langgraph.graph import MessagesState, StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver
from langchain_core.messages import HumanMessage, AIMessage

# Enhanced state with HITL control flags
class EnhancedState(MessagesState):
    # User inputs
    user_protocols: list = []

    # Agent suggestions
    suggested_protocols: list = []
    suggested_exclusions: list = []

    # Query components
    A2_block: str = ""
    B_block: str = ""
    final_query: str = ""

    # Search results
    search_results: dict = {}
    estimated_count: int = 0

    # Control flags for HITL
    needs_protocol_approval: bool = False
    needs_query_approval: bool = False
    needs_clarification: bool = False
    iteration_count: int = 0

# Build enhanced graph
builder = StateGraph(EnhancedState)

# Nodes
builder.add_node("assistant", assistant)
builder.add_node("tools", ToolNode(tools))
builder.add_node("detect_ambiguity", detect_ambiguity)
builder.add_node("request_approval", request_protocol_approval)
builder.add_node("show_query", show_query_for_approval)
builder.add_node("validate_results", validate_results)
builder.add_node("ask_next_action", ask_for_next_action)

# Flow with HITL checkpoints
builder.add_edge(START, "detect_ambiguity")
builder.add_conditional_edges("detect_ambiguity", check_needs_clarification, {
    "wait_for_clarification": END,
    "assistant": "assistant"
})

builder.add_conditional_edges("assistant", tools_condition)
builder.add_edge("tools", "request_approval")

builder.add_conditional_edges("request_approval", check_approval_needed, {
    "wait_for_human": END,
    "create_a2_block": "show_query"
})

builder.add_conditional_edges("show_query", route_after_query_generation, {
    "search_pubmed": "tools",
    "assistant": "assistant",
    "wait_for_approval": END
})

builder.add_edge("search_pubmed", "validate_results")
builder.add_conditional_edges("validate_results", route_based_on_results_feedback, {
    END: END,
    "add_exclusions": "assistant",
    "remove_filters": "assistant",
    "assistant": "assistant"
})

# Compile with memory + interrupts
memory = MemorySaver()
react_graph_hitl = builder.compile(
    checkpointer=memory,
    interrupt_before=["request_approval", "show_query", "validate_results"]  # Pause points
)
```

---

## 🎮 Usage Pattern with HITL

### Example: Full HITL Workflow

```python
config = {"configurable": {"thread_id": "user123"}}

# ============================================================
# STEP 1: Initial query with ambiguity detection
# ============================================================
state = react_graph_hitl.invoke(
    {"messages": [HumanMessage(content="IP, MS, cryo-EM")]},
    config
)
# Output: "I found some ambiguous terms..."
print(state["messages"][-1].content)

# ============================================================
# STEP 2: User clarifies
# ============================================================
state = react_graph_hitl.invoke(
    {"messages": state["messages"] + [HumanMessage(content="co-IP, mass spectrometry, cryo-EM")]},
    config
)
# Graph pauses at "request_approval"
# Output: "Based on your input, I suggest these protocols: co-IP, immunoprecipitation, ..."
print(state["messages"][-1].content)

# ============================================================
# STEP 3: User approves protocol suggestions
# ============================================================
state = react_graph_hitl.invoke(
    {"messages": state["messages"] + [HumanMessage(content="Looks good, proceed")]},
    config
)
# Graph pauses at "show_query"
# Output: "THIS IS THE SUGGESTED FINAL QUERY: ..."
print(state["messages"][-1].content)

# ============================================================
# STEP 4: User approves query and runs search
# ============================================================
state = react_graph_hitl.invoke(
    {"messages": state["messages"] + [HumanMessage(content="Run it")]},
    config
)
# Graph pauses at "validate_results"
# Output: "Found 1,234 papers. Sample PMIDs: ..."
print(state["messages"][-1].content)

# ============================================================
# STEP 5: User validates results
# ============================================================
state = react_graph_hitl.invoke(
    {"messages": state["messages"] + [HumanMessage(content="too many, narrow down")]},
    config
)
# Graph goes back to assistant to add exclusions (B_block)
# Agent asks: "What protocols would you like to exclude?"

# ... Continue iteration until satisfied
```

---

## 📊 HITL Opportunities Summary

| **Stage** | **HITL Point** | **Benefit** | **Implementation Complexity** |
|-----------|----------------|-------------|------------------------------|
| Input validation | Ambiguity detection | Prevents errors from assumptions | Low |
| Protocol suggestion | Approve/modify agent suggestions | Ensures relevance, catches misunderstandings | Medium |
| Query generation | Review final query before search | Prevents wasted API calls, educational | Low |
| Search execution | Approve estimated result count | Avoids retrieving too many/few papers | Low |
| Results validation | Confirm quality, iterate if needed | Ensures satisfaction, allows refinement | Medium |
| Iterative refinement | Multi-turn dialogue | Full user control, natural workflow | High |

---

## 🔑 Key Implementation Patterns

### Pattern 1: Interrupt-Based HITL (Recommended)

Use LangGraph's built-in interrupt mechanism:

```python
graph = builder.compile(
    checkpointer=memory,
    interrupt_before=["request_approval", "show_query", "validate_results"]
)

# Graph automatically pauses at these nodes
state = graph.invoke(input, config)
# ... user reviews ...
state = graph.invoke(updated_input, config)  # Continues from where it paused
```

**Pros:**
- Clean, stateful
- Built-in support
- Easy to resume

**Cons:**
- Requires checkpointer (memory)
- Two invocations per checkpoint

---

### Pattern 2: Conditional Routing HITL

Use conditional edges to wait for user response:

```python
def route_based_on_approval(state):
    if state.get("needs_approval"):
        return END  # Pause
    else:
        return "next_node"

builder.add_conditional_edges("node", route_based_on_approval, {
    END: END,
    "next_node": "next_node"
})
```

**Pros:**
- More flexible
- Can skip HITL based on conditions

**Cons:**
- More complex routing logic

---

### Pattern 3: Streaming HITL (Advanced)

Stream intermediate results to user:

```python
for chunk in graph.stream(input, config):
    print(chunk)  # Show progress
    if requires_user_input(chunk):
        user_input = input("Your response: ")
        # Inject user input into stream
```

**Pros:**
- Real-time feedback
- Best UX

**Cons:**
- Complex implementation
- Requires streaming support

---

## 🚀 Next Steps

### Implementation Roadmap

**Phase 1: Basic HITL** (Week 1)
- [ ] Add protocol approval checkpoint
- [ ] Add query review checkpoint
- [ ] Test with memory checkpointer

**Phase 2: Validation & Refinement** (Week 2)
- [ ] Add results validation
- [ ] Add ambiguity detection
- [ ] Implement refinement loops

**Phase 3: Advanced Features** (Week 3)
- [ ] Add B_block (exclusion) support
- [ ] Add multi-turn refinement
- [ ] Add save/resume functionality

**Phase 4: Polish** (Week 4)
- [ ] Add usage analytics
- [ ] Improve error messages
- [ ] Create user documentation

---

## 📝 Testing Checklist

- [ ] Test with clear, unambiguous input
- [ ] Test with ambiguous input (IP, MS, etc.)
- [ ] Test approval → modify → approval flow
- [ ] Test "too many results" refinement
- [ ] Test "too few results" refinement
- [ ] Test save and resume functionality
- [ ] Test multi-session memory
- [ ] Test error handling (API failures, etc.)

---

## 🔗 Related Documentation

- **Original Agent:** `/notebooks/blockA_single_Agent.py`
- **LangGraph Interrupts:** [LangGraph Documentation](https://langchain-ai.github.io/langgraph/how-tos/human_in_the_loop/breakpoints/)
- **PubMed API:** [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- **MeSH Terms:** [MeSH Browser](https://meshb.nlm.nih.gov/)

---

## 📧 Contact

For questions or feedback on this design:
- **Author:** Samuel Ahuno
- **Email:** ekwame001@gmail.com
- **Project:** proteinAgents
