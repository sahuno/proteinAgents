#!/usr/bin/env python
"""
Automated HITL Workflow Test
Simulates a complete HITL session and generates a report.
"""

import sys
import os
from datetime import datetime

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from blockA_single_Agent_HITL import (
    create_hitl_graph,
    EnhancedState,
    HumanMessage,
    AIMessage,
    MemorySaver
)


def simulate_hitl_workflow():
    """
    Simulate a complete HITL workflow with predefined inputs.

    Inputs:
    1. Initial protocols: "co immunoprecipitation, affinity purification"
    2. Protocol approval: "approve all"
    3. Query approval: "approve"
    4. Results validation: "accept"
    5. Next action: "done"
    """
    print("=" * 80)
    print("HITL WORKFLOW SIMULATION")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Create graph with memory
    builder = create_hitl_graph()
    memory = MemorySaver()
    graph = builder.compile(checkpointer=memory)
    config = {"configurable": {"thread_id": "test_session"}}

    # Track workflow steps
    workflow_log = []

    # Step 1: Initial protocol input
    print("\n" + "="*80)
    print("STEP 1: Initial Protocol Input")
    print("="*80)
    user_input_1 = "co immunoprecipitation, affinity purification"
    print(f"User: {user_input_1}")

    state = {
        "messages": [HumanMessage(content=user_input_1)],
        "user_protocols": ["co immunoprecipitation", "affinity purification"],
        "current_step": "start"
    }

    result = graph.invoke(state, config)
    workflow_log.append({
        "step": 1,
        "checkpoint": "Initial Input",
        "user_input": user_input_1,
        "state_after": result.get("current_step"),
        "protocols": result.get("user_protocols", []),
        "suggested": result.get("suggested_protocols", [])
    })

    # Display AI response
    for msg in reversed(result.get("messages", [])):
        if isinstance(msg, AIMessage):
            print(f"\nAssistant: {msg.content[:500]}...")
            break

    print(f"\nCurrent Step: {result.get('current_step')}")
    print(f"User Protocols: {result.get('user_protocols', [])}")
    print(f"Suggested Protocols: {result.get('suggested_protocols', [])}")

    # Step 2: Protocol approval
    print("\n" + "="*80)
    print("STEP 2: Protocol Approval (Checkpoint #2)")
    print("="*80)
    user_input_2 = "approve all"
    print(f"User: {user_input_2}")

    # Manually handle approval (simulate conversation loop logic)
    user_protocols = result.get("user_protocols", [])
    suggested_protocols = result.get("suggested_protocols", [])

    state = result
    state["approved_protocols"] = user_protocols + suggested_protocols
    state["needs_protocol_approval"] = False
    state["current_step"] = "create_a2_block"
    state["messages"].append(HumanMessage(content=user_input_2))

    result = graph.invoke(state, config)
    workflow_log.append({
        "step": 2,
        "checkpoint": "Protocol Approval",
        "user_input": user_input_2,
        "state_after": result.get("current_step"),
        "approved_protocols": state["approved_protocols"],
        "A2_block_created": bool(result.get("A2_block"))
    })

    # Display AI response
    for msg in reversed(result.get("messages", [])):
        if isinstance(msg, AIMessage):
            print(f"\nAssistant: {msg.content[:500]}...")
            break

    print(f"\nCurrent Step: {result.get('current_step')}")
    print(f"Approved Protocols: {state['approved_protocols']}")
    print(f"A2 Block Created: {bool(result.get('A2_block'))}")

    # Step 3: Query approval
    print("\n" + "="*80)
    print("STEP 3: Query Approval (Checkpoint #3)")
    print("="*80)
    user_input_3 = "approve"
    print(f"User: {user_input_3}")

    state = result
    state["query_approved"] = True
    state["needs_query_approval"] = False
    state["current_step"] = "search_pubmed_node"
    state["messages"].append(HumanMessage(content=user_input_3))

    result = graph.invoke(state, config)
    workflow_log.append({
        "step": 3,
        "checkpoint": "Query Approval",
        "user_input": user_input_3,
        "state_after": result.get("current_step"),
        "final_query": result.get("final_query", "")[:200] + "...",
        "search_executed": bool(result.get("search_results"))
    })

    # Display AI response
    for msg in reversed(result.get("messages", [])):
        if isinstance(msg, AIMessage):
            print(f"\nAssistant: {msg.content[:500]}...")
            break

    print(f"\nCurrent Step: {result.get('current_step')}")
    print(f"Search Results: {result.get('search_results', {}).get('count', 0)} articles found")

    # Step 4: Results validation
    print("\n" + "="*80)
    print("STEP 4: Results Validation (Checkpoint #4)")
    print("="*80)
    user_input_4 = "accept"
    print(f"User: {user_input_4}")

    state = result
    state["results_validated"] = True
    state["needs_results_validation"] = False
    state["current_step"] = "ask_next"
    state["messages"].append(HumanMessage(content=user_input_4))

    result = graph.invoke(state, config)
    workflow_log.append({
        "step": 4,
        "checkpoint": "Results Validation",
        "user_input": user_input_4,
        "state_after": result.get("current_step"),
        "results_count": result.get("search_results", {}).get("count", 0)
    })

    # Display AI response
    for msg in reversed(result.get("messages", [])):
        if isinstance(msg, AIMessage):
            print(f"\nAssistant: {msg.content[:500]}...")
            break

    print(f"\nCurrent Step: {result.get('current_step')}")

    # Step 5: Next action
    print("\n" + "="*80)
    print("STEP 5: Next Action Decision (Checkpoint #5)")
    print("="*80)
    user_input_5 = "done"
    print(f"User: {user_input_5}")

    state = result
    state["next_action"] = "done"
    state["needs_next_action"] = False
    state["messages"].append(HumanMessage(content=user_input_5))

    result = graph.invoke(state, config)
    workflow_log.append({
        "step": 5,
        "checkpoint": "Next Action",
        "user_input": user_input_5,
        "state_after": result.get("current_step"),
        "workflow_complete": True
    })

    print(f"\nWorkflow Complete!")
    print(f"Final State: {result.get('current_step')}")

    # Generate summary
    print("\n" + "="*80)
    print("WORKFLOW SUMMARY")
    print("="*80)

    summary = {
        "initial_protocols": workflow_log[0]["protocols"],
        "suggested_protocols": workflow_log[0]["suggested"],
        "approved_protocols": workflow_log[1]["approved_protocols"],
        "final_query": result.get("final_query", ""),
        "search_results_count": result.get("search_results", {}).get("count", 0),
        "sample_pmids": result.get("search_results", {}).get("ids", [])[:5],
        "total_steps": len(workflow_log),
        "checkpoints_triggered": sum(1 for log in workflow_log if "checkpoint" in log)
    }

    return workflow_log, summary, result


if __name__ == "__main__":
    try:
        workflow_log, summary, final_state = simulate_hitl_workflow()

        print("\n" + "="*80)
        print("FINAL RESULTS")
        print("="*80)
        print(f"Initial Protocols: {summary['initial_protocols']}")
        print(f"Suggested Protocols: {summary['suggested_protocols']}")
        print(f"Approved Protocols (Total): {len(summary['approved_protocols'])}")
        print(f"PubMed Results: {summary['search_results_count']:,} articles")
        print(f"Sample PMIDs: {', '.join(summary['sample_pmids'])}")
        print(f"Total Workflow Steps: {summary['total_steps']}")
        print(f"HITL Checkpoints Triggered: {summary['checkpoints_triggered']}")

        # Save workflow data for report generation
        import json
        with open("hitl_workflow_data.json", "w") as f:
            json.dump({
                "workflow_log": workflow_log,
                "summary": summary,
                "timestamp": datetime.now().isoformat()
            }, f, indent=2)

        print("\nWorkflow data saved to: hitl_workflow_data.json")

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
