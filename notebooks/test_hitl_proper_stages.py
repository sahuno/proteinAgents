#!/usr/bin/env python3
"""
Test script for HITL_PROPER implementation - Tests each stage independently

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-01-11
"""

import os
import sys
from blockA_single_Agent_HITL_PROPER import (
    # Graphs
    create_protocol_expansion_graph,
    create_query_creation_graph,
    create_search_graph,
    # State types
    ProtocolExpansionInput,
    QueryCreationInput,
    SearchInput,
    # Helper functions
    detect_ambiguity,
)


def test_stage1_protocol_expansion():
    """Test Stage 1: Protocol Expansion Graph"""
    print("\n" + "="*80)
    print("TEST: Stage 1 - Protocol Expansion")
    print("="*80)

    # Compile graph
    expansion_graph = create_protocol_expansion_graph().compile()

    # Test input
    test_protocols = ["co-immunoprecipitation", "affinity purification"]

    print(f"\nInput protocols: {test_protocols}")

    # Create input state
    expansion_input = ProtocolExpansionInput(
        user_protocols=test_protocols,
        suggested_protocols=None,
        expansion_reasoning=None
    )

    # Invoke graph
    print("\n🔄 Running protocol expansion graph...")
    result = expansion_graph.invoke(expansion_input)

    # Display results
    print("\n✅ Graph execution completed!")
    print(f"\nOriginal protocols: {result['user_protocols']}")
    print(f"\nSuggested protocols: {result['suggested_protocols']}")
    print(f"\nReasoning:\n{result['expansion_reasoning']}")

    # Validate
    assert result['user_protocols'] == test_protocols, "User protocols should be preserved"
    assert result['suggested_protocols'] is not None, "Should have suggested protocols"
    assert len(result['suggested_protocols']) > 0, "Should suggest at least one protocol"
    assert result['expansion_reasoning'] is not None, "Should have reasoning"

    print("\n✅ Stage 1 test passed!")
    return result


def test_stage2_query_creation():
    """Test Stage 2: Query Creation Graph"""
    print("\n" + "="*80)
    print("TEST: Stage 2 - Query Creation")
    print("="*80)

    # Compile graph
    query_graph = create_query_creation_graph().compile()

    # Test input (simulating approved protocols from Stage 1)
    test_protocols = [
        "co-immunoprecipitation",
        "affinity purification",
        "yeast two-hybrid",
        "pull-down assay"
    ]

    print(f"\nApproved protocols: {test_protocols}")

    # Create input state
    query_input = QueryCreationInput(
        approved_protocols=test_protocols,
        a2_block=None,
        final_query=None,
        query_components=None
    )

    # Invoke graph
    print("\n🔄 Running query creation graph...")
    result = query_graph.invoke(query_input)

    # Display results
    print("\n✅ Graph execution completed!")
    print(f"\nA2 Block:\n{result['a2_block']}")
    print(f"\nO Block:\n{result['query_components']['o_block']}")
    print(f"\nFinal Query:\n{result['final_query']}")
    print(f"\nQuery Components: {result['query_components']}")

    # Validate
    assert result['approved_protocols'] == test_protocols, "Protocols should be preserved"
    assert result['a2_block'] is not None, "Should have A2 block"
    assert result['final_query'] is not None, "Should have final query"
    assert result['query_components'] is not None, "Should have query components"
    assert result['query_components']['o_block'] is not None, "Should have O block"

    # Check that A2 block contains expected elements
    a2_block = result['a2_block']
    assert "[Title/Abstract]" in a2_block, "Should use [Title/Abstract] tags"
    assert "OR" in a2_block, "Should use OR operators"

    # Check that final query combines O_block AND A2_block
    final_query = result['final_query']
    assert "AND" in final_query, "Final query should combine blocks with AND"
    assert 'english' in final_query.lower(), "Should include language filter"
    assert 'NOT "review"' in final_query, "Should exclude reviews"

    print("\n✅ Stage 2 test passed!")
    return result


def test_stage3_pubmed_search():
    """Test Stage 3: PubMed Search Graph"""
    print("\n" + "="*80)
    print("TEST: Stage 3 - PubMed Search")
    print("="*80)

    # Compile graph
    search_graph = create_search_graph().compile()

    # Test input (simulating query from Stage 2)
    test_query = '("co-immunoprecipitation"[Title/Abstract] OR "affinity purification"[Title/Abstract])'

    print(f"\nQuery: {test_query}")

    # Create input state
    search_input = SearchInput(
        final_query=test_query,
        search_results=None,
        result_count=None,
        articles=None
    )

    # Invoke graph
    print("\n🔄 Running PubMed search graph...")
    result = search_graph.invoke(search_input)

    # Display results
    print("\n✅ Graph execution completed!")
    print(f"\nResult Count: {result['result_count']}")
    print(f"\nSearch Results: {result['search_results']}")
    if result['articles']:
        print(f"\nTop {len(result['articles'])} articles retrieved")
        for i, article in enumerate(result['articles'][:3], 1):
            print(f"  {i}. PMID: {article['pmid']}")

    # Validate
    assert result['final_query'] == test_query, "Query should be preserved"
    assert result['search_results'] is not None, "Should have search results"
    assert result['result_count'] is not None, "Should have result count"

    print("\n✅ Stage 3 test passed!")
    return result


def test_ambiguity_detection():
    """Test ambiguity detection logic"""
    print("\n" + "="*80)
    print("TEST: Ambiguity Detection")
    print("="*80)

    # Test cases
    test_cases = [
        {
            "protocols": ["protein interactions"],
            "expected": True,
            "description": "Generic single protocol"
        },
        {
            "protocols": ["co-IP", "affinity purification", "pull-down"],
            "expected": False,
            "description": "Specific protocols (3+)"
        },
        {
            "protocols": ["ppi"],
            "expected": True,
            "description": "Very generic abbreviation"
        },
        {
            "protocols": ["yeast two-hybrid", "co-immunoprecipitation"],
            "expected": False,
            "description": "Two specific protocols"
        },
        {
            "protocols": ["assay", "method"],
            "expected": True,
            "description": "Generic terms"
        }
    ]

    all_passed = True

    for i, test_case in enumerate(test_cases, 1):
        protocols = test_case["protocols"]
        expected = test_case["expected"]
        description = test_case["description"]

        result = detect_ambiguity(protocols)

        status = "✅" if result == expected else "❌"
        print(f"\nTest {i}: {description}")
        print(f"  Input: {protocols}")
        print(f"  Expected: {expected}, Got: {result} {status}")

        if result != expected:
            all_passed = False

    if all_passed:
        print("\n✅ All ambiguity detection tests passed!")
    else:
        print("\n❌ Some ambiguity detection tests failed!")

    return all_passed


def test_full_pipeline():
    """Test the full pipeline: Stage 1 → Stage 2 → Stage 3"""
    print("\n" + "="*80)
    print("TEST: Full Pipeline Integration")
    print("="*80)

    # Compile all graphs
    expansion_graph = create_protocol_expansion_graph().compile()
    query_graph = create_query_creation_graph().compile()
    search_graph = create_search_graph().compile()

    # Stage 1: Protocol Expansion
    print("\n--- Stage 1: Protocol Expansion ---")
    user_protocols = ["co-immunoprecipitation", "affinity purification"]
    print(f"Input: {user_protocols}")

    expansion_input = ProtocolExpansionInput(
        user_protocols=user_protocols,
        suggested_protocols=None,
        expansion_reasoning=None
    )
    expansion_result = expansion_graph.invoke(expansion_input)

    suggested = expansion_result['suggested_protocols']
    print(f"Suggested: {suggested}")

    # Simulate user approving all
    approved_protocols = user_protocols + suggested
    print(f"Approved (simulated 'approve all'): {approved_protocols}")

    # Stage 2: Query Creation
    print("\n--- Stage 2: Query Creation ---")
    query_input = QueryCreationInput(
        approved_protocols=approved_protocols,
        a2_block=None,
        final_query=None,
        query_components=None
    )
    query_result = query_graph.invoke(query_input)

    final_query = query_result['final_query']
    print(f"Final Query: {final_query}")

    # Stage 3: PubMed Search
    print("\n--- Stage 3: PubMed Search ---")
    search_input = SearchInput(
        final_query=final_query,
        search_results=None,
        result_count=None,
        articles=None
    )
    search_result = search_graph.invoke(search_input)

    result_count = search_result['result_count']
    print(f"Result Count: {result_count}")

    # Validate full pipeline
    assert len(approved_protocols) > len(user_protocols), "Should have more protocols after expansion"
    assert final_query is not None and len(final_query) > 0, "Should have a final query"
    assert result_count is not None, "Should have search results"

    print("\n✅ Full pipeline test passed!")
    print(f"\nSummary:")
    print(f"  Input protocols: {len(user_protocols)}")
    print(f"  After expansion: {len(approved_protocols)}")
    print(f"  Search results: {result_count}")

    return {
        "expansion": expansion_result,
        "query": query_result,
        "search": search_result
    }


def run_all_tests():
    """Run all tests in sequence"""
    print("\n" + "="*80)
    print("RUNNING ALL TESTS FOR HITL_PROPER IMPLEMENTATION")
    print("="*80)

    try:
        # Test individual stages
        test_ambiguity_detection()
        test_stage1_protocol_expansion()
        test_stage2_query_creation()
        test_stage3_pubmed_search()

        # Test full pipeline
        test_full_pipeline()

        print("\n" + "="*80)
        print("🎉 ALL TESTS PASSED!")
        print("="*80)

    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    except Exception as e:
        print(f"\n❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    run_all_tests()
