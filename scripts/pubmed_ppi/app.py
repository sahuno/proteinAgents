#!/usr/bin/env python3
"""
PPI Search Agent - Streamlit GUI
=================================
A web-based graphical interface for the PPI Search Agent.

Usage:
    streamlit run app.py

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import streamlit as st
import pandas as pd
import json
from datetime import datetime
from io import StringIO

# Import from local modules
try:
    from ppi_search_agent import (
        PPISearchAgent,
        RuleBasedExtractor,
        LLMExtractor,
        PPIQueryBuilder,
        SearchResult,
        ExtractedTerms,
        ANTHROPIC_AVAILABLE
    )
except ImportError:
    # If running from different directory
    import sys
    sys.path.insert(0, '.')
    from ppi_search_agent import (
        PPISearchAgent,
        RuleBasedExtractor,
        LLMExtractor,
        PPIQueryBuilder,
        SearchResult,
        ExtractedTerms,
        ANTHROPIC_AVAILABLE
    )


# =============================================================================
# Page Configuration
# =============================================================================

st.set_page_config(
    page_title="PPI Search Agent",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)


# =============================================================================
# Custom CSS
# =============================================================================

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #666;
        margin-bottom: 2rem;
    }
    .result-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    .metric-card {
        background-color: #ffffff;
        padding: 1rem;
        border-radius: 0.5rem;
        border: 1px solid #e0e0e0;
        text-align: center;
    }
    .stTextArea textarea {
        font-family: monospace;
    }
</style>
""", unsafe_allow_html=True)


# =============================================================================
# Session State Initialization
# =============================================================================

if 'search_history' not in st.session_state:
    st.session_state.search_history = []

if 'current_result' not in st.session_state:
    st.session_state.current_result = None


# =============================================================================
# Sidebar
# =============================================================================

with st.sidebar:
    st.header("Settings")

    # Extraction method
    use_llm = st.checkbox(
        "Use LLM (Claude) for extraction",
        value=False,
        help="Uses Claude AI for better term extraction. Requires ANTHROPIC_API_KEY environment variable.",
        disabled=not ANTHROPIC_AVAILABLE
    )

    if not ANTHROPIC_AVAILABLE:
        st.caption("Install `anthropic` package to enable LLM mode")

    # Number of results
    retmax = st.slider(
        "Max results to retrieve",
        min_value=10,
        max_value=1000,
        value=100,
        step=10,
        help="Maximum number of PMIDs to retrieve from PubMed"
    )

    st.divider()

    # Search history
    st.header("Search History")
    if st.session_state.search_history:
        for i, hist in enumerate(reversed(st.session_state.search_history[-5:])):
            with st.expander(f"{hist['query'][:30]}...", expanded=False):
                st.caption(f"Time: {hist['time']}")
                st.caption(f"Results: {hist['count']:,}")
                if st.button("Load", key=f"load_{i}"):
                    st.session_state.current_result = hist['result']
                    st.rerun()
    else:
        st.caption("No searches yet")

    st.divider()

    # About section
    st.header("About")
    st.markdown("""
    **PPI Search Agent** helps you search PubMed for protein-protein interaction literature using natural language queries.

    The agent:
    1. Extracts search terms from your query
    2. Maps terms to MeSH vocabulary
    3. Builds optimized PubMed queries
    4. Returns relevant publications

    [Documentation](./README.md)
    """)


# =============================================================================
# Main Content
# =============================================================================

# Header
st.markdown('<p class="main-header">PPI Search Agent</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Search PubMed for protein-protein interaction literature using natural language</p>', unsafe_allow_html=True)

# Search input
col1, col2 = st.columns([4, 1])

with col1:
    user_query = st.text_input(
        "What are you looking for?",
        placeholder="e.g., BRCA1 interactions with DNA repair proteins in human cells",
        help="Describe what you're looking for in natural language"
    )

with col2:
    st.write("")  # Spacer
    st.write("")  # Spacer
    search_clicked = st.button("Search", type="primary", use_container_width=True)

# Example queries
with st.expander("Example queries", expanded=False):
    example_cols = st.columns(3)
    examples = [
        "BRCA1 interactions with DNA repair proteins",
        "p53 MDM2 binding using co-IP",
        "kinase-substrate interactions in cancer from 2020-2024",
        "TP53 protein complex in human cells",
        "cryo-EM structure of protein complexes",
        "yeast two-hybrid screen for interactions"
    ]
    for i, example in enumerate(examples):
        with example_cols[i % 3]:
            if st.button(example, key=f"ex_{i}", use_container_width=True):
                user_query = example
                search_clicked = True


# =============================================================================
# Search Execution
# =============================================================================

if search_clicked and user_query:
    with st.spinner("Searching PubMed..."):
        try:
            # Initialize agent
            agent = PPISearchAgent(use_llm=use_llm)
            extraction_method = "LLM (Claude)" if agent.using_llm else "Rule-based"

            # Execute search
            result = agent.search(user_query, retmax=retmax)

            # Store in session state
            st.session_state.current_result = {
                'query': user_query,
                'result': result,
                'extraction_method': extraction_method,
                'time': datetime.now().strftime("%H:%M:%S")
            }

            # Add to history
            st.session_state.search_history.append({
                'query': user_query,
                'count': result.total_count,
                'time': datetime.now().strftime("%H:%M:%S"),
                'result': st.session_state.current_result
            })

            st.success(f"Found {result.total_count:,} matching articles!")

        except Exception as e:
            st.error(f"Search failed: {str(e)}")
            st.session_state.current_result = None


# =============================================================================
# Results Display
# =============================================================================

if st.session_state.current_result:
    current = st.session_state.current_result
    result = current['result']
    terms = result.extracted_terms

    st.divider()

    # Metrics row
    metric_cols = st.columns(4)
    with metric_cols[0]:
        st.metric("Total Matches", f"{result.total_count:,}")
    with metric_cols[1]:
        st.metric("PMIDs Retrieved", len(result.pmids))
    with metric_cols[2]:
        st.metric("Extraction Method", current['extraction_method'])
    with metric_cols[3]:
        st.metric("Proteins Found", len(terms.proteins))

    # Tabs for different views
    tab1, tab2, tab3, tab4 = st.tabs(["Extracted Terms", "PubMed Query", "Results", "Export"])

    # Tab 1: Extracted Terms
    with tab1:
        st.subheader("Extracted Search Terms")

        term_cols = st.columns(3)

        with term_cols[0]:
            st.markdown("**Proteins**")
            if terms.proteins:
                for p in terms.proteins:
                    st.code(p)
            else:
                st.caption("None detected")

            st.markdown("**Interaction Types**")
            if terms.interaction_types:
                for it in terms.interaction_types:
                    st.code(it)
            else:
                st.caption("None detected")

        with term_cols[1]:
            st.markdown("**Methods**")
            if terms.methods:
                for m in terms.methods:
                    st.code(m)
            else:
                st.caption("None detected")

            st.markdown("**Organisms**")
            if terms.organisms:
                for o in terms.organisms:
                    st.code(o)
            else:
                st.caption("None detected")

        with term_cols[2]:
            st.markdown("**Diseases**")
            if terms.diseases:
                for d in terms.diseases:
                    st.code(d)
            else:
                st.caption("None detected")

            st.markdown("**Year Range**")
            if terms.year_range:
                st.code(f"{terms.year_range[0]} - {terms.year_range[1]}")
            else:
                st.caption("Not specified")

            st.markdown("**Keywords**")
            if terms.keywords:
                for k in terms.keywords:
                    st.code(k)
            else:
                st.caption("None")

    # Tab 2: PubMed Query
    with tab2:
        st.subheader("Generated PubMed Query")
        st.caption("Copy this query to use directly in PubMed web interface")

        st.text_area(
            "Query",
            value=result.query,
            height=200,
            label_visibility="collapsed"
        )

        # Copy button hint
        st.caption("Use Ctrl+A to select all, then Ctrl+C to copy")

        # Link to PubMed
        st.markdown(f"[Open PubMed](https://pubmed.ncbi.nlm.nih.gov/)")

    # Tab 3: Results (PMIDs)
    with tab3:
        st.subheader(f"Retrieved PMIDs ({len(result.pmids)})")

        if result.pmids:
            # Display as clickable links
            pmid_cols = st.columns(5)
            for i, pmid in enumerate(result.pmids):
                with pmid_cols[i % 5]:
                    st.markdown(f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
        else:
            st.info("No PMIDs retrieved. Try broadening your search.")

    # Tab 4: Export
    with tab4:
        st.subheader("Export Results")

        export_cols = st.columns(3)

        # JSON export
        with export_cols[0]:
            st.markdown("**Full Results (JSON)**")
            json_data = {
                "timestamp": datetime.now().isoformat(),
                "user_query": current['query'],
                "extraction_method": current['extraction_method'],
                "extracted_terms": terms.to_dict(),
                "pubmed_query": result.query,
                "total_count": result.total_count,
                "pmids": result.pmids
            }
            st.download_button(
                label="Download JSON",
                data=json.dumps(json_data, indent=2),
                file_name=f"ppi_search_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json",
                use_container_width=True
            )

        # Query export
        with export_cols[1]:
            st.markdown("**PubMed Query (TXT)**")
            query_content = f"""# PPI Search Agent Query
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# User query: {current['query']}
# Total matches: {result.total_count:,}

{result.query}
"""
            st.download_button(
                label="Download Query",
                data=query_content,
                file_name=f"ppi_query_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain",
                use_container_width=True
            )

        # PMIDs export
        with export_cols[2]:
            st.markdown("**PMIDs List (TXT)**")
            pmids_content = "\n".join(result.pmids)
            st.download_button(
                label="Download PMIDs",
                data=pmids_content,
                file_name=f"pmids_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain",
                use_container_width=True
            )

        # CSV export (for spreadsheet use)
        st.markdown("---")
        st.markdown("**PMIDs Table (CSV)**")
        df = pd.DataFrame({
            'PMID': result.pmids,
            'PubMed_Link': [f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in result.pmids]
        })
        st.download_button(
            label="Download CSV",
            data=df.to_csv(index=False),
            file_name=f"pmids_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
            mime="text/csv"
        )


# =============================================================================
# Footer
# =============================================================================

st.divider()
st.caption("PPI Search Agent | Author: Samuel Ahuno (ekwame001@gmail.com) | 2025")