#!/usr/bin/env python3
"""
PubMed Base Filter Test with PPI Keywords
==========================================
Test script using the O_block (base filter) combined with user-provided
protein-protein interaction keywords.

Usage:
    python test_base_filter.py                     # Base filter only
    python test_base_filter.py "BRCA1"             # Single protein
    python test_base_filter.py "BRCA1 BRCA2"       # Multiple proteins (AND)
    python test_base_filter.py "p53" --field mesh  # Search MeSH terms
    python test_base_filter.py "kinase" --or       # Multiple terms with OR

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import argparse
from Bio import Entrez

Entrez.email = "ekwame001@gmail.com"
# Entrez.api_key = "YOUR_NCBI_API_KEY"  # uncomment for higher rate limits


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


# =============================================================================
# O: BASE FILTER ONLY
# English language, exclude unwanted publication types
# =============================================================================

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


def build_ppi_keyword_query(keywords, field="All Fields", use_or=False):
    """
    Build a PPI keyword query from user input.

    Parameters
    ----------
    keywords : list of str
        Protein or interaction keywords to search
    field : str
        PubMed field to search in:
        - "All Fields" (default)
        - "Title/Abstract" (tiab)
        - "Title" (ti)
        - "MeSH Terms" (mesh)
    use_or : bool
        If True, combine keywords with OR. If False, use AND.

    Returns
    -------
    str
        Formatted PubMed query for the keywords
    """
    if not keywords:
        return None

    # Map friendly field names to PubMed syntax
    field_map = {
        "all": "All Fields",
        "all fields": "All Fields",
        "tiab": "Title/Abstract",
        "title/abstract": "Title/Abstract",
        "ti": "Title",
        "title": "Title",
        "mesh": "MeSH Terms",
        "mesh terms": "MeSH Terms",
    }

    pubmed_field = field_map.get(field.lower(), field)

    # Build query terms
    terms = [f'"{kw}"[{pubmed_field}]' for kw in keywords]

    # Combine with AND or OR
    operator = " OR " if use_or else " AND "
    keyword_query = f"({operator.join(terms)})"

    return keyword_query


def build_full_query(keywords=None, field="All Fields", use_or=False):
    """
    Build the full query combining O_block with optional keywords.

    Parameters
    ----------
    keywords : list of str, optional
        Protein or interaction keywords
    field : str
        PubMed field to search
    use_or : bool
        Combine keywords with OR instead of AND

    Returns
    -------
    str
        Complete PubMed query
    """
    if not keywords:
        return O_block

    keyword_query = build_ppi_keyword_query(keywords, field, use_or)
    return f"({O_block}) AND ({keyword_query})"


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Search PubMed with base filter and optional PPI keywords",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Base filter only (all English research papers)
  python test_base_filter.py

  # Search for a single protein
  python test_base_filter.py "BRCA1"

  # Search for protein interaction (AND logic)
  python test_base_filter.py "BRCA1" "BRCA2"

  # Search for any of multiple proteins (OR logic)
  python test_base_filter.py "BRCA1" "BRCA2" "TP53" --or
#Would generate:
#  ("BRCA1"[All Fields] OR "BRCA2"[All Fields] OR "TP53"[All Fields])
##  If you want to be more specific, you would add --field:
  ## Search only in title/abstract
  ##python test_base_filter.py "BRCA1" "BRCA2" "TP53" --or --field tiab
  ## Search only in MeSH terms
  ##python test_base_filter.py "BRCA1" "BRCA2" "TP53" --or --field mesh

  
  # Search in title/abstract only
  python test_base_filter.py "kinase" --field tiab

  # Search MeSH terms
  python test_base_filter.py "protein interaction" --field mesh

  # Combine protein with interaction term
  python test_base_filter.py "BRCA1" "protein interaction"
        """
    )

    parser.add_argument(
        "keywords",
        nargs="*",
        help="Protein names or interaction keywords to search"
    )

    parser.add_argument(
        "--field", "-f",
        default="All Fields",
        choices=["all", "tiab", "title", "mesh"],
        help="Field to search: all (default), tiab (title/abstract), title, mesh"
    )

    parser.add_argument(
        "--or",
        dest="use_or",
        action="store_true",
        help="Combine keywords with OR instead of AND"
    )

    parser.add_argument(
        "--retmax", "-n",
        type=int,
        default=10,
        help="Number of PMIDs to retrieve (default: 10)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    print("=" * 80)
    print("PubMed Base Filter Test with PPI Keywords")
    print("=" * 80)

    # Show base filter info
    print("\nO_BLOCK (Base Filter): English language papers, excluding:")
    print("  - Meta-analyses, Reviews, Retracted publications")
    print("  - Errata, Clinical trials (all phases)")

    # Show keyword info
    if args.keywords:
        operator = "OR" if args.use_or else "AND"
        print(f"\nKEYWORDS: {args.keywords}")
        print(f"FIELD: {args.field}")
        print(f"OPERATOR: {operator}")

    # Build and show query
    full_query = build_full_query(
        keywords=args.keywords,
        field=args.field,
        use_or=args.use_or
    )

    print()
    print("-" * 80)
    print("QUERY:")
    print("-" * 80)
    print(full_query)

    # Run search
    print()
    print("-" * 80)
    print("RESULTS:")
    print("-" * 80)

    res = search_pubmed(full_query, retmax=args.retmax)
    print(f"Total matches: {res['count']:,}")
    print(f"First {len(res['ids'])} PMIDs: {res['ids']}")
    print()
    print("=" * 80)