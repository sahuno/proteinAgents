#!/usr/bin/env python3
"""
PubMed PPI Search Pipeline
==========================
Main script to search PubMed for protein-protein interaction literature,
fetch article metadata, and export results.

Usage:
    python main.py --help
    python main.py --count-only
    python main.py --fetch 100 --output results
    python main.py --save-query

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import argparse
import sys
from pathlib import Path

from query_builder import PPIQueryBuilder
from fetcher import PubMedFetcher
from export import PubMedExporter


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Search PubMed for protein-protein interaction literature",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get count of matching articles only
  python main.py --count-only

  # Search and fetch first 100 articles
  python main.py --fetch 100

  # Fetch 1000 articles and save to custom directory
  python main.py --fetch 1000 --output my_results

  # Save query to file without searching
  python main.py --save-query

  # Full pipeline with all exports
  python main.py --fetch 500 --output results --csv --json --sqlite --report
        """
    )

    parser.add_argument(
        "--email",
        default="ekwame001@gmail.com",
        help="Email for NCBI API (required by NCBI)"
    )

    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI API key for higher rate limits"
    )

    parser.add_argument(
        "--count-only",
        action="store_true",
        help="Only get count of matching articles"
    )

    parser.add_argument(
        "--fetch",
        type=int,
        metavar="N",
        help="Number of articles to fetch (default: don't fetch)"
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for fetching (default: 100)"
    )

    parser.add_argument(
        "--output",
        default="output",
        help="Output directory (default: output)"
    )

    parser.add_argument(
        "--save-query",
        action="store_true",
        help="Save query to files"
    )

    parser.add_argument(
        "--csv",
        action="store_true",
        help="Export articles to CSV"
    )

    parser.add_argument(
        "--json",
        action="store_true",
        help="Export articles to JSON"
    )

    parser.add_argument(
        "--sqlite",
        action="store_true",
        help="Export articles to SQLite database"
    )

    parser.add_argument(
        "--pmids",
        action="store_true",
        help="Save list of PMIDs to file"
    )

    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate summary report"
    )

    parser.add_argument(
        "--simple-query",
        action="store_true",
        help="Use simplified query without complex exclusion logic"
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=True,
        help="Verbose output (default: True)"
    )

    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress messages"
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    verbose = args.verbose and not args.quiet

    # Initialize components
    builder = PPIQueryBuilder()
    fetcher = PubMedFetcher(email=args.email, api_key=args.api_key)
    exporter = PubMedExporter(output_dir=args.output)

    # Build query
    if args.simple_query:
        query = builder.build_simple_query()
        if verbose:
            print("Using simplified query (no complex exclusion logic)")
    else:
        query = builder.build_query()

    if verbose:
        print("=" * 80)
        print("PubMed PPI Search Pipeline")
        print("=" * 80)
        print(f"\nQuery logic: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)")
        print(f"Output directory: {args.output}")

    # Save query if requested
    if args.save_query:
        blocks = {
            'O': builder.O_block,
            'A1': builder.A1_block,
            'A2': builder.A2_block,
            'B': builder.B_block
        }
        exporter.save_query(
            query=query,
            filename="ppi_query",
            blocks=blocks,
            metadata={
                'version': '1.0',
                'description': 'PPI literature search',
                'simple_query': args.simple_query
            }
        )

    # Count only mode
    if args.count_only:
        count = fetcher.get_article_count(query)
        print(f"\nTotal articles matching query: {count:,}")
        return 0

    # If no fetch requested, just show count and exit
    if args.fetch is None:
        if not args.save_query:
            count = fetcher.get_article_count(query)
            print(f"\nTotal articles matching query: {count:,}")
            print("\nUse --fetch N to retrieve articles, or --save-query to save the query.")
        return 0

    # Search and fetch
    if verbose:
        print(f"\nSearching PubMed...")

    search_result = fetcher.search(query, retmax=args.fetch)

    if verbose:
        print(f"Total matches: {search_result.total_count:,}")
        print(f"PMIDs to fetch: {len(search_result.pmids)}")

    # Fetch articles
    articles = fetcher.fetch_articles(
        search_result=search_result,
        batch_size=args.batch_size,
        max_articles=args.fetch,
        verbose=verbose
    )

    if not articles:
        print("No articles fetched.")
        return 1

    # Export results
    if verbose:
        print(f"\nExporting {len(articles)} articles...")

    # Default to CSV if no format specified
    if not any([args.csv, args.json, args.sqlite, args.pmids]):
        args.csv = True

    if args.csv:
        exporter.save_articles_csv(articles, filename="ppi_articles")

    if args.json:
        exporter.save_articles_json(
            articles,
            filename="ppi_articles",
            search_result=search_result
        )

    if args.sqlite:
        exporter.save_articles_sqlite(
            articles,
            filename="pubmed_ppi.db",
            search_result=search_result
        )

    if args.pmids:
        exporter.save_pmids(search_result.pmids, filename="ppi_pmids")

    if args.report:
        exporter.generate_summary_report(
            articles,
            search_result=search_result,
            filename="ppi_summary"
        )

    if verbose:
        print("\nDone!")

    return 0


if __name__ == "__main__":
    sys.exit(main())