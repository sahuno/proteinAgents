"""
PubMed Export Utilities for Protein-Protein Interactions
=========================================================
Save queries and article metadata to various formats.

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import csv
import json
import sqlite3
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional

try:
    from .fetcher import Article, SearchResult
except ImportError:
    from fetcher import Article, SearchResult


class PubMedExporter:
    """
    Export PubMed queries and article data to various formats.

    Parameters
    ----------
    output_dir : str or Path
        Base directory for output files
    """

    def __init__(self, output_dir: str = "output"):
        """Initialize exporter with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.queries_dir = self.output_dir / "queries"
        self.results_dir = self.output_dir / "results"
        self.queries_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)

    # =========================================================================
    # Query Export Methods
    # =========================================================================

    def save_query(
        self,
        query: str,
        filename: str = "ppi_query",
        blocks: Optional[Dict[str, str]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Path]:
        """
        Save query in multiple formats.

        Parameters
        ----------
        query : str
            The full PubMed query string
        filename : str
            Base filename (without extension)
        blocks : dict, optional
            Individual query blocks (O, A1, A2, B) for documentation
        metadata : dict, optional
            Additional metadata (description, version, etc.)

        Returns
        -------
        dict
            Paths to saved files
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        saved_files = {}

        # 1. Raw query for PubMed web copy-paste
        raw_path = self.queries_dir / f"{filename}_raw.txt"
        with open(raw_path, 'w') as f:
            f.write(query)
        saved_files['raw'] = raw_path

        # 2. Readable/formatted version
        readable_path = self.queries_dir / f"{filename}_readable.txt"
        with open(readable_path, 'w') as f:
            f.write(self._format_query_readable(query, blocks))
        saved_files['readable'] = readable_path

        # 3. JSON with full metadata
        json_path = self.queries_dir / f"{filename}_metadata.json"
        json_data = {
            'query': query,
            'timestamp': timestamp,
            'date': datetime.now().isoformat(),
            'blocks': blocks or {},
            'metadata': metadata or {},
            'logic': '(O AND (A1 OR A2)) NOT (A1 AND B NOT A2)',
            'description': 'PubMed query for protein-protein interaction literature'
        }
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        saved_files['json'] = json_path

        print(f"Query saved to:")
        for fmt, path in saved_files.items():
            print(f"  [{fmt}] {path}")

        return saved_files

    def _format_query_readable(
        self,
        query: str,
        blocks: Optional[Dict[str, str]] = None
    ) -> str:
        """Format query for human readability."""
        lines = [
            "=" * 80,
            "PubMed Query for Protein-Protein Interactions",
            "=" * 80,
            "",
            "SEARCH STRATEGY: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)",
            "",
            "Where:",
            "  O  = Base filter (English, exclude reviews/trials/etc.)",
            "  A1 = Inclusion concepts (protein complexes, interactions)",
            "  A2 = Inclusion methods (co-IP, affinity purification, etc.)",
            "  B  = Exclusion methods (crosslinking, two-hybrid, etc.)",
            "",
            "-" * 80,
        ]

        if blocks:
            for block_name, block_query in blocks.items():
                lines.extend([
                    f"\n{block_name} BLOCK:",
                    "-" * 40,
                    block_query,
                    ""
                ])

        lines.extend([
            "-" * 80,
            "FULL QUERY (copy-paste to PubMed):",
            "-" * 80,
            "",
            query,
            "",
            "=" * 80,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 80
        ])

        return '\n'.join(lines)

    # =========================================================================
    # Article Export Methods
    # =========================================================================

    def save_articles_csv(
        self,
        articles: List[Article],
        filename: str = "articles",
        include_timestamp: bool = True
    ) -> Path:
        """
        Save articles to CSV format.

        Parameters
        ----------
        articles : list of Article
            Articles to save
        filename : str
            Base filename (without extension)
        include_timestamp : bool
            Add timestamp to filename

        Returns
        -------
        Path
            Path to saved CSV file
        """
        if include_timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filepath = self.results_dir / f"{filename}_{timestamp}.csv"
        else:
            filepath = self.results_dir / f"{filename}.csv"

        fieldnames = [
            'pmid', 'title', 'abstract', 'authors', 'journal',
            'year', 'month', 'doi', 'mesh_terms', 'keywords', 'publication_types'
        ]

        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for article in articles:
                writer.writerow(article.to_dict())

        print(f"Saved {len(articles)} articles to {filepath}")
        return filepath

    def save_articles_json(
        self,
        articles: List[Article],
        filename: str = "articles",
        search_result: Optional[SearchResult] = None,
        include_timestamp: bool = True
    ) -> Path:
        """
        Save articles to JSON format with metadata.

        Parameters
        ----------
        articles : list of Article
            Articles to save
        filename : str
            Base filename
        search_result : SearchResult, optional
            Include search metadata
        include_timestamp : bool
            Add timestamp to filename

        Returns
        -------
        Path
            Path to saved JSON file
        """
        if include_timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filepath = self.results_dir / f"{filename}_{timestamp}.json"
        else:
            filepath = self.results_dir / f"{filename}.json"

        data = {
            'metadata': {
                'exported_at': datetime.now().isoformat(),
                'total_articles': len(articles),
            },
            'articles': [article.to_dict() for article in articles]
        }

        if search_result:
            data['metadata']['query'] = search_result.query
            data['metadata']['total_matches'] = search_result.total_count
            data['metadata']['pmids_retrieved'] = len(search_result.pmids)

        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        print(f"Saved {len(articles)} articles to {filepath}")
        return filepath

    def save_articles_sqlite(
        self,
        articles: List[Article],
        filename: str = "pubmed_ppi.db",
        search_result: Optional[SearchResult] = None
    ) -> Path:
        """
        Save articles to SQLite database.

        Parameters
        ----------
        articles : list of Article
            Articles to save
        filename : str
            Database filename
        search_result : SearchResult, optional
            Include search metadata

        Returns
        -------
        Path
            Path to SQLite database
        """
        filepath = self.results_dir / filename

        conn = sqlite3.connect(filepath)
        cursor = conn.cursor()

        # Create tables
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS articles (
                pmid TEXT PRIMARY KEY,
                title TEXT,
                abstract TEXT,
                authors TEXT,
                journal TEXT,
                year TEXT,
                month TEXT,
                doi TEXT,
                mesh_terms TEXT,
                keywords TEXT,
                publication_types TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS searches (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                query TEXT,
                total_count INTEGER,
                pmids_retrieved INTEGER,
                search_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')

        # Insert search metadata
        if search_result:
            cursor.execute(
                'INSERT INTO searches (query, total_count, pmids_retrieved) VALUES (?, ?, ?)',
                (search_result.query, search_result.total_count, len(search_result.pmids))
            )

        # Insert articles (upsert)
        for article in articles:
            data = article.to_dict()
            cursor.execute('''
                INSERT OR REPLACE INTO articles
                (pmid, title, abstract, authors, journal, year, month, doi,
                 mesh_terms, keywords, publication_types)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                data['pmid'], data['title'], data['abstract'], data['authors'],
                data['journal'], data['year'], data['month'], data['doi'],
                data['mesh_terms'], data['keywords'], data['publication_types']
            ))

        conn.commit()
        conn.close()

        print(f"Saved {len(articles)} articles to {filepath}")
        return filepath

    def save_pmids(
        self,
        pmids: List[str],
        filename: str = "pmids",
        include_timestamp: bool = True
    ) -> Path:
        """
        Save PMIDs to a text file (one per line).

        Parameters
        ----------
        pmids : list of str
            List of PMIDs
        filename : str
            Base filename
        include_timestamp : bool
            Add timestamp to filename

        Returns
        -------
        Path
            Path to saved file
        """
        if include_timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filepath = self.results_dir / f"{filename}_{timestamp}.txt"
        else:
            filepath = self.results_dir / f"{filename}.txt"

        with open(filepath, 'w') as f:
            f.write('\n'.join(pmids))

        print(f"Saved {len(pmids)} PMIDs to {filepath}")
        return filepath

    # =========================================================================
    # Summary/Report Methods
    # =========================================================================

    def generate_summary_report(
        self,
        articles: List[Article],
        search_result: Optional[SearchResult] = None,
        filename: str = "summary_report"
    ) -> Path:
        """
        Generate a summary report of the search results.

        Parameters
        ----------
        articles : list of Article
            Fetched articles
        search_result : SearchResult, optional
            Search metadata
        filename : str
            Report filename

        Returns
        -------
        Path
            Path to report file
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filepath = self.results_dir / f"{filename}_{timestamp}.txt"

        # Compute statistics
        years = [a.year for a in articles if a.year]
        journals = {}
        all_mesh = []
        for article in articles:
            if article.journal:
                journals[article.journal] = journals.get(article.journal, 0) + 1
            all_mesh.extend(article.mesh_terms)

        mesh_counts = {}
        for term in all_mesh:
            mesh_counts[term] = mesh_counts.get(term, 0) + 1

        top_journals = sorted(journals.items(), key=lambda x: x[1], reverse=True)[:10]
        top_mesh = sorted(mesh_counts.items(), key=lambda x: x[1], reverse=True)[:20]

        lines = [
            "=" * 80,
            "PubMed PPI Search Summary Report",
            "=" * 80,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
        ]

        if search_result:
            lines.extend([
                "SEARCH STATISTICS",
                "-" * 40,
                f"Total matches in PubMed: {search_result.total_count:,}",
                f"PMIDs retrieved: {len(search_result.pmids):,}",
                f"Articles fetched: {len(articles):,}",
                "",
            ])

        lines.extend([
            "ARTICLE STATISTICS",
            "-" * 40,
            f"Total articles analyzed: {len(articles)}",
            f"Year range: {min(years) if years else 'N/A'} - {max(years) if years else 'N/A'}",
            f"Unique journals: {len(journals)}",
            "",
            "TOP 10 JOURNALS",
            "-" * 40,
        ])

        for journal, count in top_journals:
            lines.append(f"  {count:4d} | {journal[:60]}")

        lines.extend([
            "",
            "TOP 20 MeSH TERMS",
            "-" * 40,
        ])

        for term, count in top_mesh:
            lines.append(f"  {count:4d} | {term}")

        lines.extend([
            "",
            "=" * 80
        ])

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

        print(f"Summary report saved to {filepath}")
        return filepath


if __name__ == "__main__":
    # Example usage
    exporter = PubMedExporter(output_dir="output")

    # Example query save
    from query_builder import PPIQueryBuilder
    builder = PPIQueryBuilder()
    query = builder.build_query()

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
        metadata={'version': '1.0', 'description': 'PPI literature search'}
    )