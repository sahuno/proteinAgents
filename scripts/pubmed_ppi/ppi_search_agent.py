#!/usr/bin/env python3
"""
PPI Search Agent
================
An AI-powered agent that helps users search PubMed for protein-protein
interaction literature by understanding natural language queries.

The agent:
1. Asks the user what they're looking for
2. Extracts relevant search terms (proteins, methods, organisms, etc.)
3. Builds appropriate PubMed queries combining O_block with extracted terms
4. Executes the search and returns results

Usage:
    python ppi_search_agent.py                    # Interactive mode
    python ppi_search_agent.py --query "..."      # Direct query mode

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import argparse
import json
import os
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Any
from Bio import Entrez

# Optional: Import anthropic for Claude API
try:
    import anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False

Entrez.email = "ekwame001@gmail.com"


# =============================================================================
# O: BASE FILTER
# =============================================================================

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
# Data Classes
# =============================================================================

@dataclass
class ExtractedTerms:
    """Container for extracted search terms from user query."""
    proteins: List[str] = field(default_factory=list)
    interaction_types: List[str] = field(default_factory=list)
    methods: List[str] = field(default_factory=list)
    organisms: List[str] = field(default_factory=list)
    diseases: List[str] = field(default_factory=list)
    keywords: List[str] = field(default_factory=list)
    year_range: Optional[tuple] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            'proteins': self.proteins,
            'interaction_types': self.interaction_types,
            'methods': self.methods,
            'organisms': self.organisms,
            'diseases': self.diseases,
            'keywords': self.keywords,
            'year_range': self.year_range
        }

    def is_empty(self) -> bool:
        return not any([
            self.proteins, self.interaction_types, self.methods,
            self.organisms, self.diseases, self.keywords
        ])


@dataclass
class SearchResult:
    """Container for PubMed search results."""
    query: str
    total_count: int
    pmids: List[str]
    extracted_terms: ExtractedTerms


# =============================================================================
# Query Builder
# =============================================================================

class PPIQueryBuilder:
    """Build PubMed queries from extracted terms."""

    # Common MeSH terms for methods
    METHOD_MESH_MAP = {
        'co-ip': '"Immunoprecipitation"[MeSH Terms]',
        'coip': '"Immunoprecipitation"[MeSH Terms]',
        'immunoprecipitation': '"Immunoprecipitation"[MeSH Terms]',
        'pulldown': '"pulldown"[All Fields]',
        'pull-down': '"pulldown"[All Fields]',
        'affinity purification': '"Chromatography, Affinity"[MeSH Terms]',
        'mass spectrometry': '"Mass Spectrometry"[MeSH Terms]',
        'yeast two-hybrid': '"Two-Hybrid System Techniques"[MeSH Terms]',
        'y2h': '"Two-Hybrid System Techniques"[MeSH Terms]',
        'cryo-em': '"Cryoelectron Microscopy"[MeSH Terms]',
        'x-ray': '"crystallography, x ray"[MeSH Terms]',
        'crystallography': '"crystallography, x ray"[MeSH Terms]',
        'nmr': '"Nuclear Magnetic Resonance, Biomolecular"[MeSH Terms]',
        'spr': '"Surface Plasmon Resonance"[MeSH Terms]',
        'fret': '"Fluorescence Resonance Energy Transfer"[MeSH Terms]',
    }

    # Common organism MeSH terms
    ORGANISM_MESH_MAP = {
        'human': '"Humans"[MeSH Terms]',
        'humans': '"Humans"[MeSH Terms]',
        'mouse': '"Mice"[MeSH Terms]',
        'mice': '"Mice"[MeSH Terms]',
        'rat': '"Rats"[MeSH Terms]',
        'rats': '"Rats"[MeSH Terms]',
        'yeast': '"Saccharomyces cerevisiae"[MeSH Terms]',
        'e. coli': '"Escherichia coli"[MeSH Terms]',
        'drosophila': '"Drosophila"[MeSH Terms]',
        'zebrafish': '"Zebrafish"[MeSH Terms]',
        'c. elegans': '"Caenorhabditis elegans"[MeSH Terms]',
    }

    def build_query(self, terms: ExtractedTerms) -> str:
        """
        Build a PubMed query from extracted terms.

        Parameters
        ----------
        terms : ExtractedTerms
            Extracted search terms

        Returns
        -------
        str
            Complete PubMed query
        """
        query_parts = []

        # Add protein terms
        if terms.proteins:
            protein_terms = [f'"{p}"[All Fields]' for p in terms.proteins]
            query_parts.append(f"({' AND '.join(protein_terms)})")

        # Add interaction type terms
        if terms.interaction_types:
            interaction_terms = [f'"{it}"[All Fields]' for it in terms.interaction_types]
            query_parts.append(f"({' OR '.join(interaction_terms)})")

        # Add method terms (use MeSH when available)
        if terms.methods:
            method_terms = []
            for method in terms.methods:
                mesh_term = self.METHOD_MESH_MAP.get(method.lower())
                if mesh_term:
                    method_terms.append(mesh_term)
                else:
                    method_terms.append(f'"{method}"[All Fields]')
            query_parts.append(f"({' OR '.join(method_terms)})")

        # Add organism terms (use MeSH when available)
        if terms.organisms:
            org_terms = []
            for org in terms.organisms:
                mesh_term = self.ORGANISM_MESH_MAP.get(org.lower())
                if mesh_term:
                    org_terms.append(mesh_term)
                else:
                    org_terms.append(f'"{org}"[MeSH Terms]')
            query_parts.append(f"({' OR '.join(org_terms)})")

        # Add disease terms
        if terms.diseases:
            disease_terms = [f'"{d}"[MeSH Terms]' for d in terms.diseases]
            query_parts.append(f"({' OR '.join(disease_terms)})")

        # Add general keywords
        if terms.keywords:
            kw_terms = [f'"{kw}"[All Fields]' for kw in terms.keywords]
            query_parts.append(f"({' AND '.join(kw_terms)})")

        # Add year range
        if terms.year_range:
            start, end = terms.year_range
            query_parts.append(f'("{start}"[Date - Publication] : "{end}"[Date - Publication])')

        # Combine with O_block
        if query_parts:
            user_query = " AND ".join(query_parts)
            return f"({O_BLOCK}) AND ({user_query})"
        else:
            return O_BLOCK


# =============================================================================
# Term Extractor (Rule-based fallback)
# =============================================================================

class RuleBasedExtractor:
    """
    Simple rule-based term extractor as fallback when LLM is not available.
    """

    # Common interaction-related keywords
    INTERACTION_KEYWORDS = [
        'interaction', 'interactions', 'interacts', 'interacting',
        'binding', 'binds', 'bound', 'complex', 'complexes',
        'association', 'associates', 'associated',
        'phosphorylation', 'ubiquitination', 'acetylation',
        'activation', 'inhibition', 'regulation'
    ]

    # Common method keywords
    METHOD_KEYWORDS = [
        'co-ip', 'coip', 'immunoprecipitation', 'pulldown', 'pull-down',
        'mass spectrometry', 'ms/ms', 'lc-ms', 'yeast two-hybrid', 'y2h',
        'cryo-em', 'x-ray', 'crystallography', 'nmr', 'spr', 'fret',
        'affinity purification', 'ap-ms'
    ]

    # Common organism keywords
    ORGANISM_KEYWORDS = [
        'human', 'humans', 'mouse', 'mice', 'rat', 'rats',
        'yeast', 'e. coli', 'drosophila', 'zebrafish', 'c. elegans'
    ]

    def extract(self, user_input: str) -> ExtractedTerms:
        """
        Extract search terms from user input using simple rules.

        Parameters
        ----------
        user_input : str
            User's natural language query

        Returns
        -------
        ExtractedTerms
            Extracted terms
        """
        terms = ExtractedTerms()
        input_lower = user_input.lower()
        words = user_input.split()

        # Extract potential protein names (capitalized words, common patterns)
        for word in words:
            clean_word = word.strip('.,;:!?()"\'')
            # Likely protein names: all caps, or mixed case with numbers
            if (clean_word.isupper() and len(clean_word) >= 2 and
                clean_word not in ['AND', 'OR', 'NOT', 'IN', 'WITH', 'THE', 'FOR']):
                terms.proteins.append(clean_word)
            elif any(c.isupper() for c in clean_word) and any(c.isdigit() for c in clean_word):
                terms.proteins.append(clean_word)

        # Extract interaction types
        for kw in self.INTERACTION_KEYWORDS:
            if kw in input_lower:
                terms.interaction_types.append(kw)

        # Extract methods
        for method in self.METHOD_KEYWORDS:
            if method in input_lower:
                terms.methods.append(method)

        # Extract organisms
        for org in self.ORGANISM_KEYWORDS:
            if org in input_lower:
                terms.organisms.append(org)

        # Extract year range (simple pattern: YYYY-YYYY or "from YYYY to YYYY")
        import re
        year_pattern = r'(\d{4})\s*[-to]+\s*(\d{4})'
        year_match = re.search(year_pattern, user_input)
        if year_match:
            terms.year_range = (year_match.group(1), year_match.group(2))

        # If no proteins found, add remaining capitalized phrases as keywords
        if not terms.proteins:
            # Look for quoted terms
            quoted = re.findall(r'"([^"]+)"', user_input)
            terms.keywords.extend(quoted)

        return terms


# =============================================================================
# LLM-based Term Extractor
# =============================================================================

class LLMExtractor:
    """
    Extract search terms using Claude API.
    """

    SYSTEM_PROMPT = """You are a scientific literature search assistant specializing in protein-protein interactions (PPI).

Your task is to extract search terms from user queries to build PubMed searches.

Extract the following categories:
- proteins: Gene/protein names (e.g., BRCA1, p53, MDM2, TP53)
- interaction_types: Types of interactions (e.g., binding, phosphorylation, complex formation)
- methods: Experimental methods (e.g., co-IP, mass spectrometry, yeast two-hybrid, cryo-EM)
- organisms: Species (e.g., human, mouse, yeast)
- diseases: Disease associations (e.g., cancer, Alzheimer's)
- keywords: Other relevant search terms
- year_range: Publication year range if mentioned (as [start, end])

Respond ONLY with a JSON object containing these fields. Use empty lists [] for categories with no matches.
Be conservative - only extract terms explicitly mentioned or clearly implied."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize with Anthropic API key."""
        self.api_key = api_key or os.environ.get('ANTHROPIC_API_KEY')
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found")
        self.client = anthropic.Anthropic(api_key=self.api_key)

    def extract(self, user_input: str) -> ExtractedTerms:
        """
        Extract search terms using Claude.

        Parameters
        ----------
        user_input : str
            User's natural language query

        Returns
        -------
        ExtractedTerms
            Extracted terms
        """
        message = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1024,
            system=self.SYSTEM_PROMPT,
            messages=[
                {"role": "user", "content": f"Extract search terms from: {user_input}"}
            ]
        )

        # Parse JSON response
        response_text = message.content[0].text

        # Try to extract JSON from response
        try:
            # Handle potential markdown code blocks
            if '```json' in response_text:
                json_str = response_text.split('```json')[1].split('```')[0]
            elif '```' in response_text:
                json_str = response_text.split('```')[1].split('```')[0]
            else:
                json_str = response_text

            data = json.loads(json_str.strip())

            year_range = None
            if data.get('year_range') and len(data['year_range']) == 2:
                year_range = tuple(data['year_range'])

            return ExtractedTerms(
                proteins=data.get('proteins', []),
                interaction_types=data.get('interaction_types', []),
                methods=data.get('methods', []),
                organisms=data.get('organisms', []),
                diseases=data.get('diseases', []),
                keywords=data.get('keywords', []),
                year_range=year_range
            )
        except (json.JSONDecodeError, KeyError, IndexError) as e:
            print(f"Warning: Could not parse LLM response, falling back to rule-based extraction")
            print(f"Response was: {response_text}")
            return RuleBasedExtractor().extract(user_input)


# =============================================================================
# PubMed Search
# =============================================================================

def search_pubmed(query: str, retmax: int = 100) -> Dict[str, Any]:
    """Execute PubMed search."""
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
# Search Result Saver
# =============================================================================

class SearchResultSaver:
    """Save search queries and results to files."""

    def __init__(self, output_dir: str = "output/searches"):
        """Initialize with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def save(
        self,
        user_query: str,
        result: SearchResult,
        extraction_method: str
    ) -> Dict[str, Path]:
        """
        Save search query and results to files.

        Parameters
        ----------
        user_query : str
            Original user query
        result : SearchResult
            Search results
        extraction_method : str
            Method used for extraction (LLM or Rule-based)

        Returns
        -------
        dict
            Paths to saved files
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Create a sanitized filename from query
        safe_query = "".join(c if c.isalnum() or c in " -_" else "" for c in user_query)
        safe_query = safe_query[:50].strip().replace(" ", "_")

        saved_files = {}

        # 1. Save full results as JSON
        json_path = self.output_dir / f"search_{timestamp}_{safe_query}.json"
        json_data = {
            "timestamp": timestamp,
            "date": datetime.now().isoformat(),
            "user_query": user_query,
            "extraction_method": extraction_method,
            "extracted_terms": result.extracted_terms.to_dict(),
            "pubmed_query": result.query,
            "total_count": result.total_count,
            "pmids_retrieved": len(result.pmids),
            "pmids": result.pmids
        }
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2)
        saved_files['json'] = json_path

        # 2. Save raw PubMed query for copy-paste
        query_path = self.output_dir / f"query_{timestamp}_{safe_query}.txt"
        with open(query_path, 'w', encoding='utf-8') as f:
            f.write(f"# PPI Search Agent Query\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# User query: {user_query}\n")
            f.write(f"# Extraction method: {extraction_method}\n")
            f.write(f"# Total matches: {result.total_count:,}\n")
            f.write(f"#\n")
            f.write(f"# PubMed Query (copy-paste to PubMed web):\n")
            f.write(f"# {'=' * 70}\n\n")
            f.write(result.query)
        saved_files['query'] = query_path

        # 3. Save PMIDs as simple list
        pmids_path = self.output_dir / f"pmids_{timestamp}_{safe_query}.txt"
        with open(pmids_path, 'w') as f:
            f.write('\n'.join(result.pmids))
        saved_files['pmids'] = pmids_path

        return saved_files


# =============================================================================
# PPI Search Agent
# =============================================================================

class PPISearchAgent:
    """
    AI-powered agent for searching PubMed for PPI literature.
    """

    def __init__(
        self,
        use_llm: bool = True,
        api_key: Optional[str] = None,
        output_dir: str = "output/searches"
    ):
        """
        Initialize the agent.

        Parameters
        ----------
        use_llm : bool
            Use LLM for term extraction (requires ANTHROPIC_API_KEY)
        api_key : str, optional
            Anthropic API key (or use ANTHROPIC_API_KEY env var)
        output_dir : str
            Directory to save search results
        """
        self.query_builder = PPIQueryBuilder()
        self.saver = SearchResultSaver(output_dir)

        # Initialize extractor
        if use_llm and ANTHROPIC_AVAILABLE:
            try:
                self.extractor = LLMExtractor(api_key)
                self.using_llm = True
            except ValueError:
                print("Warning: ANTHROPIC_API_KEY not found, using rule-based extraction")
                self.extractor = RuleBasedExtractor()
                self.using_llm = False
        else:
            self.extractor = RuleBasedExtractor()
            self.using_llm = False

    def search(self, user_query: str, retmax: int = 100) -> SearchResult:
        """
        Process user query and search PubMed.

        Parameters
        ----------
        user_query : str
            Natural language description of what user wants to find
        retmax : int
            Maximum number of results to retrieve

        Returns
        -------
        SearchResult
            Search results with query and extracted terms
        """
        # Extract terms
        terms = self.extractor.extract(user_query)

        # Build query
        pubmed_query = self.query_builder.build_query(terms)

        # Execute search
        result = search_pubmed(pubmed_query, retmax)

        return SearchResult(
            query=pubmed_query,
            total_count=result["count"],
            pmids=result["ids"],
            extracted_terms=terms
        )

    def interactive_session(self):
        """Run an interactive search session."""
        print("=" * 80)
        print("PPI Search Agent - Interactive Mode")
        print("=" * 80)
        print(f"\nExtraction method: {'LLM (Claude)' if self.using_llm else 'Rule-based'}")
        print("\nDescribe what you're looking for in natural language.")
        print("Examples:")
        print("  - 'BRCA1 interactions with DNA repair proteins in human cells'")
        print("  - 'p53-MDM2 binding studies using co-IP or crystallography'")
        print("  - 'kinase-substrate interactions in cancer from 2020-2024'")
        print("\nType 'quit' or 'exit' to end the session.\n")

        while True:
            try:
                user_input = input("\nWhat are you looking for? > ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\n\nGoodbye!")
                break

            if not user_input:
                continue

            if user_input.lower() in ['quit', 'exit', 'q']:
                print("\nGoodbye!")
                break

            # Process query
            print("\nProcessing your query...")
            result = self.search(user_input)

            # Display results
            print("\n" + "-" * 80)
            print("EXTRACTED TERMS:")
            print("-" * 80)
            terms = result.extracted_terms
            if terms.proteins:
                print(f"  Proteins: {', '.join(terms.proteins)}")
            if terms.interaction_types:
                print(f"  Interaction types: {', '.join(terms.interaction_types)}")
            if terms.methods:
                print(f"  Methods: {', '.join(terms.methods)}")
            if terms.organisms:
                print(f"  Organisms: {', '.join(terms.organisms)}")
            if terms.diseases:
                print(f"  Diseases: {', '.join(terms.diseases)}")
            if terms.keywords:
                print(f"  Keywords: {', '.join(terms.keywords)}")
            if terms.year_range:
                print(f"  Year range: {terms.year_range[0]}-{terms.year_range[1]}")

            print("\n" + "-" * 80)
            print("PUBMED QUERY:")
            print("-" * 80)
            print(result.query)

            print("\n" + "-" * 80)
            print("RESULTS:")
            print("-" * 80)
            print(f"Total matches: {result.total_count:,}")
            if result.pmids:
                print(f"First {len(result.pmids)} PMIDs: {result.pmids[:10]}")
                if len(result.pmids) > 10:
                    print(f"  ... and {len(result.pmids) - 10} more")

            print("\n" + "=" * 80)


# =============================================================================
# Main
# =============================================================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="AI-powered PubMed search agent for protein-protein interactions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python ppi_search_agent.py

  # Direct query
  python ppi_search_agent.py --query "BRCA1 interactions with DNA repair proteins"

  # Use rule-based extraction (no LLM)
  python ppi_search_agent.py --no-llm --query "p53 MDM2 binding"

  # Get more results
  python ppi_search_agent.py --query "kinase interactions" --retmax 500
        """
    )

    parser.add_argument(
        "--query", "-q",
        help="Direct query (skips interactive mode)"
    )

    parser.add_argument(
        "--no-llm",
        action="store_true",
        help="Use rule-based extraction instead of LLM"
    )

    parser.add_argument(
        "--retmax", "-n",
        type=int,
        default=100,
        help="Maximum number of results (default: 100)"
    )

    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON"
    )

    parser.add_argument(
        "--output", "-o",
        default="output/searches",
        help="Output directory for saved results (default: output/searches)"
    )

    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Don't save results to files"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Initialize agent
    agent = PPISearchAgent(use_llm=not args.no_llm, output_dir=args.output)

    if args.query:
        # Direct query mode
        extraction_method = "LLM (Claude)" if agent.using_llm else "Rule-based"

        if not args.json:
            print(f"\n[Extraction method: {extraction_method}]")

        result = agent.search(args.query, retmax=args.retmax)

        # Auto-save results (unless --no-save is specified)
        if not args.no_save:
            saved_files = agent.saver.save(args.query, result, extraction_method)
            if not args.json:
                print(f"\n[Results saved to:]")
                for file_type, path in saved_files.items():
                    print(f"  - {file_type}: {path}")

        if args.json:
            output = {
                'extraction_method': extraction_method,
                'query': result.query,
                'total_count': result.total_count,
                'pmids': result.pmids,
                'extracted_terms': result.extracted_terms.to_dict(),
                'saved_files': {k: str(v) for k, v in saved_files.items()} if not args.no_save else None
            }
            print(json.dumps(output, indent=2))
        else:
            print(f"\nExtracted terms: {result.extracted_terms.to_dict()}")
            print(f"\nPubMed query:\n{result.query}")
            print(f"\nTotal matches: {result.total_count:,}")
            print(f"First PMIDs: {result.pmids[:10]}")
    else:
        # Interactive mode
        agent.interactive_session()


if __name__ == "__main__":
    main()