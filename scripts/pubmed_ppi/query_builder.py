"""
PubMed Query Builder for Protein-Protein Interactions
======================================================
Constructs PubMed advanced queries using a modular block-based approach.

Search Strategy: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

from typing import Dict, Optional


class PPIQueryBuilder:
    """
    Build PubMed queries for protein-protein interaction literature search.

    The query follows the logic: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)

    Where:
        O  = Base filter (language, publication types to exclude)
        A1 = Inclusion concepts (protein complexes, interactions)
        A2 = Inclusion methods (experimental techniques we want)
        B  = Exclusion methods (techniques we don't want)
    """

    def __init__(self):
        """Initialize query blocks with default PPI search terms."""
        self._blocks = self._initialize_default_blocks()

    def _initialize_default_blocks(self) -> Dict[str, str]:
        """Initialize the default query blocks based on the PPI search strategy."""

        # O: BASE FILTER - English language, exclude unwanted publication types
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

        # A1: INCLUSION CRITERIA - CONCEPTS (protein complexes, interactions)
        A1_block = (
            '("nucleoproteins"[MeSH Terms] '
            'OR "protein interaction mapping"[MeSH Terms] '
            'OR (("nucleoprotein"[All Fields] OR "nucleoproteins"[All Fields] '
            'OR "multiprotein"[All Fields] OR "multiproteins"[All Fields] '
            'OR "proteins"[MeSH Terms] OR "protein"[All Fields] '
            'OR "proteins"[All Fields] OR "enzyme"[All Fields]) '
            'AND ("interact"[All Fields] OR "interacted"[All Fields] '
            'OR "interacting"[All Fields] OR "interaction"[All Fields] '
            'OR "interactions"[All Fields] OR "interactivity"[All Fields] '
            'OR "interacts"[All Fields])) '
            'OR "protein interaction"[All Fields] '
            'OR "protein interactions"[All Fields] '
            'OR "interacting protein"[All Fields] '
            'OR "interacting proteins"[All Fields] '
            'OR "multiprotein complexes"[MeSH Terms] '
            'OR (("nucleoprotein"[All Fields] OR "nucleoproteins"[All Fields] '
            'OR "multiprotein"[All Fields] OR "multiproteins"[All Fields] '
            'OR "proteins"[MeSH Terms] OR "protein"[All Fields] '
            'OR "proteins"[All Fields] OR "enzyme"[All Fields]) '
            'AND ("complex"[All Fields] OR "complexes"[All Fields] '
            'OR "heteromer"[All Fields] OR "heteromers"[All Fields] '
            'OR "homomer"[All Fields] OR "homomers"[All Fields] '
            'OR "heteromeric"[All Fields] OR "homomeric"[All Fields] '
            'OR "subunit"[All Fields] OR "subunits"[All Fields])) '
            'OR "protein complex"[All Fields] '
            'OR "protein complexes"[All Fields] '
            'OR ("protein"[All Fields] '
            'AND ("RNA"[All Fields] OR "DNA"[All Fields] '
            'OR "ribonucleic"[All Fields] OR "deoxyribonucleic"[All Fields]) '
            'AND ("interaction"[All Fields] OR "interactions"[All Fields])))'
        )

        # A2: INCLUSION CRITERIA - METHODS (experimental techniques we want)
        A2_block = (
            '("Immunoprecipitation"[MeSH Terms] '
            'OR "coimmunoprecipitation"[All Fields] '
            'OR ("co"[All Fields] AND "immunoprecipitation"[All Fields]) '
            'OR ("RNA"[All Fields] AND "immunoprecipitation"[All Fields]) '
            'OR "co immunoprecipitation"[All Fields] '
            'OR "coIP"[All Fields] '
            'OR "Chromatography, Affinity"[MeSH Terms] '
            'OR "affinity purification"[All Fields] '
            'OR "affinity isolation"[All Fields] '
            'OR "affinity chromatography"[All Fields] '
            'OR ("affinity"[All Fields] AND ("purification"[All Fields] '
            'OR "isolation"[All Fields] OR "chromatography"[All Fields])) '
            'OR "pulldown"[All Fields] '
            'OR ("pull"[All Fields] AND ("down"[All Fields] OR "downs"[All Fields])) '
            'OR "crystallography, x ray"[MeSH Terms] '
            'OR "Nuclear Magnetic Resonance, Biomolecular"[MeSH Terms] '
            'OR "Cryoelectron Microscopy"[MeSH Terms] '
            'OR "Protein Array Analysis"[MeSH Terms] '
            'OR "electrophoretic mobility shift assay"[MeSH Terms] '
            'OR "Surface Plasmon Resonance"[MeSH Terms])'
        )

        # B: EXCLUSION CRITERIA - Methods we do NOT want
        B_block = (
            '("Epitope Mapping"[MeSH Terms] '
            'OR "Precipitin Tests"[MeSH Terms] '
            'OR "Two-Hybrid System Techniques"[MeSH Terms] '
            'OR "blotting, far western"[MeSH Terms] '
            'OR "Radioimmunoprecipitation Assay"[MeSH Terms] '
            'OR "Autoantibodies"[MeSH Terms] '
            'OR "Chromatin Immunoprecipitation"[MeSH Terms] '
            'OR "Cross-Linking Reagents"[MeSH Terms] '
            'OR "Formaldehyde"[MeSH Terms] '
            'OR "Microscopy"[MeSH Terms] '
            'OR ("crosslink"[All Fields] AND "reagent"[All Fields]) '
            'OR "Formaldehyde"[All Fields] '
            'OR "DNA Footprinting"[MeSH Terms] '
            'OR "Nuclease Protection Assays"[MeSH Terms] '
            'OR "Blotting, Southwestern"[MeSH Terms] '
            'OR "Fluorescence Resonance Energy Transfer"[MeSH Terms] '
            'OR "PAR CLIP"[All Fields] '
            'OR "AlphaFold"[All Fields])'
        )

        return {
            'O': O_block,
            'A1': A1_block,
            'A2': A2_block,
            'B': B_block
        }

    @property
    def O_block(self) -> str:
        """Get the base filter block (O)."""
        return self._blocks['O']

    @property
    def A1_block(self) -> str:
        """Get the inclusion concepts block (A1)."""
        return self._blocks['A1']

    @property
    def A2_block(self) -> str:
        """Get the inclusion methods block (A2)."""
        return self._blocks['A2']

    @property
    def B_block(self) -> str:
        """Get the exclusion methods block (B)."""
        return self._blocks['B']

    def set_block(self, block_name: str, query: str) -> None:
        """
        Override a specific query block.

        Parameters
        ----------
        block_name : str
            One of 'O', 'A1', 'A2', or 'B'
        query : str
            The new query string for this block
        """
        if block_name not in self._blocks:
            raise ValueError(f"Invalid block name: {block_name}. Must be one of: O, A1, A2, B")
        self._blocks[block_name] = query

    def build_query(self) -> str:
        """
        Build the full PubMed query using the logic:
        (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)

        Returns
        -------
        str
            The complete PubMed query string
        """
        O = self._blocks['O']
        A1 = self._blocks['A1']
        A2 = self._blocks['A2']
        B = self._blocks['B']

        full_query = f"(({O}) AND (({A1}) OR ({A2}))) NOT (({A1}) AND ({B}) NOT ({A2}))"
        return full_query

    def build_simple_query(self, include_exclusions: bool = True) -> str:
        """
        Build a simplified query without the complex exclusion logic.

        Parameters
        ----------
        include_exclusions : bool
            If True, use simple NOT B logic. If False, only use inclusion criteria.

        Returns
        -------
        str
            Simplified PubMed query string
        """
        O = self._blocks['O']
        A1 = self._blocks['A1']
        A2 = self._blocks['A2']
        B = self._blocks['B']

        if include_exclusions:
            return f"(({O}) AND (({A1}) OR ({A2}))) NOT ({B})"
        else:
            return f"(({O}) AND (({A1}) OR ({A2})))"

    def get_block_descriptions(self) -> Dict[str, str]:
        """
        Get human-readable descriptions of each query block.

        Returns
        -------
        dict
            Block names mapped to their descriptions
        """
        return {
            'O': 'Base filter: English language, excludes reviews, meta-analyses, '
                 'retracted papers, errata, and clinical trials',
            'A1': 'Inclusion concepts: Protein complexes, nucleoproteins, '
                  'protein interactions, heteromers/homomers, protein-DNA/RNA interactions',
            'A2': 'Inclusion methods: co-IP, immunoprecipitation, affinity purification, '
                  'pulldown, X-ray crystallography, NMR, Cryo-EM, SPR, EMSA, protein arrays',
            'B': 'Exclusion methods: Crosslinking, formaldehyde, two-hybrid, FRET, '
                 'epitope mapping, far western, nuclease protection, ChIP, AlphaFold'
        }

    def get_query_logic_explanation(self) -> str:
        """
        Get an explanation of the query logic.

        Returns
        -------
        str
            Human-readable explanation of the search strategy
        """
        return """
PubMed Query Logic: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)

This query retrieves:
- All research papers matching base filters (O)
- AND either about protein interaction concepts (A1) OR using desired methods (A2)
- EXCLUDING papers that are about concepts (A1) with unwanted methods (B)
  but don't use any of the wanted methods (A2)

In plain English:
"Find all English research papers about protein interactions or using
interaction detection methods, but exclude papers that only use methods
we're not interested in (like crosslinking or two-hybrid assays)."
"""


if __name__ == "__main__":
    # Example usage
    builder = PPIQueryBuilder()

    print("=" * 80)
    print("QUERY BLOCK DESCRIPTIONS")
    print("=" * 80)
    for block, desc in builder.get_block_descriptions().items():
        print(f"\n{block}: {desc}")

    print("\n" + "=" * 80)
    print("QUERY LOGIC")
    print("=" * 80)
    print(builder.get_query_logic_explanation())

    print("=" * 80)
    print("FULL QUERY")
    print("=" * 80)
    print(builder.build_query())