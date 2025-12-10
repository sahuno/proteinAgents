"""
PubMed PPI Query Module
=======================
A modular toolkit for searching PubMed for protein-protein interaction literature.

Search Strategy: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)

- O  = Base filter (English, NOT reviews/meta-analysis/retracted/clinical trials)
- A1 = Inclusion criteria - Concepts (protein complexes, interactions)
- A2 = Inclusion criteria - Methods (co-IP, affinity purification, pulldown, etc.)
- B  = Exclusion criteria (crosslinking, FRET, two-hybrid, etc.)

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

from .query_builder import PPIQueryBuilder
from .fetcher import PubMedFetcher
from .export import PubMedExporter

__all__ = ['PPIQueryBuilder', 'PubMedFetcher', 'PubMedExporter']
__version__ = '1.0.0'