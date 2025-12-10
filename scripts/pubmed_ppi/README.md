# PubMed PPI Search Agent

An AI-powered agent for searching PubMed for protein-protein interaction (PPI) literature using natural language queries.

**Author:** Samuel Ahuno (ekwame001@gmail.com)
**Date:** 2025-12-08

---

## Overview

The PPI Search Agent allows users to describe what they're looking for in natural language, then automatically:

1. **Extracts relevant search terms** (proteins, methods, organisms, etc.)
2. **Maps terms to PubMed syntax** (including MeSH terms)
3. **Builds a structured query** combining base filters with user terms
4. **Executes the search** and returns results

---

## Architecture

```
┌─────────────────┐
│   User Query    │  "BRCA1 interactions with DNA repair in human cells"
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Term Extractor  │  LLM (Claude) or Rule-based
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ ExtractedTerms  │
│  - proteins     │  ["BRCA1", "DNA"]
│  - interactions │  ["interaction"]
│  - methods      │  []
│  - organisms    │  ["human"]
│  - diseases     │  []
│  - year_range   │  None
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Query Builder  │  Maps to MeSH terms, builds PubMed syntax
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    O_BLOCK      │  Base filter (English, no reviews/trials)
│       +         │
│  User Terms     │  Combined with AND/OR logic
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  PubMed Search  │  Via Entrez API
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    Results      │  PMIDs, count, query used
└─────────────────┘
```

---

## Module Structure

```
scripts/pubmed_ppi/
├── __init__.py              # Package initialization
├── query_builder.py         # Query blocks (O, A1, A2, B) and construction
├── fetcher.py               # PubMed search and article fetch functions
├── export.py                # Save/export utilities (CSV, JSON, SQLite)
├── main.py                  # CLI for full pipeline
├── test_base_filter.py      # Test script with O_block + user keywords
├── ppi_search_agent.py      # AI-powered search agent
├── app.py                   # Streamlit web GUI
├── output/searches/         # Auto-saved search results
└── README.md                # This documentation
```

---

## Key Components

### 1. Base Filter (O_BLOCK)

The foundation query that filters for:
- English language papers
- Excludes: meta-analyses, reviews, retracted papers, errata, clinical trials

```python
O_BLOCK = (
    '("english"[Language] '
    'NOT "meta-analysis"[Publication Type] '
    'NOT "review"[Publication Type] '
    'NOT "retracted publication"[Publication Type] '
    # ... etc
)
```

### 2. ExtractedTerms Data Class

Container for parsed search terms:

```python
@dataclass
class ExtractedTerms:
    proteins: List[str]           # Gene/protein names
    interaction_types: List[str]  # binding, phosphorylation, etc.
    methods: List[str]            # co-IP, mass spec, cryo-EM
    organisms: List[str]          # human, mouse, yeast
    diseases: List[str]           # cancer, Alzheimer's
    keywords: List[str]           # Other terms
    year_range: Optional[tuple]   # (start_year, end_year)
```

### 3. Term Extractors

#### Rule-Based Extractor (No API required)

Simple pattern matching:
- **Proteins**: Capitalized words (BRCA1, TP53, MDM2)
- **Methods**: Keyword matching (co-ip, pulldown, cryo-em)
- **Organisms**: Keyword matching (human, mouse, yeast)
- **Year range**: Regex pattern `\d{4}-\d{4}`

```python
class RuleBasedExtractor:
    INTERACTION_KEYWORDS = ['interaction', 'binding', 'complex', ...]
    METHOD_KEYWORDS = ['co-ip', 'mass spectrometry', 'cryo-em', ...]
    ORGANISM_KEYWORDS = ['human', 'mouse', 'yeast', ...]
```

#### LLM Extractor (Requires ANTHROPIC_API_KEY)

Uses Claude to understand natural language:

```python
class LLMExtractor:
    SYSTEM_PROMPT = """You are a scientific literature search assistant...
    Extract: proteins, interaction_types, methods, organisms, diseases, keywords, year_range
    Respond ONLY with JSON."""
```

### 4. Query Builder

Maps extracted terms to PubMed syntax:

```python
class PPIQueryBuilder:
    # Maps common method names to MeSH terms
    METHOD_MESH_MAP = {
        'co-ip': '"Immunoprecipitation"[MeSH Terms]',
        'cryo-em': '"Cryoelectron Microscopy"[MeSH Terms]',
        ...
    }

    # Maps organisms to MeSH terms
    ORGANISM_MESH_MAP = {
        'human': '"Humans"[MeSH Terms]',
        'mouse': '"Mice"[MeSH Terms]',
        ...
    }
```

---

## Usage

### Interactive Mode

```bash
python ppi_search_agent.py
```

Prompts user for natural language queries in a loop.

### Direct Query Mode

```bash
# With LLM extraction (requires ANTHROPIC_API_KEY)
python ppi_search_agent.py --query "BRCA1 interactions with DNA repair proteins"

# With rule-based extraction (no API needed)
python ppi_search_agent.py --no-llm --query "p53 MDM2 binding"
```

### Output Options

```bash
# JSON output for programmatic use
python ppi_search_agent.py --query "..." --json

# Control number of results
python ppi_search_agent.py --query "..." --retmax 500
```

---

## Examples

### Example 1: Simple protein search

```bash
python ppi_search_agent.py --no-llm --query "BRCA1 protein interaction"
```

**Extracted:**
- proteins: ["BRCA1"]
- interaction_types: ["interaction"]

**Results:** ~258 papers

### Example 2: Protein-protein interaction with method

```bash
python ppi_search_agent.py --no-llm --query "p53 MDM2 binding using co-IP"
```

**Extracted:**
- proteins: ["MDM2"] (p53 lowercase, not detected)
- interaction_types: ["binding"]
- methods: ["co-ip"] → mapped to `"Immunoprecipitation"[MeSH Terms]`

### Example 3: Full query with organism and year range

```bash
python ppi_search_agent.py --no-llm \
  --query "p53 MDM2 binding using co-IP in mouse from 2020-2024" \
  --json
```

**Output:**
```json
{
  "extracted_terms": {
    "proteins": ["MDM2"],
    "interaction_types": ["binding"],
    "methods": ["co-ip"],
    "organisms": ["mouse"],
    "year_range": ["2020", "2024"]
  },
  "total_count": 2,
  "pmids": ["33675124", "32820269"]
}
```

---

## How to Extend

### Add new method mappings

Edit `PPIQueryBuilder.METHOD_MESH_MAP`:

```python
METHOD_MESH_MAP = {
    'biolayer interferometry': '"Interferometry"[MeSH Terms]',
    'crosslinking': '"Cross-Linking Reagents"[MeSH Terms]',
    # Add more...
}
```

### Add new organism mappings

Edit `PPIQueryBuilder.ORGANISM_MESH_MAP`:

```python
ORGANISM_MESH_MAP = {
    'arabidopsis': '"Arabidopsis"[MeSH Terms]',
    'xenopus': '"Xenopus"[MeSH Terms]',
    # Add more...
}
```

### Improve rule-based extraction

Edit `RuleBasedExtractor` keyword lists:

```python
INTERACTION_KEYWORDS = [
    'interaction', 'binding', 'complex',
    'dimerization', 'oligomerization',  # Add more
]
```

---

## Dependencies

- **Required:** `biopython` (for Entrez API)
- **Optional:** `anthropic` (for LLM extraction)
- **Optional:** `streamlit`, `pandas` (for GUI)

```bash
pip install biopython
pip install anthropic  # Optional, for LLM mode
pip install streamlit pandas  # Optional, for GUI
```

---

## Web GUI (Streamlit)

The agent includes an optional web-based graphical interface built with Streamlit.

### Running the GUI

```bash
# Install dependencies
pip install streamlit pandas

# Run the app
streamlit run app.py
```

This will open a browser window with the GUI at `http://localhost:8501`.

### GUI Features

- **Natural language search input** with example queries
- **Settings sidebar** for LLM mode and result count
- **Extracted terms display** showing proteins, methods, organisms detected
- **Full PubMed query** ready to copy-paste
- **Clickable PMID links** to PubMed articles
- **Export buttons** for JSON, TXT (query), TXT (PMIDs), and CSV formats
- **Search history** to reload previous searches

### Screenshot

```
+---------------------------+--------------------------------+
| Settings                  |  PPI Search Agent              |
|  [x] Use LLM              |  [Search box..................] |
|  Max results: 100         |  [Search button]               |
|                           |                                |
| Search History            |  Extracted Terms | Query | ... |
|  - BRCA1 interactions     |  Proteins: BRCA1               |
|  - p53 MDM2 binding       |  Methods: co-ip                |
+---------------------------+--------------------------------+
```

---

## Environment Variables

| Variable | Description | Required |
|----------|-------------|----------|
| `ANTHROPIC_API_KEY` | Claude API key for LLM extraction | Only if using LLM mode |

---

## Related Files

| File | Purpose |
|------|---------|
| `test_base_filter.py` | Simple test with O_block + keywords |
| `main.py` | Full pipeline (search, fetch, export) |
| `query_builder.py` | Complete query blocks (O, A1, A2, B) |
| `fetcher.py` | Fetch article metadata (title, abstract, authors) |
| `export.py` | Export to CSV, JSON, SQLite |

---

## Search Strategy Reference

The full PPI search strategy from the slides:

```
(O AND (A1 OR A2)) NOT (A1 AND B NOT A2)
```

| Block | Description |
|-------|-------------|
| **O** | Base filter (English, no reviews/trials) |
| **A1** | Inclusion concepts (protein complexes, interactions) |
| **A2** | Inclusion methods (co-IP, affinity purification, etc.) |
| **B** | Exclusion methods (crosslinking, two-hybrid, FRET, etc.) |

The agent uses **O_block** as the foundation and adds user-specified terms on top.