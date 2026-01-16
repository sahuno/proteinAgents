# O_Block (Block 0) Reference Guide

**Date:** 2025-01-16
**Author:** Samuel Ahuno (ekwame001@gmail.com)

---

## What is the O_Block?

The **O_block** (also called "Block 0") is the **base filter** in PubMed queries that excludes unwanted publication types and ensures English-language results. It acts as a quality control filter to remove junk from search results.

---

## O_Block Definition

```python
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
```

---

## What Does O_Block Exclude?

### ✅ **Includes:**
- English-language articles
- Original research articles
- Case reports
- Letters
- Brief communications

### ❌ **Excludes:**
1. **Reviews & Meta-analyses** - Secondary sources
2. **Clinical Trials** - All phases and types (human studies)
3. **Retracted Publications** - Invalid research
4. **Published Errata** - Corrections to other articles
5. **Non-English Articles** - Language barrier

---

## Why Do We Need O_Block?

### Problem Without O_Block
If you search for:
```
("co-immunoprecipitation"[Title/Abstract] OR "affinity purification"[Title/Abstract])
```

You get:
- ✅ Original research articles
- ❌ Review articles summarizing other work
- ❌ Clinical trial protocols (not relevant for basic research)
- ❌ Meta-analyses (statistical summaries)
- ❌ Retracted/corrected articles (invalid data)
- ❌ Non-English articles

**Result:** ~38,887 articles (includes junk)

### Solution With O_Block
With O_block:
```
(O_block AND (A2_block))
```

You get:
- ✅ Original research articles only
- ✅ English-language
- ✅ Valid (not retracted)
- ✅ Primary sources

**Result:** ~23,360 articles (higher quality, focused results)

---

## Query Structure

### Complete PubMed Query Format

```
(O_block AND (A2_block))
```

Where:
- **O_block** = Base filters (language + publication type exclusions)
- **A2_block** = Methods/techniques (e.g., co-IP, Y2H, etc.)

### Example Final Query

```
(("english"[Language]
  NOT "meta-analysis"[Publication Type]
  NOT "review"[Publication Type]
  NOT "retracted publication"[Publication Type]
  NOT "retraction of publication"[Publication Type]
  NOT "published erratum"[Publication Type]
  NOT "controlled clinical trial"[Publication Type]
  NOT "clinical study"[Publication Type]
  NOT "clinical trial"[Publication Type]
  NOT "clinical trial protocol"[Publication Type]
  NOT "clinical trial, phase i"[Publication Type]
  NOT "clinical trial, phase ii"[Publication Type]
  NOT "clinical trial, phase iii"[Publication Type]
  NOT "clinical trial, phase iv"[Publication Type]
  NOT "clinical trial, veterinary"[Publication Type])
AND
  (("co-immunoprecipitation"[Title/Abstract]
    OR "affinity purification"[Title/Abstract]
    OR "yeast two-hybrid"[Title/Abstract]
    OR "pull-down assay"[Title/Abstract])))
```

---

## Impact on Results

### Test Case: co-IP + affinity purification

| Query | Result Count | Difference |
|-------|--------------|------------|
| A2_block only | 38,887 | Baseline |
| O_block AND A2_block | 23,360 | -40% (filtered out junk) |

The O_block removes approximately **40% of results** by filtering out reviews, clinical trials, and other non-primary research articles.

---

## Implementation in HITL_PROPER

### Location
- **File:** `blockA_single_Agent_HITL_PROPER.py`
- **Line:** 55-71 (constant definition)
- **Line:** 287-292 (query combination)

### How It Works

**Step 1: Define O_BLOCK constant**
```python
# blockA_single_Agent_HITL_PROPER.py:55-71
O_BLOCK = (
    '("english"[Language] '
    'NOT "meta-analysis"[Publication Type] '
    # ... rest of filters
)
```

**Step 2: Combine with A2_block in query creation**
```python
# blockA_single_Agent_HITL_PROPER.py:287-292
def create_final_query_node(state: QueryCreationInput) -> QueryCreationInput:
    a2_block = state["a2_block"]

    # Combine O_block (base filters) AND A2_block (methods)
    if a2_block:
        final_query = f"({O_BLOCK} AND ({a2_block}))"
    else:
        final_query = O_BLOCK

    return {...}
```

**Step 3: Display to user**
```python
# blockA_single_Agent_HITL_PROPER.py:439-453
def display_query(query: str, components: Dict[str, str]):
    print("\nO Block (Base Filters - excludes reviews, clinical trials, non-English):")
    print(f"  {components['o_block']}")

    print("\nA2 Block (Methods):")
    print(f"  {components['a2_block']}")

    print("\nFinal Query (O_block AND A2_block):")
    print(f"  {query}")
```

---

## Customizing O_Block

### Option 1: Modify the Constant
Edit the `O_BLOCK` constant at the top of the file to add/remove filters:

```python
# Add organism filter
O_BLOCK = (
    '("english"[Language] '
    'AND "Homo sapiens"[Organism] '  # ADD THIS
    'NOT "meta-analysis"[Publication Type] '
    # ... rest
)
```

### Option 2: Make it User-Configurable (Future Enhancement)
Add an HITL checkpoint to let users customize O_block:

```python
def run_hitl_workflow():
    # ... existing code ...

    # HITL Checkpoint: O_block customization
    print("\nDefault O_block filters:")
    print("- English language")
    print("- Excludes reviews, clinical trials, retracted")

    customize = input("Customize O_block? (yes/no): ")

    if customize == "yes":
        # Let user add organism, date range, etc.
        organism = input("Add organism filter (e.g., 'Homo sapiens'): ")
        # ... build custom O_block
```

---

## Common Customizations

### 1. Add Organism Filter
```python
O_BLOCK = (
    '("english"[Language] '
    'AND "Homo sapiens"[Organism] '  # Human only
    'NOT "meta-analysis"[Publication Type] '
    # ...
)
```

### 2. Add Date Range
```python
O_BLOCK = (
    '("english"[Language] '
    'AND "2020/01/01"[Date - Publication] : "2024/12/31"[Date - Publication] '
    'NOT "meta-analysis"[Publication Type] '
    # ...
)
```

### 3. Add Journal Quality Filter
```python
O_BLOCK = (
    '("english"[Language] '
    'AND ("Nature"[Journal] OR "Science"[Journal] OR "Cell"[Journal]) '
    'NOT "meta-analysis"[Publication Type] '
    # ...
)
```

---

## FAQ

**Q: Why exclude reviews?**
A: Reviews summarize existing research and don't report original experiments. For finding primary protein interaction studies, we want original research articles.

**Q: Why exclude clinical trials?**
A: Clinical trials focus on human therapeutic interventions, which are typically not relevant for basic PPI research using biochemical methods.

**Q: Can I include reviews?**
A: Yes! Remove the `'NOT "review"[Publication Type] '` line from O_BLOCK if you want to include review articles.

**Q: Why English only?**
A: Most researchers work in English, and it ensures you can read the papers. Remove `'("english"[Language] '` if you want all languages.

**Q: Does O_block affect result quality?**
A: Yes - it improves quality by filtering out secondary sources and invalid research, giving you more focused, actionable results.

---

## References

- **Original Implementation:** `blockA_single_Agent_HITL.py:163-181`
- **Proper Implementation:** `blockA_single_Agent_HITL_PROPER.py:55-71, 287-292`
- **PubMed Search Tags:** https://pubmed.ncbi.nlm.nih.gov/help/#search-tags

---

**Author:** Samuel Ahuno (ekwame001@gmail.com)
**Date:** 2025-01-16
**Version:** 1.0
