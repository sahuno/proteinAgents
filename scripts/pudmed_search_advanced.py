from Bio import Entrez

Entrez.email = "ekwame001@gmail.com"
# Entrez.api_key = "YOUR_NCBI_API_KEY"  # uncomment and add your key for higher rate limits

def search_ppi(query, retmax=200):
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax,          # number of IDs to actually return now
        usehistory="y"          # store results on NCBI history server
    )
    record = Entrez.read(handle)
    handle.close()
    
    count = int(record["Count"])    # total matches
    id_list = record["IdList"]      # first up-to-retmax PMIDs
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    
    return {
        "count": count,
        "ids": id_list,
        "webenv": webenv,
        "query_key": query_key,
    }

# =============================================================================
# PubMed Advanced Query for Protein-Protein Interactions
# Search Strategy: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)
#
# O  = Base filter (English, NOT reviews/meta-analysis/retracted/clinical trials)
# A1 = Inclusion criteria - Concepts (protein complexes, interactions)
# A2 = Inclusion criteria - Methods (co-IP, affinity purification, pulldown, etc.)
# B  = Exclusion criteria (crosslinking, FRET, two-hybrid, etc.)
# =============================================================================

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

# =============================================================================
# Final Query Logic: (O AND (A1 OR A2)) NOT (A1 AND B NOT A2)
#
# Meaning: All research papers on protein complexes (O AND (A1 OR A2)),
# excluding those that contain only methodologies we are not interested in
# (papers in A1 AND B but NOT in A2)
# =============================================================================

full_query = f"(({O_block}) AND (({A1_block}) OR ({A2_block}))) NOT (({A1_block}) AND ({B_block}) NOT ({A2_block}))"

print("=" * 80)
print("FULL PUBMED QUERY")
print("=" * 80)
print(full_query)


res = search_ppi(full_query, retmax=50)
print("Total matches:", res["count"])
print("First 10 PMIDs:", res["ids"][:10])

print("First 10 PMIDs:", res["ids"][:10])
