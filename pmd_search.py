from Bio import Entrez

Entrez.email = "ekwame001@gmail.com"
Entrez.api_key = "YOUR_NCBI_API_KEY"  # comment out if you don’t have one

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

# Example usage
ppi_block = """
(
  "nucleoproteins"[MeSH Terms]
  OR "protein interaction mapping"[MeSH Terms]
  OR (
       "nucleoprotein"[All Fields]
       OR "nucleoproteins"[All Fields]
       OR "multiprotein"[All Fields]
       OR "multiproteins"[All Fields]
       OR "proteins"[MeSH Terms]
       OR "protein"[All Fields]
       OR "proteins"[All Fields]
       OR "enzyme"[All Fields]
     )
)
"""

filter_block = """
(
  "english"[Language]
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
  NOT "clinical trial, veterinary"[Publication Type]
)
"""

full_query = f"({ppi_block}) AND ({filter_block})"
print(full_query)


res = search_ppi(full_query, retmax=50)
print("Total matches:", res["count"])
print("First 10 PMIDs:", res["ids"][:10])

print("First 10 PMIDs:", res["ids"][:10])
