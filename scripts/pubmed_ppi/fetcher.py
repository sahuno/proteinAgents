"""
PubMed Fetcher for Protein-Protein Interactions
================================================
Search PubMed and fetch article metadata using the Entrez API.

Author: Samuel Ahuno (ekwame001@gmail.com)
Date: 2025-12-08
"""

import time
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from Bio import Entrez
from xml.etree import ElementTree as ET


@dataclass
class SearchResult:
    """Container for PubMed search results."""
    query: str
    total_count: int
    pmids: List[str]
    webenv: str
    query_key: str
    retmax: int
    search_time: str = ""


@dataclass
class Article:
    """Container for article metadata."""
    pmid: str
    title: str = ""
    abstract: str = ""
    authors: List[str] = field(default_factory=list)
    journal: str = ""
    year: str = ""
    month: str = ""
    doi: str = ""
    mesh_terms: List[str] = field(default_factory=list)
    keywords: List[str] = field(default_factory=list)
    publication_types: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert article to dictionary."""
        return {
            'pmid': self.pmid,
            'title': self.title,
            'abstract': self.abstract,
            'authors': '; '.join(self.authors),
            'journal': self.journal,
            'year': self.year,
            'month': self.month,
            'doi': self.doi,
            'mesh_terms': '; '.join(self.mesh_terms),
            'keywords': '; '.join(self.keywords),
            'publication_types': '; '.join(self.publication_types)
        }


class PubMedFetcher:
    """
    Fetch articles from PubMed using the NCBI Entrez API.

    Parameters
    ----------
    email : str
        Email address for NCBI API (required by NCBI)
    api_key : str, optional
        NCBI API key for higher rate limits (10 requests/sec vs 3/sec)
    """

    def __init__(self, email: str, api_key: Optional[str] = None):
        """Initialize the fetcher with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        # Rate limiting: 3 requests/sec without key, 10 with key
        self._request_delay = 0.1 if api_key else 0.34

    def search(self, query: str, retmax: int = 10000) -> SearchResult:
        """
        Search PubMed with a query string.

        Parameters
        ----------
        query : str
            PubMed query string
        retmax : int
            Maximum number of PMIDs to retrieve (default: 10000)

        Returns
        -------
        SearchResult
            Container with search results and history server info
        """
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=retmax,
            usehistory="y"
        )
        record = Entrez.read(handle)
        handle.close()

        return SearchResult(
            query=query,
            total_count=int(record["Count"]),
            pmids=record["IdList"],
            webenv=record["WebEnv"],
            query_key=record["QueryKey"],
            retmax=retmax
        )

    def fetch_articles(
        self,
        search_result: Optional[SearchResult] = None,
        pmids: Optional[List[str]] = None,
        batch_size: int = 100,
        max_articles: Optional[int] = None,
        verbose: bool = True
    ) -> List[Article]:
        """
        Fetch article metadata for PMIDs.

        Parameters
        ----------
        search_result : SearchResult, optional
            Result from search() to use history server (more efficient)
        pmids : list of str, optional
            List of PMIDs to fetch (if not using search_result)
        batch_size : int
            Number of articles to fetch per request (default: 100)
        max_articles : int, optional
            Maximum number of articles to fetch (default: all)
        verbose : bool
            Print progress messages

        Returns
        -------
        list of Article
            List of Article objects with metadata
        """
        if search_result is None and pmids is None:
            raise ValueError("Must provide either search_result or pmids")

        articles = []

        if search_result is not None:
            # Use history server for efficient batch fetching
            total = min(search_result.total_count, max_articles or search_result.total_count)
            total = min(total, len(search_result.pmids)) if search_result.pmids else total

            if verbose:
                print(f"Fetching {total} articles using history server...")

            for start in range(0, total, batch_size):
                end = min(start + batch_size, total)
                if verbose:
                    print(f"  Fetching articles {start + 1} to {end}...")

                try:
                    handle = Entrez.efetch(
                        db="pubmed",
                        rettype="xml",
                        retmode="xml",
                        retstart=start,
                        retmax=batch_size,
                        webenv=search_result.webenv,
                        query_key=search_result.query_key
                    )
                    batch_articles = self._parse_xml(handle.read())
                    articles.extend(batch_articles)
                    handle.close()
                except Exception as e:
                    print(f"  Error fetching batch {start}-{end}: {e}")

                time.sleep(self._request_delay)

        else:
            # Fetch by PMID list
            pmids_to_fetch = pmids[:max_articles] if max_articles else pmids

            if verbose:
                print(f"Fetching {len(pmids_to_fetch)} articles by PMID...")

            for i in range(0, len(pmids_to_fetch), batch_size):
                batch = pmids_to_fetch[i:i + batch_size]
                if verbose:
                    print(f"  Fetching batch {i // batch_size + 1}...")

                try:
                    handle = Entrez.efetch(
                        db="pubmed",
                        id=",".join(batch),
                        rettype="xml",
                        retmode="xml"
                    )
                    batch_articles = self._parse_xml(handle.read())
                    articles.extend(batch_articles)
                    handle.close()
                except Exception as e:
                    print(f"  Error fetching batch: {e}")

                time.sleep(self._request_delay)

        if verbose:
            print(f"Successfully fetched {len(articles)} articles")

        return articles

    def _parse_xml(self, xml_data: bytes) -> List[Article]:
        """
        Parse PubMed XML response into Article objects.

        Parameters
        ----------
        xml_data : bytes
            Raw XML response from Entrez.efetch

        Returns
        -------
        list of Article
            Parsed article metadata
        """
        articles = []

        try:
            root = ET.fromstring(xml_data)

            for article_elem in root.findall('.//PubmedArticle'):
                article = self._parse_article_element(article_elem)
                if article:
                    articles.append(article)

        except ET.ParseError as e:
            print(f"XML parse error: {e}")

        return articles

    def _parse_article_element(self, elem: ET.Element) -> Optional[Article]:
        """Parse a single PubmedArticle XML element."""
        try:
            # Get PMID
            pmid_elem = elem.find('.//PMID')
            if pmid_elem is None:
                return None
            pmid = pmid_elem.text

            # Get title
            title_elem = elem.find('.//ArticleTitle')
            title = title_elem.text if title_elem is not None and title_elem.text else ""

            # Get abstract
            abstract_parts = []
            for abs_elem in elem.findall('.//AbstractText'):
                label = abs_elem.get('Label', '')
                text = ''.join(abs_elem.itertext()) if abs_elem.text or len(abs_elem) > 0 else ''
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
            abstract = ' '.join(abstract_parts)

            # Get authors
            authors = []
            for author in elem.findall('.//Author'):
                lastname = author.find('LastName')
                forename = author.find('ForeName')
                if lastname is not None:
                    name = lastname.text or ""
                    if forename is not None and forename.text:
                        name = f"{name} {forename.text}"
                    authors.append(name)

            # Get journal
            journal_elem = elem.find('.//Journal/Title')
            journal = journal_elem.text if journal_elem is not None else ""

            # Get publication date
            year = ""
            month = ""
            pub_date = elem.find('.//PubDate')
            if pub_date is not None:
                year_elem = pub_date.find('Year')
                month_elem = pub_date.find('Month')
                year = year_elem.text if year_elem is not None else ""
                month = month_elem.text if month_elem is not None else ""

            # Get DOI
            doi = ""
            for article_id in elem.findall('.//ArticleId'):
                if article_id.get('IdType') == 'doi':
                    doi = article_id.text or ""
                    break

            # Get MeSH terms
            mesh_terms = []
            for mesh in elem.findall('.//MeshHeading/DescriptorName'):
                if mesh.text:
                    mesh_terms.append(mesh.text)

            # Get keywords
            keywords = []
            for kw in elem.findall('.//Keyword'):
                if kw.text:
                    keywords.append(kw.text)

            # Get publication types
            pub_types = []
            for pt in elem.findall('.//PublicationType'):
                if pt.text:
                    pub_types.append(pt.text)

            return Article(
                pmid=pmid,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=journal,
                year=year,
                month=month,
                doi=doi,
                mesh_terms=mesh_terms,
                keywords=keywords,
                publication_types=pub_types
            )

        except Exception as e:
            print(f"Error parsing article: {e}")
            return None

    def get_article_count(self, query: str) -> int:
        """
        Get the total count of articles matching a query without fetching PMIDs.

        Parameters
        ----------
        query : str
            PubMed query string

        Returns
        -------
        int
            Total number of matching articles
        """
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=0  # Only get count, no IDs
        )
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])


if __name__ == "__main__":
    # Example usage
    from query_builder import PPIQueryBuilder

    # Build query
    builder = PPIQueryBuilder()
    query = builder.build_query()

    # Initialize fetcher (replace with your email)
    fetcher = PubMedFetcher(email="ekwame001@gmail.com")

    # Get count only
    count = fetcher.get_article_count(query)
    print(f"Total articles matching query: {count:,}")

    # Search and fetch a small sample
    print("\nSearching for first 10 articles...")
    result = fetcher.search(query, retmax=10)
    print(f"Found {result.total_count:,} total matches")
    print(f"Retrieved {len(result.pmids)} PMIDs")

    # Fetch article details
    articles = fetcher.fetch_articles(search_result=result, max_articles=5)
    for article in articles:
        print(f"\n[{article.pmid}] {article.title[:80]}...")
        print(f"  Authors: {', '.join(article.authors[:3])}...")
        print(f"  Journal: {article.journal}, {article.year}")