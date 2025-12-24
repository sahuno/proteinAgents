#!/usr/bin/env python3
"""
Parse NLM MeSH desc2025.xml -> per-descriptor text -> embeddings.

Outputs:
- mesh_descriptors.jsonl : one record per DescriptorUI
- mesh_embeddings.npy    : embeddings matrix (N x D)
- mesh_faiss.index       : optional FAISS index

Example:
  python /Users/ahunos/myWork/proteinAgents/notebooks/scripts/mesh_desc_to_embeddings.py \
    --desc_xml /Users/ahunos/myWork/proteinAgents/data/mesh_2025/desc2025 \
    --out_dir mesh_out \
    --backend openai \
    --openai_model text-embedding-3-small \
    --batch_size 128 \
    --max_terms 0

Notes:
- desc2025.xml is large; we stream-parse using iterparse.
- OpenAI backend requires env var OPENAI_API_KEY.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import xml.etree.ElementTree as ET


# -----------------------------
# Data model
# -----------------------------
@dataclass
class MeshDescriptor:
    descriptor_ui: str
    descriptor_name: str
    scope_note: str
    entry_terms: List[str]
    tree_numbers: List[str]

    def to_embedding_text(
        self,
        include_scope_note: bool = True,
        include_entry_terms: bool = True,
        include_tree_numbers: bool = True,
        max_entry_terms: int = 30,
    ) -> str:
        parts: List[str] = []
        parts.append(f"MeSH Descriptor: {self.descriptor_name}")
        parts.append(f"DescriptorUI: {self.descriptor_ui}")

        if include_scope_note and self.scope_note:
            # Keep scope note relatively prominent; it tends to be highly informative.
            parts.append(f"Definition: {self.scope_note.strip()}")

        if include_entry_terms and self.entry_terms:
            terms = self.entry_terms[:max_entry_terms]
            parts.append("Synonyms/Entry terms: " + "; ".join(terms))

        if include_tree_numbers and self.tree_numbers:
            parts.append("Tree numbers: " + "; ".join(self.tree_numbers))

        return "\n".join(parts).strip()


# -----------------------------
# XML parsing (streaming)
# -----------------------------
def _text_or_empty(elem: Optional[ET.Element]) -> str:
    return (elem.text or "").strip() if elem is not None else ""


def parse_desc_xml(desc_xml_path: str, max_terms: int = 0) -> Iterable[MeshDescriptor]:
    """
    Stream-parse descXXXX.xml and yield MeshDescriptor records.
    max_terms=0 means no limit.
    """
    # Events: start/end; we only act on end of DescriptorRecord.
    context = ET.iterparse(desc_xml_path, events=("end",))
    count = 0

    for event, elem in context:
        if elem.tag != "DescriptorRecord":
            continue

        # DescriptorUI
        descriptor_ui = _text_or_empty(elem.find("DescriptorUI"))

        # DescriptorName/String
        dn = elem.find("DescriptorName/String")
        descriptor_name = _text_or_empty(dn)

        # ScopeNote (can be absent)
        scope_note = _text_or_empty(elem.find("ConceptList/Concept/ScopeNote"))

        # Entry terms: ConceptList/Concept/TermList/Term/String
        entry_terms: List[str] = []
        for term_str in elem.findall("ConceptList/Concept/TermList/Term/String"):
            s = _text_or_empty(term_str)
            if s:
                entry_terms.append(s)

        # Deduplicate while preserving order
        seen = set()
        entry_terms = [t for t in entry_terms if not (t in seen or seen.add(t))]

        # Tree numbers
        tree_numbers = [_text_or_empty(tn) for tn in elem.findall("TreeNumberList/TreeNumber")]
        tree_numbers = [t for t in tree_numbers if t]

        yield MeshDescriptor(
            descriptor_ui=descriptor_ui,
            descriptor_name=descriptor_name,
            scope_note=scope_note,
            entry_terms=entry_terms,
            tree_numbers=tree_numbers,
        )

        # Free memory by clearing processed subtree
        elem.clear()

        count += 1
        if max_terms and count >= max_terms:
            break


# -----------------------------
# Embedding backends
# -----------------------------
def embed_texts_openai(texts: List[str], model: str, batch_size: int) -> np.ndarray:
    """
    Uses OpenAI embeddings. Requires:
      pip install openai
      export OPENAI_API_KEY=...
    """
    try:
        from openai import OpenAI
    except Exception as e:
        raise RuntimeError("OpenAI package not installed. Run: pip install openai") from e

    client = OpenAI()

    vectors: List[List[float]] = []
    for i in range(0, len(texts), batch_size):
        batch = texts[i : i + batch_size]
        resp = client.embeddings.create(model=model, input=batch)
        # Preserve order
        vectors.extend([d.embedding for d in resp.data])

    arr = np.array(vectors, dtype=np.float32)
    return arr


def embed_texts_sentence_transformers(texts: List[str], model: str, batch_size: int) -> np.ndarray:
    """
    Offline embeddings with sentence-transformers.
      pip install sentence-transformers
    """
    try:
        from sentence_transformers import SentenceTransformer
    except Exception as e:
        raise RuntimeError(
            "sentence-transformers not installed. Run: pip install sentence-transformers"
        ) from e

    st = SentenceTransformer(model)
    emb = st.encode(
        texts,
        batch_size=batch_size,
        show_progress_bar=True,
        convert_to_numpy=True,
        normalize_embeddings=True,
    )
    return emb.astype(np.float32)


def maybe_build_faiss_index(embeddings: np.ndarray, out_path: Path) -> None:
    """
    Build a FAISS index (Inner Product). With normalized vectors this approximates cosine similarity.
      pip install faiss-cpu
    """
    try:
        import faiss  # type: ignore
    except Exception as e:
        raise RuntimeError("faiss not installed. Run: pip install faiss-cpu") from e

    d = embeddings.shape[1]
    index = faiss.IndexFlatIP(d)
    index.add(embeddings)
    faiss.write_index(index, str(out_path))


# -----------------------------
# Main
# -----------------------------
def main() -> int:
    p = argparse.ArgumentParser(description="Parse MeSH descXXXX.xml to embeddings")
    p.add_argument("--desc_xml", required=True, help="Path to desc2025.xml")
    p.add_argument("--out_dir", required=True, help="Output directory")

    p.add_argument(
        "--backend",
        choices=["openai", "st"],
        default="openai",
        help="Embedding backend: openai or st (sentence-transformers)",
    )
    p.add_argument(
        "--openai_model",
        default="text-embedding-3-small",
        help="OpenAI embedding model name",
    )
    p.add_argument(
        "--st_model",
        default="sentence-transformers/all-MiniLM-L6-v2",
        help="SentenceTransformers model name",
    )

    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--max_terms", type=int, default=0, help="0 = no limit (for testing)")
    p.add_argument("--build_faiss", action="store_true", help="Also write a FAISS index")

    p.add_argument("--no_scope_note", action="store_true")
    p.add_argument("--no_entry_terms", action="store_true")
    p.add_argument("--no_tree_numbers", action="store_true")
    p.add_argument("--max_entry_terms", type=int, default=30)

    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    jsonl_path = out_dir / "mesh_descriptors.jsonl"
    emb_path = out_dir / "mesh_embeddings.npy"
    faiss_path = out_dir / "mesh_faiss.index"

    # 1) Parse and write JSONL + collect embedding texts
    records: List[MeshDescriptor] = []
    texts: List[str] = []

    with jsonl_path.open("w", encoding="utf-8") as f:
        for rec in parse_desc_xml(args.desc_xml, max_terms=args.max_terms):
            records.append(rec)
            text = rec.to_embedding_text(
                include_scope_note=not args.no_scope_note,
                include_entry_terms=not args.no_entry_terms,
                include_tree_numbers=not args.no_tree_numbers,
                max_entry_terms=args.max_entry_terms,
            )
            texts.append(text)

            # Write JSONL record (includes the exact embedding_text used)
            row = asdict(rec)
            row["embedding_text"] = text
            f.write(json.dumps(row, ensure_ascii=False) + "\n")

    if not records:
        print("No records parsed. Check your XML path.", file=sys.stderr)
        return 2

    # 2) Embed
    if args.backend == "openai":
        embeddings = embed_texts_openai(texts, model=args.openai_model, batch_size=args.batch_size)
    else:
        embeddings = embed_texts_sentence_transformers(texts, model=args.st_model, batch_size=args.batch_size)

    # 3) Save embeddings
    np.save(emb_path, embeddings)

    # 4) Optional FAISS
    if args.build_faiss:
        # If you didn’t normalize (OpenAI embeddings aren’t normalized by default), do it for cosine/IP usage:
        norms = np.linalg.norm(embeddings, axis=1, keepdims=True) + 1e-12
        emb_norm = embeddings / norms
        maybe_build_faiss_index(emb_norm.astype(np.float32), faiss_path)

    print(f"Wrote: {jsonl_path}")
    print(f"Wrote: {emb_path}  shape={embeddings.shape} dtype={embeddings.dtype}")
    if args.build_faiss:
        print(f"Wrote: {faiss_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
