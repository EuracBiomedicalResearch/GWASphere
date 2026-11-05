#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


if len(sys.argv) != 4:
    raise SystemExit("Usage: generate_report_insights.py <evidence.json> <assessment.json> <status.json>")

EVIDENCE_PATH   = Path(sys.argv[1]).resolve()
ASSESSMENT_PATH = Path(sys.argv[2]).resolve()
STATUS_PATH     = Path(sys.argv[3]).resolve()
MODEL           = os.getenv("OLLAMA_MODEL", "qwen3:8b")
OLLAMA_URL      = os.getenv("OLLAMA_URL", "http://127.0.0.1:11434/api/generate")
NUM_CTX         = int(os.getenv("OLLAMA_NUM_CTX", "16384"))
NUM_PREDICT     = int(os.getenv("OLLAMA_NUM_PREDICT", "2600"))
NUM_THREADS     = int(os.getenv("OLLAMA_NUM_THREADS", "8"))


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f"{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def write_status(status: str, error: str | None = None) -> None:
    atomic_json(STATUS_PATH, {
        "status":     status,
        "job_id":     os.getenv("SLURM_JOB_ID"),
        "model":      MODEL,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "error":      error,
    })


def ask_ollama(prompt: str) -> dict[str, Any]:
    request = {
        "model":  MODEL,
        "prompt": prompt,
        "stream": False,
        "think":  False,
        "format": "json",
        "options": {
            "num_ctx":        NUM_CTX,
            "num_predict":    NUM_PREDICT,
            "num_thread":     NUM_THREADS,
            "temperature":    0.0,
            "top_p":          0.8,
            "repeat_penalty": 1.05,
            "seed":           42,
        },
    }
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            response = requests.post(OLLAMA_URL, json=request, timeout=(10, 1800))
            response.raise_for_status()
            result = json.loads(response.json()["response"])
            if not isinstance(result, dict):
                raise ValueError("Model response is not a JSON object")
            return result
        except Exception as exc:
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Ollama request failed: {last_error}")


def make_prompt(evidence: dict[str, Any]) -> str:
    trait    = evidence.get("trait", "")
    datasets = evidence.get("datasets", [])
    
    # Build a structured summary of loci for the prompt
    loci_summary = []
    for ds in datasets:
        label = ds.get("dataset_label", "")
        novelty_check_performed = ds.get("novelty_check_performed", True)
        # Build loci summary directly from novelty rows — each row is a lead SNP
        # with GENE_NAME, NOVEL, KNOWN_* already populated by novel_detection.py
        for nov in ds.get("novelty", [])[:10]:
                if not isinstance(nov, dict):
                    continue

                snp_id   = str(nov.get("SNPID", "")).strip()
                rsid     = str(nov.get("rsID",  "")).strip()
                gene     = str(nov.get("GENE_NAME", "")).strip()
                log10p   = nov.get("LOG10P", "")
                chrom    = nov.get("CHR", "")
                is_novel = nov.get("NOVEL", None)

                # Novel status
                if not novelty_check_performed:
                    novel_str    = "Not assessed"
                    known_detail = "Novelty classification was not performed (no trait string or EFO mapping provided). Cross-trait GWAS Catalog associations were assessed."
                elif is_novel is True or str(is_novel).lower() in ("true", "1"):
                    novel_str    = "Novel"
                    known_detail = "Novel: no association reported within 500 kb for this trait in the GWAS Catalog"
                elif is_novel is False or str(is_novel).lower() in ("false", "0"):
                    novel_str = "Known"
                    author  = str(nov.get("KNOWN_AUTHOR", "")).strip()
                    pubmed  = str(nov.get("KNOWN_PUBMED_ID", "")).strip()
                    dist    = nov.get("DISTANCE_TO_KNOWN", "")
                    known_detail = "Already reported in the GWAS Catalog for this trait"
                    if author and author != "nan":
                        known_detail += f" (first reported by {author}"
                        if pubmed and pubmed not in ("nan", "None", ""):
                            known_detail += f", PMID: {pubmed}"
                        known_detail += ")"
                    if dist and str(dist) not in ("nan", "None", ""):
                        known_detail += f"; nearest known association is {int(float(dist))} bp away"
                else:
                    novel_str    = "Unknown"
                    known_detail = ""

                # PubMed hits for this gene from the pubmed section
                pubmed_hits = []
                if gene and gene != "nan":
                    for pub in ds.get("pubmed", [])[:50]:
                        if isinstance(pub, dict):
                            pub_gene = str(pub.get("gene", "")).strip().lower()
                            if gene.lower() in pub_gene or pub_gene in gene.lower():
                                pubmed_hits.append({
                                    "title": pub.get("title", ""),
                                    "year":  pub.get("year", ""),
                                    "pmid":  pub.get("pmid", ""),
                                })

                loci_summary.append({
                    "dataset":      label,
                    "lead_snp":     snp_id,
                    "lead_rsid":    rsid,
                    "chr":          chrom,
                    "log10p":       log10p,
                    "nearest_gene": gene if gene and gene != "nan" else "not annotated",
                    "novel":        novel_str,
                    "known_detail": known_detail,
                    "pubmed_hits":  pubmed_hits[:3],
                })

    # Top MAGMA gene-sets (top 10, flag significance)
    top_genesets = []
    for ds in datasets:
        label = ds.get("dataset_label", "")
        for gs in ds.get("magma_genesets", [])[:10]:
            if isinstance(gs, dict):
                p_val = gs.get("P", None)
                try:
                    p_float = float(p_val)
                    significant = p_float < 0.05
                except (TypeError, ValueError):
                    p_float = None
                    significant = False
                top_genesets.append({
                    "dataset":     label,
                    "geneset":     gs.get("FULL_NAME", gs.get("VARIABLE", "")),
                    "p":           p_val,
                    "significant": significant,
                })

    # Top MAGMA tissues (top 10, flag significance)
    top_tissues = []
    for ds in datasets:
        label = ds.get("dataset_label", "")
        for tis in ds.get("magma_tissues", [])[:10]:
            if isinstance(tis, dict):
                p_val = tis.get("P", None)
                try:
                    p_float = float(p_val)
                    significant = p_float < 0.05
                except (TypeError, ValueError):
                    p_float = None
                    significant = False
                top_tissues.append({
                    "dataset":     label,
                    "tissue":      tis.get("VARIABLE", ""),
                    "p":           p_val,
                    "significant": significant,
                })

    # Top MAGMA genes (top 10, flag significance)
    top_genes = []
    for ds in datasets:
        label = ds.get("dataset_label", "")
        for g in ds.get("magma_genes", [])[:10]:
            if isinstance(g, dict):
                p_val = g.get("P", None)
                try:
                    p_float = float(p_val)
                    significant = p_float < 0.05
                except (TypeError, ValueError):
                    p_float = None
                    significant = False
                symbol = str(g.get("SYMBOL", g.get("GENE", ""))).strip()
                top_genes.append({
                    "dataset":     label,
                    "gene":        symbol,
                    "p":           p_val,
                    "significant": significant,
                })

    # Compute locus-gene overlap: loci whose nearest gene is also a top MAGMA gene
    top_gene_names = {g["gene"].lower() for g in top_genes if g["gene"]}
    for locus in loci_summary:
        gene = locus.get("nearest_gene", "").lower()
        locus["magma_gene_overlap"] = gene in top_gene_names and gene not in ("", "not annotated")

    # Compute counts for the prompt context
    n_loci = len(loci_summary)
    n_novel = sum(1 for l in loci_summary if l["novel"] == "Novel")
    n_known = sum(1 for l in loci_summary if l["novel"] == "Known")
    n_with_pubmed = sum(1 for l in loci_summary if l["pubmed_hits"])
    n_sig_genesets = sum(1 for g in top_genesets if g["significant"])
    n_sig_tissues = sum(1 for t in top_tissues if t["significant"])
    n_overlap = sum(1 for l in loci_summary if l.get("magma_gene_overlap"))

    return f"""
You are writing a concise, prioritized scientific summary of GWAS results for trait "{trait}".
Use ONLY the evidence provided below. Do NOT use prior knowledge or external information.
Do NOT invent SNP IDs, gene names, p-values, PMIDs, or any other values.
If a field is empty or missing, omit it rather than speculating.

CONTEXT:
- Total loci: {n_loci} ({n_known} known, {n_novel} novel)
- Novelty classification: {"NOT PERFORMED — no trait string or EFO mapping was provided. Do not mention novel/known status. Cross-trait GWAS Catalog comparison was performed." if any(l["novel"] == "Not assessed" for l in loci_summary) else "performed"}
- Loci with PubMed support for the nearest gene: {n_with_pubmed}
- Loci whose nearest gene is also a top MAGMA gene (convergent evidence): {n_overlap}
- Significant MAGMA gene-sets (P < 0.05): {n_sig_genesets}
- Significant MAGMA tissues (P < 0.05): {n_sig_tissues}

YOUR TASK — follow these rules strictly:

SUMMARY (2-3 sentences):
- State how many loci were found
- If novelty classification was performed: state how many are novel vs known. If most loci are novel and lack PubMed support, say so explicitly — do not speculate about their function
- If novelty classification was NOT performed: do NOT mention novel/known status. Instead note that cross-trait GWAS Catalog associations were assessed for each locus
- Highlight only loci where pubmed_hits is non-empty as having biological support

LOCI (prioritize as follows, in this order):
1. Known loci with pubmed_hits AND magma_gene_overlap=true — strongest convergent evidence
2. Novel loci with pubmed_hits AND magma_gene_overlap=true — top novel candidates
3. Known loci with pubmed_hits — mention PMID
4. Novel loci with pubmed_hits — mention gene and PMID
5. Any locus with magma_gene_overlap=true but no pubmed_hits — note the MAGMA convergence
6. Known loci without pubmed_hits — replications only, brief
7. Novel loci without pubmed_hits — list gene name only, write "No supporting literature found"
DO NOT write speculative interpretations like "may influence expression of X" when pubmed_hits is empty.
IMPORTANT: Known loci always have GWAS Catalog support — NEVER write "No supporting literature found" for a Known locus.
For Known loci, "pubmed_support" must always mention the GWAS Catalog reference from known_detail.
Note that the nearest_gene in this study may differ from the gene in the original GWAS Catalog entry — this is expected and should be stated explicitly.
If magma_gene_overlap is true, note: "Gene also identified in MAGMA gene-based analysis — convergent evidence."

PATHWAY ENRICHMENT:
- Report the top gene-sets with their p-values
- Comment on whether they are biologically plausible for the trait "{trait}"
- If none are an obvious fit for the trait, explicitly state: "The top pathways do not show an obvious connection to {trait} — this may reflect limited statistical power"

TISSUE ENRICHMENT:
- Report the top tissues
- Comment on biological plausibility for "{trait}"
- If a tissue is unexpected (e.g. Testis for a metabolic trait), explicitly note: "unexpected for this trait — interpret with caution"

KEY TAKEAWAY:
- One sentence on the single most actionable or interesting finding
- Prioritize: (a) a novel locus with PubMed support, or (b) a known locus replication, or (c) an enriched tissue that makes biological sense
- If nothing stands out clearly, write: "Results require replication in a larger sample to draw firm conclusions"

LOCI DATA (magma_gene_overlap=true means the nearest gene is also a top MAGMA gene):
{json.dumps(loci_summary, ensure_ascii=False, indent=2)}

TOP MAGMA GENES (gene-based analysis):
{json.dumps(top_genes, ensure_ascii=False, indent=2)}

TOP MAGMA GENE-SETS:
{json.dumps(top_genesets, ensure_ascii=False, indent=2)}

TOP MAGMA TISSUES:
{json.dumps(top_tissues, ensure_ascii=False, indent=2)}

Return exactly this JSON structure:
{{
  "summary": "2-3 sentence overall summary, honest about strength of evidence",
  "loci": [
    {{
      "locus_id": "locus identifier",
      "lead_snp": "SNP ID",
      "nearest_gene": "gene name",
      "log10p": "value",
      "novelty": "Novel | Known | Unknown",
      "novelty_detail": "copy the known_detail field exactly from the locus data",
      "pubmed_support": "Follow exactly one of these four cases: (1) Known + pubmed_hits non-empty: 'Locus previously reported in GWAS Catalog ([known_detail]). Additional PubMed literature for nearest gene: [title, PMID from pubmed_hits].' (2) Known + pubmed_hits empty: 'Locus previously reported in GWAS Catalog ([known_detail]). No additional PubMed literature found for nearest gene [nearest_gene] in {trait}.' (3) Novel + pubmed_hits non-empty: 'Novel locus. Supporting literature for nearest gene [nearest_gene]: [title, PMID from pubmed_hits].' (4) Novel + pubmed_hits empty: 'Novel locus. No supporting literature found for [nearest_gene] in {trait}.'",
      "magma_convergence": "'Gene also identified in MAGMA gene-based analysis' if magma_gene_overlap is true, else ''",
      "interpretation": "one sentence grounded only in provided data. Case 1 (Known+pubmed): state locus replicates known association and note nearest gene has independent literature support. Case 2 (Known+no pubmed): state locus replicates known association; note nearest gene [nearest_gene] may differ from originally reported gene. Case 3 (Novel+pubmed): briefly describe gene biology from pubmed_hits. Case 4 (Novel+no pubmed): write 'Novel locus; no supporting literature found for [nearest_gene] in {trait}'."
    }}
  ],
  "pathway_enrichment": "top gene-sets with p-values and biological plausibility comment, or empty string",
  "tissue_enrichment": "top tissues with plausibility comment; flag unexpected tissues explicitly, or empty string",
  "key_takeaway": "one sentence most actionable finding, or statement that replication is needed"
}}
"""


def main() -> None:
    try:
        write_status("running")
        with EVIDENCE_PATH.open(encoding="utf-8") as handle:
            evidence = json.load(handle)
        raw        = ask_ollama(make_prompt(evidence))
        assessment = {
            "schema_version": "2.0",
            "model":          MODEL,
            "trait":          evidence.get("trait", ""),
            "generated_at":   datetime.now(timezone.utc).isoformat(),
            "summary":        raw.get("summary", ""),
            "loci":           raw.get("loci", []),
            "pathway_enrichment": raw.get("pathway_enrichment", ""),
            "tissue_enrichment":  raw.get("tissue_enrichment", ""),
            "key_takeaway":   raw.get("key_takeaway", ""),
        }
        atomic_json(ASSESSMENT_PATH, assessment)
        write_status("completed")
    except Exception as exc:
        write_status("failed", str(exc))
        raise


if __name__ == "__main__":
    main()
