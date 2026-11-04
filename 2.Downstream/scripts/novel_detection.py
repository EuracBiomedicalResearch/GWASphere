import gwaslab as gl
import pandas as pd
from pathlib import Path
import numpy as np
import os
import sys

# GWAS Catalog novelty check:
# Lead SNPs identified through clumping were assessed for novelty using GWASLab, by querying the GWAS Catalog (full associations download) filtered to genome-wide significant associations (p ≤ 5×10⁻⁸). For each trait, catalog entries were split into same-trait (matched by EFO identifier) and other-trait subsets; 
# novelty was defined as no known association within 500 kb in the same-trait subset, and cross-trait associations were summarised by aggregating all catalog traits with at least one entry within 500 kb of each lead SNP.

# =========================================================
# CONFIG
# =========================================================
OUT_DIR    = Path(sys.argv[1]) / "novelty_check"
INPUT_DIR  = Path(sys.argv[2])
LOCI_TSV   = Path(sys.argv[1]) / "loci" / "ALL_TRAITS.loci_table.tsv"

TRAIT_MAP_XLSX  =  Path(sys.argv[3])
CATALOG_TSV     =  Path(sys.argv[4])

DEFAULT_TRAIT_STRING = sys.argv[5]

PATTERN         = "*.regenie.gz"
OUT_DIR.mkdir(parents=True, exist_ok=True)

WINDOW_KB  = 500
WINDOW_BP  = WINDOW_KB * 1000
SIG_LEVEL  = 5e-8

# =========================================================
# HELPER FUNCTIONS
# =========================================================
def uniq_sorted(vals):
    vals = [
        str(v).strip()
        for v in vals
        if pd.notna(v)
        and str(v).strip() != ""
        and str(v).strip().lower() != "nan"
    ]
    return sorted(set(vals))


def collapse_for_trait(df_trait):
    trait    = df_trait[TRAIT_COL].iloc[0]
    rsids    = uniq_sorted(df_trait[RS_COL].tolist()) if RS_COL else []
    pubmeds  = uniq_sorted(df_trait["PUBMEDID"].tolist()) if "PUBMEDID" in df_trait.columns else []
    pos_list = uniq_sorted(df_trait["KNOWN_ID"].tolist())
    minp     = df_trait["P_VALUE"].min()
    return pd.Series({
        "trait":        trait,
        "minp":         minp,
        "rsids_str":    "|".join(rsids),
        "pubmeds_str":  "|".join(pubmeds),
        "pos_str":      "|".join(pos_list),
    })


def aggregate_one_lead(df_lead):
    per_trait = (
        df_lead
        .dropna(subset=[TRAIT_COL])
        .groupby(TRAIT_COL, as_index=False)
        .apply(collapse_for_trait)
    )
    per_trait = per_trait.sort_values("minp", ascending=True, na_position="last")
    return pd.Series({
        "N_OTHER_TRAITS":    per_trait.shape[0],
        "OTHER_TRAITS":      ",".join(per_trait["trait"].astype(str)),
        "OTHER_RSID_LIST":   ",".join(per_trait["rsids_str"].astype(str)),
        "OTHER_PUBMED_LIST": ",".join(per_trait["pubmeds_str"].astype(str)),
        "OTHER_POS_LIST":    ",".join(per_trait["pos_str"].astype(str)),
        "MIN_P_OTHER":       per_trait["minp"].min(),
    })


def add_regional_signal_metrics(
    leads, all_hits,
    trait_col="trait", chr_col="CHROM", pos_col="GENPOS",
    score_col="LOG10P", id_col="ID", rsid_col="RSID",
    window_bp=WINDOW_BP
):
    leads    = leads.copy()
    all_hits = all_hits.copy()
    leads[pos_col]    = pd.to_numeric(leads[pos_col],    errors="coerce")
    leads[score_col]  = pd.to_numeric(leads[score_col],  errors="coerce")
    all_hits[pos_col]   = pd.to_numeric(all_hits[pos_col],   errors="coerce")
    all_hits[score_col] = pd.to_numeric(all_hits[score_col], errors="coerce")

    out = []
    for _, lead in leads.iterrows():
        trait      = lead[trait_col]
        chrom      = lead[chr_col]
        pos        = int(lead[pos_col])
        lead_score = float(lead[score_col])

        region = all_hits[
            (all_hits[trait_col].astype(str) == str(trait)) &
            (all_hits[chr_col].astype(str)   == str(chrom)) &
            (all_hits[pos_col] >= pos - window_bp) &
            (all_hits[pos_col] <= pos + window_bp)
        ].dropna(subset=[pos_col, score_col]).sort_values(score_col, ascending=False)

        n_gws = len(region)

        if n_gws >= 2:
            second           = region.iloc[1]
            second_score     = float(second[score_col])
            second_pos       = int(second[pos_col])
            second_id        = second.get(id_col,   "NA")
            second_rsid      = second.get(rsid_col, "NA")
            dist_to_second   = abs(pos - second_pos)
            log10p_gap       = lead_score - second_score
        else:
            second_score = second_pos = dist_to_second = log10p_gap = None
            second_id = second_rsid = "NA"

        if   n_gws >= 10 and log10p_gap is not None and log10p_gap <= 1.5:
            peak_support = "strong_genome_wide_regional_peak"
        elif n_gws >= 3  and log10p_gap is not None and log10p_gap <= 2.0:
            peak_support = "genome_wide_regional_peak"
        elif n_gws >= 2:
            peak_support = "weak_genome_wide_regional_support"
        else:
            peak_support = "isolated_genome_wide_lead_variant"

        lead = lead.copy()
        lead["REGIONAL_WINDOW_BP"]                    = window_bp
        lead["REGIONAL_VARIANTS_GWS"]                 = n_gws
        lead["REGIONAL_SECOND_BEST_ID"]               = second_id
        lead["REGIONAL_SECOND_BEST_RSID"]             = second_rsid
        lead["REGIONAL_SECOND_BEST_LOG10P"]           = second_score
        lead["REGIONAL_SECOND_BEST_POS"]              = second_pos
        lead["REGIONAL_DISTANCE_TO_SECOND_BEST_BP"]   = dist_to_second
        lead["REGIONAL_LEAD_TO_SECOND_LOG10P_GAP"]    = log10p_gap
        lead["REGIONAL_PEAK_SUPPORT"]                 = peak_support
        lead["REGIONAL_LEAD_IS_ISOLATED"]             = peak_support == "isolated_genome_wide_lead_variant"
        out.append(lead)

    return pd.DataFrame(out)


def mask_by_efo(df, efo_id):
    if not efo_id or pd.isna(efo_id) or uri_col is None:
        return pd.Series(False, index=df.index)
    return df[uri_col].astype(str).str.contains(str(efo_id), na=False)


def mask_by_name(df, trait_name):
    trait_name = (trait_name or "").strip()
    if not trait_name:
        return pd.Series(False, index=df.index)
    col = "TRAIT_MAPPED" if "TRAIT_MAPPED" in df.columns else "TRAIT_REPORTED"
    return df[col].astype(str).str.contains(trait_name, case=False, na=False)


def run_novelty_and_other(sig_ss, catalog_same, catalog_other, base, out_dir, window_kb, window_bp, sig_level, subset=None):
    """
    Run GWASLab novelty check for same-trait and other-trait catalog subsets.
    Saves all output files and returns nothing.
    """
    # ---- same trait ----
    if catalog_same is not None and len(catalog_same) > 0:
        res_same = sig_ss.get_novel(
            known=catalog_same,
            only_novel=False,
            windowsizekb_for_novel=window_kb,
            sig_level=sig_level,
            if_get_lead=False,
            use_p=True,
            anno=False,
            wc_correction=False,
            use_cache=True,
            cache_dir="./",
            build="38",
            source="ensembl",
            gwascatalog_source="EBI",
            output_known=True,
            verbose=True,
        )
        novel_same, meta_same = res_same[0], res_same[1]

        if novel_same is not None and len(novel_same) > 0:
            novel_same.to_csv(out_dir / f"{base}_novel.tsv", sep="\t", index=False)

        if meta_same is not None and len(meta_same) > 0:
            meta_same.to_csv(out_dir / f"{base}_novel_metadata.tsv", sep="\t", index=False)
    else:
        print(f"  [same trait] no same-trait GWAS Catalog subset for {base} — writing all loci as NOVEL")
        # Write a minimal _novel.tsv with all lead SNPs marked as novel
        # so PubMed.R can still run a generic functional literature search
        try:
            src = subset.copy() if subset is not None else pd.DataFrame()
            novel_fallback = pd.DataFrame({
                "SNPID":     src.get("ID",        pd.Series(dtype=str)),
                "rsID":      src.get("RSID",      pd.Series(dtype=str)),
                "CHR":       src.get("CHROM",     pd.Series(dtype=str)),
                "POS":       src.get("GENPOS",    pd.Series(dtype=str)),
                "GENE_NAME": src.get("GENE_NAME", pd.Series([""] * len(src))),
                "NOVEL":     True,
            })
            novel_fallback.to_csv(out_dir / f"{base}_novel.tsv", sep="\t", index=False)
            print(f"  [same trait] wrote fallback _novel.tsv with {len(novel_fallback)} loci (all NOVEL=TRUE)")
        except Exception as e:
            print(f"  [same trait] could not write fallback _novel.tsv: {e}")

    # ---- other traits ----
    res_other = sig_ss.get_novel(
        known=catalog_other,
        only_novel=False,
        windowsizekb_for_novel=window_kb,
        sig_level=sig_level,
        if_get_lead=False,
        use_p=True,
        anno=False,
        wc_correction=False,
        use_cache=True,
        cache_dir="./",
        build="38",
        source="ensembl",
        gwascatalog_source="EBI",
        output_known=True,
        verbose=True,
    )
    other_novel, other_meta = res_other[0], res_other[1]

    if other_novel is None or len(other_novel) == 0:
        print(f"  [other traits] no associations for {base}")
        return

    other_novel.to_csv(out_dir / f"{base}_association.tsv", sep="\t", index=False)

    # filter to window
    df = other_novel.copy()
    df["DISTANCE_TO_KNOWN"] = pd.to_numeric(df["DISTANCE_TO_KNOWN"], errors="coerce")
    df_inwin = df[df["DISTANCE_TO_KNOWN"].abs() <= window_bp].copy()
    print(f"  [other traits] total={len(df)}, in window={len(df_inwin)}")

    if df_inwin.empty or other_meta is None or len(other_meta) == 0:
        return

    keep_ids   = set(df_inwin["KNOWN_ID"].dropna().astype(str))
    meta_inwin = other_meta[other_meta["SNPID"].astype(str).isin(keep_ids)].copy()

    joined = df_inwin.merge(
        meta_inwin[["SNPID", "TRAIT_MAPPED", "TRAIT_REPORTED", "P_VALUE", "PUBMEDID", "AUTHOR", "SNPS"]],
        left_on="KNOWN_ID",
        right_on="SNPID",
        how="left",
        suffixes=("", "_CAT"),
    )

    global TRAIT_COL, RS_COL
    TRAIT_COL = "TRAIT_MAPPED" if "TRAIT_MAPPED" in joined.columns else "TRAIT_REPORTED"
    RS_COL    = "SNPS" if "SNPS" in joined.columns else None

    joined["P_VALUE"] = pd.to_numeric(joined["P_VALUE"], errors="coerce")

    out = (
        joined
        .groupby(["SNPID", "rsID", "CHR", "POS"], as_index=False)
        .apply(aggregate_one_lead)
        .reset_index(drop=True)
    )
    out.to_csv(out_dir / f"{base}_othertraits_in_window.tsv", sep="\t", index=False)
    print(f"  Wrote {base}_othertraits_in_window.tsv -> {out.shape}")


# =========================================================
# 0) Load GWAS Catalog
# =========================================================
catalog = pd.read_csv(CATALOG_TSV, sep="\t", dtype=str, low_memory=False)
catalog = catalog.rename(columns={
    "CHR_ID": "CHR",
    "CHR_POS": "POS",
    "MAPPED_TRAIT": "TRAIT_MAPPED",
    "DISEASE/TRAIT": "TRAIT_REPORTED",
    "PUBMEDID": "PUBMEDID",
    "FIRST AUTHOR": "AUTHOR",
    "STUDY ACCESSION": "STUDY_ACCESSION",
    "P-VALUE": "P_VALUE",
    "SNPS": "SNPS",
    "SNP_ID_CURRENT": "SNP_ID_CURRENT"
})
catalog["CHR"] = catalog["CHR"].astype(str).str.replace("chr", "", regex=False)
catalog["CHR"] = catalog["CHR"].replace({"X": "23", "Y": "24", "MT": "25", "M": "25"})
catalog["CHR"] = pd.to_numeric(catalog["CHR"], errors="coerce")
catalog["POS"] = pd.to_numeric(catalog["POS"], errors="coerce")
catalog = catalog.dropna(subset=["CHR", "POS"]).copy()
catalog["CHR"] = catalog["CHR"].astype("Int64")
catalog["POS"] = catalog["POS"].astype(int)
catalog["P_VALUE"] = pd.to_numeric(catalog["P_VALUE"], errors="coerce")
catalog = catalog[(catalog["P_VALUE"].isna()) | (catalog["P_VALUE"] <= SIG_LEVEL)].copy()
print("Catalog loaded:", catalog.shape)

URI_CANDIDATES = ["TRAIT_URI", "MAPPED_TRAIT_URI", "MAPPED_TRAIT_URIS", "EFO_ID", "EFO_URI"]
uri_col = next((c for c in URI_CANDIDATES if c in catalog.columns), None)
print("URI column:", uri_col)

# =========================================================
# 1) Load trait -> EFO mapping
# =========================================================
trait_to_efo = {}
use_trait_map = False

trait_map_path = Path(TRAIT_MAP_XLSX)

if trait_map_path.exists():

    suffix = trait_map_path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        trait_map_df = pd.read_excel(trait_map_path, dtype=str)
    elif suffix == ".csv":
        trait_map_df = pd.read_csv(trait_map_path, dtype=str, sep=",")
    elif suffix in (".tsv", ".txt"):
        trait_map_df = pd.read_csv(trait_map_path, dtype=str, sep="\t")
    else:
        # Try tab-separated as fallback
        try:
            trait_map_df = pd.read_csv(trait_map_path, dtype=str, sep="\t")
        except Exception:
            trait_map_df = pd.read_csv(trait_map_path, dtype=str, sep=",")

    trait_map_df.columns = [c.strip() for c in trait_map_df.columns]

    if "trait" not in trait_map_df.columns or "EFO_ID" not in trait_map_df.columns:
        raise ValueError("Trait map file must contain columns: 'trait' and 'EFO_ID' (accepts .xlsx, .xls, .csv, .tsv)")

    trait_map_df["trait"] = trait_map_df["trait"].astype(str).str.strip()
    trait_map_df["EFO_ID"] = trait_map_df["EFO_ID"].astype(str).str.strip()

    trait_map_df.loc[
        trait_map_df["EFO_ID"].isin(["", "nan", "None", "NA"]),
        "EFO_ID"
    ] = np.nan

    trait_to_efo = dict(zip(trait_map_df["trait"], trait_map_df["EFO_ID"]))

    use_trait_map = True

    print(f"Loaded trait-EFO mapping: {trait_map_path}")

else:

    print(f"[WARNING] Trait-EFO mapping file not found: {trait_map_path}")
    print(f"Fallback: using trait string '{DEFAULT_TRAIT_STRING}'")

    trait_to_efo = {}

# =========================================================
# 2) Load lead SNPs from clumping output
# =========================================================
loci = pd.read_csv(LOCI_TSV, sep="\t", dtype=str)
loci["trait"]    = loci["trait"].str.strip()
loci["leadSNPs"] = loci["leadSNPs"].astype(str)

# explode comma-separated lead SNPs -> one row per lead SNP
loci_exploded = (
    loci.assign(leadSNPs=loci["leadSNPs"].str.split(","))
        .explode("leadSNPs")
)
loci_exploded["leadSNPs"] = loci_exploded["leadSNPs"].str.strip()
loci_exploded = loci_exploded[loci_exploded["leadSNPs"].ne("") & loci_exploded["leadSNPs"].notna()]

# =========================================================
# 3) Process each trait
# =========================================================
TRAIT_COL = None  # will be set inside run_novelty_and_other
RS_COL    = None

for f in sorted(INPUT_DIR.glob(PATTERN)):
    base = f.name.removesuffix(".regenie.gz")
    print(f"\n{'='*40}\n[TRAIT] {base}")

    # --- lead SNPs for this trait ---
    sub_loci = loci_exploded[loci_exploded["trait"] == base].copy()
    if sub_loci.empty:
        print(f"  [SKIP] no lead SNPs in loci table for {base}")
        continue
    lead_set = set(sub_loci["leadSNPs"].dropna().astype(str))

    # --- load regenie, filter to lead SNPs ---
    table = pd.read_csv(f, sep="\t", compression="gzip", dtype=str)
    subset = table[table["ID"].astype(str).isin(lead_set)].copy()
    if subset.empty:
        print(f"  [SKIP] no lead SNPs found in regenie file for {base}")
        continue

    subset["LOG10P"] = pd.to_numeric(subset.get("LOG10P"), errors="coerce")
    subset["P"]      = 10 ** (-subset["LOG10P"])

    # --- add_regional_signal_metrics ---
    # we need all GWS hits for the trait to compute regional metrics
    all_hits_trait = table.copy()
    all_hits_trait["LOG10P"] = pd.to_numeric(all_hits_trait["LOG10P"], errors="coerce")
    all_hits_trait["P"]      = 10 ** (-all_hits_trait["LOG10P"])
    all_hits_trait["trait"]  = base
    all_hits_trait = all_hits_trait[all_hits_trait["P"] <= SIG_LEVEL].copy()

    subset["trait"] = base
    subset_with_metrics = add_regional_signal_metrics(
        leads=subset,
        all_hits=all_hits_trait,
        trait_col="trait",
        chr_col="CHROM",
        pos_col="GENPOS",
        score_col="LOG10P",
        id_col="ID",
        rsid_col="RSID",
        window_bp=WINDOW_BP,
    )
    subset_with_metrics.to_csv(
        OUT_DIR / f"{base}_lead_snps_regional_metrics.tsv",
        sep="\t", index=False
    )

    # --- build GWASLab Sumstats ---
    sig_ss = gl.Sumstats(
        subset,
        snpid="ID",
        rsid="RSID" if "RSID" in subset.columns else None,
        chrom="CHROM",
        pos="GENPOS",
        ea="ALLELE1",
        nea="ALLELE0",
        neaf="A1FREQ",
        beta="BETA",
        se="SE",
        p="P",
        n="N",
        sep="\t",
        build="38",
    )
    
    # --- EFO / trait matching ---
    efo_id = trait_to_efo.get(base)

    if isinstance(efo_id, str) and not efo_id.strip():
        efo_id = None

    print(f"  EFO: {efo_id}")

    if efo_id and uri_col is not None:

        # preferred: EFO-based matching
        mask_same = mask_by_efo(catalog, efo_id)

    else:

        # fallback: trait-name matching
        if not DEFAULT_TRAIT_STRING or not DEFAULT_TRAIT_STRING.strip():
            print("[WARNING] No trait string or EFO mapping provided for this trait.")
            print("[WARNING] Novelty classification (known/novel) will be skipped.")
            print("[WARNING] Other-trait associations in window will still be reported.")
            mask_same = pd.Series([False] * len(catalog), index=catalog.index)
        else:
            mask_same = mask_by_name(catalog, DEFAULT_TRAIT_STRING)

    catalog_same  = catalog[mask_same].copy()
    catalog_other = catalog[~mask_same].copy()

    run_novelty_and_other(
        sig_ss=sig_ss,
        catalog_same=catalog_same,
        catalog_other=catalog_other,
        base=base,
        out_dir=OUT_DIR,
        window_kb=WINDOW_KB,
        window_bp=WINDOW_BP,
        sig_level=SIG_LEVEL,
        subset=subset,
    )

print("\nDone.")