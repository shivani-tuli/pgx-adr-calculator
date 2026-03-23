#!/usr/bin/env python3
"""
Update risk CSV files with meta-analysis-derived risk ratios.

Replaces cherry-picked single-study risk ratios with pooled estimates
from published meta-analyses, adding citation provenance (PMID, CI, etc.)
"""
import csv
import os
from datetime import datetime, timezone

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# ============================================================================
# Meta-analysis-derived risk ratios
# ============================================================================
# Keys: (gene_or_variant, drug_name, phenotype_or_allele)
# Values: dict with updated risk_ratio, ci_lower, ci_upper, meta_analysis_id, pmid

RSID_UPDATES = {
    # CYP2C9/Warfarin: Defined 2013 meta-analysis, PMID 23800783
    ("CYP2C9", "Warfarin", "*2"):   {"rr": 2.05, "ci_lo": 1.36, "ci_hi": 3.10, "ma_id": "MA008", "pmid": "23800783"},
    ("CYP2C9", "Warfarin", "*3"):   {"rr": 3.12, "ci_lo": 1.38, "ci_hi": 17.14,"ma_id": "MA009", "pmid": "23800783"},

    # VKORC1/Warfarin: Defined 2013, PMID 23800783
    ("VKORC1", "Warfarin", "-1639G>A"): {"rr": 2.14, "ci_lo": 1.75, "ci_hi": 2.62, "ma_id": "MA011", "pmid": "23800783"},

    # SLCO1B1/Statins: SEARCH 2008 NEJM (PMID 18650507) + Xiang 2018 meta
    ("SLCO1B1", "Simvastatin", "*5"):  {"rr": 4.50, "ci_lo": 2.60, "ci_hi": 7.70, "ma_id": "MA004", "pmid": "18650507"},
    ("SLCO1B1", "Atorvastatin", "*5"): {"rr": 2.00, "ci_lo": 1.11, "ci_hi": 3.52, "ma_id": "MA007", "pmid": "31636355"},
    ("SLCO1B1", "Rosuvastatin", "*5"): {"rr": 1.60, "ci_lo": 1.20, "ci_hi": 2.16, "ma_id": "MA007", "pmid": "31636355"},

    # CYP2C19/Clopidogrel: Mega 2010 JAMA (PMID 20801498)
    ("CYP2C19", "Clopidogrel", "*2"):  {"rr": 1.57, "ci_lo": 1.13, "ci_hi": 2.16, "ma_id": "MA001", "pmid": "20801498"},
    ("CYP2C19", "Clopidogrel", "*3"):  {"rr": 1.57, "ci_lo": 1.13, "ci_hi": 2.16, "ma_id": "MA001", "pmid": "20801498"},
    ("CYP2C19", "Clopidogrel", "*17"): {"rr": 1.52, "ci_lo": None, "ci_hi": None, "ma_id": "MA003", "pmid": "36450637"},

    # CYP2D6/Codeine: CPIC guideline evidence (PMID 28520350)
    ("CYP2D6", "Codeine", "*4"):     {"rr": 2.00, "ci_lo": None, "ci_hi": None, "ma_id": None, "pmid": "28520350"},
    ("CYP2D6", "Codeine", "*6"):     {"rr": 2.00, "ci_lo": None, "ci_hi": None, "ma_id": None, "pmid": "28520350"},

    # CYP2D6/Tamoxifen: PMID 39540781 meta (2025)
    ("CYP2D6", "Tamoxifen", "*4"):   {"rr": 1.34, "ci_lo": 1.10, "ci_hi": 1.63, "ma_id": "MA021", "pmid": "39540781"},
    ("CYP2D6", "Tamoxifen", "*10"):  {"rr": 1.34, "ci_lo": 1.10, "ci_hi": 1.63, "ma_id": "MA021", "pmid": "39540781"},

    # CYP2B6/Efavirenz: PMID 31636355
    ("CYP2B6", "Efavirenz", "*6"):   {"rr": 1.67, "ci_lo": 1.15, "ci_hi": 2.44, "ma_id": "MA019", "pmid": "31636355"},

    # DPYD/Fluorouracil: PMID 27296574
    ("DPYD", "Fluorouracil", "*2A"):  {"rr": 2.54, "ci_lo": 2.15, "ci_hi": 3.00, "ma_id": "MA012", "pmid": "27296574"},

    # CYP3A5/Tacrolimus: PMID 25201288
    ("CYP3A5", "Tacrolimus", "*3"):  {"rr": 1.80, "ci_lo": None, "ci_hi": None, "ma_id": "MA022", "pmid": "25201288"},
}

# Diplotype-level updates: (gene, phenotype_contains, drug)
DIPLOTYPE_UPDATES = {
    # CYP2C19/Clopidogrel
    ("CYP2C19", "Poor Metabolizer", "Clopidogrel"):          {"rr": 2.81, "ci_lo": 1.81, "ci_hi": 4.37, "ma_id": "MA002", "pmid": "20801498"},
    ("CYP2C19", "Intermediate Metabolizer", "Clopidogrel"):  {"rr": 1.57, "ci_lo": 1.13, "ci_hi": 2.16, "ma_id": "MA001", "pmid": "20801498"},
    ("CYP2C19", "Ultrarapid Metabolizer", "Clopidogrel"):    {"rr": 1.52, "ci_lo": None, "ci_hi": None, "ma_id": "MA003", "pmid": "36450637"},
    ("CYP2C19", "Rapid Metabolizer", "Clopidogrel"):         {"rr": 1.30, "ci_lo": None, "ci_hi": None, "ma_id": "MA003", "pmid": "36450637"},

    # CYP2C9/Warfarin
    ("CYP2C9", "Intermediate Metabolizer", "Warfarin"):      {"rr": 2.05, "ci_lo": 1.36, "ci_hi": 3.10, "ma_id": "MA008", "pmid": "23800783"},
    ("CYP2C9", "Poor Metabolizer", "Warfarin"):              {"rr": 4.87, "ci_lo": 1.38, "ci_hi": 17.14,"ma_id": "MA009", "pmid": "23800783"},

    # SLCO1B1/Simvastatin
    ("SLCO1B1", "Intermediate Function", "Simvastatin"):     {"rr": 2.90, "ci_lo": 1.59, "ci_hi": 5.34, "ma_id": "MA006", "pmid": "31636355"},
    ("SLCO1B1", "Poor Function", "Simvastatin"):             {"rr": 16.90,"ci_lo": 4.70, "ci_hi": 61.10,"ma_id": "MA005", "pmid": "18650507"},

    # DPYD/Fluorouracil + Capecitabine
    ("DPYD", "Intermediate Metabolizer", "Fluorouracil"):    {"rr": 2.54, "ci_lo": 2.15, "ci_hi": 3.00, "ma_id": "MA012", "pmid": "27296574"},
    ("DPYD", "Poor Metabolizer", "Fluorouracil"):            {"rr": 25.60,"ci_lo":12.10, "ci_hi": 53.90,"ma_id": "MA013", "pmid": "27296574"},
    ("DPYD", "Intermediate Metabolizer", "Capecitabine"):    {"rr": 2.54, "ci_lo": 2.15, "ci_hi": 3.00, "ma_id": "MA012", "pmid": "27296574"},
    ("DPYD", "Poor Metabolizer", "Capecitabine"):            {"rr": 25.60,"ci_lo":12.10, "ci_hi": 53.90,"ma_id": "MA013", "pmid": "27296574"},

    # TPMT/Thiopurines
    ("TPMT", "Intermediate Metabolizer", "Mercaptopurine"):  {"rr": 3.90, "ci_lo": 2.50, "ci_hi": 6.10, "ma_id": "MA014", "pmid": "31342537"},
    ("TPMT", "Poor Metabolizer", "Mercaptopurine"):          {"rr": 6.67, "ci_lo": 3.88, "ci_hi": 11.47,"ma_id": "MA015", "pmid": "25201288"},
    ("TPMT", "Intermediate Metabolizer", "Azathioprine"):    {"rr": 3.90, "ci_lo": 2.50, "ci_hi": 6.10, "ma_id": "MA014", "pmid": "31342537"},
    ("TPMT", "Poor Metabolizer", "Azathioprine"):            {"rr": 6.67, "ci_lo": 3.88, "ci_hi": 11.47,"ma_id": "MA015", "pmid": "25201288"},
    ("TPMT", "Intermediate Metabolizer", "Thioguanine"):     {"rr": 3.90, "ci_lo": 2.50, "ci_hi": 6.10, "ma_id": "MA014", "pmid": "31342537"},

    # UGT1A1/Irinotecan
    ("UGT1A1", "Intermediate Metabolizer", "Irinotecan"):    {"rr": 1.90, "ci_lo": 1.44, "ci_hi": 2.51, "ma_id": "MA018", "pmid": "23529007"},
    ("UGT1A1", "Poor Metabolizer", "Irinotecan"):            {"rr": 3.44, "ci_lo": 2.45, "ci_hi": 4.82, "ma_id": "MA016", "pmid": "23529007"},

    # CYP2D6/Tamoxifen
    ("CYP2D6", "Poor Metabolizer", "Tamoxifen"):             {"rr": 1.34, "ci_lo": 1.10, "ci_hi": 1.63, "ma_id": "MA021", "pmid": "39540781"},
    ("CYP2D6", "Intermediate Metabolizer", "Tamoxifen"):     {"rr": 1.34, "ci_lo": 1.10, "ci_hi": 1.63, "ma_id": "MA021", "pmid": "39540781"},

    # CYP2B6/Efavirenz
    ("CYP2B6", "Intermediate Metabolizer", "Efavirenz"):     {"rr": 1.47, "ci_lo": 1.10, "ci_hi": 1.96, "ma_id": "MA020", "pmid": "31636355"},
    ("CYP2B6", "Poor Metabolizer", "Efavirenz"):             {"rr": 1.67, "ci_lo": 1.15, "ci_hi": 2.44, "ma_id": "MA019", "pmid": "31636355"},

    # CYP3A5/Tacrolimus
    ("CYP3A5", "Poor Metabolizer", "Tacrolimus"):            {"rr": 1.80, "ci_lo": None, "ci_hi": None, "ma_id": "MA022", "pmid": "25201288"},
    ("CYP3A5", "Intermediate Metabolizer", "Tacrolimus"):    {"rr": 1.30, "ci_lo": None, "ci_hi": None, "ma_id": "MA022", "pmid": "25201288"},

    # CYP2D6/Codeine
    ("CYP2D6", "Poor Metabolizer", "Codeine"):               {"rr": 2.50, "ci_lo": None, "ci_hi": None, "ma_id": None, "pmid": "28520350"},
    ("CYP2D6", "Intermediate Metabolizer", "Codeine"):       {"rr": 1.50, "ci_lo": None, "ci_hi": None, "ma_id": None, "pmid": "28520350"},
    ("CYP2D6", "Ultrarapid Metabolizer", "Codeine"):         {"rr": 3.50, "ci_lo": None, "ci_hi": None, "ma_id": None, "pmid": "28520350"},
}

def update_rsid_csv():
    """Update pgx_adr_associations.csv with meta-analysis risk ratios."""
    path = os.path.join(DATA_DIR, "pgx_adr_associations.csv")
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    # Add provenance columns if not present
    new_cols = ["meta_analysis_id", "pmid", "ci_lower", "ci_upper", "rr_source"]
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    updated = 0
    now = datetime.now(timezone.utc).isoformat()

    for row in rows:
        gene = row.get("gene", "")
        drug = row.get("drug_name", "")
        allele = row.get("allele_name", "")
        key = (gene, drug, allele)

        if key in RSID_UPDATES:
            u = RSID_UPDATES[key]
            old_rr = row.get("risk_ratio", "")
            row["risk_ratio"] = str(u["rr"])
            row["meta_analysis_id"] = u["ma_id"] or ""
            row["pmid"] = u["pmid"] or ""
            row["ci_lower"] = str(u["ci_lo"]) if u["ci_lo"] else ""
            row["ci_upper"] = str(u["ci_hi"]) if u["ci_hi"] else ""
            row["rr_source"] = "meta-analysis"
            row["source"] = "Meta-analysis"
            row["last_updated"] = now
            updated += 1
            print(f"  [{gene}/{drug}/{allele}] {old_rr} → {u['rr']} (PMID:{u['pmid']})")
        else:
            row.setdefault("meta_analysis_id", "")
            row.setdefault("pmid", "")
            row.setdefault("ci_lower", "")
            row.setdefault("ci_upper", "")
            row.setdefault("rr_source", "curated")

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Updated {updated}/{len(rows)} rsID rows with meta-analysis data\n")


def update_diplotype_csv():
    """Update diplotype_risk_associations.csv with meta-analysis risk ratios."""
    path = os.path.join(DATA_DIR, "diplotype_risk_associations.csv")
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    new_cols = ["meta_analysis_id", "pmid", "ci_lower", "ci_upper", "rr_source"]
    for col in new_cols:
        if col not in fieldnames:
            fieldnames.append(col)

    updated = 0
    now = datetime.now(timezone.utc).isoformat()

    for row in rows:
        gene = row.get("gene", "")
        drug = row.get("drug_name", "")
        pheno = row.get("phenotype", "")
        key = (gene, pheno, drug)

        if key in DIPLOTYPE_UPDATES:
            u = DIPLOTYPE_UPDATES[key]
            old_rr = row.get("risk_ratio", "")
            row["risk_ratio"] = str(u["rr"])
            row["meta_analysis_id"] = u["ma_id"] or ""
            row["pmid"] = u["pmid"] or ""
            row["ci_lower"] = str(u["ci_lo"]) if u["ci_lo"] else ""
            row["ci_upper"] = str(u["ci_hi"]) if u["ci_hi"] else ""
            row["rr_source"] = "meta-analysis"
            row["source"] = "Meta-analysis"
            row["last_updated"] = now
            updated += 1
            print(f"  [{gene}/{pheno}/{drug}] {old_rr} → {u['rr']} (PMID:{u['pmid']})")
        else:
            row.setdefault("meta_analysis_id", "")
            row.setdefault("pmid", "")
            row.setdefault("ci_lower", "")
            row.setdefault("ci_upper", "")
            row.setdefault("rr_source", "curated")

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Updated {updated}/{len(rows)} diplotype rows with meta-analysis data\n")


if __name__ == "__main__":
    print("=== Meta-Analysis Risk Ratio Updater ===\n")

    print("--- Updating rsID associations ---")
    update_rsid_csv()

    print("--- Updating diplotype associations ---")
    update_diplotype_csv()

    print("=== Done ===")
