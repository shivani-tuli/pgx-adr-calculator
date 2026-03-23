#!/usr/bin/env python3
"""
FAERS Validation Module
=======================
Validates the PGx-ADR Calculator's predicted risk ratios against real-world
adverse event data from the FDA Adverse Event Reporting System (FAERS) via
the public openFDA API.

Computes Proportional Reporting Ratios (PRRs) and correlates them with our
meta-analysis-derived risk ratios using Spearman rank correlation.

Usage:
    python3 scripts/validate_faers.py
"""
import csv
import json
import math
import os
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime, timezone

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
OPENFDA_BASE = "https://api.fda.gov/drug/event.json"

# Rate limit: openFDA allows 240 requests per minute without API key
RATE_LIMIT_DELAY = 0.35  # seconds between requests

# ============================================================================
# ADR Term Mapping: our primary_adr → MedDRA preferred terms for openFDA
# ============================================================================
# openFDA uses MedDRA preferred terms (patient.reaction.reactionmeddrapt)
# We map our database ADR labels to the closest MedDRA terms
ADR_TO_MEDDRA = {
    # Warfarin
    "Major Bleeding": ["HAEMORRHAGE", "GASTROINTESTINAL HAEMORRHAGE",
                        "CEREBROVASCULAR ACCIDENT", "HAEMATURIA"],
    "Dose Sensitivity (Bleeding if unadjusted)": ["INTERNATIONAL NORMALISED RATIO INCREASED",
                                                    "HAEMORRHAGE"],
    "Dose Sensitivity": ["INTERNATIONAL NORMALISED RATIO INCREASED"],
    "Dose Sensitivity (Mild)": ["INTERNATIONAL NORMALISED RATIO INCREASED"],
    "Dose Sensitivity (High Risk)": ["INTERNATIONAL NORMALISED RATIO INCREASED",
                                      "HAEMORRHAGE"],
    "Major Bleeding Risk": ["HAEMORRHAGE", "GASTROINTESTINAL HAEMORRHAGE"],
    "Dose Adjustment Needed": ["INTERNATIONAL NORMALISED RATIO INCREASED"],

    # Statins
    "Myopathy/Rhabdomyolysis": ["RHABDOMYOLYSIS", "MYOPATHY", "MYALGIA",
                                 "BLOOD CREATINE PHOSPHOKINASE INCREASED"],
    "Myopathy": ["MYOPATHY", "MYALGIA", "MUSCULAR WEAKNESS"],
    "Rhabdomyolysis": ["RHABDOMYOLYSIS"],

    # Clopidogrel
    "Major Cardiovascular Events": ["MYOCARDIAL INFARCTION", "STENT THROMBOSIS",
                                     "CEREBROVASCULAR ACCIDENT"],
    "Increased MACE Risk": ["MYOCARDIAL INFARCTION", "STENT THROMBOSIS"],
    "Bleeding": ["HAEMORRHAGE", "GASTROINTESTINAL HAEMORRHAGE"],
    "Increased Bleeding": ["HAEMORRHAGE"],

    # SSRIs
    "Treatment Failure": ["DRUG INEFFECTIVE", "DEPRESSION",
                           "THERAPEUTIC RESPONSE DECREASED"],
    "Treatment Failure/QT Prolongation": ["DRUG INEFFECTIVE",
                                           "ELECTROCARDIOGRAM QT PROLONGED"],
    "QT Prolongation/Toxicity": ["ELECTROCARDIOGRAM QT PROLONGED"],
    "Subtherapeutic Levels": ["DRUG INEFFECTIVE", "DRUG LEVEL DECREASED"],

    # Opioids
    "Treatment Failure (Analgesic)": ["DRUG INEFFECTIVE", "PAIN"],
    "Reduced Analgesic Effect": ["DRUG INEFFECTIVE", "PAIN"],
    "Respiratory Depression/Toxicity": ["RESPIRATORY DEPRESSION",
                                         "OVERDOSE"],
    "Ultra-rapid Metabolism Toxicity": ["RESPIRATORY DEPRESSION", "OVERDOSE"],
    "Altered Response": ["DRUG INEFFECTIVE"],

    # Efavirenz
    "CNS Toxicity": ["DIZZINESS", "INSOMNIA", "HALLUCINATION",
                      "ABNORMAL DREAMS", "DEPRESSION"],
    "Severe CNS Toxicity": ["DIZZINESS", "HALLUCINATION", "PSYCHOTIC DISORDER"],

    # Tacrolimus
    "Nephrotoxicity": ["RENAL IMPAIRMENT", "BLOOD CREATININE INCREASED",
                        "RENAL FAILURE"],
    "Nephrotoxicity (Overexposure)": ["RENAL IMPAIRMENT",
                                       "DRUG LEVEL INCREASED"],
    "Standard Dosing Needed": ["DRUG LEVEL INCREASED"],

    # Fluoropyrimidines
    "Severe Toxicity": ["NEUTROPENIA", "DIARRHOEA", "STOMATITIS",
                         "MUCOSITIS ORAL"],
    "Fatal Toxicity": ["DEATH", "NEUTROPENIA", "SEPSIS"],

    # Thiopurines
    "Myelosuppression": ["NEUTROPENIA", "LEUKOPENIA", "PANCYTOPENIA",
                          "BONE MARROW FAILURE"],
    "Severe Myelosuppression/Death": ["PANCYTOPENIA", "FEBRILE NEUTROPENIA",
                                       "DEATH"],
    "Severe Myelosuppression": ["PANCYTOPENIA", "FEBRILE NEUTROPENIA"],

    # Irinotecan
    "Neutropenia/Diarrhea": ["NEUTROPENIA", "DIARRHOEA"],
    "Severe Neutropenia/Diarrhea": ["FEBRILE NEUTROPENIA", "DIARRHOEA"],

    # Voriconazole
    "Hepatotoxicity": ["HEPATOTOXICITY", "HEPATIC FUNCTION ABNORMAL",
                        "ALANINE AMINOTRANSFERASE INCREASED"],
    "Hepatotoxicity/Toxicity": ["HEPATOTOXICITY", "DRUG LEVEL INCREASED"],
    "Reduced Efficacy": ["DRUG INEFFECTIVE"],
    "Reduced Efficacy/Recurrence": ["DRUG INEFFECTIVE", "NEOPLASM RECURRENCE"],

    # Tamoxifen
    "Reduced Efficacy/Recurrence": ["DRUG INEFFECTIVE",
                                     "BREAST CANCER RECURRENT"],

    # Atomoxetine
    "Increased Side Effects": ["NAUSEA", "INSOMNIA", "DIZZINESS",
                                "TACHYCARDIA"],
    "Increased Side Effects/Toxicity": ["NAUSEA", "TACHYCARDIA", "OVERDOSE"],

    # Methadone
    "QT Prolongation": ["ELECTROCARDIOGRAM QT PROLONGED"],
    "QT Prolongation Risk": ["ELECTROCARDIOGRAM QT PROLONGED"],

    # Phenytoin
    "Toxicity/SJS Risk": ["STEVENS-JOHNSON SYNDROME", "DRUG TOXICITY"],
    "Severe Toxicity/SJS Risk": ["STEVENS-JOHNSON SYNDROME", "DRUG TOXICITY"],

    # Celecoxib
    "GI/CV Toxicity Risk": ["GASTROINTESTINAL HAEMORRHAGE",
                             "MYOCARDIAL INFARCTION"],
    "Severe GI/CV Toxicity": ["GASTROINTESTINAL HAEMORRHAGE",
                               "MYOCARDIAL INFARCTION"],

    # Atazanavir
    "Hyperbilirubinemia": ["HYPERBILIRUBINAEMIA", "JAUNDICE"],
    "Severe Jaundice": ["JAUNDICE"],

    # G6PD
    "Hemolytic Anemia Risk": ["HAEMOLYTIC ANAEMIA"],
    "Severe Hemolytic Anemia": ["HAEMOLYTIC ANAEMIA"],

    # IFNL3
    "Reduced Response": ["DRUG INEFFECTIVE"],

    # Catch-all
    "Reduced Efficacy (Beneficial)": ["DRUG INEFFECTIVE"],
    "Increased Exposure (Beneficial)": ["DRUG LEVEL INCREASED"],
}


def query_openfda(drug_name, meddra_terms=None, count_field=None):
    """Query the openFDA FAERS API."""
    params = {}

    # Build search query
    search_parts = [f'patient.drug.openfda.generic_name:"{drug_name.lower()}"']
    if meddra_terms:
        term_queries = " OR ".join(
            [f'patient.reaction.reactionmeddrapt.exact:"{t}"' for t in meddra_terms]
        )
        search_parts.append(f"({term_queries})")

    params["search"] = " AND ".join(search_parts)

    if count_field:
        params["count"] = count_field
    else:
        params["limit"] = "1"

    url = f"{OPENFDA_BASE}?{urllib.parse.urlencode(params)}"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PGx-ADR-Calculator/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
            return data
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None  # No results
        print(f"    HTTP {e.code} for {drug_name}: {e.reason}")
        return None
    except Exception as e:
        print(f"    Error for {drug_name}: {e}")
        return None


def get_total_reports(drug_name):
    """Get total FAERS reports for a drug."""
    data = query_openfda(drug_name)
    if data and "meta" in data and "results" in data["meta"]:
        return data["meta"]["results"]["total"]
    return 0


def get_adr_reports(drug_name, meddra_terms):
    """Get count of reports mentioning specific ADR terms for a drug."""
    data = query_openfda(drug_name, meddra_terms)
    if data and "meta" in data and "results" in data["meta"]:
        return data["meta"]["results"]["total"]
    return 0


def get_global_adr_count(meddra_terms):
    """Get total reports with these ADR terms across ALL drugs."""
    term_queries = " OR ".join(
        [f'patient.reaction.reactionmeddrapt.exact:"{t}"' for t in meddra_terms]
    )
    params = {
        "search": term_queries,
        "limit": "1"
    }
    url = f"{OPENFDA_BASE}?{urllib.parse.urlencode(params)}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PGx-ADR-Calculator/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
            if data and "meta" in data and "results" in data["meta"]:
                return data["meta"]["results"]["total"]
    except Exception:
        pass
    return 0


def compute_prr(a, b, c, d):
    """
    Compute Proportional Reporting Ratio.

    a = reports of drug with ADR
    b = reports of drug without ADR
    c = reports of all other drugs with ADR
    d = reports of all other drugs without ADR

    PRR = (a / (a+b)) / (c / (c+d))
    """
    if (a + b) == 0 or (c + d) == 0 or c == 0:
        return None
    prr = (a / (a + b)) / (c / (c + d))
    return round(prr, 3)


def run_validation():
    """Main validation pipeline."""
    print("=== FAERS Validation Pipeline ===\n")

    # Load diplotype risk data (unique drug-ADR pairs)
    dip_path = os.path.join(DATA_DIR, "diplotype_risk_associations.csv")
    with open(dip_path, "r") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    # Extract unique drug-ADR-RR combinations
    drug_adr_pairs = {}
    for row in all_rows:
        drug = row["drug_name"]
        adr = row["primary_adr"]
        rr = float(row["risk_ratio"])
        key = (drug, adr)
        if key not in drug_adr_pairs:
            drug_adr_pairs[key] = {
                "drug_name": drug,
                "drug_class": row.get("drug_class", ""),
                "primary_adr": adr,
                "max_rr": rr,
                "gene": row.get("gene", ""),
            }
        else:
            drug_adr_pairs[key]["max_rr"] = max(drug_adr_pairs[key]["max_rr"], rr)

    print(f"Found {len(drug_adr_pairs)} unique drug-ADR pairs to validate\n")

    # Get total FAERS database size (approximate)
    total_db_resp = query_openfda("aspirin")  # Common drug as proxy
    time.sleep(RATE_LIMIT_DELAY)

    results = []
    skipped = 0

    for (drug, adr), info in sorted(drug_adr_pairs.items()):
        meddra_terms = ADR_TO_MEDDRA.get(adr)
        if not meddra_terms:
            print(f"  SKIP: No MedDRA mapping for '{adr}'")
            skipped += 1
            continue

        print(f"  Querying: {drug} / {adr[:40]}...", end=" ", flush=True)

        # Get drug-specific counts
        time.sleep(RATE_LIMIT_DELAY)
        total_drug = get_total_reports(drug)

        if total_drug == 0:
            print(f"→ No FAERS data for {drug}")
            skipped += 1
            continue

        time.sleep(RATE_LIMIT_DELAY)
        adr_drug = get_adr_reports(drug, meddra_terms)

        time.sleep(RATE_LIMIT_DELAY)
        adr_global = get_global_adr_count(meddra_terms)

        # Approximate total database size: ~23M reports
        total_db = 23_000_000

        # PRR calculation
        a = adr_drug
        b = total_drug - adr_drug
        c = adr_global - adr_drug
        d = total_db - total_drug - c

        prr = compute_prr(a, b, c, d)

        result = {
            "drug_name": drug,
            "drug_class": info["drug_class"],
            "primary_adr": adr,
            "gene": info["gene"],
            "predicted_rr": info["max_rr"],
            "faers_total_reports": total_drug,
            "faers_adr_reports": adr_drug,
            "faers_prr": prr,
            "meddra_terms": meddra_terms,
            "log_rr": round(math.log(info["max_rr"]), 3) if info["max_rr"] > 0 else None,
            "log_prr": round(math.log(prr), 3) if prr and prr > 0 else None,
        }
        results.append(result)

        print(f"→ {adr_drug}/{total_drug} reports, PRR={prr}")

    print(f"\n  Validated: {len(results)}, Skipped: {skipped}\n")

    # Compute Spearman rank correlation
    valid = [(r["predicted_rr"], r["faers_prr"])
             for r in results if r["faers_prr"] is not None and r["predicted_rr"] > 0]

    if len(valid) >= 5:
        rho = spearman_correlation(valid)
        print(f"  Spearman ρ (predicted RR vs FAERS PRR): {rho:.3f}")
        print(f"  N pairs: {len(valid)}")
    else:
        rho = None
        print(f"  Too few valid pairs ({len(valid)}) for correlation")

    # Build report
    report = {
        "validation_date": datetime.now(timezone.utc).isoformat(),
        "data_source": "FDA FAERS via openFDA API",
        "api_endpoint": OPENFDA_BASE,
        "methodology": "Proportional Reporting Ratio (PRR) vs predicted Risk Ratio",
        "n_drug_adr_pairs_tested": len(results),
        "n_skipped": skipped,
        "spearman_rho": round(rho, 4) if rho else None,
        "n_valid_pairs": len(valid),
        "interpretation": interpret_correlation(rho) if rho else "Insufficient data",
        "results": results,
    }

    out_path = os.path.join(DATA_DIR, "faers_validation_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\n  Report saved to: {out_path}")
    print(f"\n=== Validation Complete ===")

    return report


def spearman_correlation(pairs):
    """Compute Spearman rank correlation coefficient."""
    n = len(pairs)
    if n < 2:
        return 0

    # Rank x and y
    x_vals = [p[0] for p in pairs]
    y_vals = [p[1] for p in pairs]

    def rank(values):
        indexed = sorted(enumerate(values), key=lambda x: x[1])
        ranks = [0] * len(values)
        for rank_pos, (orig_idx, _) in enumerate(indexed):
            ranks[orig_idx] = rank_pos + 1
        # Handle ties (average rank)
        i = 0
        while i < len(indexed):
            j = i
            while j < len(indexed) and indexed[j][1] == indexed[i][1]:
                j += 1
            avg_rank = (i + j + 1) / 2  # 1-indexed average
            for k in range(i, j):
                ranks[indexed[k][0]] = avg_rank
            i = j
        return ranks

    rx = rank(x_vals)
    ry = rank(y_vals)

    # Spearman = 1 - 6*sum(d^2) / (n*(n^2-1))
    d_sq = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    rho = 1 - (6 * d_sq) / (n * (n ** 2 - 1))
    return rho


def interpret_correlation(rho):
    """Interpret Spearman ρ value."""
    if rho is None:
        return "Insufficient data"
    abs_rho = abs(rho)
    if abs_rho >= 0.7:
        strength = "strong"
    elif abs_rho >= 0.4:
        strength = "moderate"
    elif abs_rho >= 0.2:
        strength = "weak"
    else:
        strength = "negligible"

    direction = "positive" if rho > 0 else "negative"
    return (f"{strength.capitalize()} {direction} correlation (ρ={rho:.3f}). "
            f"{'Predicted risk ratios align well with' if rho > 0.3 else 'Limited alignment between predicted RR and'} "
            f"real-world FAERS adverse event reporting patterns.")


if __name__ == "__main__":
    run_validation()
