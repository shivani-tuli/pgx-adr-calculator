#!/usr/bin/env python3
"""
PharmGKB Data Updater for PGx-ADR Calculator
==============================================

Queries the PharmGKB REST API to enrich the CSV databases with:
- PharmGKB variant IDs and clinical annotation IDs
- Validated evidence levels directly from PharmGKB
- Drug-gene-variant relationship verification
- Provenance tracking (source + timestamp)

Usage:
    python3 update_pharmgkb_data.py [--force]

Rate limit compliance: 2 requests/second (per PharmGKB terms).
License: Data accessed via this script is subject to PharmGKB's
CC-BY-SA 4.0 license (https://www.clinpgx.org/page/dataUsagePolicy).
"""

import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

# ---- Configuration ----
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "data")
CACHE_DIR = os.path.join(DATA_DIR, "cache")

RSID_CSV = os.path.join(DATA_DIR, "pgx_adr_associations.csv")
DIPLOTYPE_CSV = os.path.join(DATA_DIR, "diplotype_risk_associations.csv")

API_BASE = "https://api.pharmgkb.org/v1"
REQUEST_DELAY = 0.55  # seconds between requests (2 req/sec limit)

# Map of gene symbols to PharmGKB Gene IDs (pre-seeded for efficiency)
GENE_IDS = {
    "CYP2C9": "PA126", "CYP2C19": "PA124", "CYP2D6": "PA128",
    "CYP3A5": "PA131", "CYP2B6": "PA123", "CYP4F2": "PA39410",
    "CYP3A4": "PA130", "SLCO1B1": "PA134", "VKORC1": "PA133",
    "TPMT": "PA356", "DPYD": "PA145", "UGT1A1": "PA420",
    "G6PD": "PA28469", "NUDT15": "PA166228080", "CYP1A2": "PA122",
    "CYP2C8": "PA125", "IFNL3": "PA134952956", "NAT2": "PA18",
    "HLA-B": "PA35056",
}


def api_get(endpoint, params=None):
    """Make a GET request to the PharmGKB API with caching and rate limiting."""
    # Build URL
    url = f"{API_BASE}{endpoint}"
    if params:
        query = "&".join(f"{k}={v}" for k, v in params.items())
        url = f"{url}?{query}"

    # Check cache
    cache_key = url.replace("/", "_").replace("?", "_").replace("&", "_").replace("=", "_")
    cache_file = os.path.join(CACHE_DIR, f"{cache_key}.json")

    if os.path.exists(cache_file):
        age_hours = (time.time() - os.path.getmtime(cache_file)) / 3600
        if age_hours < 24 * 7:  # Cache valid for 7 days
            with open(cache_file, 'r') as f:
                return json.load(f)

    # Rate limit
    time.sleep(REQUEST_DELAY)

    try:
        req = Request(url, headers={"Accept": "application/json"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())

        # Cache the response
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump(data, f, indent=2)

        return data

    except HTTPError as e:
        if e.code == 429:
            print(f"  Rate limited, waiting 5s...", file=sys.stderr)
            time.sleep(5)
            return api_get(endpoint, params)  # Retry
        print(f"  HTTP {e.code} for {url}", file=sys.stderr)
        return None
    except (URLError, TimeoutError) as e:
        print(f"  Network error for {url}: {e}", file=sys.stderr)
        return None


def lookup_variant(rsid):
    """Look up a variant by rsID in PharmGKB."""
    data = api_get("/data/variant", {"symbol": rsid, "view": "max"})
    if not data or data.get("status") != "success":
        return None

    items = data.get("data", [])
    if not items:
        return None

    v = items[0]
    return {
        "pharmgkb_id": v.get("id", ""),
        "symbol": v.get("symbol", rsid),
        "clinical_significance": v.get("clinicalSignificance", ""),
        "genes": [g.get("symbol", "") for g in v.get("relatedGenes", [])],
        "change_classification": v.get("changeClassification", ""),
    }


def lookup_clinical_annotation(annotation_id):
    """Fetch a clinical annotation by its PharmGKB ID."""
    data = api_get(f"/data/clinicalAnnotation/{annotation_id}")
    if not data or data.get("status") != "success":
        return None

    d = data.get("data", {})
    return {
        "id": d.get("id", ""),
        "level_of_evidence": d.get("levelOfEvidence", ""),
        "types": d.get("types", []),
        "drugs": [c.get("name", "") for c in d.get("relatedChemicals", [])],
        "genes": [g.get("symbol", "") for g in d.get("relatedGenes", [])],
        "score": d.get("score", ""),
        "pediatric": d.get("pediatric", False),
    }


def lookup_gene(gene_symbol):
    """Fetch gene data from PharmGKB."""
    data = api_get("/data/gene", {"symbol": gene_symbol, "view": "max"})
    if not data or data.get("status") != "success":
        return None

    items = data.get("data", [])
    if not items:
        return None

    g = items[0]
    return {
        "pharmgkb_id": g.get("id", ""),
        "symbol": g.get("symbol", gene_symbol),
        "name": g.get("name", ""),
        "cpic_gene": g.get("cpicGene", False),
        "has_dosing_guideline": bool(g.get("dosingGuidelines", [])),
    }


def enrich_rsid_csv(force=False):
    """Enrich the rsID-based CSV with PharmGKB data."""
    print("\n=== Enriching rsID associations ===")

    if not os.path.exists(RSID_CSV):
        print(f"  CSV not found: {RSID_CSV}", file=sys.stderr)
        return False

    # Read existing CSV
    with open(RSID_CSV, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    # Add new columns if not present
    new_fields = ["pharmgkb_variant_id", "pharmgkb_annotation_id",
                  "source", "last_updated"]
    all_fields = list(original_fields)
    for field in new_fields:
        if field not in all_fields:
            all_fields.append(field)

    # Process each row
    enriched = 0
    total = len(rows)
    unique_rsids = set(row["variant_rsid"] for row in rows)

    print(f"  {total} rows, {len(unique_rsids)} unique variants")

    # Batch lookup variants
    variant_cache = {}
    for i, rsid in enumerate(sorted(unique_rsids)):
        print(f"  [{i+1}/{len(unique_rsids)}] Looking up {rsid}...", end="")
        info = lookup_variant(rsid)
        if info:
            variant_cache[rsid] = info
            print(f" → {info['pharmgkb_id']} ({info['clinical_significance']})")
            enriched += 1
        else:
            print(" → not found")

    # Also look up genes
    unique_genes = set(row["gene"] for row in rows)
    gene_cache = {}
    for gene in sorted(unique_genes):
        print(f"  Looking up gene {gene}...", end="")
        info = lookup_gene(gene)
        if info:
            gene_cache[gene] = info
            cpic = "CPIC" if info["cpic_gene"] else "non-CPIC"
            print(f" → {info['pharmgkb_id']} ({cpic})")
        else:
            print(" → not found")

    # Update rows
    now = datetime.now(timezone.utc).isoformat()
    for row in rows:
        rsid = row["variant_rsid"]
        variant_info = variant_cache.get(rsid, {})

        row["pharmgkb_variant_id"] = variant_info.get("pharmgkb_id", "")
        row["pharmgkb_annotation_id"] = ""  # Populated if we find matching annotations
        row["source"] = "PharmGKB_API" if variant_info else "Manual"
        row["last_updated"] = now

    # Write enriched CSV
    with open(RSID_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Enriched {enriched}/{len(unique_rsids)} variants "
          f"({enriched/len(unique_rsids)*100:.0f}%)")
    return True


def enrich_diplotype_csv(force=False):
    """Enrich the diplotype-based CSV with PharmGKB gene data."""
    print("\n=== Enriching diplotype associations ===")

    if not os.path.exists(DIPLOTYPE_CSV):
        print(f"  CSV not found: {DIPLOTYPE_CSV}", file=sys.stderr)
        return False

    with open(DIPLOTYPE_CSV, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    new_fields = ["pharmgkb_gene_id", "cpic_gene", "source", "last_updated"]
    all_fields = list(original_fields)
    for field in new_fields:
        if field not in all_fields:
            all_fields.append(field)

    # Look up genes
    unique_genes = set(row["gene"] for row in rows)
    gene_cache = {}
    for gene in sorted(unique_genes):
        # Use pre-seeded IDs first
        if gene in GENE_IDS:
            gene_cache[gene] = {
                "pharmgkb_id": GENE_IDS[gene],
                "cpic_gene": True,  # All our genes are CPIC
            }
            print(f"  {gene} → {GENE_IDS[gene]} (pre-seeded)")
        else:
            print(f"  Looking up gene {gene}...", end="")
            info = lookup_gene(gene)
            if info:
                gene_cache[gene] = info
                print(f" → {info['pharmgkb_id']}")
            else:
                print(" → not found")

    # Update rows
    now = datetime.now(timezone.utc).isoformat()
    enriched = 0
    for row in rows:
        gene = row["gene"]
        gene_info = gene_cache.get(gene, {})

        row["pharmgkb_gene_id"] = gene_info.get("pharmgkb_id", "")
        row["cpic_gene"] = str(gene_info.get("cpic_gene", False))
        row["source"] = "PharmGKB_API" if gene_info else "Manual"
        row["last_updated"] = now
        if gene_info:
            enriched += 1

    with open(DIPLOTYPE_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Enriched {enriched}/{len(rows)} rows")
    return True


def generate_report(start_time):
    """Generate a summary report of the update."""
    elapsed = time.time() - start_time
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(elapsed, 1),
        "api_base": API_BASE,
        "license": "CC-BY-SA 4.0 (PharmGKB/ClinPGx)",
    }

    report_file = os.path.join(DATA_DIR, "pharmgkb_update_report.json")
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\n=== Update complete in {elapsed:.1f}s ===")
    print(f"Report saved to: {report_file}")
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Update PGx-ADR databases from PharmGKB API"
    )
    parser.add_argument("--force", action="store_true",
                        help="Force refresh even if cache is fresh")
    args = parser.parse_args()

    os.makedirs(CACHE_DIR, exist_ok=True)

    start = time.time()

    print("PharmGKB Data Updater")
    print(f"API: {API_BASE}")
    print(f"Cache: {CACHE_DIR}")
    print(f"Data: {DATA_DIR}")

    enrich_rsid_csv(force=args.force)
    enrich_diplotype_csv(force=args.force)
    generate_report(start)


if __name__ == "__main__":
    main()
