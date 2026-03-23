#!/usr/bin/env python3
"""
PharmCAT Wrapper for PGx-ADR Calculator
=========================================

Runs the PharmCAT pipeline on a VCF file and parses the output into
a standardized TSV format for consumption by the R Shiny app.

Usage:
    python3 pharmcat_wrapper.py <input_vcf> <output_tsv> [--pharmcat-dir <path>]

Output TSV columns:
    gene, diplotype, phenotype, activity_score, lookup_key, variants, source
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import glob


# Default PharmCAT installation directory (relative to this script)
DEFAULT_PHARMCAT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "pharmcat"
)

# Java 17 path (Homebrew)
JAVA_HOME = "/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home"


def find_pharmcat(pharmcat_dir):
    """Locate PharmCAT pipeline executable and JAR."""
    pipeline = os.path.join(pharmcat_dir, "pharmcat_pipeline")
    jar = os.path.join(pharmcat_dir, "pharmcat.jar")
    venv_python = os.path.join(pharmcat_dir, ".venv", "bin", "python3")

    if not os.path.exists(jar):
        raise FileNotFoundError(f"PharmCAT JAR not found at {jar}")
    if not os.path.exists(pipeline):
        raise FileNotFoundError(f"PharmCAT pipeline not found at {pipeline}")

    return pipeline, jar, venv_python


def run_pharmcat(input_vcf, output_dir, pharmcat_dir):
    """Run PharmCAT pipeline on a VCF file."""
    pipeline, jar, venv_python = find_pharmcat(pharmcat_dir)

    env = os.environ.copy()
    env["JAVA_HOME"] = JAVA_HOME
    env["PATH"] = f"{JAVA_HOME}/bin:{env.get('PATH', '')}"

    # Activate venv if it exists
    if os.path.exists(venv_python):
        env["PATH"] = f"{os.path.dirname(venv_python)}:{env['PATH']}"
        env["VIRTUAL_ENV"] = os.path.join(pharmcat_dir, ".venv")

    cmd = [
        pipeline,
        input_vcf,
        "-o", output_dir,
    ]

    print(f"Running PharmCAT: {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env,
        cwd=pharmcat_dir,
        timeout=300  # 5 minute timeout
    )

    if result.returncode != 0:
        print(f"PharmCAT STDERR:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"PharmCAT failed with exit code {result.returncode}")

    print(f"PharmCAT completed successfully", file=sys.stderr)
    return result


def parse_pharmcat_json(json_path):
    """Parse PharmCAT JSON output into structured records."""
    with open(json_path, 'r') as f:
        data = json.load(f)

    records = []

    # Extract gene calls from the "geneResults" section
    gene_results = data.get("results", [])

    for result in gene_results:
        gene = result.get("gene", "")
        diplotypes = result.get("diplotypes", [])

        if not diplotypes:
            continue

        for diplo in diplotypes:
            diplotype_str = diplo.get("name", "")
            phenotype = diplo.get("phenotype", "")
            activity_score = diplo.get("activityScore", "")

            # Build a lookup key for the risk database
            # Use phenotype-based lookup (e.g., "CYP2D6:Poor Metabolizer")
            lookup_key = f"{gene}:{phenotype}" if phenotype else f"{gene}:{diplotype_str}"

            # Get variant details
            variant_list = []
            for allele in [diplo.get("allele1"), diplo.get("allele2")]:
                if allele and isinstance(allele, dict):
                    for var in allele.get("variants", []):
                        rsid = var.get("rsid", "")
                        if rsid:
                            variant_list.append(rsid)

            records.append({
                "gene": gene,
                "diplotype": diplotype_str,
                "phenotype": phenotype if phenotype else "Unknown",
                "activity_score": str(activity_score) if activity_score else "N/A",
                "lookup_key": lookup_key,
                "variants": ";".join(variant_list) if variant_list else "N/A",
                "source": "PharmCAT"
            })

    return records


def find_json_output(output_dir, input_vcf):
    """Find the PharmCAT JSON output file."""
    basename = os.path.splitext(os.path.basename(input_vcf))[0]

    # PharmCAT outputs files like: <basename>.pharmcat.json
    # or in preprocessed mode: <basename>.preprocessed.pharmcat.json
    patterns = [
        os.path.join(output_dir, f"*pharmcat.json"),
        os.path.join(output_dir, f"{basename}*.json"),
        os.path.join(output_dir, "*.json"),
    ]

    for pattern in patterns:
        matches = glob.glob(pattern)
        if matches:
            return matches[0]

    raise FileNotFoundError(
        f"No PharmCAT JSON output found in {output_dir}. "
        f"Files present: {os.listdir(output_dir)}"
    )


def write_tsv(records, output_path):
    """Write parsed records to TSV."""
    columns = ["gene", "diplotype", "phenotype", "activity_score",
                "lookup_key", "variants", "source"]

    with open(output_path, 'w') as f:
        f.write("\t".join(columns) + "\n")
        for rec in records:
            f.write("\t".join(rec.get(c, "") for c in columns) + "\n")

    print(f"Wrote {len(records)} gene calls to {output_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Run PharmCAT and parse output for PGx-ADR Calculator"
    )
    parser.add_argument("input_vcf", help="Input VCF file (GRCh38-aligned)")
    parser.add_argument("output_tsv", help="Output TSV file path")
    parser.add_argument(
        "--pharmcat-dir",
        default=DEFAULT_PHARMCAT_DIR,
        help=f"PharmCAT installation directory (default: {DEFAULT_PHARMCAT_DIR})"
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_vcf):
        print(f"Error: Input VCF not found: {args.input_vcf}", file=sys.stderr)
        sys.exit(1)

    # Create temp output directory for PharmCAT
    with tempfile.TemporaryDirectory(prefix="pharmcat_") as temp_dir:
        # Run PharmCAT
        run_pharmcat(args.input_vcf, temp_dir, args.pharmcat_dir)

        # Find and parse the JSON output
        json_path = find_json_output(temp_dir, args.input_vcf)
        print(f"Parsing PharmCAT output: {json_path}", file=sys.stderr)

        records = parse_pharmcat_json(json_path)

        if not records:
            print("Warning: No gene calls found in PharmCAT output", file=sys.stderr)

        # Write TSV
        write_tsv(records, args.output_tsv)

    print("PharmCAT wrapper completed successfully", file=sys.stderr)


if __name__ == "__main__":
    main()
