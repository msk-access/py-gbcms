/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE_SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Collect per-sample MAF filtering statistics and produce a run summary.

    Input:  All *.filter_stats.tsv files from FILTER_MAF (collected)
    Output: pipeline_summary.tsv with one row per sample

    Runs only when filter_by_sample is enabled.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PIPELINE_SUMMARY {
    label 'process_low'

    publishDir "${params.outdir}/gbcms", mode: params.publish_dir_mode

    input:
    path stats_files  // collected *.filter_stats.tsv files from FILTER_MAF

    output:
    path "pipeline_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env python3
    \"\"\"Aggregate per-sample filtering stats into a pipeline summary.\"\"\"
    import glob
    import sys

    stats = sorted(glob.glob("*.filter_stats.tsv"))
    rows = []
    for f in stats:
        with open(f) as fh:
            next(fh)  # skip header
            for line in fh:
                parts = line.strip().split("\\t")
                if len(parts) >= 7:
                    rows.append(parts)

    # Write consolidated TSV
    with open("pipeline_summary.tsv", "w") as out:
        out.write(
            "sample\\tstatus\\tvariants_matched\\tvariants_total"
            "\\tduplicates_removed\\tmatched_tsbs\\tmode\\n"
        )
        for r in sorted(rows, key=lambda x: x[0]):
            out.write("\\t".join(r) + "\\n")

    # Print formatted table to stdout (visible in Nextflow log)
    print()
    print("=" * 90)
    print("GBCMS PIPELINE SUMMARY — MAF Filtering")
    print("=" * 90)
    print(
        f"{'Sample':<30} {'Status':<10} {'Variants':<12} {'Dups':<6} {'Matched TSBs'}"
    )
    print("-" * 90)
    for r in sorted(rows, key=lambda x: x[0]):
        sample, status, matched, total, dups, tsbs, mode = r[:7]
        print(f"{sample:<30} {status:<10} {matched}/{total:<8} {dups:<6} {tsbs}")
    print("-" * 90)
    processed = sum(1 for r in rows if r[1] == "processed")
    skipped_count = sum(1 for r in rows if r[1] == "skipped")
    print(
        f"Total: {len(rows)} sample(s), {processed} processed, {skipped_count} skipped"
    )
    print("=" * 90)
    print()
    """
}
