#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/.." && pwd)}"; export GMO_ROOT
# ============================================================================
# GMO loci pipeline — driver.   usage:  run_pipeline.sh <stage> [round] [options]
#
#   registry              01  build loci registry from all_blast_stats.tsv          (login node, seconds)
#   fetch [--only FILE]   02  download contigs + assembly metadata                  (login node, NCBI-bound, ~1-2 h for 2000 contigs)
#   guides                03  build T1/T2(/T3) guide proteins                       (login node, minutes)
#   annotate <round> [--only FILE|--local]
#                         04  miniprot + AUGUSTUS per locus, LSF job array          (cluster, ~minutes wall time)
#   validate <round>      05  family HMMs, hmmsearch, status per locus              (login node, minutes)
#   diamond <round>       06a collect neighbours/targets -> diamond vs nr on LSF    (cluster, 20-60 min, one job)
#   classify <round>      06b parse diamond, taxonomy, per-locus context            (login node, NCBI taxonomy lookups)
#   summarize <round>     07  loci_table.tsv, per-family models, T3 guides          (login node)
#   all <round>           runs 01..07 sequentially, waiting for LSF jobs (-K)
#
# Round 2 = same registry/contigs; rerun 'guides' (picks up T3 from round 1), then annotate/validate/diamond/classify/summarize 2.
# Set NCBI_API_KEY in the environment for 10 req/s at NCBI.
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GMO_CONFIG="${GMO_CONFIG:-$HERE/config.yaml}"
S="$HERE/scripts"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$(python3 -c "import yaml;print(yaml.safe_load(open('$GMO_CONFIG'))['conda_env'])" 2>/dev/null || echo annot)"
cfg() { python3 - "$1" <<'EOF'
import sys,yaml,os; c=yaml.safe_load(open(os.environ['GMO_CONFIG'])); v=c
for k in sys.argv[1].split('.'): v=v[k]
print(os.path.expandvars(str(v)))
EOF
}
OUT="$(cfg project_dir)/$(cfg out_dir)"; mkdir -p "$OUT/logs"
Q=$(cfg lsf.queue)
stage="${1:-}"; shift || true

case "$stage" in
  registry) python3 "$S/01_build_registry.py" ;;
  fetch)    python3 "$S/02_fetch_contigs.py" "$@" ;;
  guides)   python3 "$S/03_build_guides.py" ;;
  annotate)
    R="$1"; shift; RD="$OUT/round$R"; mkdir -p "$RD/loci"
    LIST="$RD/loci_list.txt"
    if [[ "${1:-}" == "--only" ]]; then cp "$2" "$LIST"; shift 2
    else python3 - <<EOF > "$LIST"
import csv,os
c={r['contig'] for r in csv.DictReader(open('$OUT/contigs/contigs.tsv'),delimiter='\t')}
for r in csv.DictReader(open('$OUT/registry/loci.tsv'),delimiter='\t'):
    if r['contig'] in c: print(r['locus_id'])
EOF
    fi
    N=$(wc -l < "$LIST"); echo "loci to annotate: $N (list: $LIST)"
    if [[ "${1:-}" == "--local" ]]; then          # small tests on the login node
      for i in $(seq 1 "$N"); do python3 "$S/04_annotate_locus.py" "$R" "$i"; done
    else
      # LSF job arrays do not accept -K and are capped by MAX_JOB_ARRAY_SIZE (2000 here): submit chunks with an index offset,
      # then wait for every chunk with bwait (fallback: poll bjobs).
      rm -f "$OUT"/logs/annot_r$R.*.out "$OUT"/logs/annot_r$R.*.err
      MAXA=$(bparams -a 2>/dev/null | awk '/MAX_JOB_ARRAY_SIZE/{print $3}'); MAXA=${MAXA:-1000}
      JOBIDS=(); OFF=0
      while (( OFF < N )); do
        CH=$(( N - OFF )); (( CH > MAXA )) && CH=$MAXA
        JID=$(bsub -q "$Q" -J "annot_r${R}_$OFF[1-$CH]%$(cfg lsf.annotate_concurrency)" -n 2 -R "span[hosts=1]" -R "rusage[mem=$(cfg lsf.mem_per_slot_mb)]" \
             -o "$OUT/logs/annot_r$R.%J.%I.out" -e "$OUT/logs/annot_r$R.%J.%I.err" \
             "source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot; export GMO_CONFIG=$GMO_CONFIG; python3 $S/04_annotate_locus.py $R \$(( $OFF + \$LSB_JOBINDEX ))" \
             | grep -oP 'Job <\K[0-9]+')
        echo "submitted array $JID: indices $((OFF+1))-$((OFF+CH))"; JOBIDS+=("$JID"); OFF=$(( OFF + CH ))
      done
      sleep 5
      for JID in "${JOBIDS[@]}"; do
        if ! bwait -w "ended($JID)" 2>/dev/null; then
          seen=0
          while true; do
            out=$(bjobs -noheader "$JID" 2>/dev/null || true)
            if [[ -z "$out" ]]; then (( seen )) && break || { sleep 10; continue; }; fi
            seen=1; echo "$out" | grep -q -E 'PEND|RUN' || break; sleep 30
          done
        fi
      done
      JOBID="*"
      cat "$OUT"/logs/annot_r$R.*.out 2>/dev/null | grep -P '^\S+\tbest=' > "$RD/annotate_summary.tsv" || true
      echo "annotated: $(wc -l < "$RD/annotate_summary.tsv") ; failures: $(grep -l -E 'Traceback|Error' "$OUT"/logs/annot_r$R.*.err 2>/dev/null | wc -l)"
    fi ;;
  validate) python3 "$S/05_validate_models.py" "$1" ;;
  diamond)
    R="$1"; RD="$OUT/round$R"
    python3 "$S/06_classify_neighbors.py" "$R" collect
    T=$(cfg diamond.threads); B=$(cfg diamond.block_size); C=$(cfg diamond.index_chunks); SENS=$(cfg diamond.sensitivity)
    bsub -K -q "$Q" -J dmnd_r$R -n "$T" -R "span[hosts=1]" -R "rusage[mem=$(cfg lsf.diamond_mem_per_slot_mb)]" \
         -o "$OUT/logs/diamond_r$R.%J.out" -e "$OUT/logs/diamond_r$R.%J.err" \
         "source $(conda info --base)/etc/profile.d/conda.sh; conda activate annot; diamond blastp -q $RD/diamond_query.faa -d $(cfg diamond.db) -o $RD/diamond_nr.tsv --threads $T -b $B -c $C $SENS --max-target-seqs $(cfg diamond.max_target_seqs) --evalue $(cfg diamond.evalue) --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore stitle"
    echo "diamond hits: $(wc -l < "$RD/diamond_nr.tsv")" ;;
  classify)  python3 "$S/06_classify_neighbors.py" "$1" parse ;;
  summarize) python3 "$S/07_summarize.py" "$1" ;;
  all)
    R="$1"
    "$0" registry; "$0" fetch; "$0" guides; "$0" annotate "$R"; "$0" validate "$R"; "$0" diamond "$R"; "$0" classify "$R"; "$0" summarize "$R" ;;
  *) sed -n '2,20p' "$0"; exit 1 ;;
esac
