#!/usr/bin/env bash
# =============================================================================
# ONT mitochondrial genome pipeline
#   1. FASTQ QC (length + mean Q)
#   2. Flye assembly (optional replicates, pick best vs reference)
#   3. Racon polish (2 rounds)
#   4. Medaka polish
# =============================================================================
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: mt_pipeline.sh [options]

Required:
  -i, --input FASTQ          Raw ONT FASTQ
  -o, --outdir DIR           Output directory

Optional:
  -r, --reference FASTA      Reference for picking best Flye replicate
  -t, --threads INT          Threads (default: 16)
  -g, --genome-size STR      Flye genome size (default: 15k)
  --min-length INT           Min read length (default: 1000)
  --min-mean-q FLOAT         Min mean Q (default: 7)
  --flye-runs INT            Flye replicates (default: 1)
  --racon-rounds INT         Racon rounds (default: 2)
  --medaka-model STR         Medaka model (default: r1041_e82_400bps_sup_v5.0.0)
  --skip-qc                  Skip filtering if filtered FASTQ exists
  --skip-flye                Skip Flye (use existing draft)
  --skip-racon               Skip Racon
  --skip-medaka              Skip Medaka
  -h, --help                 Show help

Environment (override tool paths if needed):
  FLYE, MINIMAP2, RACON, MEDAKA_CONSENSUS, PYTHON, MUMMER

Example:
  ./mt_pipeline.sh \\
    -i reads.fastq \\
    -o results \\
    -r PX108868.fasta \\
    --flye-runs 6 \\
    -t 16
EOF
}

# Defaults
INPUT=""
OUTDIR=""
REFERENCE=""
THREADS=16
GENOME_SIZE="15k"
MIN_LENGTH=1000
MIN_MEAN_Q=7
FLYE_RUNS=1
RACON_ROUNDS=2
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"
SKIP_QC=0
SKIP_FLYE=0
SKIP_RACON=0
SKIP_MEDAKA=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
FLYE="${FLYE:-flye}"
MINIMAP2="${MINIMAP2:-minimap2}"
RACON="${RACON:-racon}"
MEDAKA_CONSENSUS="${MEDAKA_CONSENSUS:-medaka_consensus}"
MUMMER="${MUMMER:-mummer}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -r|--reference) REFERENCE="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -g|--genome-size) GENOME_SIZE="$2"; shift 2 ;;
    --min-length) MIN_LENGTH="$2"; shift 2 ;;
    --min-mean-q) MIN_MEAN_Q="$2"; shift 2 ;;
    --flye-runs) FLYE_RUNS="$2"; shift 2 ;;
    --racon-rounds) RACON_ROUNDS="$2"; shift 2 ;;
    --medaka-model) MEDAKA_MODEL="$2"; shift 2 ;;
    --skip-qc) SKIP_QC=1; shift ;;
    --skip-flye) SKIP_FLYE=1; shift ;;
    --skip-racon) SKIP_RACON=1; shift ;;
    --skip-medaka) SKIP_MEDAKA=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT" && -n "$OUTDIR" ]] || { usage; exit 1; }
[[ -f "$INPUT" ]] || { echo "Input FASTQ not found: $INPUT" >&2; exit 1; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"
LOG="$OUTDIR/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "ONT mt assembly pipeline"
echo "Started: $(date)"
echo "Input:   $INPUT"
echo "Outdir:  $OUTDIR"
echo "============================================================"

# -----------------------------------------------------------------------------
# Step 1: QC
# -----------------------------------------------------------------------------
FILTERED="$OUTDIR/filter_q${MIN_MEAN_Q}_${MIN_LENGTH}.fastq"
if [[ "$SKIP_QC" -eq 0 ]]; then
  if [[ -s "$FILTERED" ]]; then
    echo "[QC] Reusing existing $FILTERED"
  else
    echo "[QC] Filtering reads (length>=${MIN_LENGTH}, mean Q>=${MIN_MEAN_Q}) ..."
    "$PYTHON" "$SCRIPT_DIR/filter_fastq.py" \
      -i "$INPUT" \
      -o "$FILTERED" \
      -l "$MIN_LENGTH" \
      -q "$MIN_MEAN_Q"
  fi
else
  FILTERED="$INPUT"
fi
READS="$FILTERED"

# -----------------------------------------------------------------------------
# Step 2: Flye
# -----------------------------------------------------------------------------
DRAFT="$OUTDIR/draft.fasta"
if [[ "$SKIP_FLYE" -eq 0 ]]; then
  echo "[Flye] Running $FLYE_RUNS replicate(s) ..."
  SCORES="$OUTDIR/flye_scores.tsv"
  echo -e "run\tassembly\tindels\ttotal_variants" > "$SCORES"
  BEST_INDELS=999999
  BEST_RUN=""

  for i in $(seq 1 "$FLYE_RUNS"); do
    RUN_DIR="$OUTDIR/flye_run_${i}"
    ASM="$RUN_DIR/assembly.fasta"
    if [[ ! -s "$ASM" ]]; then
      echo "  [Flye] run_${i} ..."
      "$FLYE" \
        --nano-raw "$READS" \
        --out-dir "$RUN_DIR" \
        --threads "$THREADS" \
        --genome-size "$GENOME_SIZE" \
        --iterations 1 \
        --min-overlap "$MIN_LENGTH" \
        --meta
    else
      echo "  [Flye] run_${i} exists, skip"
    fi

    if [[ -n "$REFERENCE" && -f "$REFERENCE" ]]; then
      PREFIX="$OUTDIR/stats/flye_run_${i}"
      mkdir -p "$OUTDIR/stats"
      SCORE=$("$PYTHON" "$SCRIPT_DIR/score_vs_ref.py" \
        --ref "$REFERENCE" \
        --assembly "$ASM" \
        --mummer "$MUMMER" \
        --prefix "$PREFIX" | tail -1)
      INDELS=$(echo "$SCORE" | cut -f4)
      TOTAL=$(echo "$SCORE" | cut -f8)
      echo -e "run_${i}\t$ASM\t$INDELS\t$TOTAL" >> "$SCORES"
      if [[ "$INDELS" -lt "$BEST_INDELS" ]]; then
        BEST_INDELS=$INDELS
        BEST_RUN=$i
        cp "$ASM" "$DRAFT"
      fi
    elif [[ "$FLYE_RUNS" -eq 1 ]]; then
      cp "$ASM" "$DRAFT"
      BEST_RUN=1
    fi
  done

  if [[ -n "$REFERENCE" && -n "$BEST_RUN" ]]; then
    echo "[Flye] Best replicate: run_${BEST_RUN} (indels=$BEST_INDELS vs reference)"
    echo "run_${BEST_RUN}" > "$OUTDIR/best_flye_run.txt"
  elif [[ "$FLYE_RUNS" -gt 1 && -z "$REFERENCE" ]]; then
    echo "[Flye] Warning: multiple runs but no --reference; using run_1" >&2
    cp "$OUTDIR/flye_run_1/assembly.fasta" "$DRAFT"
  fi
else
  [[ -s "$DRAFT" ]] || { echo "draft.fasta not found; remove --skip-flye" >&2; exit 1; }
fi

# Reformat draft for medaka (60-char lines)
FMT_DRAFT="$OUTDIR/draft_fmt.fasta"
"$PYTHON" - <<PY
from pathlib import Path
src = Path("$DRAFT")
dst = Path("$FMT_DRAFT")
lines = src.read_text().splitlines()
hdr = lines[0]
seq = "".join(l.strip() for l in lines[1:])
with dst.open("w") as f:
    f.write(hdr + "\\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\\n")
print(f"Formatted draft: {len(seq)} bp")
PY
DRAFT="$FMT_DRAFT"

# -----------------------------------------------------------------------------
# Step 3: Racon
# -----------------------------------------------------------------------------
POLISH_IN="$DRAFT"
if [[ "$SKIP_RACON" -eq 0 ]]; then
  echo "[Racon] $RACON_ROUNDS round(s) ..."
  for r in $(seq 1 "$RACON_ROUNDS"); do
    OUT_FA="$OUTDIR/racon${r}.fasta"
    SAM="$OUTDIR/racon${r}.sam"
    if [[ -s "$OUT_FA" ]]; then
      echo "  [Racon] round ${r} exists, skip"
    else
      echo "  [Racon] round ${r}: mapping ..."
      "$MINIMAP2" -ax map-ont -t "$THREADS" "$POLISH_IN" "$READS" > "$SAM"
      echo "  [Racon] round ${r}: polishing ..."
      "$RACON" -t "$THREADS" "$READS" "$SAM" "$POLISH_IN" > "$OUT_FA"
    fi
    POLISH_IN="$OUT_FA"
  done
  cp "$POLISH_IN" "$OUTDIR/racon.fasta"
else
  POLISH_IN="${POLISH_IN:-$DRAFT}"
fi

# -----------------------------------------------------------------------------
# Step 4: Medaka
# -----------------------------------------------------------------------------
if [[ "$SKIP_MEDAKA" -eq 0 ]]; then
  echo "[Medaka] Polishing ..."
  MEDAKA_DIR="$OUTDIR/medaka"
  if [[ -s "$MEDAKA_DIR/consensus.fasta" ]]; then
    echo "  [Medaka] exists, skip"
  else
    "$MEDAKA_CONSENSUS" \
      -i "$READS" \
      -d "$POLISH_IN" \
      -o "$MEDAKA_DIR" \
      -m "$MEDAKA_MODEL" \
      -t "$THREADS"
  fi
  cp "$MEDAKA_DIR/consensus.fasta" "$OUTDIR/medaka.fasta"
  FINAL="$OUTDIR/medaka.fasta"
else
  FINAL="$OUTDIR/racon.fasta"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "============================================================"
echo "Pipeline finished: $(date)"
echo "Filtered reads:  $READS"
echo "Draft:           $DRAFT"
[[ -f "$OUTDIR/racon.fasta" ]] && echo "Racon:           $OUTDIR/racon.fasta"
[[ -f "$OUTDIR/medaka.fasta" ]] && echo "Medaka:          $OUTDIR/medaka.fasta"
echo "Final assembly:  $FINAL"
if [[ -n "$REFERENCE" && -f "$REFERENCE" ]]; then
  echo "[Stats] Final assembly vs reference:"
  "$PYTHON" "$SCRIPT_DIR/score_vs_ref.py" \
    --ref "$REFERENCE" \
    --assembly "$FINAL" \
    --mummer "$MUMMER" \
    --prefix "$OUTDIR/stats/final_vs_ref"
fi
echo "Log: $LOG"
echo "============================================================"
