#!/bin/bash
#SBATCH --job-name=kh2013_map
#SBATCH --output=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/kh2013_map_%j.out
#SBATCH --error=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/kh2013_map_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dtjaernberg+unity@umass.edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Add KH2013 gene IDs to the KY21 dictionary.
#
#   ssh unity 'cd /work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/code && sbatch --parsable kh2013_map.sbatch'
#
# Bulk lives in the crispytargets scratch workspace: staged inputs, the BLAST
# database and the BLAST table. Only results and provenance stay under /work.
# BLAST direction matches the other species columns: KY21 proteins are the
# query, the species proteome is the database. The mapper runs on a copy of
# the existing dictionary so the incremental backbone keeps every column.

set -eo pipefail

set +u
source /etc/profile
source /etc/profile.d/*.sh 2>/dev/null || true
module purge
module load blast-plus/2.16.0
module load conda/latest
eval "$(conda shell.bash hook)"
conda activate mycrispytargets

PROJ=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets
RES="$PROJ/results/kh2013_mapping"
QA="$HOME/projects/src/crispytargets/tools/qa_kh2013_mapping.py"

# Workspace resolved by name at run time, never hardcoded.
SCRATCH=$(ws_list -v | awk '/^id: crispytargets$/{found=1; next} \
    found && /workspace directory/{print $NF; exit}')
if [ -z "$SCRATCH" ]; then
    echo "no crispytargets workspace; run: ws_allocate crispytargets 14"
    exit 1
fi
DATA="$SCRATCH/kh2013"
INPUTS="$DATA/inputs"
TMP="$DATA/work"

KY21_PROT="$INPUTS/HT.KY21Gene.protein.2.fasta"
KY21_GFF="$INPUTS/HT.KY21Gene.2.gff3"
KH_PROT="$INPUTS/KH.KHGene.2013.Longest.protein.clean.fasta"
DICT="$INPUTS/KY21_dictionary.tsv"
GHOST="$INPUTS/KHtoKY21.txt.gz"

echo "=== kh2013_map on $(hostname) at $(date -Is)"
echo "job ${SLURM_JOB_ID}, ${SLURM_CPUS_PER_TASK} cpus"
echo "scratch bundle: $DATA"
blastp -version
python -V

cd "$INPUTS"
sha256sum -c PROVENANCE.sha256
for input in "$KY21_PROT" "$KY21_GFF" "$KH_PROT" "$DICT" "$GHOST" "$QA"; do
    test -s "$input" || { echo "missing input: $input"; exit 1; }
done

mkdir -p "$RES" "$TMP"

echo "=== makeblastdb"
makeblastdb -in "$KH_PROT" -dbtype prot -out "$TMP/kh2013_db"

echo "=== blastp"
blastp -query "$KY21_PROT" \
       -db "$TMP/kh2013_db" \
       -out "$TMP/kh2013_vs_ky21.tsv" \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 20 \
       -num_threads "$SLURM_CPUS_PER_TASK"

echo "=== map (incremental on a copy of the dictionary)"
cp "$DICT" "$RES/KY21_dictionary.extended.tsv"
crispytargets map \
    --gff "$KY21_GFF" \
    --blast-file "$TMP/kh2013_vs_ky21.tsv" \
    --species KH2013 \
    --evalue 1e-50 \
    --output "$RES/KY21_dictionary.extended.tsv"

echo "=== clean-mapping"
crispytargets clean-mapping \
    -i "$RES/KY21_dictionary.extended.tsv" \
    -o "$RES/KY21_dictionary.extended.cleaned.tsv" \
    -c KH2013_primary,KH2013_secondary
mv "$RES/KY21_dictionary.extended.cleaned.tsv" "$RES/KY21_dictionary.extended.tsv"

echo "=== qa"
python "$QA" \
    --base "$DICT" \
    --extended "$RES/KY21_dictionary.extended.tsv" \
    --ghost "$GHOST" \
    --projection "$RES/KY21_to_KH2013.tsv" \
    --report "$RES/qa_report.tsv"

echo "=== pointer for /work"
cat > "$RES/SCRATCH_BUNDLE.txt" <<EOF
Bulk for this run lives in the crispytargets scratch workspace:
  $DATA

inputs/  staged inputs and PROVENANCE (sha256sum -c PROVENANCE.sha256)
work/    kh2013_db.* and kh2013_vs_ky21.tsv, regenerable by this sbatch

Scratch has no snapshots. Everything here is either re-downloadable or
reproducible from tools/kh2013_map.sbatch.
EOF

echo "=== checksums"
(cd "$INPUTS" && sha256sum "$(basename "$KY21_PROT")" "$(basename "$KY21_GFF")" \
    "$(basename "$KH_PROT")" "$(basename "$DICT")" "$(basename "$GHOST")")
(cd "$RES" && sha256sum KY21_dictionary.extended.tsv KY21_to_KH2013.tsv qa_report.tsv)

echo "=== done at $(date -Is)"
