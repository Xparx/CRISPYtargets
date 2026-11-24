#!/bin/bash
#SBATCH --job-name=map_species
#SBATCH --output=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/map_species_%j.out
#SBATCH --error=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/map_species_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dtjaernberg+unity@umass.edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Add one species' gene IDs to the KY21 dictionary.
#
#   ssh unity 'cd /work/.../projects/crispytargets/code && \
#       SPECIES=ENSEMBL PROTEOME=... RESULT_NAME=ensembl_mapping sbatch map_species.sbatch'
#
# Defaults reproduce the KH2013 run (job 1792936). Bulk lives in the
# crispytargets scratch workspace under $SCRATCH/<species lower case>: staged
# inputs, the BLAST database and the BLAST table. Only results and provenance
# stay under /work. BLAST direction matches the other species columns: KY21
# proteins are the query, the species proteome is the database. The mapper
# runs on a copy of the existing dictionary so the incremental backbone keeps
# every column.
#
# Environment:
#   SPECIES      column label, e.g. KH2013 or ENSEMBL
#   PROTEOME     staged proteome filename in $INPUTS
#   SUBJECT_MAP  optional staged TSV (protein_id, gene_id) for proteomes
#                whose BLAST subjects are not gene identifiers
#   REFERENCE    staged two-column TSV (gene_id, SPECIES) for the QA
#   RESULT_NAME  results subdirectory under /work
#   PROJECTION   standalone projection filename

set -eo pipefail

set +u
source /etc/profile
source /etc/profile.d/*.sh 2>/dev/null || true
module purge
module load blast-plus/2.16.0
module load conda/latest
eval "$(conda shell.bash hook)"
conda activate mycrispytargets

SPECIES=${SPECIES:-KH2013}
PROTEOME=${PROTEOME:-KH.KHGene.2013.Longest.protein.clean.fasta}
SUBJECT_MAP=${SUBJECT_MAP:-}
REFERENCE=${REFERENCE:-KH2013.reference.tsv}
RESULT_NAME=${RESULT_NAME:-kh2013_mapping}
PROJECTION=${PROJECTION:-KY21_to_KH2013.tsv}
SP_LOWER=$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')

PROJ=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets
RES="$PROJ/results/$RESULT_NAME"
QA="$HOME/projects/src/crispytargets/tools/qa_species_mapping.py"

# Workspace resolved by name at run time, never hardcoded.
SCRATCH=$(ws_list -v | awk '/^id: crispytargets$/{found=1; next} \
    found && /workspace directory/{print $NF; exit}')
if [ -z "$SCRATCH" ]; then
    echo "no crispytargets workspace; run: ws_allocate crispytargets 14"
    exit 1
fi
DATA="$SCRATCH/$SP_LOWER"
INPUTS="$DATA/inputs"
TMP="$DATA/work"

KY21_PROT="$INPUTS/HT.KY21Gene.protein.2.fasta"
KY21_GFF="$INPUTS/HT.KY21Gene.2.gff3"
PROT="$INPUTS/$PROTEOME"
DICT="$INPUTS/KY21_dictionary.tsv"
REF="$INPUTS/$REFERENCE"

echo "=== map_species $SPECIES on $(hostname) at $(date -Is)"
echo "job ${SLURM_JOB_ID}, ${SLURM_CPUS_PER_TASK} cpus"
echo "scratch bundle: $DATA"
blastp -version
python -V

cd "$INPUTS"
sha256sum -c PROVENANCE.sha256
for input in "$KY21_PROT" "$KY21_GFF" "$PROT" "$DICT" "$REF" "$QA"; do
    test -s "$input" || { echo "missing input: $input"; exit 1; }
done
if [ -n "$SUBJECT_MAP" ]; then
    test -s "$INPUTS/$SUBJECT_MAP" || { echo "missing subject map: $SUBJECT_MAP"; exit 1; }
fi

mkdir -p "$RES" "$TMP"

echo "=== makeblastdb"
makeblastdb -in "$PROT" -dbtype prot -out "$TMP/${SP_LOWER}_db"

echo "=== blastp"
blastp -query "$KY21_PROT" \
       -db "$TMP/${SP_LOWER}_db" \
       -out "$TMP/${SP_LOWER}_vs_ky21.tsv" \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 20 \
       -num_threads "$SLURM_CPUS_PER_TASK"

MAP_ARGS=()
if [ -n "$SUBJECT_MAP" ]; then
    MAP_ARGS+=(--subject-map "$INPUTS/$SUBJECT_MAP")
fi

echo "=== map (incremental on a copy of the dictionary)"
cp "$DICT" "$RES/KY21_dictionary.extended.tsv"
crispytargets map \
    --gff "$KY21_GFF" \
    --blast-file "$TMP/${SP_LOWER}_vs_ky21.tsv" \
    --species "$SPECIES" \
    --evalue 1e-50 \
    --output "$RES/KY21_dictionary.extended.tsv" \
    "${MAP_ARGS[@]}"

echo "=== clean-mapping"
crispytargets clean-mapping \
    -i "$RES/KY21_dictionary.extended.tsv" \
    -o "$RES/KY21_dictionary.extended.cleaned.tsv" \
    -c "${SPECIES}_primary,${SPECIES}_secondary"
mv "$RES/KY21_dictionary.extended.cleaned.tsv" "$RES/KY21_dictionary.extended.tsv"

echo "=== qa"
python "$QA" \
    --base "$DICT" \
    --extended "$RES/KY21_dictionary.extended.tsv" \
    --species "$SPECIES" \
    --reference "$REF" \
    --projection "$RES/$PROJECTION" \
    --report "$RES/qa_report.tsv"

echo "=== pointer for /work"
cat > "$RES/SCRATCH_BUNDLE.txt" <<EOF
Bulk for this run lives in the crispytargets scratch workspace:
  $DATA

inputs/  staged inputs and PROVENANCE (sha256sum -c PROVENANCE.sha256)
work/    ${SP_LOWER}_db.* and ${SP_LOWER}_vs_ky21.tsv, regenerable by this sbatch

Scratch has no snapshots. Everything here is either re-downloadable or
reproducible from tools/map_species.sbatch.
EOF

echo "=== checksums"
(cd "$INPUTS" && sha256sum "$(basename "$KY21_PROT")" "$(basename "$KY21_GFF")" \
    "$(basename "$PROT")" "$(basename "$DICT")" "$(basename "$REF")")
(cd "$RES" && sha256sum KY21_dictionary.extended.tsv "$PROJECTION" qa_report.tsv)

echo "=== done at $(date -Is)"
