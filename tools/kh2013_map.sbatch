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
# BLAST direction matches the other species columns: KY21 proteins are the
# query, the species proteome is the database. The mapper is then run on a
# copy of the existing dictionary so the incremental backbone keeps every
# column that is already there.

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
DATA="$PROJ/data/kh2013"
RES="$PROJ/results/kh2013_mapping"
QA="$HOME/projects/src/crispytargets/tools/qa_kh2013_mapping.py"

KY21_PROT="$DATA/HT.KY21Gene.protein.2.fasta"
KY21_GFF="$DATA/HT.KY21Gene.2.gff3"
KH_PROT="$DATA/KH.KHGene.2013.Longest.protein.fasta"
DICT="$DATA/KY21_dictionary.tsv"
GHOST="$DATA/KHtoKY21.txt.gz"

echo "=== kh2013_map on $(hostname) at $(date -Is)"
echo "job ${SLURM_JOB_ID}, ${SLURM_CPUS_PER_TASK} cpus"
blastp -version
python -V

for input in "$KY21_PROT" "$KY21_GFF" "$KH_PROT" "$DICT" "$GHOST" "$QA"; do
    test -s "$input" || { echo "missing input: $input"; exit 1; }
done

mkdir -p "$RES"

echo "=== makeblastdb"
makeblastdb -in "$KH_PROT" -dbtype prot -out "$DATA/kh2013_db"

echo "=== blastp"
blastp -query "$KY21_PROT" \
       -db "$DATA/kh2013_db" \
       -out "$RES/kh2013_vs_ky21.tsv" \
       -outfmt 6 \
       -evalue 1e-5 \
       -max_target_seqs 20 \
       -num_threads "$SLURM_CPUS_PER_TASK"

echo "=== map (incremental on a copy of the dictionary)"
cp "$DICT" "$RES/KY21_dictionary.extended.tsv"
crispytargets map \
    --gff "$KY21_GFF" \
    --blast-file "$RES/kh2013_vs_ky21.tsv" \
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

echo "=== checksums"
(cd "$DATA" && sha256sum "$(basename "$KY21_PROT")" "$(basename "$KY21_GFF")" \
    "$(basename "$KH_PROT")" "$(basename "$DICT")" "$(basename "$GHOST")")
(cd "$RES" && sha256sum KY21_dictionary.extended.tsv KY21_to_KH2013.tsv qa_report.tsv)

echo "=== done at $(date -Is)"
