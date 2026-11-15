#!/bin/bash
#SBATCH --job-name=env_update
#SBATCH --output=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/env_update_%j.out
#SBATCH --error=/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/logs/env_update_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dtjaernberg+unity@umass.edu
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

# Add the external tools the design and lift paths need to mycrispytargets.
#
#   ssh unity 'cd /work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/code && sbatch --parsable env_update.sbatch'
#
# --freeze-installed keeps the existing solve intact; if that cannot
# resolve, the fallback allows changes and the list diff shows them.

set -eo pipefail

set +u
source /etc/profile
source /etc/profile.d/*.sh 2>/dev/null || true
module purge
module load conda/latest
eval "$(conda shell.bash hook)"
conda activate mycrispytargets

BEFORE=$(conda list --export)

echo "=== installing openjdk, minimap2 and tabix"
conda install -y --freeze-installed -c conda-forge -c bioconda openjdk minimap2 tabix \
    || conda install -y -c conda-forge -c bioconda openjdk minimap2 tabix

echo "=== packages changed by the install"
diff <(echo "$BEFORE") <(conda list --export) || true

echo "=== tool check"
java -version
minimap2 --version
python -V

echo "=== suite"
python -m pytest "$HOME/projects/src/crispytargets/tests" -q

echo "=== done at $(date -Is)"
