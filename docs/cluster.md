# Cluster: running CRISPYtargets on Unity

Conventions for this project on the UMass Unity cluster. Paths and job
patterns that a session would otherwise rediscover.

## Layout

```
/work/pi_ybernadskaya_umass_edu/dtjaernberg_umass_edu/projects/crispytargets/
├── code/     staged sbatch scripts
├── logs/     SLURM stdout/stderr
└── results/  deliverables only: mappings, reports, checksums
```

The code of record is the GitHub clone at
`~/projects/src/crispytargets` (installed editable into the
`mycrispytargets` conda env). Job scripts are authored in the repo under
`tools/` and copied to `code/` for submission.

## Scratch: what goes there

Scratch is the default home for anything bulk, downloaded, or
regenerable: genomes and proteomes, annotation files, BLAST databases and
BLAST tables, converted intermediates. `/work` keeps results, provenance
and non-regenerable artifacts only.

The project workspace is `crispytargets`, allocated with
`ws_allocate crispytargets 14`. Resolve its path at run time from
`ws_list -v`; never write the `/scratchN/...` path into a file. If the
workspace is missing or expired, reallocate and re-stage, do not fall
back to `/work`.

```
<workspace>/kh2013/inputs/   staged inputs + PROVENANCE.sha256
<workspace>/kh2013/work/     databases, BLAST tables, intermediates
```

Results carry a `SCRATCH_BUNDLE.txt` pointing at the scratch trees they
came from, so a run stays traceable after the workspace expires.

## Jobs

- Before submitting, run the suite: `python -m pytest tests/` locally and,
  after the clone is pulled on the cluster, the same suite in
  `mycrispytargets`. Java and minimap2 tests skip there unless those tools
  are on PATH.
- Submit with `sbatch --parsable` from the directory the script assumes;
  the script cd's to absolute paths anyway, so the submit cwd does not
  matter.
- Every script sources `/etc/profile` itself (the module trap) and keeps
  nounset off across activation.
- External tools come from modules: `blast-plus/2.16.0`, `openjdk/21`,
  `bcftools/1.19`, `samtools/1.19.2`. There is no minimap2 module; the
  lift path needs it installed in the env.
- The Python env is `mycrispytargets`:
  `module load conda/latest`, then `conda activate mycrispytargets`.
- Wall time is fixed at submit and cannot be raised; the `cpu` partition
  allows 14 days.

## Inputs

Reference downloads (GHOST, UniProt) are fetched on the local machine and
scp'd, because GHOST refuses compute-node egress. Each staged bundle
carries `PROVENANCE.txt` and `PROVENANCE.sha256`; the job verifies the
checksums before doing any work.

The FlashFry JAR is not in git. After cloning into `~/projects/src/`, copy
it in so the design tests run instead of skipping:

```
scp tools/FlashFry-assembly-1.15.jar unity:projects/src/crispytargets/tools/
```

With java, minimap2, tabix and the JAR present, the cluster suite runs
complete: 60 passed, 0 skipped (2026-11-16).

## Known input quirks

GHOST's KH2013 longest-ORF protein FASTA (`KH.KHGene.2013.Longest.protein.fasta`)
uses `+` as a stop marker, which `makeblastdb` rejects. Stage the cleaned
copy (`KH.KHGene.2013.Longest.protein.clean.fasta`, `+` stripped) instead;
1,357 lines carry it.
