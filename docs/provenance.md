# Provenance and verification log

Audit record for the *Ciona robusta* worked example: where every input came
from, its checksum, and the verification performed when the recipes in
[tutorial_snps.md](tutorial_snps.md), [tutorial_design.md](tutorial_design.md)
and [tutorial_mapping.md](tutorial_mapping.md) were established. This file is
dense by design; the tutorials are the user-facing documents.

All source files verified 2026-10-19 against their upstream archives.

## Source inputs

| File | Source | sha256 |
|------|--------|--------|
| `HT.Ref.fasta` (KY21 genome, 67 sequences) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.Ref.fasta.zip | `c72fcfc29ba51b564bcf3fa3623bb7e425bdf7eb942910390a42d2c17b2f4d75` |
| `HT.KY21Gene.2.gff3` (18,788 genes, 627,048 data lines) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.2.gff3.zip | `17f81536a5dbc33a69bc759827aac7f2f2b79c8fdf1963fd0495a0266fdff9d3` |
| `HT.KY21Gene.protein.2.fasta` (55,505 proteins, mapping query) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.KY21Gene.protein.2.fasta.zip | `af9cb597ed9a5b256375bf6db6f44c02fef5a175abb3ed9349c738458b964cd0` |
| `HT.Gene.protein.fasta` (68,338 KY-era proteins, renamed from `HT.KYGene.protein.fasta`; KY19 BLAST database) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.KYGene.protein.fasta.zip | `228aa89cd16cb1c1c8976b3d2f4ad69b1a1bdb7617a90a8068ae1f55a823a394` |
| `HT.RefwMG0.fasta` + `HT.KY21Gene.2.wMG0.sort.gtf` | https://ghost.zool.kyoto-u.ac.jp/datas/HT.RefwMG0.zip | `9604b7f3ee2b1a8161431085991577b4663a63c5ffb3a196e776ebabd1677690` / `de415fc73ef17df6d7d199b6a5f5240f0962ceb64e05d3ac466f95cd36964fe5` |
| `Ghost_SNPs.vcf.gz` + `.tbi` (3,233,449 SNPs, 1,271 contigs; upstream files dated 2016-10-04) | https://crispor.gi.ucsc.edu/genomes/ci3/ | `ef6d9e6a1e501a4a0caca388ee497c5bfdcbc1ddd241de9e734f4bf51e450ad2` / `1395a25635324eca7a046a210eafe3beb993a9bbe35aee0cd566b60fb5c201af` |
| `ci3.2bit` (old assembly, 1,272 sequences, `Kh*` names) | https://crispor.gi.ucsc.edu/genomes/ci3/ | `ca376f14e7dd0507edae43b5fd3612d31b12ca5bdeb1bb3f90a6044e0df544ee` |
| `uniprotkb_Human_AND_reviewed_true_AND_m_2026_04_21.fasta` | UniProt reviewed query recorded in the filename, retrieved 2026-04-20 (living source) | `b205ed9970fee62413e2f192c6fb82daab86e4a5897ef521f907ecbfcd7b3e63` |
| `uniprotkb_Mouse_AND_reviewed_true_AND_m_2026_04_21.fasta` | as above | `a655fdec2f9fe80ac036670ae9fd7bfb012537fabe3f818f2566967e01423343` |
| `genelist.tsv` (5-gene worked example) | authored locally, not fetchable | `65da29c5c2b18e1df9b54c46cb35c91d0b627299e30f097b84c4d92db809634a` |

The canonical copy of `Ghost_SNPs.vcf.gz` is kept in cionalib
(`cionalib/data/SNPs/`, https://gitlab.com/Xparx/cionalib) because it is the
hardest file to track down; byte-identical to the copy verified above.

## Derived inputs

| File | How it is derived | Verification |
|------|-------------------|--------------|
| `KH2013.fa` | FASTA conversion of `ci3.2bit` (`twoBitToFa ci3.2bit KH2013.fa`) | all 1,272 sequence names and lengths match the 2bit index; five windows (KhM0, KhC1, KhC2, KhC7, KhC14) base-identical |
| `Ghost_SNPs_main.vcf.gz` | main-chromosome subset of `Ghost_SNPs.vcf.gz`; exact command recorded in its own header: `bcftools view -r KhC1,...,KhC14 Ghost_SNPs.vcf.gz` (bcftools 1.20, 2026-08-03) | filtering the source by the chromosome list reproduces the 2,179,863 record count exactly |

Main chromosome list (from `main_chroms.txt`): KhC1, KhC2, KhC3, KhC4, KhC5,
KhC6, KhC7, KhC8, KhC9, KhC10, KhC11, KhC12, KhC13, KhC14.

## Deliverables and their verification

| Deliverable | sha256 | Verification |
|-------------|--------|--------------|
| `geneIDmapping.tsv` (becomes cionalib `KY21_dictionary.tsv`) | `1a276c7e2636ea3094379654287a117f8a36c5102dc5cfb295d434263420f9c2` | reproduced **byte-identically** on 2026-10-19 from the BLAST results of that run (since removed as regenerable intermediates; rebuild them per [tutorial_mapping.md](tutorial_mapping.md)), running `map` in KY19, Human, Mouse order then `clean-mapping` (18,788 rows). BLAST parameters consistent with those results: e-value ceiling exactly 1e-5; never more than 5 distinct subjects per query, so the `-max_target_seqs 20` cap never binds |
| `KY21_dictionary.extended.tsv` (13 columns, 18,788 rows; adds `KH2013_primary`, `KH2013_secondary`) | `aa05c786e9fe7d163f0885e374a153f009d98f3417825f7afac093283296b3a2` | produced 2026-11-15 by the KH2013 job (now `tools/map_species.sbatch`) on Unity (job 1792936, 57 min): `blastp` of the 55,505 KY21 proteins against GHOST's 55,017 KH2013 longest-ORF proteins (`+` stop markers stripped for `makeblastdb`), `-evalue 1e-5 -max_target_seqs 20`, then `map --species KH2013 --evalue 1e-50` and `clean-mapping`; the 11 earlier columns are unchanged. 14,119 rows (75.2%) gain a KH2013 primary; 90.9% of GHOST's 13,906 KH-to-KY21 pairs appear as primary or secondary, 89.9% of primaries agree with GHOST. Secondary hits never repeat the primary. The secondary columns were regenerated on 2026-11-16 with the fixed `clean-mapping`, without rerunning BLAST, after the mapper began deduplicating stripped accessions. Standalone projection `KY21_to_KH2013.tsv` sha `b5684115bea2ed05b6c827d97a69abb4162b4061e955b8bc1c79999eb48cb820`; method and cionalib options in [kh2013_dictionary_update.md](kh2013_dictionary_update.md) |
| `KY21_dictionary.extended.tsv` (15 columns, 18,788 rows; adds `ENSEMBL_primary`, `ENSEMBL_secondary`) | `e477bb5e4c3f9d8f08dc98032afa12a8a92bbbed3d56365d5536809c6bcf622a` | produced 2026-11-24 by `tools/map_species.sbatch` (`SPECIES=ENSEMBL`) on Unity (job 2233878, 15 min): `blastp` of the 55,505 KY21 proteins against Ensembl release 116 `pep.all` (17,302 peptides, 16,671 genes, KH assembly GCA_000224145.1), `-evalue 1e-5 -max_target_seqs 20`, then `map --species ENSEMBL --evalue 1e-50 --subject-map ensembl_protein_to_gene.tsv` (ENSCINP to ENSCING from the FASTA headers) and `clean-mapping`; the 13 earlier columns are unchanged. 12,378 rows (65.9%) gain an Ensembl primary over 11,145 distinct genes; 90.8% of the legacy KY19-chained bridge's 10,635 pairs appear as primary or secondary, and 80.1% of primaries agree where both map. Secondary hits never repeat the primary. Standalone `KY21_to_ENSEMBL.tsv` sha `7fca09eca467ab70530d4a377a7123aa87b1fbe9f502c5b5e78ce1c13c168815`; method and caveats in [ensembl_dictionary_update.md](ensembl_dictionary_update.md) |
| `Ghost_SNPs_main.bed_lifted.bed` | `916c5729c4b676f7dd967b6b4ae71a6f14cd4c54c034f7a0e5701bc5b016e642` | content-identical reproduction on 2026-10-19 via **two** verified commands, both giving sorted-content sha256 `e173234205e6a61e` over 2,021,692 lines: `--old-genome ci3.2bit` (auto-converts) and `--old-genome KH2013_main.fa`. Raw line order varies with `--threads` (minimap2 scheduling is not order-stable); compare sorted content, not file bytes |
| `Ghost_SNPs_main.vcf_lifted.vcf.gz` | `f8db9b769e6fda78398a2a3aca8dc33be2217739609d1af634e9a2d40b5eb6f2` | same runs; gzip member timestamps make byte comparison meaningless, compare decompressed content |
| `final_sgrnas.tsv` (40 columns incl. `Mouse_uniq` and `feature_index`, 426 rows) | `ddd0bd3160f7286be661a3c07f1ea397e2d22ccfede3153db5b4be98010cfbd0` | produced 2026-10-30 by the `design` command in tutorial_design.md plus the ortholog join; run log reports prepare 6.2 s, FlashFry 82.3 s, post-process 1126.5 s, total 1215.0 s; supersedes `4f2a3808…`, matching it byte-for-byte on all 39 previous columns and adding `feature_index` (`final_sgrnas.bed` unchanged, same sha) |

Liftover input accounting: 2,179,863 records in, 2,021,692 lifted BED lines
out; records that fail to align are reported per chromosome in the log.

The BLAST databases for the verified run were built 2026-04-20 from
`uniprotkb_Human_...fasta`, `uniprotkb_Mouse_...fasta` and
`HT.Gene.protein.fasta`, confirmed from the database headers before the
databases and result TSVs were removed as regenerable intermediates.

## Tool versions for the verified runs

| Tool | Version |
|------|---------|
| CRISPYtargets | commit `522a8f6` (no code change since `a25fbc0`, which fixed guide-level SNP assignment) |
| cionalib | commit `1ba45dc` at the time of verification |
| minimap2 | 2.28-r1209 |
| bcftools / tabix | 1.20 |
| BLAST+ (blastp) | 2.12.0+ |
| FlashFry | 1.15 (`FlashFry-assembly-1.15.jar`, OpenJDK 21) |
| Python | 3.10, pandas / pyranges / pyfaidx stack from `pyproject.toml` |