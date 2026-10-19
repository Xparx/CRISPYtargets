# Tutorial: lift SNPs onto your assembly

**Deliverable:** a BED file (and lifted VCF) of your SNP variants in your
assembly's coordinates, ready for `crispytargets design --snps`.

**Read this first:** this step is *optional* and it is *slow*. The design
pipeline runs perfectly well without any SNP data, and everything else in
CRISPYtargets is unaffected if you skip it. You only need this tutorial if
you want SNP-aware output (variant-rescored guides, entropy, distances). The
slowness comes from one whole-genome alignment: roughly 5 minutes for the
Ciona example below (125 Mb genome, 2.2 million SNP records), and it scales
with genome size.

## What you need

Prerequisites on your system:

- `crispytargets` installed (`pip install git+https://github.com/Xparx/CRISPYtargets.git`)
- `minimap2` and `bcftools` (with `tabix`) on your PATH, e.g.
  `conda install -c bioconda minimap2 bcftools`

Set up your workspace. Everything below assumes this one variable; nothing
has to live in any particular place:

```bash
working_dir=~/ciona_snps                    # change to taste
mkdir -p "$working_dir"/inputs "$working_dir"/outputs
```

Download the three inputs into `"$working_dir"/inputs/`:

| Input | Where it comes from |
|-------|---------------------|
| `Ghost_SNPs.vcf.gz` (SNPs, in ci3/old coordinates) | kept in [cionalib](https://gitlab.com/Xparx/cionalib) at `cionalib/data/SNPs/` because it is hard to track down anywhere else: [file link](https://gitlab.com/Xparx/cionalib/blob/main/cionalib/data/SNPs/Ghost_SNPs.vcf.gz) |
| `ci3.2bit` (the old assembly those SNPs were called on) | https://crispor.gi.ucsc.edu/genomes/ci3/ci3.2bit |
| `HT.Ref.fasta` (the assembly you design against, KY21) | https://ghost.zool.kyoto-u.ac.jp/datas/HT.Ref.fasta.zip |

For any organism other than Ciona the shape is the same: your variant calls
in the coordinates of the assembly they were called on, that old assembly,
and the assembly you want to design against.

## Step 1: keep the main chromosomes only (optional but recommended)

The Ciona SNP file also contains thousands of unplaced scaffolds. Restricting
to the 14 main chromosomes keeps the alignment short and the output clean:

```bash
bcftools view -r KhC1,KhC2,KhC3,KhC4,KhC5,KhC6,KhC7,KhC8,KhC9,KhC10,KhC11,KhC12,KhC13,KhC14 \
    "$working_dir"/inputs/Ghost_SNPs.vcf.gz \
  | bgzip > "$working_dir"/inputs/Ghost_SNPs_main.vcf.gz
tabix -p vcf "$working_dir"/inputs/Ghost_SNPs_main.vcf.gz
```

If your SNPs are already on main chromosomes only, skip this step.

## Step 2: lift

```bash
crispytargets lift \
    --old-genome "$working_dir"/inputs/ci3.2bit \
    --new-genome "$working_dir"/inputs/HT.Ref.fasta \
    --vcf "$working_dir"/inputs/Ghost_SNPs_main.vcf.gz \
    --output "$working_dir"/outputs/snps_lifted/ \
    --bed \
    --threads 8 \
    --preset asm10
```

Notes:

- `--old-genome` accepts a `.2bit` file directly and converts it for you.
- `--preset asm10` fits two assemblies of the same species that have drifted
  a little; use `asm5` for very close assemblies, `asm20` for distant ones.
- Variant positions that cannot be placed on the new assembly are dropped and
  summarised per chromosome in the log. Do not expect every record to survive:
  unplaced scaffolds and rearranged regions are the usual losses.

## What you get

In `"$working_dir"/outputs/snps_lifted/`:

| File | Use it for |
|------|-----------|
| `Ghost_SNPs_main.bed_lifted.bed` | pass to `crispytargets design --snps` (see [tutorial_design.md](tutorial_design.md)) |
| `Ghost_SNPs_main.vcf_lifted.vcf.gz` (+ `.tbi`) | the same lifted variants in VCF form, for other tools |

Naming follows the input file name, so a different input gives differently
named outputs.