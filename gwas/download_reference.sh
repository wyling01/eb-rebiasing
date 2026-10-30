#!/usr/bin/env bash
# Download chr18-22 and retain the exact phased EUR variants used in the study.
set -euo pipefail
cd "$(dirname "$0")"

for program in bcftools curl unzip; do command -v "$program" >/dev/null; done
base_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
map_url="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"
reference_dir="data/simulation/reference"
download_dir="data/downloads"
samples="data/simulation/EUR.samples.txt"
selection="data/simulation/snp_selection.tsv"
mkdir -p "$reference_dir" "$download_dir"
work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT

# Keep completed downloads; resume interrupted downloads from their .part files.
download() {
  if [[ ! -f "$2" ]]; then
    curl --fail --location --retry 3 --continue-at - --output "$2.part" "$1"
    mv "$2.part" "$2"
  fi
}

# Reproduce the original EUR sample-ID order, independent of system locale.
panel_name="integrated_call_samples_v3.20130502.ALL.panel"
download "$base_url/$panel_name" "$download_dir/$panel_name"
awk '$3 == "EUR" {print $1}' "$download_dir/$panel_name" |
  LC_ALL=C sort -u > "$work_dir/EUR.samples.txt"
if [[ $(wc -l < "$work_dir/EUR.samples.txt") -ne 503 ]]; then
  echo "Expected 503 EUR samples in the Phase 3 metadata." >&2
  exit 1
fi
mv "$work_dir/EUR.samples.txt" "$samples"

download "$map_url" "$download_dir/plink.GRCh37.map.zip"
unzip -p "$download_dir/plink.GRCh37.map.zip" README.txt \
  > "$reference_dir/map_source.txt"
chromosomes=("$@")
if [[ ${#chromosomes[@]} -eq 0 ]]; then chromosomes=(18 19 20 21 22); fi
for chr in "${chromosomes[@]}"; do
  case "$chr" in 18|19|20|21|22) ;; *) echo "Expected chromosomes 18-22." >&2; exit 1 ;; esac
  name="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  output="$reference_dir/chr${chr}.vcf.gz"
  if [[ ! -f "$output" ]]; then
    download "$base_url/$name" "$download_dir/$name"
    download "$base_url/$name.tbi" "$download_dir/$name.tbi"
    awk -F '\t' -v chr="$chr" 'NR>1 && $1==chr {print $1 "\t" $2}' \
      "$selection" > "$work_dir/sites.tsv"
    bcftools view -R "$work_dir/sites.tsv" -S "$samples" -m2 -M2 -v snps \
      -Ou "$download_dir/$name" |
      bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF |
      bcftools view -i 'INFO/MAF>0.05' -Oz -o "$work_dir/chr.vcf.gz"
    output="$work_dir/chr.vcf.gz"
  fi

  # Check sample order and variant/allele identity before simulating families.
  awk -F '\t' -v chr="$chr" 'NR>1 && $1==chr {print $1 "\t" $2 "\t" $3 "\t" $4}' \
    "$selection" > "$work_dir/expected.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$output" > "$work_dir/observed.tsv"
  bcftools query -l "$output" > "$work_dir/samples.txt"
  diff "$work_dir/expected.tsv" "$work_dir/observed.tsv"
  diff "$samples" "$work_dir/samples.txt"
  if [[ "$output" == "$work_dir/chr.vcf.gz" ]]; then
    mv "$output" "$reference_dir/chr${chr}.vcf.gz"
  fi
  unzip -p "$download_dir/plink.GRCh37.map.zip" "plink.chr${chr}.GRCh37.map" \
    > "$reference_dir/chr${chr}.map"
  echo "chr${chr}: selected EUR reference ready"
done
