###########
## bash script to parse vcf files.
## I use vcftools, whose documentation may be found here: http://vcftools.sourceforge.net/man_latest.html
## To install vcftools, simply compile it then add its path to .bashrc via
##    export PATH="/path/to/dir:$PATH"
## Then the `vcftools` function should be available for use in this file.
###########

## Output bi-allelic counts for all sites in each chromosome. I.e. exclude sites that have more than two alleles.
vcftools --vcf ../data/vcf_Guilherme/myGBSGenos_mergedSNPs_mergedTaxa_chr1.vcf --extract-FORMAT-info  --min-alleles 2 --max-alleles 2 --out ../output/allele_counts/chr1_ct
