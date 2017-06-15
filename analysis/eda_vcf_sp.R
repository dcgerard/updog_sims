
library(vcfR)
library(tidyverse)

dat <- vcfR::read.vcfR(file = "~/Data/shirasawa_etal/KDRIsweetpotatoXushu18S1LG2017.vcf")


queryMETA(dat)
unique_vals <- stringr::str_replace(queryMETA(dat), "^.+=ID=", "")
for (index in 1:length(unique_vals)) {
  print(queryMETA(dat, element = unique_vals[index]))
}

queryMETA(dat, element = "AD")
queryMETA(dat, element = "RD")
queryMETA(dat, element = "GQ")

## WT = number of samples called reference
## HET = number of samples called heterozygous-variant
## GT = Genotype

alt_mat <- apply(extract.gt(x = dat, element = "AD"), 2, as.numeric)
ref_mat <- apply(extract.gt(x = dat, element = "RD"), 2, as.numeric)
gq_mat  <- apply(extract.gt(x = dat, element = "GQ"), 2, as.numeric)


var_ave <- rowMeans(ref_mat + alt_mat, na.rm = TRUE)
order_vec <- order(var_ave, decreasing = TRUE)

var_ave[order_vec[1]]
var_ave[order_vec[length(order_vec)]]

library(updog)
par(ask = TRUE)
for (index in 500:600) {
  refvec <- ref_mat[order_vec[index], ]
  altvec <- alt_mat[order_vec[index], ]
  cat(index)
  print(plot_geno(ocounts = refvec, osize = refvec + altvec, ploidy = 6))
}
par(ask = FALSE)

ratio_ave <- rowMeans((ref_mat - alt_mat) / ref_mat, na.rm = TRUE)

index <- which.max(ratio_ave[var_ave > 100])
refvec <- ref_mat[var_ave > 100, ][index, ]
altvec <- alt_mat[var_ave > 100, ][index, ]
print(plot_geno(ocounts = refvec, osize = refvec + altvec, ploidy = 6))

## Extract the one's I want
snp_nums <- order_vec[c(506, 519, 517)] ## OD, bias, and outlier, in that order

example_ref <- t(ref_mat[snp_nums, ])
example_alt <- t(alt_mat[snp_nums, ])
example_tot <- example_ref + example_alt
colnames(example_ref) <- c("SNP1", "SNP2", "SNP3")
colnames(example_tot) <- c("SNP1", "SNP2", "SNP3")

write.csv(example_ref, file = "../output/shirasawa_snps/example_refcounts.csv", row.names = TRUE)
write.csv(example_tot, file = "../output/shirasawa_snps/example_readcounts.csv", row.names = TRUE)
