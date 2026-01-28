# This code was used to preprocessed the summary data used in the manuscript
# Writted by Yang Xiao, PKU, 2025
# @: xiaoyang9604@gmail.com

rm(list=ls())

library(stringr)
library(dplyr)
library(tidyr)
library(plinkbinr)
library(ieugwasr)
library(LDlinkR)
library(data.table)
library(rtracklayer)
chain <- import.chain("hg38ToHg19.over.chain")
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snps <-SNPlocs.Hsapiens.dbSNP155.GRCh38

# set workplace
input_path <-  ''
output_path <-  ''
IDPs_path <- ''

exposure_list <- c("")
exposure_filenames <- paste0(input_path, exposure_list,'.saige.txt.gz')

# remove LD Clumping
exposure_data <- fread(exposure_filenames, 
                       select = c("CHR", "POS", "Allele1", "Allele2", "AF_Allele2", 
                                  "BETA", "SE", "p.value", "N", "QC"))
exposure_data <- exposure_data[exposure_data$QC != "FAIL", ]
exposure_data$CHR <- sub("chr", "", exposure_data$CHR)

gr38 <- GRanges(seqnames = exposure_data$CHR,
                ranges = IRanges(start = exposure_data$POS, 
                                 end = exposure_data$POS))
snp.res <- snpsByOverlaps(snps, gr38)
snp.res <- as.data.table(snp.res)

exposure_data <- merge(exposure_data, snp.res[, .(seqnames, pos, RefSNP_id)], 
                       by.x = c("CHR", "POS"), by.y = c("seqnames", "pos"), all.x = TRUE)
exposure_data <- exposure_data %>% drop_na()

exposure_data$CHR <- paste0("chr", exposure_data$CHR)
gr38 <- GRanges(seqnames = exposure_data$CHR,
                ranges = IRanges(start = exposure_data$POS, 
                                 end = exposure_data$POS),
                id = exposure_data$RefSNP_id)

gr37 <- liftOver(gr38, chain)
gr37 <- unlist(gr37)

gr37_dt <- data.table(id = mcols(gr37)$id,
                      chromosome_GRCh37 = sub("chr", "", as.character(seqnames(gr37))),
                      pos_GRCh37 = start(gr37))

merged_data <- merge(exposure_data, gr37_dt, by.x = "RefSNP_id", by.y = "id", all.x = TRUE)                            
exposure_data <- dplyr::select(merged_data, c("RefSNP_id", 'chromosome_GRCh37','pos_GRCh37',
                                              'Allele1','Allele2',
                                              'AF_Allele2','BETA',
                                              'SE','p.value',"N"))
colnames(exposure_data) <- c('SNP','CHR','BP','A1','A2','MAF','BETA','SE','P','N')
exposure_data <- exposure_data %>% drop_na()

exposure_data <- exposure_data %>%
  dplyr::filter(!(CHR == 6 & BP >= 25000000 & BP <= 35000000))

exposure_data <- subset(exposure_data, P < 5e-8)

exposure_data <- data.frame(exposure_data)
exposure_stress <- format_data(
  exposure_data,
  type = 'exposure',
  header = T,
  phenotype_col = 'Phenotype',
  snp_col = 'SNP',
  beta_col = 'BETA',
  effect_allele_col = 'A1',
  other_allele_col = 'A2',
  pval_col = 'P',
  se_col = 'SE',
  eaf_col = 'EAF',
  samplesize_col = 'N'
)

# locally clumpping
exposure_stress$id <- exposure_stress$id.exposure
exposure_stress$rsid <- exposure_stress$SNP
exposure_stress$pval <- exposure_stress$pval.exposure

bfile <- "g1000_eas"
exposure_stress_clumped <- ld_clump(exposure_stress,
                                    plink_bin = get_plink_exe(),
                                    bfile = bfile,
                                    clump_kb = 10000, 
                                    clump_p = 5e-8,
                                    clump_r2 = 0.001
)
exposure_stress_clumped$id <- NULL
exposure_stress_clumped$rsid <- NULL
exposure_stress_clumped$pval <- NULL

print(paste0("Number of IVs: ", as.character(length(exposure_stress_clumped$SNP))))
write.table(exposure_stress_clumped, 
            str_c(output_path, 'exposure_', exposure_list, '_clumped.tsv'), 
            col.names = T, sep = "\t", row.names = F, quote = F)

