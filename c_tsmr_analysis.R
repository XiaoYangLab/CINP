# This code was used to conduct two-sample Medenlian randomization
# Writted by Yang Xiao, PKU, 2025
# @: xiaoyang9604@gmail.com

rm(list = ls())

library(TwoSampleMR)
library(stringr)
library(dplyr)
library(tidyr)
library(MRPRESSO)
library(mr.raps)
library(RadialMR)
library(ggplot2)
library(plinkbinr)
library(ieugwasr)
library(LDlinkR)
library(data.table)

# 1. set workplace
input_path <-  ' '
output_path <-  ' '
IDPs_path <- ' '

# 2. import exposure data
exposure_list <- c(" ")

exposure_filenames <- paste0(output_path, 'exposure_', exposure_list,'_clumped.tsv')
exposure_clumped <- fread(exposure_filenames)
exposure_clumped <- data.frame(exposure_clumped)

# 3. import outcome data
IDPfile <- list.files(path = IDPs_path)

res_df <- data.frame(
  IDP = character(),
  nsnp = numeric(),
  beta_mr = numeric(),
  se_mr = numeric(),
  pIVW_mr = numeric(),
  hetero_Q = numeric(),
  MR_egger_inter = numeric(),
  MR_presso = numeric(),
  MR_direct = numeric(),
  MR_loo = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(IDPfile)) {
  IDPsnp_dat <- fread(paste0(IDPs_path, IDPfile[i]),
                      select = c("rsid", "chromosome", "base_pair_location", 
                                 "effect_allele", "other_allele", 
                                 "effect_allele_frequency", "beta", "standard_error", 
                                 "p_value"))
  IDPsnp_dat$N <- 7058
  IDPsnp_dat <- data.frame(IDPsnp_dat)
  IDPsnp_dat <- IDPsnp_dat %>%
    dplyr::filter(!(chromosome == 6 & base_pair_location >= 25000000 & base_pair_location <= 35000000))
  
  snps <- intersect(exposure_clumped$SNP,IDPsnp_dat$rsid)
  if (length(snps) == 0){
    IDPsnp_outcome <- format_data(
      snps = IDPsnp_dat$SNP[1],
      IDPsnp_dat,
      type = 'outcome',
      header = TRUE,
      phenotype_col = 'Phenotype',
      snp_col = "rsid",
      beta_col = "beta",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      se_col = 'standard_error',
      eaf_col = 'effect_allele_frequency',
      samplesize_col = "N"
    )
  } else {
    IDPsnp_outcome <- format_data(
      snps = snps, # use data after clean
      IDPsnp_dat,
      type = 'outcome',
      header = TRUE,
      phenotype_col = 'Phenotype',
      snp_col = "rsid",
      beta_col = "beta",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      se_col = 'standard_error',
      eaf_col = 'effect_allele_frequency',
      samplesize_col = "N"
    )
  }
  
  # Identifying & printing exposure instruments missing from outcome GWAS
  missing_IVs <- exposure_clumped$SNP[!(exposure_clumped$SNP %in% IDPsnp_outcome$SNP)]

  # Replacing missing instruments from outcome GWAS with proxies
  if(length(missing_IVs) == 0) {
    print("All exposure IVs found in outcome GWAS.")
  } else {
    print("Some exposure IVs missing from outcome GWAS.")

    for (m in 1:length(missing_IVs)) {
      print(paste0("Processing missing_ivs [", m, "/", length(missing_IVs), "]"))
      proxies <- LDproxy(snp = missing_IVs[m],
                         pop = "EAS",
                         r2d = "r2",
                         token = "c91e1d4d9dee",
                         file = FALSE)
      proxies <- proxies[proxies$R2 > 0.8, ]
      proxies <- proxies[grepl('rs', proxies$RS_Number), ]
      proxy_present = FALSE

      if(length(proxies$RS_Number) == 0){

        print(paste0("No proxy SNP available for ", missing_IVs[m]))

      } else {

        for (n in 1:length(proxies$RS_Number)) {

          proxy_present <- proxies$RS_Number[n] %in% IDPsnp_dat$SNP

          if (proxy_present) {
            proxy_SNP = proxies$RS_Number[n]
            proxy_SNP_allele_1 = str_sub(proxies$Alleles[n], 2, 2)
            proxy_SNP_allele_2 = str_sub(proxies$Alleles[n], 4, 4)
            original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
            original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
            break
          }
        }
      }

      if(proxy_present == TRUE) {
        print(paste0("Proxy SNP found. ", missing_IVs[m], " replaced with ", proxy_SNP))
        proxy_row <- IDPsnp_outcome[1, ]
        proxy_row$SNP = missing_IVs[m]
        proxy_row$beta.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "beta"])
        proxy_row$se.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "standard_error"])
        if (IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
        if (IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
        if (IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
        if (IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
        proxy_row$pval.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "p_value"])
        proxy_row$samplesize.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "N"])

        if("N_case" %in% colnames(IDPsnp_dat)) proxy_row$ncase.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "N_case"])
        if("N_control" %in% colnames(IDPsnp_dat)) proxy_row$ncontrol.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "N_control"])

        if("effect_allele_frequency" %in% colnames(IDPsnp_dat)) proxy_row$eaf.outcome = as.numeric(IDPsnp_dat[IDPsnp_dat$rsid == proxy_SNP, "effect_allele_frequency"])
        IDPsnp_outcome <- rbind(IDPsnp_outcome, proxy_row)
      }

      if(proxy_present == FALSE) {
        print(paste0("No proxy SNP available for ", missing_IVs[m], " in outcome GWAS."))
      }
    }

    if (length(snps) == 0){
      IDPsnp_outcome <- IDPsnp_outcome[2:nrow(IDPsnp_outcome), ]
    }

  }
  
  # 4. harmonization
  harm_dat <- harmonise_data(exposure_dat = exposure_clumped,
                             outcome_dat = IDPsnp_outcome,
                             action = 2)
  harm_dat <- harm_dat[which(harm_dat$mr_keep == T), ]
  
  if (length(harm_dat$SNP)  == 0) {
    next
  }
  
  # 5. test outliers by using ivw_radial & egger_radial
  if (length(harm_dat$SNP) > 4){
    outlier1 <- ivw_radial(harm_dat, 0.05, 1, 1e-4)
    outlier2 <- egger_radial(harm_dat, 0.05, 1)
    
    if (is.data.frame(outlier1$outliers)) {
      a1 <- cbind(outlier1$outliers[1])
      a1 <- as.matrix(a1)
    } else {
      a1 <- matrix()
    }
    
    if (is.data.frame(outlier2$outliers)) {
      a2 <- cbind(outlier2$outliers[1])
      a2 <- as.matrix(a2)
    } else {
      a2 <- matrix()
    }
    
    a3 <- rbind(a1, a2)
    a3 <- na.omit(a3)
    a3 <- as.data.frame(a3)
    a3 <- unique(a3)
    harm_dat <- harm_dat[which(!(harm_dat$SNP %in% a3$SNP)), ]
  }
  
  # 6. Calculating the F-statistic
  harm_dat$rsq <- harm_dat$beta.exposure * harm_dat$beta.exposure / (harm_dat$beta.exposure * harm_dat$beta.exposure + harm_dat$se.exposure * harm_dat$se.exposure * harm_dat$samplesize.exposure)
  harm_dat$f <- sum(harm_dat$rsq) * (harm_dat$samplesize.exposure-length(harm_dat$SNP)-1) / length(harm_dat$SNP)*(1-sum(harm_dat$rsq))
  harm_dat <- harm_dat[harm_dat$f > 10, ]
  
  # 7. heterogeneity test (p>0.05)
  if (length(harm_dat$SNP) > 1){
    heterogeneity = mr_heterogeneity(harm_dat)
    if (heterogeneity$Q_pval[1] < 0.05) {
      method_list = c('mr_ivw_mre', 'mr_weighted_median', 
                      'mr_weighted_mode','mr_simple_mode','mr_egger_regression')
    } else {
      method_list = c('mr_ivw_fe', 'mr_weighted_median',
                      'mr_weighted_mode','mr_simple_mode','mr_egger_regression')
    }
  } else {
    method_list = c('mr_ivw_fe', 'mr_weighted_median',
                    'mr_weighted_mode','mr_simple_mode','mr_egger_regression')
  }
  
  # 8. Performing Mendelian randomization analyses
  if (length(harm_dat$SNP) > 1){
    tsmr <- mr(harm_dat, method_list = method_list)
    raps <- mr_raps(harm_dat$beta.exposure,
                    harm_dat$beta.outcome,
                    harm_dat$se.exposure,
                    harm_dat$se.outcome
    )
    
    tsmr[6,] = tsmr[1,]
    tsmr$method[6] = 'MR RAPS'
    tsmr$b[6] = raps$b
    tsmr$se[6] = raps$se
    tsmr$pval[6] = raps$pval
    
    # 9. sensitivity analysis
    # 9.1 horizontal pleiotropy by using MR_PRESSO global test
    if (length(harm_dat$SNP) > 3){
      Hori_pleio <- mr_presso(BetaOutcome = 'beta.outcome',
                              BetaExposure = 'beta.exposure',
                              SdOutcome = 'se.outcome',
                              SdExposure = 'se.exposure',
                              data = harm_dat,
                              OUTLIERtest = T,
                              seed = 123456)
      
      tsmr[7,] = tsmr[1,]
      tsmr$method[7] = 'PRESSO'
      tsmr$b[7] = Hori_pleio$`Main MR results`$`Causal Estimate`[1]
      tsmr$se[7] = Hori_pleio$`Main MR results`$`Sd`[1]
      tsmr$pval[7] = Hori_pleio$`Main MR results`$`P-value`[1]
    }
  } else {
    tsmr <- mr(harm_dat)
  }
  
  print(tsmr$pval)
  
  # save significant IDP
  write.table(tsmr,
              paste0(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10),'_mr.tsv'),
              col.names = T, sep = "\t", row.names = F, quote = F)
  
  if (length(harm_dat$SNP) > 3){
    # save horizontal pleiotropy
    write.table(Hori_pleio$`MR-PRESSO results`$`Global Test`,
                str_c(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10) ,'_hori_pleio.tsv'),
                col.names = T, sep = "\t", row.names = F, quote = F)
  }
  # 9.2 pleiotropy using MR-Egger intercept
  mean_pleio <- mr_pleiotropy_test(harm_dat)
  write.table(mean_pleio,
              str_c(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10) ,'_mean_pleio.tsv'),
              col.names = T, sep = "\t", row.names = F, quote = F)
  
  # 9.3 leave-one-out analysis
  tsmr_single <- mr_singlesnp(harm_dat)
  write.table(tsmr_single,
              str_c(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10) ,'_tsmr_single.tsv'),
              col.names = T, sep = "\t", row.names = F, quote = F)
  
  # 9.4 steiger directionality test
  whole_direct <- directionality_test(harm_dat)
  write.table(whole_direct,
              str_c(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10) ,'_whole_direct.tsv'),
              col.names = T, sep = "\t", row.names = F, quote = F)
  
  snp_direct <- steiger_filtering(harm_dat)
  write.table(snp_direct,
              str_c(output_path, exposure_list, '_', str_sub(IDPfile[i], 1, -10) ,'_snp_direct.tsv'),
              col.names = T, sep = "\t", row.names = F, quote = F)
  
  # save beta of MR-IVW
  res_df <- rbind(res_df,
                  data.frame(
                    IDP = str_sub(IDPfile[i], 1, -10),
                    nsnp = length(harm_dat$SNP),
                    beta_mr = tsmr$b[1],
                    se_mr = tsmr$se[1],
                    pIVW_mr = tsmr$pval[1],
                    hetero_Q = heterogeneity$Q_pval[1],
                    MR_egger_inter = mean_pleio$pval,
                    MR_presso = tsmr$b[7],
                    MR_direct = whole_direct$steiger_pval,
                    MR_loo = sum(tsmr_single$p < 0.05)
                  ))
  
}

write.table(res_df, str_c(output_path, exposure_list,'_IDPs_mr.csv'),  
            col.names = T, sep = ",", row.names = FALSE, quote = FALSE)
