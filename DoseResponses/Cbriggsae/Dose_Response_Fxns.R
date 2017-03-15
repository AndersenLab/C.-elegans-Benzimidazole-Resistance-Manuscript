## HERITABILITY
# data is data frame that contains strain and phenotype column
# indicies are used by the boot function to sample from the 'data' data.frame
H2.test.boot <- function(data, indicies){
  
  d <- data[indicies,]
  
  pheno <- as.data.frame(dplyr::select(d, phenotype))[,1]
  strain <- as.factor(d$strain)
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)  
}

# data is data frame that contains strain and phenotype column
H2.test <- function(data){
  
  pheno <- as.data.frame(dplyr::select(data, phenotype))[,1]
  strain <- as.factor(data$strain)
  
  reffMod <- lme4::lmer(pheno ~ 1 + (1|strain))
  
  Variances <- as.data.frame(lme4::VarCorr(reffMod, comp = "Variance"))
  
  Vg <- Variances$vcov[1]
  Ve <- Variances$vcov[2]
  H2 <- Vg/(Vg+Ve)
  
  # errors <- sqrt(diag(lme4::VarCorr(reffMod, comp = "Variance")$strain))
  
  return(H2)  
}

# df is data frame that contains strain and phenotype column
H2.calc <- function(df, boot = T){
  if(boot == T){
    # bootstrapping with 1000 replications 
    results <- boot(data=df, statistic=H2.test.boot, R=1000)
    # get 95% confidence interval 
    ci <- boot.ci(results, type="bca")
    H2_errors <- data.frame(H2 = ci$t0, ci_l = ci$bca[4], ci_r = ci$bca[5])
    return(H2_errors)
  } else {
    H2 <- data.frame(H2 = H2.test(data = df), ci_l = NA, ci_r = NA)
    return(H2)
  }
}