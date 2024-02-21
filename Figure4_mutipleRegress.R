library(relaimpo)
library(tcltk)  
perm_count <- 5000
u <- 1:perm_count  
k <- 1:8 
res <- array(0, dim = c(8, perm_count))
res_real <- array(0, dim = c(8, 1))
p_perm <- array(0, dim = c(8, 1))
r <- array(0, dim = c(1,perm_count ))
x <- read.csv('receptor data.csv')

################################################################################ step 1: linear fitting receptor

fit_receptor <- lm(gamlss_t ~ Serotonin + Dopamine + Histamine + Acetylcholine
                   + Cannabinoid + Glutamate + GABA+Noradrenaline ,
                  data = x)
for(j in k) { 
  res_real[j,1] <- unname(fit_receptor$coefficients[j+1])
}
write.table(fit_receptor$fitted.values,
            "data/receptor/recptorfit.csv",
            row.names=FALSE,
            col.names=TRUE,
            sep=",")

B=boot.relimp(fit_receptor, 
              bty = "perc" , 
              level = 0.95, 
              type = c( "pmvd"), 
              rela = TRUE, 
              fixed =FALSE)
booteval.relimp(B, 
                bty = "perc", 
                level = 0.95)

################################################################################ step 2: linear fitting fitted Value
corrFit <- data.frame(gamlss_t=x$gamlss_t, 
                      fit_value = fit_receptor$fitted.values)
fit_corr <- lm(fit_value ~ gamlss_t ,
           data = corrFit)

new_gamlss_t <- seq(min(corrFit$gamlss_t),
                    max(corrFit$gamlss_t),
                    0.01)
gamlss_td_wt <- data.frame(predict(fit_corr, 
                                   newdata = data.frame(gamlss_t = new_gamlss_t),
                                   interval = "confidence"), 
                           new_gamlss_t = new_gamlss_t)

write.csv(corrFit,file = "data/receptor/fit.csv",
          row.names = FALSE)
write.csv(gamlss_td_wt,file = "data/receptor/ci_df.csv",
          row.names = FALSE)


################################################################################ step 3: permutation test
t_map_perm = read.csv('data/receptor/t_map_perm.csv',
                      header = FALSE)

pb <- tkProgressBar("schedule","Completed %", 0, 100)  
for(i in u) {  
  info<- sprintf("Completed %d%%", round(i*100/max(u)))  
  setTkProgressBar(pb, 
                   i*100/max(u), 
                   sprintf("schedule (%s)", 
                           info),
                   info)  
  x[,1] <- t_map_perm[,i]
  
  fit_perm <- lm(gamlss_t ~ Serotonin + Dopamine + Histamine + Acetylcholine
                 + Cannabinoid + Glutamate + GABA+Noradrenaline ,
                 data = x)
  for(j in k) { 
    res[j,i] <- unname(fit_perm$coefficients[j+1])
  }
  corrFit <- data.frame(gamlss_t=x$gamlss_t, 
                        fit_value = fit_perm$fitted.values)
  fit_corr_perm <- lm(fit_value ~ gamlss_t ,
                  data = corrFit)
  r[i] <- summary(fit_corr_perm)$r.squared
  
}    
write.csv(res,file = 'data/receptor/t_map_perm_result.csv')
for(j in k) { 
  if(res_real[j] > 0){
    p_perm[j,1] <- sum(res[j,] > res_real[j])/perm_count
  }
  if(res_real[j] < 0){
    p_perm[j,1] <- sum(res[j,] < res_real[j])/perm_count
  }
  
}

adjusted_p_perm <- p.adjust(p_perm, method = "fdr")
                
  

