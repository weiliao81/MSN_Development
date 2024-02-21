
library("gamlss.dist")
x <- read.csv('betweenDispersion.csv')
centile <- data.frame(
  age = seq.int(8, 89,  length.out = 1000)
)
between_gamlss_t <-array(0,dim=c(1,21))
between_gamlss_p <-array(0,dim=c(1,21))
quantiles <-c(0.5,0.75,0.25)
centile_names <-c("mid","high","low")
min_age <-min(x$age)
max_age <-max(x$age)

for (i in 1:21 ) {
################################################################################ step 1: Gamlss fitting
  betweenDispersion <- x[,i+3]
  
  m0 <- gamlss(betweenDispersion ~ cs(age, df = 5)+sex+eTIV, 
               sigma.fo = ~cs(age, df = 3)+sex+eTIV,
               family = GG,
               method = mixed(30, 100),
               data = x,
               trace = FALSE)
  fn <- function(p) AIC(gamlss(betweenDispersion ~ cs(age, df = p[1])+sex+eTIV, 
                               sigma.fo = ~cs(age,p[2])+sex+eTIV,
                               family = GG,
                               method = mixed(30, 100),
                               data = x, trace = FALSE, 
                               start.from = m0),
                        k = 2)
  opAIC <- optim(par = c(5, 3), 
                 fn, 
                 method = "L-BFGS-B", 
                 lower = c(1,1), 
                 upper = c(15, 15))

  
  command_gamlss <-paste0("model <- gamlss(betweenDispersion ~ cs(age, df = opAIC$par[1])+sex+eTIV, 
                                    sigma.fo = ~cs(age, df = opAIC$par[2])+sex+eTIV, 
                                    family = GG,
                                    method = mixed(30, 100),
                                    data = x, 
                                    trace = FALSE)")
  eval(parse(text = command_gamlss))
  
  resInfo <- summary(model)
  between_gamlss_t[i] = resInfo[2,3]
  between_gamlss_p[i] = resInfo[2,4]
  
################################################################################ step 2: Get quantiles
  for (j in 1:length(quantiles)){
    
    Qua <- getQuantile(model, quantile=quantiles[j],term="age")
    out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=FALSE,col="red",n = 1000) 
    centile_name = paste(centile_names[j], i, sep = "_") 
    centile[centile_name] <- out$y
  }
################################################################################ step 3: Plot Figure
  command_plot <-paste0(" p",i," <-ggplot(data = x)+
                              geom_point(aes(x = age,y = x[,",i+3,"]),shape=21, alpha = 0.7,size = 1.5)+
                              geom_point(centile, mapping = aes(x = age,y = low_",i,") )+
                              geom_point(centile, mapping = aes(x = age,y = mid_",i,")) +
                              geom_point(centile, mapping = aes(x = age,y = high_",i,")) +
                              theme_minimal() +
                              theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
                   ")
  eval(parse(text = command_plot))
  
}
################################################################################ step 4: Save Data
write.csv(centile, file = 'betweenDispersion_centile.csv', row.names = FALSE)
write.csv(between_gamlss_t, file = 'betweenDisTvalueGamlss.csv', row.names = FALSE)
write.csv(between_gamlss_p, file = 'betweenDisPvalueGamlss.csv', row.names = FALSE)

