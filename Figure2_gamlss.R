library("gamlss")
library("mgcv")
x <- read.csv('global dispersion.csv')
centile <- data.frame(
  age = seq.int(8, 89,  length.out = 1000)
)
quantiles <-c(0.25,0.5,0.75)
centile_names <-c("low","mid","high")
min_age <-min(x$age)
max_age <-max(x$age)
################################################################################ step 1: Gamlss fitting
globalDispersion <- x[,4]

m0 <- gamlss(globalDispersion ~ cs(age, df = 5)+sex+eTIV, 
             sigma.fo = ~cs(age, df = 3)+sex+eTIV,
             family = GG,
             method = mixed(30, 100),
             data = x,
             trace = FALSE)
fn <- function(p) AIC(gamlss(globalDispersion ~ cs(age, df = p[1])+sex+eTIV, 
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


command_gamlss <-paste0("model <- gamlss(globalDispersion ~ cs(age, df = opAIC$par[1])+sex+eTIV, 
                                    sigma.fo = ~cs(age, df = opAIC$par[2])+sex+eTIV, 
                                    family = GG,
                                    method = mixed(30, 100),
                                    data = x, 
                                    trace = FALSE)")
eval(parse(text = command_gamlss))

resInfo <- summary(model,type = "qr")
global_gamlss_t = resInfo[2,3]
global_gamlss_p = resInfo[2,4]

################################################################################ step 2: Get quantiles
for (j in 1:length(quantiles)){
  
  Qua <- getQuantile(model, quantile=quantiles[j],term="age")
  out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=FALSE,col="red",n = 1000) 
  centile_name = centile_names[j] 
  centile[centile_name] <- out$y
}
################################################################################ step 3: Plot Figure
command_plot <-paste0(" p <-ggplot(data = x)+
                              geom_point(aes(x = age,y = x[,4]),shape=21, alpha = 0.7,size = 1.5)+
                              geom_point(centile, mapping = aes(x = age,y = low) )+
                              geom_point(centile, mapping = aes(x = age,y = mid)) +
                              geom_point(centile, mapping = aes(x = age,y = high)) +
                              theme_minimal() +
                              theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
                   ")
eval(parse(text = command_plot))

################################################################################ step 4: Save Data
write.csv(centile, file = 'filename.csv', row.names = FALSE)


################################################################################ Figure2B
x <- read.csv('whole brain regions' dispersion.csv')
res_t <- array(0, dim = c(1533, 1))
res_p <- array(0, dim = c(1533, 1))
u <- 4:1536
age <-x$age
sex <- x$gender
eTIV <- x$eTIV
df <- data.frame(regionGamlss = x[, 4], age = age, sex  = sex ,eTIV = eTIV)
################################################################################ step 1: Gamlss fitting
for(i in u) {  
  df[,1] <- x[,i]
  m0 <- gamlss(regionGamlss ~ cs(age, df = 5 )+sex+eTIV, 
               sigma.fo = ~cs(age, df = 3)+sex+eTIV,
               family = GG,
               method = mixed(30, 100),
               data = df, 
               trace = FALSE)
  fn <- function(p) AIC(gamlss(regionGamlss ~ cs(age, df = p[1])+sex+eTIV, 
                               sigma.fo = ~cs(age, p[2])+sex+eTIV, 
                               data = df,
                               family = GG,
                               method = mixed(30, 100), 
                               trace = FALSE, 
                               start.from = m0), 
                        k = 2)
  
  opAIC <- optim(par = c(5, 3), 
                 fn, 
                 method = "L-BFGS-B", 
                 lower = c(1,1), 
                 upper = c(15, 15))
  
  m1 <- gamlss(regionGamlss ~ cs(age, df = opAIC$par[1])+sex+eTIV, 
               sigma.fo = ~cs(age, df = opAIC$par[2])+sex+eTIV,
               family = GG,
               method = mixed(30, 100),
               data = df, 
               trace = FALSE)
  resInfo <- summary(m1,type = "qr")
################################################################################ step 2: Save Data
  res_t[i-3, 1] <-resInfo[2,3]
  res_p[i-3, 1] <- resInfo[2,4]
} 
write.csv(res_t, file = 'whole brain regions' gamlssTvalues.csv', row.names = FALSE)     






