library("lavaan")
library("readxl")
data <- read_excel("EF.xls")
X <- data$age
M1 <- data$priMotor
M2 <- data$Asso1
M3 <- data$Asso2
M4 <- data$secSensory
M5 <- data$priSensory
M6 <- data$limbic
M7 <- data$insular
Y <- data$EFC

Data <- data.frame(X = X, Y = Y, M1 = M1, M2 = M2,M3 = M3, M4 = M4,M5 = M5, M6 = M6,M7 = M7)
################################################################################ step 1: sem fitting
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M1 ~ a1*X
             M2 ~ a2*X
             M3 ~ a3*X
             M4 ~ a4*X
             M5 ~ a5*X
             M6 ~ a6*X
             M7 ~ a7*X
             Y ~ b1*M1 + b2*M2 + b3*M3+ b4*M4 + b5*M5 + b6*M6 + b7*M7
           # indirect effect (a*b)
             ab := a1*b1 + a2*b2+ a3*b3 + a4*b4 + a5*b5 + a6*b6 + a7*b7
             ab1 := a1*b1
             ab2 := a2*b2
             ab3 := a3*b3
             ab4 := a4*b4
             ab5 := a5*b5
             ab6 := a6*b6
             ab7 := a7*b7
             
           # total effect 
             total := c + (a1*b1 + a2*b2+ a3*b3 + a4*b4 + a5*b5 + a6*b6 + a7*b7)
         '
fit <- sem(model, data=Data,se = "bootstrap")
summary(fit)



