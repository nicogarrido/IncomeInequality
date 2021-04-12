library(readxl)
library(dplyr)
#library(tidyverse)
library(foreign)



Datos <- read_excel("SAM 2016.xlsx", range="SAM 2016_final_FE Simu2!C4:AQ44", col_names=FALSE)
D <- as.matrix(Datos)
D[is.na(D)]<-0


#Ordenar matriz: Sacar flia y columna de gobierno

Dnew = D * 0 

totCtas = 41
ctas = 1:totCtas
cExo = c(33,41)
cEndo = setdiff(ctas, cExo)
nCtasEx = length(cExo)
nCtasEn = length(cEndo)

Dnew[1:nCtasEn,1:nCtasEn] = D[cEndo,cEndo]
for (i in 1:nCtasEx){
  Dnew[(nCtasEn+i),1:nCtasEn] = D[cExo[i], cEndo]
  Dnew[(nCtasEn+i),((nCtasEn+1):totCtas)] = D[cExo[i], cExo]
  Dnew[1:nCtasEn,(nCtasEn+i)] = D[cEndo,cExo[i]]
  Dnew[((nCtasEn+1):totCtas),(nCtasEn+i)] = D[cExo,cExo[i]]
}

xf <- colSums(Dnew)

AnularColumna = which(xf==0) #para que no se indetermina la columna de margenes, la dejo en 0 directamente.


#Matriz A de coef t?cnicos

A <- Dnew %*% diag(as.vector(1/xf)) #Esto es un vector diagonalizado 1/xf, en donde el vector xf es un vector de los valores de la suma de cada fila o an?logamente de cada columna.
A[,AnularColumna] = 0

Ann = as.matrix(A[1:nCtasEn,1:nCtasEn])
Ank = as.matrix(A[1:nCtasEn,(nCtasEn+1):totCtas])
Akn = as.matrix(A[(nCtasEn+1):totCtas,1:nCtasEn])
Akk = as.matrix(A[(nCtasEn+1):totCtas,(nCtasEn+1):totCtas])
yn = as.vector(xf[1:nCtasEn])
yk = as.vector(xf[(nCtasEn+1):totCtas])


#Efectos totales

if (abs(sum(yn - (Ann %*% yn + Ank %*% yk)))<0.000001){
  print("Test de consistencia 1 OK")
}

Ma = solve(diag(nCtasEn)-Ann)
x = Ank %*% yk

if (abs(sum(Ma %*% x - (Ann %*% yn + Ank %*% yk)))<0.000001){
  print("Test de consistencia 2 OK")
}

##################################################################
Ec2<-Ma%*%x #Ecuaci?n (2), ecuaci?n principal del paper


#Matriz redistributiva

e = matrix(1, nCtasEn, 1)

zn = yn / as.vector((t(e) %*% yn))

R1 = 1 / as.vector((t(e) %*% yn))
R2 = diag(nCtasEn) - zn%*%t(e)

R = R1 * R2 %*% Ma


#Impacto


#Efecto redistributivo-> Cambio en la posici?n relativa
dx1 = c(1,0)
dx = Ank %*% dx1
impacto = R %*% dx


#Efectos de difusi?n del aumento en una unidad de la cuenta ex?gena
dYn <- Ma%*%dx


#Diffusion_effects<-data.frame(Diffusion_effects)
#Como exportar a stata: write.dta(Diffusion_effects, "Deff.dta")

dYn<-data.frame(dYn)
write.dta(dYn, "efecto difusion.dta")

impacto<-data.frame(impacto)
#write.dta(impacto, "impacto.dta")

barplot(R[,2])

############################     Graficos     ######################################################















