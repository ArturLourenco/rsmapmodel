###############################################################################################################################
################################################# Clear Console and Memory ####################################################
###############################################################################################################################

rm(list = ls())
cat("\014")

#========================================================================================
#=====     SMAP - Soil Moisture Accounting Procedure Model (Monthly Version)  ===========
#========================================================================================

#===========================================================================================================
#====  Copyright (c) 2013, Artur Gonçalves All rights reserved.=============================================
#====  Redistribution and use in source and binary forms, with or without ==================================
#====  modification, are permitted =========================================================================
#===========================================================================================================
#====   Reference: LOPES J.E.G., BRAGA B.P.F., CONEJO J.G.L. (1982), ======================================= 
#====   SMAP - A Simplified Hydrological Model, Applied Modelling in Catchment Hydrology, ed. V.P.Singh, ===
#====   Water Resourses Publications.   ==================================================================== 
#===========================================================================================================

#===========================================================================================================           
#=========================== Definindo diretorio de trabalho ===============================================
#===========================================================================================================

setwd("C:/Users/Artur/Documents/Universidade/Mestrado/Hidrologia II/SMAP/Trab_Final")

#===========================================================================================================           
#============================== Carregando Pacotes R.matlab & hydroGOF=======================================
#===========================================================================================================           

#obs.:Conjunto de dados está em um arquivo ".mat", mas poderia estar em ".txt",".csv" ou ".xls", por exemplo.

require(R.matlab)
require(hydroGOF)
require(hydromad)
require(mco)
require(XLConnect)
require(mopsocd)

#===========================================================================================================           
#================================ Carregando Conjunto de Dados==============================================
#===========================================================================================================

Dados <- readMat(file.path("C:/Users/Artur/Documents/Universidade/Mestrado/Hidrologia II/SMAP/Trab_final","dados_pianco.mat"))

#===========================================================================================================           
#================================ Definindos Dados de Entrada ==============================================
#===========================================================================================================           

P<- as.matrix(Dados[['Pm.4.postos.1963.1991']][140:324]) #<---- Serie de Chuvas Mensais 1974-1989
Dategen<-t(t(format(seq(as.Date("1974-08-01"), as.Date("1989-12-31"), 'months'), format="%Y-%m-%d", tz="UTC")));
#P1<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.8.postos.1963.1988']] #<---- Serie de Chuvas Mensais
#P<- as.matrix(Dados[['Pm.8.postos.1963.1988']][7:312]) #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.7.postos.1963.1988']] #<---- Serie de Chuvas Mensais
#Ep1<- as.matrix(Dados[['Ep.pianco.1963.1991']][133:324]) #<---- Serie de Evapotranspiração Potencial Mensal
Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][140:324]) #<---- Serie de Evapotranspiração Potencial Mensal 1974-1989
#Ep1<- Dados[['Ep.pianco.1963.1991']]#<---- Serie de Evapotranspiração Potencial Mensal
#Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][7:(348-36)]) #<---- Serie de Evapotranspiração Potencial Mensal 1963 - 1988
#Ep<- as.matrix(Dados[['Ep.boqueirao.1963.1991']][7:(348-36)]) #<---- Serie de Evapotranspiração Potencial Mensal 1963 - 1988
#Ep<- as.matrix(Dados[['Ep.boqueirao.1963.1991']][140:324]) #<---- Serie de Evapotranspiração Potencial Mensal 1974-1989
#Ep2<- Dados[['Ep.boqueirao.1963.1991']][1:(348-36)] #<---- Serie de Evapotranspiração Potencial Mensal 1963-1988
#Qob<- as.matrix(Dados[['Qob.1963.1991']][1:(348-36)]) #<--- Série de Vazões Médias Mensais 1963 - 1988
Qob<- as.matrix(Dados[['Qob.1963.1991']][140:324]) #<--- Série de Vazões Médias Mensais 1974 - 1989
#Qob<- as.matrix(Dados[['Qob.1963.1991']][7:(348-36)]) #<--- Série de Vazões Médias Mensais 1963 - 1988
#Qob<- Dados[['Qob.1963.1991']] #<--- Série de Vazões Médias Mensais 1963 - 1988
Ad<- Dados[['Ad']] #<--- Área de Drenagem da Bacia


#===========================================================================================================           
#=========================== Inicialização dos Parametros de Calibração Manual==============================
#=========================================================================================================== 

#obs.:Caso queria definir seus proprios parametros comentar os definidos e descomentar os manuais.

#print('Defina: 500< str < 4000 ='); str=scan()
#print('Defina: 0.01 < crec <  0.64 ='); crec=scan()
#print('Defina: 0.5 < kk < 0.9 ='); kk=scan()
#print('Defina: 1< pes < 10 ='); pes=scan()
#kk<- Dados[['kk']] #<--- Constante de recessão do escoamento básico (mês^-1) = 1 < kk < 6


#===========================================================================================================           
#============================= Inicialização dos Parametros de Estado=======================================
#===========================================================================================================            

tuin<- 0.2 #<--- Teor de umidade inicial (ad.)
ebin<- 0.002 #<--- Vazão básica inicial (m3/s)      

#===========================================================================================================           
#====================== Inicialização dos Parametros de Calibração Automatica===============================
#=========================================================================================================== 

kk<- 0.5 #<--- Constante de recessão do escoamento básico (mês^-1) = 1 < kk < 6 

source('Fun_Mult_Obj_Lopes_1982.R')

varcount <- 3
fncount <- 2
lbound <- c(500,1,0.01)
ubound <- c(4000,10,0.64)
optmin <- 1

Resultoptm <- mopsocd(Nash,varcnt=varcount,fncnt=fncount,
               lowerbound=lbound,upperbound=ubound,opt=optmin,popsize=100,maxgen=25,archivesize=500)

# Resultoptm <- nsga2(Nash, 3, 2,
#           generations=150, popsize=100,
#           cprob=0.7, cdist=20,
#           mprob=0.2, mdist=20,
#           lower.bounds=lbound,
#           upper.bounds=ubound)

str<- Resultoptm$paramvalues[Resultoptm$numsols,1]
pes<- Resultoptm$paramvalues[Resultoptm$numsols,2]
crec<- Resultoptm$paramvalues[Resultoptm$numsols,3]


#===========================================================================================================           
#================================== Inicializando os Reservatorios==========================================
#=========================================================================================================== 


Rsoloin = tuin * str #<--- reservatório do solo (zona aerada)

Rsubin = (ebin * 2630 / (1-0.5^(1/kk))) / Ad #<--- reservatório subterrâneo (zona saturada)

#===========================================================================================================           
#======================================== Modelagem no SMAP=================================================
#===========================================================================================================           
Qcal<-0           

for(i in 1:length(P))
  
{
  
  Es = P[i]*tuin^pes; #<--- calcula o escoamento superficial
  
  Er = tuin * Ep[i]; #<--- calcula a evapotranspiraçao real
  
  Rec = crec * tuin^4 * Rsoloin; #<--- calcula a recarga subterrânea
  
  Eb = (1-0.5^(1/kk)) * Rsubin; #<--- calcula o escoamento básico
  
  Qcal[i] = (Es + Eb) * Ad / 2630; #<--- Vazão calculada
  
  Qbas = Eb * Ad / 2630; 
  
  Rsoloin = Rsoloin + P[i] - Es - Er - Rec; #<--- Atualiza o reservatório do solo
  
  Rsubin = Rsubin + Rec - Eb; #<--- Atualiza o reservatório subterrâneo
  
  tuin = Rsoloin / str; #<--- Atualiza o teor de umidade
  
}

#===========================================================================================================           
#=========================================== Resultados=====================================================
#===========================================================================================================

#Gráfico para da vazão calculada e vazão observada
plot(1:length(P),Qcal, ann = FALSE, type = "l", col = "black")
lines(1:length(P),Qob, ann = FALSE, type = "l", col = "red")
title(main = "Resultados para o Posto Piancó - Piancó Ano 1963/2008",
      xlab = "Tempo (Meses)",ylab = "Vazão (m^3/s)",
      col.main = "blue",cex.main = 1.2, cex.lab = 1.0, 
      font.main = 12, font.lab = 7)
legend(12,320,c("Vazão Observada","Vazão Calculada"),col=c("red","black"),lty=1)

#Gráfico de Regressão
plot(Qob,Qcal, ann = FALSE, type = "p", col = "black")
reg1 <- lm(Qob~Qcal)
abline(reg1)
title(main = "Gráfico de Regressão",
      xlab = "Vazão (m^3/s)",ylab = "Vazão (m^3/s)",
      col.main = "blue",cex.main = 1.2, cex.lab = 1.0, 
      font.main = 12, font.lab = 7)

#===========================================================================================================           
#=================================Indice Estatisticos e Eficiencia==========================================
#===========================================================================================================

r<- cor(Qob,Qcal) #<--- Coeficiente de Correlação de Pearson
Bias<- ((mean(Qcal)-mean(Qob))/mean(Qob))*100 #<--- Viés Relativo
r2.lm = lm(Qob ~ Qcal) 
r2<- summary(r2.lm)$r.squared #<--- Coeficiente de Determinaçao
Nash<- 1 -(sum(rowSums(Qcal - Qob)^2)/sum(rowSums(Qob - mean(Qob))^2)) #<--- Eficiencia de Nash-Sutclif

#===========================================================================================================           
#======================================== Goodness-of-fit ==================================================
#===========================================================================================================

Qcal<-as.numeric(Qcal)
Qobs<-as.numeric(Qob)
gof(Qcal,Qob)
Eff<-gof(Qcal,Qob)                     

#===========================================================================================================           
#======================================== Goodness-of-fit II================================================
#===========================================================================================================

# Calculando a eficiencia substituindo os valores de Qcal por zero onde a Qob é zero

for(i in 1:length(Qcal))
  
{
  
  if  (Qobs[i]==0)
    
    Qcal[i]<-Qobs[i] else
      
      Qcal[i]<-Qcal[i]
  
}

Eff2<-gof(Qcal,Qobs)


#===========================================================================================================           
#======================================== Export to Excel ==================================================
#===========================================================================================================

#output1 = writeWorksheetToFile("Smap_Lopes.xlsx",data = list( Qcal = Qcal, Qob = Qob), 
#                               sheet = "Sheet1", startRow = 1, 
#                               startCol = c(2,3), header = FALSE) 

#a<- readWorksheetFromFile("Trab-Final-Smap_mensal_ultimo_preenchida_seriecompletavazao.xlsx", sheet = "Sheet1", 
#                          header = FALSE, startCol = 5, 
#                          startRow = 2, endCol = 5, 
#                          endRow = 349)

