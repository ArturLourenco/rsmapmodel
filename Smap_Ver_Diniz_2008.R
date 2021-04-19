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
#====   Reference: Diniz (2008),============================================================================ 
#====   SMAP - Versão modificada para regiões semiaridas brasileiras========================================
#====   UFRG================================================================================================ 
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
require(XLConnect)
require(hydroPSO)
#===========================================================================================================           
#================================ Carregando Conjunto de Dados==============================================
#===========================================================================================================

Dados <- readMat(file.path("C:/Users/Artur/Documents/Universidade/Mestrado/Hidrologia II/SMAP/Trab_final","dados_pianco.mat"))

#===========================================================================================================           
#================================ Definindos Dados de Entrada ==============================================
#===========================================================================================================           

P<- as.matrix(Dados[['Pm.4.postos.1963.1991']][140:324]) #<---- Serie de Chuvas Mensais 1974-1989
Dategen<-t(t(format(seq(as.Date("1974-08-01"), as.Date("1989-12-31"), 'months'), format="%Y-%m-%d", tz="UTC")));
#P<- as.matrix(Dados[['Pm.8.postos.1963.1988']][7:312]) #<---- Serie de Chuvas Mensais
#P1<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P2<- Dados[['Pm.8.postos.1963.1988']] #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.7.postos.1963.1988']] #<---- Serie de Chuvas Mensais
#Ep1<- as.matrix(Dados[['Ep.pianco.1963.1991']][133:324]) #<---- Serie de Evapotranspiração Potencial Mensal
#Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][7:(348-36)]) #<---- Serie de Evapotranspiração Potencial Mensal 1963 - 1988
Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][140:324]) #<---- Serie de Evapotranspiração Potencial Mensal 1974-1989
#Ep<- as.matrix(Dados[['Ep.boqueirao.1963.1991']][7:(348-36)]) #<---- Serie de Evapotranspiração Potencial Mensal 1963 - 1988
#Ep1<- Dados[['Ep.pianco.1963.1991']]#<---- Serie de Evapotranspiração Potencial Mensal
#Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][1:(348-36)]) #<---- Serie de Evapotranspiração Potencial Mensal 1963 - 1988
#Ep2<- Dados[['Ep.boqueirao.1963.1991']] #<---- Serie de Evapotranspiração Potencial Mensal
#Ep2<- Dados[['Ep.boqueirao.1963.1991']][1:(348-36)] #<---- Serie de Evapotranspiração Potencial Mensal 1963-1988
#Qob<- as.matrix(Dados[['Qob.1963.1991']][7:(348-36)]) #<--- Série de Vazões Médias Mensais 1963 - 1988
Qob<- as.matrix(Dados[['Qob.1963.1991']][140:324]) #<--- Série de Vazões Médias Mensais 1974 - 1989
#Qob<- Dados[['Qob.1963.1991']] #<--- Série de Vazões Médias Mensais 1963 - 1988
Ad<- Dados[['Ad']] #<--- Área de Drenagem da Bacia


#===========================================================================================================           
#=========================== Inicialização dos Parametros de Calibração Manual==============================
#=========================================================================================================== 

#obs.:Caso queria definir seus proprios parametros comentar os definidos e descomentar os manuais.

#print('Defina: 40< str < 8000 ='); str=scan()
#print('Defina: 0 < e1 < 20 ='); pes=scan()
#print('Defina: 0< e2 < 15 ='); pes=scan()
#print('Defina: 0 < cinf < 0.99 ='); pes=scan()


#===========================================================================================================           
#============================= Inicialização dos Parametros de Estado=======================================
#===========================================================================================================            

tuin<- 0.2 #<--- Teor de umidade inicial (ad.)
ebin<- 0.002 #<--- Vazão básica inicial (m3/s)      

#===========================================================================================================           
#====================== Inicialização dos Parametros de Calibração Automatica===============================
#=========================================================================================================== 


source('Fun_Obj_Diniz_1994.R')
Resultoptm <- SCEoptim(Nash,c(2000,2,2,2),lower = c(40,0.1,0,0),upper = c(8000,20,15,0.99),
                      control = list(trace = 1,fnscale = -1, ncomplex = 5, elitism = 2))

# Resultoptm2<-hydroPSO(fn=Nash, lower=c(40,0.1,0,0), upper=c(8000,20,15,0.99),
#                       control=list(topology="vonNeumann", reltol=1E-20, normalise=TRUE,REPORT=50, write2disk=FALSE,MinMax=('max')));

str<- Resultoptm[[3]][1]
e2<- Resultoptm[[3]][2]
e1<- Resultoptm[[3]][3]
cinf<- Resultoptm[[3]][4]

# str<- Resultoptm2[[1]][1]
# e2<- Resultoptm2[[1]][2]
# e1<- Resultoptm2[[1]][3]
# cinf<- Resultoptm2[[1]][4]


#===========================================================================================================           
#================================== Inicializando os Reservatorios==========================================
#=========================================================================================================== 


Rsoloin = tuin * str #<--- reservatório do solo (zona aerada)

#===========================================================================================================           
#======================================== Modelagem no SMAP=================================================
#===========================================================================================================           
Qcal<-0           

for(i in 1:length(P))
  
{
  
  Es = P[i]*tuin^e2; #<--- calcula o escoamento superficial
  
  Er = tuin * Ep[i]^e1; #<--- calcula a evapotranspiraçao real
    
  Qcal[i] = (Es) * Ad / 2630; #<--- Vazão calculada
  
  Rsoloin = Rsoloin + P[i] - Es - Er; #<--- Atualiza o reservatório do solo
    
  tuin = (Rsoloin + cinf*P[i]) / str; #<--- Atualiza o teor de umidade
  
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
legend(45,3500,c("Vazão Observada","Vazão Calculada"),col=c("black","red"),lty=1)

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
Qob<-as.numeric(Qob)
gof(Qcal,Qob)
Eff<-gof(Qcal,Qob)

#===========================================================================================================           
#======================================== Goodness-of-fit II================================================
#===========================================================================================================

# Calculando a eficiencia substituindo os valores de Qcal por zero onde a Qob é zero

for(i in 1:length(Qcal))
  
{
  
  if  (Qob[i]==0)
    
    Qcal[i]<-Qob[i] else
      
      Qcal[i]<-Qcal[i]
  
}

Eff2<-gof(Qcal,Qob)


#===========================================================================================================           
#======================================== Export to Excel ==================================================
#===========================================================================================================

#output1 = writeWorksheetToFile("Smap_Diniz.xlsx",data = list( Qcal = Qcal, Qob = Qob), 
#                               sheet = "Sheet1", startRow = 1, 
#                               startCol = c(2,3), header = FALSE)

sim<-zoo(Qcal,Dategen);
obs<-zoo(Qob,Dategen);

ggof(sim=Qcal, obs=Qob, dates = as.character(Dategen), main = "Vazão Observada vs. Simulada",xlab = "Tempo", ylab=c("Q, [m3/s]"))

#write.table(data,file="data.csv",dec=",",sep=";",row.names=F)