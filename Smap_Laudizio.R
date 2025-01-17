#========================================================================================
#=====     SMAP - Soil Moisture Accounting Procedure Model (Monthly Version)  ===========
#========================================================================================

#===========================================================================================================
#====  Copyright (c) 2013, Artur Gon�alves All rights reserved.=============================================
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
#============================== Carregando Pacote R.matlab & hydroGOF=======================================
#===========================================================================================================           

#obs.:Conjunto de dados est� em um arquivo do .mat, mas poderia estar em ".txt",".csv" ou ".xls", por exemplo.

require(R.matlab)
require(hydroGOF)
require(hydromad)
require(XLConnect)

#===========================================================================================================           
#================================ Carregando Conjunto de Dados==============================================
#===========================================================================================================

Dados <- readMat(file.path("C:/Users/Artur/Documents/Universidade/Mestrado/Hidrologia II/SMAP/Trab_final","dados_pianco.mat"))

#===========================================================================================================           
#================================ Definindos Dados de Entrada ==============================================
#===========================================================================================================           

#P<- as.matrix(Dados[['Pm.4.postos.1963.1991']][133:324]) #<---- Serie de Chuvas Mensais 1974-1989
#P1<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P<- Dados[['Pm.4.postos.1963.1991']] #<---- Serie de Chuvas Mensais
#P2<- Dados[['Pm.8.postos.1963.1988']] #<---- Serie de Chuvas Mensais
P<- Dados[['Pm.7.postos.1963.1988']] #<---- Serie de Chuvas Mensais
#Ep1<- as.matrix(Dados[['Ep.pianco.1963.1991']][133:324]) #<---- Serie de Evapotranspira��o Potencial Mensal
#Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][133:324]) #<---- Serie de Evapotranspira��o Potencial Mensal
#Ep1<- Dados[['Ep.pianco.1963.1991']]#<---- Serie de Evapotranspira��o Potencial Mensal
Ep<- as.matrix(Dados[['Ep.pianco.1963.1991']][1:(348-36)]) #<---- Serie de Evapotranspira��o Potencial Mensal 1963 - 1988
#Ep2<- Dados[['Ep.boqueirao.1963.1991']] #<---- Serie de Evapotranspira��o Potencial Mensal
#Ep2<- Dados[['Ep.boqueirao.1963.1991']][1:(348-36)] #<---- Serie de Evapotranspira��o Potencial Mensal 1963-1988
Qob<- as.matrix(Dados[['Qob.1963.1991']][1:(348-36)]) #<--- S�rie de Vaz�es M�dias Mensais 1963 - 1988
#Qob<- as.matrix(Dados[['Qob.1963.1991']][133:324]) #<--- S�rie de Vaz�es M�dias Mensais 1974 - 1989
#Qob<- Dados[['Qob.1963.1991']] #<--- S�rie de Vaz�es M�dias Mensais 1963 - 1988
Ad<- Dados[['Ad']] #<--- �rea de Drenagem da Bacia


#===========================================================================================================           
#=========================== Inicializa��o dos Parametros de Calibra��o Manual==============================
#=========================================================================================================== 

#obs.:Caso queria definir seus proprios parametros comentar os definidos e descomentar os manuais.

#print('Defina: 500< str < 4000 ='); str=scan()
#print('Defina: 0.01 < crec <  0.64 ='); crec=scan()
#print('Defina: 1 < kk < 6 ='); kk=scan()
#print('Defina: 1< pes < 10 ='); pes=scan()
crec<- 0
kk<- Dados[['kk']] #<--- Constante de recess�o do escoamento b�sico (m�s^-1) = 1 < kk < 6

#===========================================================================================================           
#============================= Inicializa��o dos Parametros de Estado=======================================
#===========================================================================================================            

tuin<- 0 #<--- Teor de umidade inicial (ad.)
ebin<- 0 #<--- Vaz�o b�sica inicial (m3/s)      

#===========================================================================================================           
#====================== Inicializa��o dos Parametros de Calibra��o Automatica===============================
#=========================================================================================================== 

source('smapfun2.R')
Resultoptm <- SCEoptim(Nash,c(1000,5,0.3,0.3),lower = c(400,0.1,0,0),upper = c(5000,10,Inf,Inf),
                       control = list(trace = 1,fnscale = -1, ncomplex = 5, elitism = 2))

str<- Resultoptm[[3]][1]
e2<- Resultoptm[[3]][2]
e1<- Resultoptm[[3]][3]
cinf<- Resultoptm[[3]][4]


#===========================================================================================================           
#================================== Inicializando os Reservatorios==========================================
#=========================================================================================================== 


Rsoloin = tuin * str #<--- reservat�rio do solo (zona aerada)

Rsubin = (ebin * 2630 / (1-0.5^(1/kk))) / Ad #<--- reservat�rio subterr�neo (zona saturada)

#===========================================================================================================           
#======================================== Modelagem no SMAP=================================================
#===========================================================================================================           
Qcal<-0           

for(i in 1:length(P))
  
{
  
  Es = P[i]*tuin^e2; #<--- calcula o escoamento superficial
  
  Er = tuin * Ep[i]^e1; #<--- calcula a evapotranspira�ao real
  
  Rec = crec * tuin^4 * Rsoloin; #<--- calcula a recarga subterr�nea
  
  Eb = (1-0.5^(1/kk)) * Rsubin; #<--- calcula o escoamento b�sico
  
  Qcal[i] = (Es + Eb) * Ad / 2630; #<--- Vaz�o calculada
  
  Qbas = Eb * Ad / 2630; 
  
  Rsoloin = Rsoloin + P[i] - Es - Er - Rec; #<--- Atualiza o reservat�rio do solo
  
  Rsubin = Rsubin + Rec - Eb; #<--- Atualiza o reservat�rio subterr�neo
  
  tuin = (Rsoloin + cinf*P[i]) / str; #<--- Atualiza o teor de umidade
  
}

#===========================================================================================================           
#=========================================== Resultados=====================================================
#===========================================================================================================

#Gr�fico para da vaz�o calculada e vaz�o observada
plot(1:length(P),Qcal, ann = FALSE, type = "l", col = "black")
lines(1:length(P),Qob, ann = FALSE, type = "l", col = "red")
title(main = "Resultados para o Posto Pianc� - Pianc� Ano 1963/2008",
      xlab = "Tempo (Meses)",ylab = "Vaz�o (m^3/s)",
      col.main = "blue",cex.main = 1.2, cex.lab = 1.0, 
      font.main = 12, font.lab = 7)
legend(45,3500,c("Vaz�o Observada","Vaz�o Calculada"),col=c("black","red"),lty=1)

#Gr�fico de Regress�o
plot(Qob,Qcal, ann = FALSE, type = "p", col = "black")
reg1 <- lm(Qob~Qcal)
abline(reg1)
title(main = "Gr�fico de Regress�o",
      xlab = "Vaz�o (m^3/s)",ylab = "Vaz�o (m^3/s)",
      col.main = "blue",cex.main = 1.2, cex.lab = 1.0, 
      font.main = 12, font.lab = 7)

#===========================================================================================================           
#=================================Indice Estatisticos e Eficiencia==========================================
#===========================================================================================================

r<- cor(Qob,Qcal) #<--- Coeficiente de Correla��o de Pearson
Bias<- ((mean(Qcal)-mean(Qob))/mean(Qob))*100 #<--- Vi�s Relativo
r2.lm = lm(Qob ~ Qcal) 
r2<- summary(r2.lm)$r.squared #<--- Coeficiente de Determina�ao
Nash<- 1 -(sum(rowSums(Qcal - Qob)^2)/sum(rowSums(Qob - mean(Qob))^2)) #<--- Eficiencia de Nash-Sutclif

#===========================================================================================================           
#======================================== Goodness-of-fit ==================================================
#===========================================================================================================

Qcal<-as.numeric(Qcal)
Qob<-as.numeric(Qob)
gof(Qcal,Qob)
Eff<-gof(Qcal,Qob)                     

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

output1 = writeWorksheetToFile("plots.xlsx",data = list( Qcal = Qcal, Qob = Qob), 
                               sheet = "Sheet1", startRow = 1, 
                               startCol = c(2,3), header = FALSE)
