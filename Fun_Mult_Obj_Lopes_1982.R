Nash = function(x){
  
  
  str <- x[1]
  pes <- x[2]
  crec <- x[3]
  
  Rsoloin = tuin * str #<--- reservatório do solo (zona aerada)
  Rsubin = (ebin * 2630 / (1-0.5^(1/kk))) / Ad #<--- reservatório subterrâneo (zona saturada)
  
  
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
  #fn2<- ((length(P)*sum(rowSums(Qcal * Qob)))-(sum(Qob) * sum(Qcal)))/((sqrt(length(P)*(sum(Qob^2)) - (sum(Qob)^2)))*(sqrt(length(P)*(sum(Qcal^2)) - (sum(Qcal)^2))))
  fn2<- mse(as.numeric(Qcal),as.numeric(Qob))*-1
  #fn1<- (1 -(sum(rowSums(Qcal - Qob)^2)/sum(rowSums(Qob - mean(Qob))^2))) #<--- Eficiencia de Nash-Sutclif
  fn1<- mNSE(as.numeric(Qcal),as.numeric(Qob))
  return(c(fn1,fn2))  
} 