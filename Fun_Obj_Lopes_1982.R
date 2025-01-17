Nash = function(x){
  
  str <- x[1]
  pes <- x[2]
  crec <- x[3]
    
  Rsoloin = tuin * str #<--- reservat�rio do solo (zona aerada)
  Rsubin = (ebin * 2630 / (1-0.5^(1/kk))) / Ad #<--- reservat�rio subterr�neo (zona saturada)
  
  
  Qcal<-0           
  
  for(i in 1:length(P))
    
  {
    
    Es = P[i]*tuin^pes; #<--- calcula o escoamento superficial
    
    Er = tuin * Ep[i]; #<--- calcula a evapotranspira�ao real
    
    Rec = crec * tuin^4 * Rsoloin; #<--- calcula a recarga subterr�nea
    
    Eb = (1-0.5^(1/kk)) * Rsubin; #<--- calcula o escoamento b�sico
    
    Qcal[i] = (Es + Eb) * Ad / 2630; #<--- Vaz�o calculada
    
    Qbas = Eb * Ad / 2630; 
    
    Rsoloin = Rsoloin + P[i] - Es - Er - Rec; #<--- Atualiza o reservat�rio do solo
    
    Rsubin = Rsubin + Rec - Eb; #<--- Atualiza o reservat�rio subterr�neo
    
    tuin = Rsoloin / str; #<--- Atualiza o teor de umidade
    
  }
  
  Nash<- (1 -(sum(rowSums(Qcal - Qob)^2)/sum(rowSums(Qob - mean(Qob))^2))) #<--- Eficiencia de Nash-Sutclif
  
} 