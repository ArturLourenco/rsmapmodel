Nash = function(x){
  
  str <- x[1]
  e2 <- x[2]
  e1 <- x[3]
  cinf <- x[4]
  
  Rsoloin = tuin * str #<--- reservatório do solo (zona aerada)
    
  Qcal<-0           
  
  for(i in 1:length(P))
    
  {
    
    Es = P[i]*tuin^e2; #<--- calcula o escoamento superficial
    
    Er = tuin * Ep[i]^e1; #<--- calcula a evapotranspiraçao real
      
    Qcal[i] = (Es) * Ad / 2630; #<--- Vazão calculada
      
    Rsoloin = Rsoloin + P[i] - Es - Er; #<--- Atualiza o reservatório do solo
    
    tuin = (Rsoloin + cinf*P[i]) / str; #<--- Atualiza o teor de umidade
    
  }
  
  fn2<- (((length(P)*sum(rowSums(Qcal * Qob)))-(sum(Qob) * sum(Qcal)))/((sqrt(length(P)*(sum(Qob^2)) - 
                      (sum(Qob)^2)))*(sqrt(length(P)*(sum(Qcal^2)) - (sum(Qcal)^2)))))
  fn1<- (1 -(sum(rowSums(Qcal - Qob)^2)/sum(rowSums(Qob - mean(Qob))^2))) #<--- Eficiencia de Nash-Sutclif
    return(c(fn1,fn2))
  
} 