
writeLines("

model{
  
  # this controls for studies with one outcome not reported by setting the correlation equal to 
  # zero
  
  for (k in 1:(2*Ns-N2h)){
  
  control[k] <- step(9999-var1[k])*step(9999-var2[k])
  
  }
  
  # two-arm studies
  
  for(i in 1:N2h){
    
    s[i,1,1] <- var1[i]+psi1.sq
    
    s[i,2,2] <- var2[i]+psi2.sq
    
    s[i,1,2] <- control[i]*rho1*sqrt(var1[i]+psi1.sq)*sqrt(var2[i]+psi2.sq)
    
    s[i,2,1] <- s[i,1,2]
    
    prec2A[i,1:2,1:2] <- inverse(s[i,,])
    
    y[(2*i-1):(2*i)] ~ dmnorm(mean[(2*i-1):(2*i)],prec2A[i,,])
    
    }
  
  # three-arm studies
  
  for (i in 1:(Ns-N2h)){
    S[i,1,1]<- var1[N2h+2*i-1]+psi1.sq

    S[i,2,2] <- var2[N2h+2*i-1]+psi2.sq

    S[i,3,3] <- var1[N2h+2*i]+psi1.sq

    S[i,4,4] <- var2[N2h+2*i]+psi2.sq

    S[i,1,2] <- control[i]*rho1*sqrt(S[i,1,1])*sqrt(S[i,2,2])

    S[i,2,1] <- S[i,1,2]

    S[i,1,3] <- control[i]*sqrt(S[i,1,1])*sqrt(S[i,3,3])/2

    S[i,3,1] <- S[i,1,3]

    S[i,1,4] <- control[i]*rho1*sqrt(S[i,1,1])*sqrt(S[i,4,4])/2

    S[i,4,1] <- S[i,1,4]

    S[i,2,3] <- control[i]*rho1*sqrt(S[i,2,2])*sqrt(S[i,3,3])/2

    S[i,3,2] <- S[i,2,3]

    S[i,2,4] <- control[i]*sqrt(S[i,2,2])*sqrt(S[i,4,4])/2

    S[i,4,2] <- S[i,2,4]

    S[i,4,3] <- control[i]*rho1*sqrt(S[i,4,4])*sqrt(S[i,3,3])

    S[i,3,4] <- S[i,4,3]
  }

  for (k in 1:(Ns-N2h)){

    prec3A[k,1:4,1:4] <- inverse(S[k,,])

    y[(2*N2h+4*k-3):(2*N2h+4*k)] ~ dmnorm(mean[(2*N2h+4*k-3):(2*N2h+4*k)],prec3A[k,,])

    }
  
  # Parameterization of the means
  
  for(i in 1:N2h) { 
    
    mean[2*i-1] <- d1[T2[i]] - d1[T1[i]]
    
    mean[2*i] <- d2[T2[i]] -  d2[T1[i]]
    
    }
  
  for(i in 1:(Ns-N2h)) {
    
    mean[2*N2h+4*i-3] <- d1[T2[N2h+i]] - d1[T1[N2h+i]]

    mean[2*N2h+4*i-2] <- d2[T2[N2h+i]] - d2[T1[N2h+i]]
   
    mean[2*N2h+4*i-1] <- d1[T3[N2h+i]] - d1[T1[N2h+i]]

    mean[2*N2h+4*i] <- d2[T3[N2h+i]] - d2[T1[N2h+i]]
  #  
   }
  
  # Priors		

  for(k in 1:(ref-1)) {	
  
  d1[k] ~ dnorm(0,1e-03)
  
  }
  
  for(k in (ref+1):NT) {
  
  d1[k] ~ dnorm(0,1e-03)
  
  }		
  
  for(k in 1:(ref-1)) {
  
  d2[k] ~ dnorm(0,1e-03)
  
  }
  
  for(k in (ref+1):NT) {
  
  d2[k] ~ dnorm(0,1e-03)
  
  }		
  
  psi1.sq <- psi1*psi1
  
  #psi1 ~ dunif(0,1)
  
  psi1 ~ dnorm(0,1)T(0,)
  
  psi2.sq <- psi2*psi2
  
  psi2 ~ dnorm(0,1)T(0,)
  
  #psi2 ~ dunif(0,1)
  
  rho1 ~ dunif(lb.rho,ub.rho)
  
  #Estimated  Effect Sizes
  
  d1[ref]<- 0
  
  for (c in 1:(ref-1)) { 
  
  res.ref1[c] <- d1[c] - d1[ref] 
  
  } 
  
  for (c in (ref+1):NT) { 
  
  res.ref1[c] <- d1[c] - d1[ref]	
  
  } 
  
  for (c in 1:(NT-1)) {
  
    for (k in (c+1):NT) { 
      
     res.1[c,k] <- d1[k] - d1[c]
    }
      }
  
  d2[ref]<- 0
  
  for (c in 1:(ref-1)) { 
  
  res.ref2[c] <- d2[c] - d2[ref]
  
  } 
  
  for (c in (ref+1):NT) { 
  
  res.ref2[c] <- d2[c] - d2[ref]	
  
  }
  
  for (c in 1:(NT-1)) {
  
    for (k in (c+1):NT) {
    
      res.2[c,k] <- d2[k] - d2[c] 
    }
      }	
  
  # # SUCRA rankings
  # 
  # # Ranking of treatments for response. This part of the code is adjusted for the acute mania
  # # dataset 
  # for (k in 1:13){
  # dd1[k] <- d1[k]
  # }
  # 
  # for(k in 1:13) {
  # 
  # order1[k]<- 13- rank(dd1[],k)
  # 
  # most.effective1[k]<-equals(order1[k],1)
  #   
  #   for(j in 1: 13) {
  #     
  #     effectiveness1[k,j]<- equals(order1[k],j)
  #     
  #     cumeffectiveness1[k,j]<- sum(effectiveness1[k,1:j])}}		
  # 
  # for(k in 1:13) {
  # SUCRA2[k] <- sum(cumeffectiveness1[k,1:(13-1)]) /(13-1)
  # }
  # 
  # #Ranking of treatments for dropout
  # 
  # for(k in 1:NT) {
  #   
  #   order2[k]<- rank(d2[],k)		
  #   
  #   most.effective2[k]<-equals(order2[k],1)
  #   
  #   for(j in 1: NT) {
  #     
  #     effectiveness2[k,j]<- equals(order2[k],j)
  #     
  #     cumeffectiveness2[k,j]<- sum(effectiveness2[k,1:j])}}	
  # 
  # for(k in 1:NT) {
  #   
  #   SUCRA2[k]< - sum(cumeffectiveness2[k,1:( NT-1)]) /(NT-1)}}
    
    }",
   con="mvNMA_2out.txt"
    )


modfile_2out="mvNMA_2out.txt"
