#### caluclo de Rogan gladen con intervalo de agresti-coull ajustado por Se y Sp conocidas

rg.agresti <- function(x, n, se, sp, conf = 0.95){
  z <- qnorm(1-(1-conf)/2)
  np = n+z^2
  ap = x/n
  AP = (n*ap+(z^2)/2)/(np)
  P_rg <- (AP+sp-1)/(se+sp-1)
  if(P_rg<0){
    P_rg <- 0
  }
  else(P_rg = P_rg)
  ee <- z*1/(se+sp-1)*sqrt(AP*(1-AP)/np)
  if(P_rg == 0){
    lower <- 0
    upper <- P_rg+ee
  }
  else(lower <-  P_rg-ee)
  upper <-  P_rg+ee
  r <- tibble(method = "Lang",x = x, n = n, mean = P_rg, lower = lower, upper = upper)
  return(r)
}


