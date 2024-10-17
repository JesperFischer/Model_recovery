
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)

logit<-function(delta,spread){
  beta=2 * log( (1 - delta) / delta ) / spread
  return(log(beta))
}

probit<-function(delta,spread){
  beta=2 * sqrt(2) * erfinv(1 - (2 * delta)) / spread
  return(log(beta))
}

gumbel<-function(delta,spread){
  beta=log10( log(delta) / log(1 - delta) ) / spread
  return(log(beta))
}

hyperbolic<-function(delta,spread){
  beta=2 * log(tan((pi/2) * (1 - delta)) / tan((pi /2) * delta)) / (pi * spread)
  return(log(beta))
}

delta=0.025
spread=c(1,200)

sprintf('Prior for the slope of the logistic: N(%.1f,%.1f)',(logit(delta,spread[1])+logit(delta,spread[2]))/2,(logit(delta,spread[1])-logit(delta,spread[2]))/4)
sprintf('Prior for the slope of the cumulative normal: N(%.1f,%.1f)',(probit(delta,spread[1])+probit(delta,spread[2]))/2,(probit(delta,spread[1])-probit(delta,spread[2]))/4)
sprintf('Prior for the slope of the Gumbel: N(%.1f,%.1f)',(gumbel(delta,spread[1])+gumbel(delta,spread[2]))/2,(gumbel(delta,spread[1])-gumbel(delta,spread[2]))/4)
sprintf('Prior for the slope of the hyperbolic sequent: N(%.1f,%.1f)',(hyperbolic(delta,spread[1])+hyperbolic(delta,spread[2]))/2,(hyperbolic(delta,spread[1])-hyperbolic(delta,spread[2]))/4)


