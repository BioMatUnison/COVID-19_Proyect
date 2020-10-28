ModelFunction<-function(Time,State,Pars)
{
  with(as.list(c(State, Pars)),{
  
    if(Time<StartTimePropP)
    {
      valpha_ps<-0
      valpha_pa<-0
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0
    }
    if(Time>=StartTimePropP & Time<TimeGetPropP)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega20<-0
    }
    if(Time>=TimeGetPropP & Time<ReturnStartTime)
    {
      valpha_pls<-0
      valpha_pla<-0
      vomega10<-0
      vomega20<-0  
    }
    if(Time>=ReturnStartTime & Time<TimeGetPropPr)
    {
      vomega10<-0
    }
    if(Time>TimeGetPropPr)
    {
      vomega10<-0
      vomega20<-0
    }
    N  <- S + E + Ia + Is + R + P + Pl
    dS <- -valpha_s*S*Is/N - valpha_a*S*Ia/N - vomega10 * S
    dE <- valpha_s*S*Is/N + valpha_a*S*Ia/N + valpha_ps*P*Is/N + valpha_pa*P*Ia/N + valpha_pls*Pl*Is/N + valpha_pla*Pl*Ia/N- vdelta*E
    dIa<- (1-vtheta)*vdelta*E - veta_a*Ia
    dIs<- vtheta*vdelta*E - vgamma_s*Is
    dH <- vbeta*vgamma_s*Is - vmu*H + vtau*vpsi*Q
    dD <- vpi*vmu*H
    dQ <- (1-vbeta)*vgamma_s*Is - vpsi*Q
    dR <- veta_a*Ia + (1-vpi)*vmu*H + (1-vtau)*vpsi*Q
    dP <- vomega10 * S + - valpha_ps*P*Is/N - valpha_pa*P*Ia/N - vomega20*P
    dPl <- vomega20*P - valpha_pls*Pl*Is/N - valpha_pla*Pl*Ia/N
    dIc <- vtheta*vdelta*E
    dHc <- vbeta*vgamma_s*Is 
    dIt <- vdelta*E
    dQc<-(1-vbeta)*vgamma_s*Is
    list(c(dS,dE,dIa,dIs,dH,dD,dQ,dR,dP,dPl,dIc,dHc,dIt,dQc))
  })
}