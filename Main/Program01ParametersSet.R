rm(list=ls())

#Calling some libraries:
library(deSolve)
library(dplyr)
library(extraDistr) 


#External function used in this program
source("ModelFunction.R",local=FALSE)

#######################################
#                Data                 #
#######################################
DataBase<-read.csv("DatosHillo200719.csv",header=TRUE)[,-1] #Reading Hermosillo database
n<-82
DatIN<-DataBase[1:n,2]
DatHN<-DataBase[1:n,3]
DatQN<-DataBase[1:n,4]
DatDN<-DataBase[1:n,5]



####################################################################
####################################################################
####################################################################
###
###           Empirical restriction on epidemic curves
###
####################################################################
####################################################################
####################################################################

# Initial conditions for the model
ValS<-930668      
ValE<-0          
ValIa<-0
ValIs<-1
ValH<-0
ValD<-0
ValQ<-0
ValR<-0
ValP<-0
ValPr<-0

Valt0<-0
ValT0<-200

valuem<-1000 #m = 1000 parameter sets for each scenario; used later to calculate daily incidence of: symptomatic
             #infections, hospitalized cases, ambulatory cases, and deaths.
Matrixparam<-matrix(NA,valuem,21)

colnames(Matrixparam)<-c("valpha_a","valpha_s","valpha_pa","valpha_ps","valpha_pla","valpha_pls",
                         "vbeta","vdelta","vetaA","vgamma_s","vmu","vpi","vpsi","vtau",
                         "vtheta","valomega10","valomega20","sumssquarederrorsI","sumssquarederrorsH",
                         "sumssquarederrorsD","sumssquarederrorsQ")
while(valuem>0){
  
  #Model parameter distributions
  valpha_a<-rtnorm(n = 1, mean = 1.245815,sd=0.0115464, 1.22, 1.27)
  valpha_s<-rtnorm(n = 1, mean = 1.207603,sd=0.04644592, 1.10, 1.30)
  valpha_pa<-rtnorm(n = 1, mean = 0.007402952 ,sd=0.01143118, 0,.04)
  valpha_ps<-rtnorm(n = 1, mean = 1.201071,sd=0.00428834, 1.19, 1.21)
  valpha_pla<-runif(1,0.05,.25)
  valpha_pls<-rtnorm(n = 1, mean =0.5839841,sd=0.05096, .46, .66)
  
  vbeta<-runif(1,0.14,.25)
  
  vdelta<-rtnorm(n = 1, mean = 0.2923478,sd=0.02117284, .25, .35)
  
  veta_a<-rtnorm(n = 1, mean = 0.1503479,sd=0.02530067, .08, .2)
  
  vgamma_s<-runif(1,.5,2)
  
  vmu<-runif(1,0.1,.8)
  
  vpi<-runif(1,0.14,.3)
  
  vpsi<-runif(1,.01,1)
  
  vtau<-runif(1,0,.2)
  
  vtheta<-runif(1,.01,.2)
  
  
  #On-and-off periods of social distancing
  
  ValStartTimePropP<-5         
  ValPeriodGetPropP<-30 #T_U1 - T_L1 
  ValReturnStartTime<-50    
  ValPeriodGetPropPr<-15 #T_U2 - T_L2
  
  
  valomega10<- runif(1,.04,.08)
  valomega20<- runif(1,0.003,0.015)
  ValTimeGetPropP<-ValStartTimePropP+ValPeriodGetPropP
  ValTimeGetPropPr<-ValReturnStartTime+ValPeriodGetPropPr
  
  #Solving systems of differential equations
  yini<-c(S=ValS,E=ValE,Ia=ValIa,Is=ValIs,H=ValH,D=ValD,Q=ValQ,R=ValR,P=ValP,
          Pl=ValPr,Ic=ValIs,Hc=ValH,It=ValIs+ValIa,Qc=ValQ)
  TimeGrid<-seq(from=Valt0,to=ValT0,by=0.1) 
  vecpar<-c(valpha_s=valpha_s,valpha_a=valpha_a,valpha_ps=valpha_ps,valpha_pa=valpha_pa,
            valpha_pls=valpha_pls,valpha_pla=valpha_pla,vdelta=vdelta,vtheta=vtheta,
            veta_a=veta_a,vgamma_s=vgamma_s,vbeta=vbeta,vmu=vmu,vpi=vpi,
            vpsi=vpsi,vtau=vtau,vomega10=valomega10,vomega20=valomega20,
            StartTimePropP=ValStartTimePropP,TimeGetPropP= ValTimeGetPropP,
            ReturnStartTime=ValReturnStartTime,TimeGetPropPr= ValTimeGetPropPr)
  sol<-ode(yini,TimeGrid,ModelFunction,vecpar)
  vecobsFinal<-Valt0:ValT0
  Accum<-sol[,c(7,12,13,15)][sol[,1] %in% vecobsFinal,]
  NewCases_temp<-Accum[2:(length(vecobsFinal)),]-Accum[1:(length(vecobsFinal)-1),] 
  NewCases<-rbind(Accum[1,],NewCases_temp)
  NewCasesDc<-NewCases[,1]
  NewCasesIc<-NewCases[,2]
  NewCasesHc<-NewCases[,3]
  NewCasesQc<-NewCases[,4]
  
  sumssquarederrorsI<-sum((DatIN-NewCasesIc[1:n])^2)
  sumssquarederrorsH<-sum((DatHN-NewCasesHc[1:n])^2)
  sumssquarederrorsD<-sum((DatDN-NewCasesDc[1:n])^2)
  sumssquarederrorsQ<-sum((DatQN-NewCasesQc[1:n])^2)
  
  Infaccumulated<-sol[ValT0*10+1,14]
  InfSaccumulated<-sum(NewCasesIc)
  Haccumulated<-sum(NewCasesHc)
  Daccumulated<-sum(NewCasesDc)
  Qaccumulated<-sum(NewCasesQc)
  
  #Considering an empirical constraint for prevalence
  if(Infaccumulated<=.216*ValS){
    Matrixparam[valuem,]<-c(valpha_a,valpha_s,valpha_pa,valpha_ps,valpha_pla,valpha_pls,vbeta,
                            vdelta,veta_a,vgamma_s,vmu,vpi,vpsi,vtau,vtheta,valomega10,valomega20,
                            sumssquarederrorsI,sumssquarederrorsH,sumssquarederrorsD,sumssquarederrorsQ)
    valuem<-valuem-1
  }
}


# Calculation of some specific upper bounds
UpperBoundI<-quantile(Matrixparam[,18],c(.25))
UpperBoundQ<-quantile(Matrixparam[,21],c(.25))
UpperBoundH<-quantile(Matrixparam[,19],c(.25))
UpperBoundD<-quantile(Matrixparam[,20],c(.25))

####################################################################
####################################################################
####################################################################
###
###          Selection of 5000 solutions
###
####################################################################
####################################################################
####################################################################

# Initial conditions for the model
ValS<-930668      
ValE<-0         
ValIa<-0
ValIs<-1
ValH<-0
ValD<-0
ValQ<-0
ValR<-0
ValP<-0
ValPr<-0


Valt0<-0
ValT0<-200

valuej<-5000 # To select 5000 solutions from the system of differential equations
Matrixparam<-matrix(NA,valuej,21)

colnames(Matrixparam)<-c("valpha_a","valpha_s","valpha_pa","valpha_ps","valpha_pla","valpha_pls",
                         "vbeta","vdelta","vetaA","vgamma_s","vmu","vpi","vpsi","vtau",
                         "vtheta","valomega10","valomega20","sumssquarederrorsI","sumssquarederrorsH",
                         "sumssquarederrorsD","sumssquarederrorsQ")
while(valuej>0){
  
  # Model parameter distributions
  valpha_a<-rtnorm(n = 1, mean = 1.245815,sd=0.0115464, 1.22, 1.27)
  valpha_s<-rtnorm(n = 1, mean = 1.207603,sd=0.04644592, 1.10, 1.30)
  valpha_pa<-rtnorm(n = 1, mean = 0.007402952 ,sd=0.01143118, 0,.04)
  valpha_ps<-rtnorm(n = 1, mean = 1.201071,sd=0.00428834, 1.19, 1.21)
  valpha_pla<-runif(1,0.05,.25)
  valpha_pls<-rtnorm(n = 1, mean =0.5839841,sd=0.05096, .46, .66)
  
  vbeta<-runif(1,0.14,.25)
  
  vdelta<-rtnorm(n = 1, mean = 0.2923478,sd=0.02117284, .25, .35)
  
  veta_a<-rtnorm(n = 1, mean = 0.1503479,sd=0.02530067, .08, .2)
  
  vgamma_s<-runif(1,.5,2)
  
  vmu<-runif(1,0.1,.8)
  
  vpi<-runif(1,0.14,.3)
  
  vpsi<-runif(1,.01,1)
  
  vtau<-runif(1,0,.2)
  
  vtheta<-runif(1,.01,.2)
  
  
  #On-and-off periods of social distancing
  ValStartTimePropP<-5         
  ValPeriodGetPropP<-30 #T_U1 - T_L1
  ValReturnStartTime<-50    
  ValPeriodGetPropPr<-15 #T_U2 - T_L2
  
  
  valomega10<- runif(1,.04,.08)
  valomega20<- runif(1,0.003,0.015)
  ValTimeGetPropP<-ValStartTimePropP+ValPeriodGetPropP
  ValTimeGetPropPr<-ValReturnStartTime+ValPeriodGetPropPr
  
  #Solving systems of differential equations
  yini<-c(S=ValS,E=ValE,Ia=ValIa,Is=ValIs,H=ValH,D=ValD,Q=ValQ,R=ValR,P=ValP,
          Pl=ValPr,Ic=ValIs,Hc=ValH,It=ValIs+ValIa,Qc=ValQ)
  TimeGrid<-seq(from=Valt0,to=ValT0,by=0.1) 
  vecpar<-c(valpha_s=valpha_s,valpha_a=valpha_a,valpha_ps=valpha_ps,valpha_pa=valpha_pa,
            valpha_pls=valpha_pls,valpha_pla=valpha_pla,vdelta=vdelta,vtheta=vtheta,
            veta_a=veta_a,vgamma_s=vgamma_s,vbeta=vbeta,vmu=vmu,vpi=vpi,
            vpsi=vpsi,vtau=vtau,vomega10=valomega10,vomega20=valomega20,
            StartTimePropP=ValStartTimePropP,TimeGetPropP= ValTimeGetPropP,
            ReturnStartTime=ValReturnStartTime,TimeGetPropPr= ValTimeGetPropPr)
  sol<-ode(yini,TimeGrid,ModelFunction,vecpar)
  vecobsFinal<-Valt0:ValT0
  Accum<-sol[,c(7,12,13,15)][sol[,1] %in% vecobsFinal,]
  NewCases_temp<-Accum[2:(length(vecobsFinal)),]-Accum[1:(length(vecobsFinal)-1),] 
  NewCases<-rbind(Accum[1,],NewCases_temp)
  NewCasesDc<-NewCases[,1]
  NewCasesIc<-NewCases[,2]
  NewCasesHc<-NewCases[,3]
  NewCasesQc<-NewCases[,4]
  
  sumssquarederrorsI<-sum((DatIN-NewCasesIc[1:n])^2)
  sumssquarederrorsH<-sum((DatHN-NewCasesHc[1:n])^2)
  sumssquarederrorsD<-sum((DatDN-NewCasesDc[1:n])^2)
  sumssquarederrorsQ<-sum((DatQN-NewCasesQc[1:n])^2)
  if(sumssquarederrorsI<UpperBoundI && sumssquarederrorsQ<UpperBoundQ && sumssquarederrorsH< UpperBoundH  && sumssquarederrorsD< UpperBoundD){
    Infaccumulated<-sol[ValT0*10+1,14]
    InfSaccumulated<-sum(NewCasesIc)
    Haccumulated<-sum(NewCasesHc)
    Daccumulated<-sum(NewCasesDc)
    Qaccumulated<-sum(NewCasesQc)
    
    #Considering an empirical constraint for prevalence
    if(Infaccumulated<=.216*ValS){
      Matrixparam[valuej,]<-c(valpha_a,valpha_s,valpha_pa,valpha_ps,valpha_pla,valpha_pls,vbeta,
                                vdelta,veta_a,vgamma_s,vmu,vpi,vpsi,vtau,vtheta,valomega10,valomega20,
                                sumssquarederrorsI,sumssquarederrorsH,sumssquarederrorsD,sumssquarederrorsQ)
      valuej<-valuej-1
    }
  }
}

write.csv(Matrixparam,file="Matrixparam.csv")
