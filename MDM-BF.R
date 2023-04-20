## LOAD PACKAGES
library(bnlearn)
library(pcalg)
library(gRain)
library(catnet)
library(fpc)

#--STEP 1-- ESTIMATE NETWORK (TABU METHOD)
NET12 <- tabu(DATA)
plot(NET12)

## FINDING MARKOV EQUIVALENT NETWORKS 
cpdag.P12 <- cpdag(NET12)
plot(cpdag.P12)

# MATRIZES ADJECENTES
AMATP12 <- amat(NET12)
AMcpdag.P12 <- amat(cpdag.P12)

# ARCS (DIRECTED & NOT)
directed.arcs(NET12)
undirected.arcs(NET12)
directed.arcs(cpdag.P12)
undirected.arcs(cpdag.P12)

#--STEP 2-- CALCULATE PARTIAL CONDIT. DAG
# For instance, each subject has 23 optodes, then 46 nodes in total
############### INDICATING ARCS DIRECTION ##############
######################## EQV777 ########################
Kk11= set.arc(cpdag.P12, from = "n26", to = "n24")      # ADD
Kk22= set.arc(Kk11, from = "n24", to = "n25")           # ADD
Kk33= set.arc(Kk22, from = "n26", to = "n25")           # ADD
EQV777= set.arc(Kk33, from = "n24", to = "n21")         # ADD

directed.arcs(EQV777)
undirected.arcs(EQV777)

########################### EQV888 #####################
Ll11= set.arc(cpdag.P12, from = "n26", to = "n24") 
Ll22= set.arc(Ll11, from = "n24", to = "n25")           
Ll33= set.arc(Ll22, from = "n26", to = "n25") 
EQV888= set.arc(Ll33, from = "n21", to = "n24") 

directed.arcs(EQV888)
undirected.arcs(EQV888)

########################### EQV1000 ####################
nn11= set.arc(cpdag.P12, from = "n26", to = "n25") 
nn22= set.arc(nn11, from = "n26", to = "n24")           
nn33= set.arc(nn22, from = "n25", to = "n24") 
EQV1000= set.arc(nn33, from = "n24", to = "n21") 

directed.arcs(EQV1000)
undirected.arcs(EQV1000)

########################### EQV1111 ####################
Oo11= set.arc(cpdag.P12, from = "n26", to = "n25") 
Oo22= set.arc(Oo11, from = "n26", to = "n24")           
Oo33= set.arc(Oo22, from = "n25", to = "n24") 
EQV1111= set.arc(Oo33, from = "n21", to = "n24") 

directed.arcs(EQV1111)
undirected.arcs(EQV1111)

########################### EQV1222 ####################
Pp11= set.arc(cpdag.P12, from = "n26", to = "n25") 
Pp22= set.arc(Pp11, from = "n24", to = "n25")           
Pp33= set.arc(Pp22, from = "n24", to = "n26") 
EQV1222= set.arc(Pp33, from = "n24", to = "n21") 

directed.arcs(EQV1222)
undirected.arcs(EQV1222)

########################### EQV1333 ####################
Qq11= set.arc(cpdag.P12, from = "n26", to = "n25") 
Qq22= set.arc(Qq11, from = "n24", to = "n25")     
Qq33= set.arc(Qq22, from = "n24", to = "n26") 
EQV1333= set.arc(Qq33, from = "n21", to = "n24") 

directed.arcs(EQV1333)
undirected.arcs(EQV1333)

########################### EQV1444 ####################
Rr11= set.arc(cpdag.P12, from = "n24", to = "n25") 
Rr22= set.arc(Rr11, from = "n24", to = "n26")     
Rr33= set.arc(Rr22, from = "n25", to = "n26") 
EQV1444= set.arc(Rr33, from = "n24", to = "n21") 

directed.arcs(EQV1444)
undirected.arcs(EQV1444)

########################### EQV1555 ####################
Ss11= set.arc(cpdag.P12, from = "n24", to = "n25") 
Ss22= set.arc(Ss11, from = "n24", to = "n26")     
Ss33= set.arc(Ss22, from = "n25", to = "n26") 
EQV1555= set.arc(Ss33, from = "n21", to = "n24") 

directed.arcs(EQV1555)
undirected.arcs(EQV1555)

########################### EQV1666 ####################
Tt11= set.arc(cpdag.P12, from = "n25", to = "n24") 
Tt22= set.arc(Tt11, from = "n25", to = "n26")     
Tt33= set.arc(Tt22, from = "n24", to = "n26") 
EQV1666= set.arc(Tt33, from = "n24", to = "n21") 

directed.arcs(EQV1666)
undirected.arcs(EQV1666)

########################### EQV1777 ####################
Uu11= set.arc(cpdag.P12, from = "n25", to = "n24") 
Uu22= set.arc(Uu11, from = "n25", to = "n26")     
Uu33= set.arc(Uu22, from = "n24", to = "n26") 
EQV1777= set.arc(Uu33, from = "n21", to = "n24") 

directed.arcs(EQV1777)
undirected.arcs(EQV1777)

########################### EQV1888 ####################
Vv11= set.arc(cpdag.P12, from = "n25", to = "n24") 
Vv22= set.arc(Vv11, from = "n25", to = "n26")     
Vv33= set.arc(Vv22, from = "n26", to = "n24") 
EQV1888= set.arc(Vv33, from = "n24", to = "n21") 

directed.arcs(EQV1888)
undirected.arcs(EQV1888)

########################### EQV1999 ####################
ww11= set.arc(cpdag.P12, from = "n25", to = "n26") 
ww22= set.arc(ww11, from = "n26", to = "n24")     
ww33= set.arc(ww22, from = "n25", to = "n24") 
EQV1999= set.arc(ww33, from = "n21", to = "n24") 

directed.arcs(EQV1999)
undirected.arcs(EQV1999)

################ BN characteristic ###############
AMP12<- list(amat(EQV777), amat(EQV888), amat(EQV1000), amat(EQV1111),
             amat(EQV1222), amat(EQV1333), amat(EQV1444), amat(EQV1555), 
             amat(EQV1666), amat(EQV1777), amat(EQV1888), amat(EQV1999))

AMP12[[1]]


#--STEP 3-- ESTIMATE MDM-BF
############################################################################################
# Creating function for DLM with Filtering for unknown observational and state variances
############################################################################################

# Input:
# Yt = the vector of observed time series with length T
# Ft = the matrix of covariates with dimension: number of thetas (p) X sample size (T)  
# delta = discount factor | Wt=Ctx(1-delta)/delta
# Gt = the matrix of state equation with dimension: p X p X T. The default is identity matrix block.
# m0 = the vector of prior mean at time t=0 with length p. The default is non-informative prior, with zero mean.
# CS0 = the squared matrix of prior variance - C*0 | C*0Vt = C0, with length p. The default is non-informative prior, with prior variance equal to 3 times the observed variance.
# n0 and d0 = the prior hypermarameters of precision phi ~ G(n0/2; d0/2). The default is non-informative prior, with value of 0.001. n0 has to be higher than 0.

# output: 
# mt = the matrix of posterior mean with dimension p X T
# Ct = the squared matrix of posterior variance with dimension p X p X T 
# Rt = the squared matrix of prior variance with dimension p X p X T
# nt and dt = the vector of prior hypermarameters of precision phi with length T
# ft = the vector of one-step forecast mean with length T
# Qt = the vector of one-step forecast variance with length T
# ets = the vector of standardised residuals with length T
# lpl = Log Predictive Likelihood with length T

dlm_filt <- function(Yt, Ft, delta, Gt = array(diag(nrow(Ft)), dim=c(nrow(Ft),nrow(Ft),length(Yt))), m0 = rep(0,nrow(Ft)), CS0 = 3*diag(nrow(Ft)), n0 = 0.001, d0 = 0.001) {
  
  # defining objects
  p = nrow(Ft) # the number of thetas
  Nt = length(Yt)+1 # the sample size + t=0
  if (n0 == 0){
    n0 = 0.001
    warning("n0 is 0.001")
  }
  Y = rep(0, Nt)
  Y[2:Nt] = Yt
  F = array(0, dim=c(p,Nt))
  F[,2:Nt] = Ft
  G = array(0, dim=c(p,p,Nt))
  G[,,2:Nt] = Gt
  mt = array(0, dim=c(p,Nt))
  mt[,1] = m0
  Ct = array(0, dim=c(p,p,Nt)) 
  Ct[,,1] = CS0*d0/n0
  Rt = array(0, dim=c(p,p,Nt))
  nt = rep(0, Nt)
  nt[1] = n0
  dt = rep(0, Nt)
  dt[1] = d0
  ft = rep(0, Nt)
  Qt = rep(0, Nt)
  ets = rep(0, Nt)
  lpl = rep(0, Nt)
  
  for (i in 2:Nt){
    
    # Posterior at {t-1}: (theta_{t-1}|y_{t-1}) ~ t_{n_{t-1}}[m_{t-1}, C_{t-1} = C*_{t-1}xd_{t-1}/n_{t-1}]
    # Prior at {t}: (theta_{t}|y_{t-1}) ~ t_{n_{t-1}}[a_{t}, R_{t}]
    at = G[,,i] %*% mt[,(i-1)]
    RSt = G[,,i] %*% (Ct[,,(i-1)]*nt[(i-1)]/dt[(i-1)]) %*% t(G[,,i]) / delta
    Rt[,,i] = RSt * dt[(i-1)] / nt[(i-1)]
    
    # One-step forecast: (Y_{t}|y_{t-1}) ~ t_{n_{t-1}}[f_{t}, Q_{t}]
    ft[i] = t(F[,i]) %*% at
    QSt = t(F[,i]) %*% RSt %*% F[,i] + 1
    Qt[i] = QSt * dt[(i-1)] / nt[(i-1)]
    et = Y[i] - ft[i]
    ets[i] = et / sqrt(Qt[i])
    
    # Posterior at t: (theta_{t}|y_{t}) ~ t_{n_{t}}[m_{t}, C_{t}]
    At = Rt[,,i] %*% F[,i] / Qt[i]
    mt[,i] = at + At * et
    nt[i] = nt[(i-1)] + 1
    dt[i] = dt[(i-1)] + (et^2) / QSt
    CSt = RSt - (At %*% t(At)) * QSt[1]
    Ct[,,i] = CSt * dt[i] / nt[i]
    
    # Log Predictive Likelihood
    lpl[i] <- lgamma((nt[i]+1)/2)-lgamma(nt[i]/2)-0.5*log(pi*nt[i]*Qt[i])-((nt[i]+1)/2)*log(1+(1/nt[i])*et^2/Qt[i])         
    
  }
  
  result <- list(mt=mt[,2:Nt], Ct=Ct[,,2:Nt], Rt=Rt[,,2:Nt], nt=nt[2:Nt], dt=dt[2:Nt], ft=ft[2:Nt], Qt=Qt[2:Nt], ets=ets[2:Nt], lpl=lpl[2:Nt])
  
  return(result)
}


############################################################################################
# Creating function for DLM with Smoothing for unknown observational and state variances
############################################################################################

# Input: all objects are resulted from "dlm_filt", except Gt
# mt = the matrix of posterior mean with dimension p X T
# Ct = the squared matrix of posterior variance with dimension p X p X T 
# Rt = the squared matrix of prior variance with dimension p X p X T
# nt and dt = the vector of prior hypermarameters of precision phi with length T
# Gt = the matrix of state equation with dimension: p X p X T. The default is identity matrix block.

# output: 
# smt = the matrix of smoothing posterior mean with dimension p X T
# sCt = the squared matrix of smoothing posterior variance with dimension p X p X T 

dlm_smoo <- function(mt, Ct, Rt, nt, dt, Gt = 0) {
  
  # defining objects
  if (is.vector(mt)){
    mt = array(mt, dim=c(1,length(mt)))
    Ct = array(Ct, dim=c(1,1,length(mt)))
    Rt = array(Rt, dim=c(1,1,length(Rt)))     
  }
  if (Gt == 0){Gt = array(diag(nrow(mt)), dim=c(nrow(mt),nrow(mt),ncol(mt)))}
  p = nrow(mt) # the number of thetas
  Nt = ncol(mt) # the sample size
  smt = array(0, dim=c(p,Nt))
  sCt = array(0, dim=c(p,p,Nt)) 
  
  # in the last time point
  smt[,Nt] = mt[,Nt]
  sCt[,,Nt] = Ct[,,Nt] 
  
  # for other time points
  for (i in (Nt-1):1){
    RSt = Rt[,,(i+1)]*nt[i]/dt[i]
    CSt = Ct[,,i]*nt[i]/dt[i]
    inv.sR = solvecov(RSt, cmax = 1e+10)$inv
    B = CSt %*% t(Gt[,,(i+1)]) %*% inv.sR
    smt[,i] = mt[, i] + B %*% (smt[,(i+1)] - Gt[,,(i+1)] %*% mt[,i])
    sCS = CSt + B %*% (sCt[,,(i+1)]*nt[Nt]/dt[Nt] - RSt) %*% t(B)     
    sCt[,,i] = sCS * dt[Nt] / nt[Nt]
  }
  
  result <- list(smt=smt, sCt=sCt)
  return(result)
}

###################################################################
### Choosing the delta
###################################################################

# Input:
#  dts = the matrix with dataset; Number of timepoints X Number of nodes
#  m_ad = # Square Matrix Adjacent with dimension = Number of nodes # 1 if edge exists; 0 otherwise
#  nbf => the Log Predictive Likelihood will be calculate from this time point. It has to be a positive integer number. The default is 15.
#  delta = the vector with the sequence of all discount factors. The default is seq(from=0.5, to=1.0, by=0.01).

# Output: 
# lpldet = LPL for each value of delta and for each node; length(delta) X Number of nodes;
# DF_hat = vector with delta that maximizes the LPL for each node with dimension = Number of nodes.

CDELT <- function(dts,m_ad,nbf=15,delta=seq(from=0.5, to=1.0, by=0.01)) {
  nd = length(delta)
  Nn = ncol(dts)
  Nt = nrow(dts)
  lpldet = array(0,dim=c(nd,Nn))
  for (k in 1:nd){
    for (i in 1:Nn){
      # Initials:
      p = sum(m_ad[,i])
      if (m_ad[i,i] == 0) {p = p + 1}  	
      Ft = array(1, dim=c(Nt,p))
      aux = c(1:Nn)[m_ad[,i]>0]
      aux2 = aux[aux!=i] 
      if (length(aux2)>0){ Ft[,2:(length(aux2)+1)] = dts[,aux2]}
      Yt = dts[,i]
      # DLM
      a=dlm_filt(Yt, t(Ft), delta=delta[k])
      lpldet[k,i]=sum(a$lpl[nbf:Nt])
    }
  }
  #DF_hat=delta[max.col(t(lpldet))] # with some deltas provide NA in lpl, it means that this delta is not good for this particular dataset so we have to use:
  DF_hat = rep(0,Nn)
  for (i in 1:Nn){
    DF_hat[i] = na.omit(delta[lpldet[,i]==max(lpldet[,i],na.rm=TRUE)])[1]
  }
  result <- list(lpldet = lpldet, DF_hat = DF_hat)
  return(result)
}

############################################################################################
# Creating function for MDM with Filtering for unknown observational and state variances
############################################################################################

# Input:
# dts = the matrix with dataset; Number of timepoints X Number of nodes
# m_ad = # Square Matrix Adjacent with dimension = Number of nodes # 1 if edge exists; 0 otherwise
# DF_hat = vector with delta that maximizes the LPL for each node with dimension = Number of nodes.

# output: 
# mt = list with the matrix of posterior mean with dimension p X T
# Ct = list with the squared matrix of posterior variance with dimension p X p X T 
# Rt = list with the squared matrix of prior variance with dimension p X p X T
# nt and dt = list with the vector of prior hypermarameters of precision phi with length T
# ft = list with the vector of one-step forecast mean with length T
# Qt = list with the vector of one-step forecast variance with length T
# ets = list with the vector of standardised residuals with length T
# lpl = list with Log Predictive Likelihood with length T

mdm_filt <- function(dts, m_ad, DF_hat) {
  Nn = ncol(dts)
  Nt = nrow(dts)
  mt = vector(Nn, mode = "list")
  Ct = vector(Nn, mode = "list")
  Rt = vector(Nn, mode = "list")
  nt = vector(Nn, mode = "list")
  dt = vector(Nn, mode = "list")
  ft = vector(Nn, mode = "list")
  Qt = vector(Nn, mode = "list")
  ets = vector(Nn, mode = "list")
  lpl = vector(Nn, mode = "list")
  for (i in 1:Nn){
    # Initials:
    p = sum(m_ad[,i])
    if (m_ad[i,i] == 0) {p = p + 1}   	         
    Ft = array(1, dim=c(Nt,p))
    aux = c(1:Nn)[m_ad[,i]>0]
    aux2 = aux[aux!=i] 
    if (length(aux2)>0){ Ft[,2:(length(aux2)+1)] = dts[,aux2]}
    Yt = dts[,i]
    # DLM
    a=dlm_filt(Yt, t(Ft), delta=DF_hat[i])
    mt[[i]] = a$mt
    Ct[[i]] = a$Ct
    Rt[[i]] = a$Rt
    nt[[i]] = a$nt
    dt[[i]] = a$dt
    ft[[i]] = a$ft
    Qt[[i]] = a$Qt
    ets[[i]] = a$ets
    lpl[[i]] = a$lpl         
  }
  result <- list(mt=mt, Ct=Ct, Rt=Rt, nt=nt, dt=dt, ft=ft, Qt=Qt, ets=ets, lpl=lpl)
  return(result)
}

############################################################################################
# Creating function for MDM with Smoothing for unknown observational and state variances
############################################################################################

# Input: all objects are resulted from "mdm_filt"
# mt = list with the matrix of posterior mean with dimension p X T
# Ct = list with the squared matrix of posterior variance with dimension p X p X T 
# Rt = list with the squared matrix of prior variance with dimension p X p X T
# nt and dt = list with the vector of prior hypermarameters of precision phi with length T

# output: 
# smt = list with the matrix of smoothing posterior mean with dimension p X T
# sCt = list with the squared matrix of smoothing posterior variance with dimension p X p X T 

mdm_smoo <- function(mt, Ct, Rt, nt, dt) {
  Nn = length(mt) # the number of nodes
  smt = vector(Nn, mode = "list")
  sCt = vector(Nn, mode = "list")
  for (i in 1:Nn){
    a=dlm_smoo(mt=mt[[i]], Ct=Ct[[i]], Rt=Rt[[i]], nt=nt[[i]], dt=dt[[i]])
    smt[[i]] = a$smt
    sCt[[i]] = a$sCt
  }
  result <- list(smt=smt, sCt=sCt)
  return(result)
}


###########################################
## LPL  
### Choosing DF for node and model
###########################################

P12=as.matrix(P12m)

DF0 = CDELT(P12,AMATP12) #DBN
amat0lpl<-sum(apply(DF0$lpldet,2,max,na.rm=TRUE)) #LPL

DF1 = CDELT(P12,amat(EQV777)) 
amat1lpl<-sum(apply(DF1$lpldet,2,max,na.rm=TRUE))

DF2 = CDELT(P12,amat(EQV888)) 
amat2lpl<-sum(apply(DF2$lpldet,2,max,na.rm=TRUE))

DF3 = CDELT(P12,amat(EQV1000))
amat3lpl<-sum(apply(DF3$lpldet,2,max,na.rm=TRUE))

DF4 = CDELT(P12,amat(EQV1111))
amat4lpl<-sum(apply(DF4$lpldet,2,max,na.rm=TRUE))

DF5 = CDELT(P12,amat(EQV1222))
amat5lpl<-sum(apply(DF5$lpldet,2,max,na.rm=TRUE))

DF6 = CDELT(P12,amat(EQV1333))
amat6lpl<-sum(apply(DF6$lpldet,2,max,na.rm=TRUE))

DF7 = CDELT(P12,amat(EQV1444))
amat7lpl<-sum(apply(DF7$lpldet,2,max,na.rm=TRUE))

DF8 = CDELT(P12,amat(EQV1555))
amat8lpl<-sum(apply(DF8$lpldet,2,max,na.rm=TRUE))

DF9 = CDELT(P12,amat(EQV1666))
amat9lpl<-sum(apply(DF9$lpldet,2,max,na.rm=TRUE))

DF10 = CDELT(P12,amat(EQV1777))
amat10lpl<-sum(apply(DF10$lpldet,2,max,na.rm=TRUE))

DF11 = CDELT(P12,amat(EQV1888))
amat11lpl<-sum(apply(DF11$lpldet,2,max,na.rm=TRUE))

DF12 = CDELT(P12,amat(EQV1999))
amat12lpl<-sum(apply(DF12$lpldet,2,max,na.rm=TRUE))

## LPL values
LPL_P12<-c(amat0lpl,amat1lpl,amat2lpl,amat3lpl,amat4lpl,amat5lpl,
           amat6lpl,amat7lpl,amat8lpl,amat9lpl,amat10lpl,amat11lpl,amat12lpl)

## Analyzing the models' LPL
max(LPL_P12)
min(LPL_P12)
plot(LPL_P12,type="b")

LPLMAX_P12<-amat1lpl  #MAX LPL

## APPLYING FILTER & SMOOTHING BASED ON THE BEST MODEL 
mdmfP12 = mdm_filt(P12, amat(EQV777), DF1$DF_hat)          ##Filtered
mdmsP12 = mdm_smoo(mdmfP12$mt, mdmfP12$Ct, mdmfP12$Rt, mdmfP12$nt, mdmfP12$dt) ##Smoothed

MTFP12<-mdmfP12$mt             #Filtered evolution parameter    
MTSP12<-mdmsP12$smt            #Smoothed evolution parameter


####################### THETAs MEAN #################
# FILTER
FTHETA_P12<-lapply(MTFP12,rowMeans)

MATP12=amat(EQV777)
## CREATING THE ADJACENTE MATRIZ BASED ON THETA MEAN
for(i in 1:length(FTHETA_P12)){
  MATP12[MATP12[,i]==1,i]<-FTHETA_P12[[i]][-1]
}

colnames(MATP12)<-c("P1_1",	"P1_2",	"P1_3",	"P1_4",	"P1_5",	"P1_6",	"P1_7",	"P1_8",	"P1_9",	"P1_10",	"P1_11",	"P1_12",	"P1_13",	"P1_14",	"P1_15",	"P1_16",	"P1_17",	"P1_18",	"P1_19",	"P1_20",	"P1_21",	"P1_22",	"P1_23","P2_1",	"P2_2",	"P2_3",	"P2_4",	"P2_5",	"P2_6",	"P2_7",	"P2_8",	"P2_9",	"P2_10",	"P2_11",	"P2_12",	"P2_13",	"P2_14",	"P2_15",	"P2_16",	"P2_17",	"P2_18",	"P2_19",	"P2_20",	"P2_21",	"P2_22",	"P2_23")
rownames(MATP12)<-c("P1_1",	"P1_2",	"P1_3",	"P1_4",	"P1_5",	"P1_6",	"P1_7",	"P1_8",	"P1_9",	"P1_10",	"P1_11",	"P1_12",	"P1_13",	"P1_14",	"P1_15",	"P1_16",	"P1_17",	"P1_18",	"P1_19",	"P1_20",	"P1_21",	"P1_22",	"P1_23","P2_1",	"P2_2",	"P2_3",	"P2_4",	"P2_5",	"P2_6",	"P2_7",	"P2_8",	"P2_9",	"P2_10",	"P2_11",	"P2_12",	"P2_13",	"P2_14",	"P2_15",	"P2_16",	"P2_17",	"P2_18",	"P2_19",	"P2_20",	"P2_21",	"P2_22",	"P2_23")

# SMOOTH
STHETA_P12<-lapply(MTSP12,rowMeans)

SMATP12=amat(EQV777)
## CREATING THE ADJACENTE MATRIZ BASED ON THETA MEAN
for(i in 1:length(STHETA_P12)){
  SMATP12[SMATP12[,i]==1,i]<-STHETA_P12[[i]][-1]
}

colnames(SMATP12)<-c("P1_1",	"P1_2",	"P1_3",	"P1_4",	"P1_5",	"P1_6",	"P1_7",	"P1_8",	"P1_9",	"P1_10",	"P1_11",	"P1_12",	"P1_13",	"P1_14",	"P1_15",	"P1_16",	"P1_17",	"P1_18",	"P1_19",	"P1_20",	"P1_21",	"P1_22",	"P1_23","P2_1",	"P2_2",	"P2_3",	"P2_4",	"P2_5",	"P2_6",	"P2_7",	"P2_8",	"P2_9",	"P2_10",	"P2_11",	"P2_12",	"P2_13",	"P2_14",	"P2_15",	"P2_16",	"P2_17",	"P2_18",	"P2_19",	"P2_20",	"P2_21",	"P2_22",	"P2_23")
rownames(SMATP12)<-c("P1_1",	"P1_2",	"P1_3",	"P1_4",	"P1_5",	"P1_6",	"P1_7",	"P1_8",	"P1_9",	"P1_10",	"P1_11",	"P1_12",	"P1_13",	"P1_14",	"P1_15",	"P1_16",	"P1_17",	"P1_18",	"P1_19",	"P1_20",	"P1_21",	"P1_22",	"P1_23","P2_1",	"P2_2",	"P2_3",	"P2_4",	"P2_5",	"P2_6",	"P2_7",	"P2_8",	"P2_9",	"P2_10",	"P2_11",	"P2_12",	"P2_13",	"P2_14",	"P2_15",	"P2_16",	"P2_17",	"P2_18",	"P2_19",	"P2_20",	"P2_21",	"P2_22",	"P2_23")

#--STEP 4-- DYNAMIC ANALYSIS
# Time-varying Parameter (Filt$mt[[node]][param,])
# For instance, node=1 param=4(beta)
tt = qt(0.05, df=mdmfP12$nt[[1]])
plot(mdmfP12$mt[[1]][4,], type="l", xlab="Time", col="blue",ylim=c(-0.3,0.3), ylab="Posterior")
     lines((mdmfP12$mt[[1]][4,]-tt*sqrt(mdmfP12$Ct[[1]][4,4,])), lty = "dashed", col="blue")
     lines((mdmfP12$mt[[1]][4,]+tt*sqrt(mdmfP12$Ct[[1]][4,4,])), lty = "dashed", col="blue")
     lines(mdmsP12$smt[[1]][4,], lty = "solid", col="green")
     lines((mdmsP12$smt[[1]][4,]-tt[length(tt)]*sqrt(mdmsP12$sCt[[1]][4,4,])), lty = "dashed", col="green")
     lines((mdmsP12$smt[[1]][4,]+tt[length(tt)]*sqrt(mdmsP12$sCt[[1]][4,4,])), lty = "dashed", col="green")
     legend("top", legend = c("Smooth Post. Mean","Filter Post. Mean","Credible Interval"), lty = c("solid", "solid","dashed"), col = c("green","blue", "blue"), bty = "n")

# ----- Define a function for CAUSAL plotting a matrix ----- #
myImagePlot_b <- function(x, cl=1, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  if (cl==1){
    # from black to white
    ColorRamp <- gray(seq(0,1,length=256))} 
  else {
    ColorRamp <- gray(seq(1,0,length=256))
  }
  
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #


#################################################################################################################
############## FILTER VISUALIZATION #####################
myImagePlot_b(MATP12,cl=0,zlim=c(0,0.088),title="Indivíduos 1 e 2 conjuntamente")  

############## SMOOTH VISUALIZATION #####################
myImagePlot_b(SMATP12,cl=0,zlim=c(0,max(SMATP12)),title="Indivíduos 1 e 2 conjuntamente")  

############## CORRELATION #####################
MCORP12<-cor(P12m)
myImagePlot_b(MCORP12,cl=0,zlim=c(0,1.0),title="Indivíduos 1 e 2 conjuntamente")    

