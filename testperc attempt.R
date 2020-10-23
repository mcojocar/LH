###### 1. Load package. #####
require("deSolve")
library(readxl)
library("openxlsx")

###### 2. Setup #####
Mobile <- read_excel("C:/Users/roief/Documents/School/Research/Group Vaccination Allocation Project/Mobile.xlsx")
Ontario6 <- read_excel("C:/Users/roief/Documents/School/Research/Group Vaccination Allocation Project/Ontario6.xlsx")
SVEILR.model<- function(t, x, vparameters){ #function of time(t), x (your compartments), and your parameters
  
  #Each compartment is set as a matrix to account for the multiple age groups in each compartment
  S1    <- x[1] # susceptible 
  V1    <- x[2] # vaccinated
  E1    <- x[3] # exposed. Not contagious. Asymptomatic. 
  I1    <- x[4] # infected. Contagious. Asymptomatic.
  A1    <- x[5] # asymptomatic. Contagious. Asymptomatic
  Y1    <- x[6] # symptomatic. Contagious. Symptomatic
  L1    <- x[7] #isoLated. isolates to prevent all spread
  R1    <- x[8] #recovered
  H1    <- x[9] #hospitalized from population 1
  IR1   <- x[10] #total infected (E+I+A+Y)
  
  S2    <- x[11] # susceptible 
  V2    <- x[12] # vaccinated
  E2    <- x[13] # exposed. Not contagious. Asymptomatic. 
  I2    <- x[14] # infected. Contagious. Asymptomatic.
  A2    <- x[15] # asymptomatic. Contagious. Asymptomatic
  Y2    <- x[16] # symptomatic. Contagious. Symptomatic
  L2    <- x[17] #isoLated. isolates to prevent all spread
  R2    <- x[18] #recovered
  H2    <- x[19] #hospitalized from population 2
  IR2   <- x[20] #total infected (E+I+A+Y)
  
  S3    <- x[21] # susceptible 
  V3    <- x[22] # vaccinated
  E3    <- x[23] # exposed. Not contagious. Asymptomatic. 
  I3    <- x[24] # infected. Contagious. Asymptomatic.
  A3    <- x[25] # asymptomatic. Contagious. Asymptomatic
  Y3    <- x[26] # symptomatic. Contagious. Symptomatic
  L3    <- x[27] #isoLated. isolates to prevent all spread
  R3    <- x[28] #recovered
  H3    <- x[29] #hospitalized from population 3
  IR3   <- x[30] 
  
  IRT   <- x[31]
  
  
  with(as.list(vparameters),{ #how compartments interact with each other
    N1 <- S1+V1+E1+I1+A1+Y1+L1+R1
    N2 <- S2+V2+E2+I2+A2+Y2+L2+R2
    N3 <- S3+V3+E3+I3+A3+Y3+L3+R3
    N_tot <- N1+N2+N3
    
    
    # Group 1 (younger)
    dS1 = (-beta11*S1*(I1+A1+Y1))/N1 - (beta12*S1*(I2+A2+Y2))/N2 - (beta13*S1*(I3+A3+Y3))/N3 - theta1*N_tot
    dV1 = theta1*N_tot
    dE1 = (beta11*S1*(I1+A1+Y1))/N1 + (beta12*S1*(I2+A2+Y2))/N2 + (beta13*S1*(I3+A3+Y3))/N3 - sigma*E1 - omega3*E1
    dI1 = sigma*E1 - psi*I1 - omega3*I1
    dA1 = alpha1*psi*I1 - gamma*A1 - omega3*E1
    dY1 = (1 - alpha1)*psi*I1 - gamma*Y1 - epsilon*Y1 - omega3*Y1
    dL1 = epsilon*Y1 - gamma*L1
    dR1 = gamma*(L1+A1+Y1)
    dH1 = omega1*(L1+Y1)-H1
    #dIR1 = alpha1*E1 - IR1
    dIR1 = alpha1*psi*I1 - IR1 
    
    
    # Group 2 (older)
    dS2 = (-beta22*S2*(I2+A2+Y2))/N2- (beta21*S2*(I1+A1+Y1))/N1 - (beta23*S2*(I3+A3+Y3))/N3 - theta2*N_tot
    dV2 = theta2*N_tot
    dE2 = beta22*S2*(I2+A2+Y2)/N2 + (beta21*S2*(I1+A1+Y1))/N1 + (beta23*S2*(I3+A3+Y3))/N3 - sigma*E2 - omega3*E2
    dI2 = sigma*E2 - psi*I2 - omega3*I2
    dA2 = alpha2*psi*I2 - gamma*A2 - omega3*A2
    dY2 = (1 - alpha2)*psi*I2 - gamma*Y2 - epsilon*Y2 - omega3*Y2
    dL2 = epsilon*Y2 - gamma*L2
    dR2 = gamma*(L2+A2+Y2)
    dH2 = omega2*(L2+Y2)-H2
    #dIR2 = alpha2*E2 - IR2
    dIR2 = alpha2*psi*I2 - IR2
    
    # Group 3 (oldest)
    dS3 = (-beta33*S3*(I3+A3+Y3))/N3 - (beta32*S3*(I1+A1+Y1))/N1 - (beta31*S3*(I2+A2+Y2))/N2 - theta3*N_tot
    dV3 = theta3*N_tot
    dE3 = beta33*S3*(I3+A3+Y3)/N3 + (beta32*S3*(I1+A1+Y1))/N1 + (beta31*S3*(I2+A2+Y2))/N2 - sigma*E3 - omega3*E3
    dI3 = sigma*E3 - psi*I3 - omega3*I3
    dA3 = alpha3*psi*I3 - gamma*A3 - omega3*A3
    dY3 = (1 - alpha3)*psi*I3 - gamma*Y3 - epsilon*Y3 - omega3*Y3
    dL3 = epsilon*Y3 - gamma*L3
    dR3 = gamma*(L3+A3+Y3)
    dH3 = omega2*(L3+Y3)-H3
    #dIR3 = alpha3*E3 - IR3
    dIR3 = alpha3*psi*I3 - IR3
    
    dIRT = dIR1 + dIR2 + dIR3
    
    dx  <- c(dS1,dV1,dE1,dI1,dA1,dY1,dL1,dR1,dH1,dIR1,dS2,dV2,dE2,dI2,dA2,dY2,dL2,dR2,dH2,dIR2,dS3,dV3,dE3,dI3,dA3,dY3,dL3,dR3,dH3,dIR3,dIRT) #dx is the grouping of the individual equations
    list(dx) #link equations together as a list so it runs them all at the same time
  })
}



###### 3. Assign value: time, parameters and initial count of different compartments #####
N_tot <- 14566547  #total starting population
N1 <- 3141693 # Group 3 population - ages 0-19 in Ontario
N2 <- 7977131 # Group 1 population - ages 20-59 in Ontario
N3 <- 3447723 # Group 2 population - ages 60+ in Ontario




I1_0  <- 0
I2_0  <- 1
I3_0  <- 0
S1_0  <- N1-I1_0 #starting conditions is N-the amount of infected people
S2_0  <- N2-I2_0
S3_0  <- N3-I3_0
V1_0  <- 0
V2_0  <- 0
V3_0  <- 0
E1_0  <- 0
E2_0  <- 0
E3_0  <- 0
R1_0  <- 0
R2_0  <- 0
R3_0  <- 0
A1_0  <- 0
A2_0  <- 0
A3_0  <- 0
Y1_0  <- 0
Y2_0  <- 0
Y3_0  <- 0
L1_0  <- 0
L2_0  <- 0
L3_0  <- 0
H1_0  <- 0
H2_0  <- 0
H3_0  <- 0
IR1_0  <- 0
IR2_0  <- 0
IR3_0  <- 0
IRT_0  <- 0

beta <- 0.104 # beta scalar

delayer <- 2


TInter1 <- 30
TInter2 <- TInter1+7
TInter3 <- TInter2+7
TInter4 <- TInter3+7
TInter5 <- TInter4+7
TInter6 <- TInter5+7
TInter7 <- TInter6+7
TInter8 <- TInter7+7
TInter9 <- TInter8+7
TInter10 <- TInter9+7
TInter11 <- TInter10+7
TInter12 <- TInter11+7
TInter13 <- TInter12+7

TInter = c(0,TInter1,TInter2,TInter3,TInter4,TInter5,TInter6,TInter7,TInter8,TInter9,TInter10,TInter11,TInter12,TInter13)


Tmax <- 365*1.5 #integrate for __ days
tau <- 1/4 #size of time step (in days)


#Effective Contact Rates
beta11 <- beta*8.565645854
beta12 <- beta*4.661272358 
beta13 <- beta*0.304014985
beta21 <- beta*2.987996842
beta22 <- beta*11.49063721 
beta23 <- beta*0.522285468
beta31 <- beta*1.248771202
beta32 <- beta*3.686037351
beta33 <- beta*1.982952539

phi <- 1/Tmax # vaccine production rate (per day)
theta1 <- 0 # vaccination allocation to group 1 (WE CONTROL - INDEPENDENT VARIABLE)
theta2 <- 0 # vaccination allocation to group 2 (WE CONTROL - INDEPENDENT VARIABLE)
theta3 <- 0 # vaccination allocation to group 3 (WE CONTROL - INDEPENDENT VARIABLE)
sigma <-  1/2.5 # time from exposure until contagious
psi <- 1/3.5 # time from contagious to symptoms
alpha1 <- 0.5 # proportion of aysmptomatic cases in group 1
alpha2 <- 0.5 # proportion of aysmptomatic cases in group 2
alpha3 <- 0.5 # proportion of aysmptomatic cases in group 3
epsilon <- 0.95 # isolation compliance rate
gamma <- 1/7 # recovery rate once syptoms are shown
omega1 <- 0.03166 # hospitalization rate per day, based on 20.1 percent (340/1686) total 
omega2 <- 0.07295 # hospitalization rate per day, based on 41.1 percent (314/763) total 
omega3 <- 0.0
testperc <- 0.1





times <- seq(TInter1,TInter2 + delayer,by=tau) #function seq returns a sequence
params <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc) #gather parameters and initial conditions
params2 <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc)
inits <- c(S1=S1_0,V1=V1_0,E1=E1_0,I1=I1_0,A1=A1_0,Y1=Y1_0,L1=L1_0,R1=R1_0,H1=H1_0,IR1=IR1_0,S2=S2_0,V2=V2_0,E2=E2_0,I2=I2_0,A2=A2_0,Y2=Y2_0,L2=L2_0,R2=R2_0,H2=H2_0,IR2=IR2_0,S3=S3_0,V3=V3_0,E3=E3_0,I3=I3_0,A3=A3_0,Y3=Y3_0,L3=L3_0,R3=R3_0,H3=H3_0,IR3=IR3_0,IRT=IRT_0)



##### 4. Solve the model #####


out2 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out3 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out4 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out5 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out6 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out7 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out8 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out9 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out10 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out11 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out12 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
out13 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))

times <- seq(0,TInter1 + delayer,by=tau) #function seq returns a sequence
params <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc) #gather parameters and initial conditions
inits <- c(S1=S1_0,V1=V1_0,E1=E1_0,I1=I1_0,A1=A1_0,Y1=Y1_0,L1=L1_0,R1=R1_0,H1=H1_0,IR1=IR1_0,S2=S2_0,V2=V2_0,E2=E2_0,I2=I2_0,A2=A2_0,Y2=Y2_0,L2=L2_0,R2=R2_0,H2=H2_0,IR2=IR2_0,S3=S3_0,V3=V3_0,E3=E3_0,I3=I3_0,A3=A3_0,Y3=Y3_0,L3=L3_0,R3=R3_0,H3=H3_0,IR3=IR3_0,IRT=IRT_0)
out1 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))



#### Use optimal betas from before (q)

beta1 <- 0.04508497
beta2 <- 0.0130156
beta3 <- 0.01137626
beta4 <- 0.0130156
beta5 <- 0.0125625
beta6 <- 0.009456877
beta7 <- 0.006351257
beta8 <- 0.008270636
beta9 <- 0.008270636
beta10 <- 0.01064312
beta11 <- 0.007537498
beta12 <- 0.009456877
beta13 <- 0.004431878


beta <- 1

output <- list(out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13)
betas <- c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10,beta11,beta12,beta13)




times <- seq(0,TInter1 + delayer,by=tau) #function seq returns a sequence
params <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc) #gather parameters and initial conditions
inits <- c(S1=S1_0,V1=V1_0,E1=E1_0,I1=I1_0,A1=A1_0,Y1=Y1_0,L1=L1_0,R1=R1_0,H1=H1_0,IR1=IR1_0,S2=S2_0,V2=V2_0,E2=E2_0,I2=I2_0,A2=A2_0,Y2=Y2_0,L2=L2_0,R2=R2_0,H2=H2_0,IR2=IR2_0,S3=S3_0,V3=V3_0,E3=E3_0,I3=I3_0,A3=A3_0,Y3=Y3_0,L3=L3_0,R3=R3_0,H3=H3_0,IR3=IR3_0,IRT=IRT_0)
out1 <- as.data.frame(lsoda(inits, times, SVEILR.model, params))



beta <- betas[2]
beta11 <- beta*8.565645854
beta12 <- beta*4.661272358 #
beta13 <- beta*0.304014985
beta21 <- beta*2.987996842
beta22 <- beta*11.49063721
beta23 <- beta*0.522285468
beta31 <- beta*1.248771202
beta32 <- beta*3.686037351
beta33 <- beta*1.982952539

times <- seq(TInter[2],TInter[3] + delayer,by=tau) #function seq returns a sequence
params <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc) #gather parameters and initial conditions
inits <- c(S1=S1_0,V1=V1_0,E1=E1_0,I1=I1_0,A1=A1_0,Y1=Y1_0,L1=L1_0,R1=R1_0,H1=H1_0,IR1=IR1_0,S2=S2_0,V2=V2_0,E2=E2_0,I2=I2_0,A2=A2_0,Y2=Y2_0,L2=L2_0,R2=R2_0,H2=H2_0,IR2=IR2_0,S3=S3_0,V3=V3_0,E3=E3_0,I3=I3_0,A3=A3_0,Y3=Y3_0,L3=L3_0,R3=R3_0,H3=H3_0,IR3=IR3_0,IRT=IRT_0)

x <- as.data.frame(lsoda(inits, times, SVEILR.model, params))
output[2] <- x

for (i in 2:3) { 
  beta <- betas[i]
  beta11 <- beta*8.565645854
  beta12 <- beta*4.661272358 #
  beta13 <- beta*0.304014985
  beta21 <- beta*2.987996842
  beta22 <- beta*11.49063721
  beta23 <- beta*0.522285468
  beta31 <- beta*1.248771202
  beta32 <- beta*3.686037351
  beta33 <- beta*1.982952539
  
  times <- seq(TInter[i],TInter[i+1] + delayer,by=tau) #function seq returns a sequence
  params <- c(beta11=beta11,beta12=beta12,beta13=beta13,beta21=beta21,beta22=beta22,beta23=beta23,beta31=beta31,beta32=beta32,beta33=beta33,phi=phi,theta1=theta1,theta2=theta2,theta3=theta3,sigma=sigma,psi=psi,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,epsilon=epsilon,gamma=gamma,omega1=omega1,omega2=omega2,omega3=omega3,testperc=testperc) #gather parameters and initial conditions
  inits <- c(S1=S1_0,V1=V1_0,E1=E1_0,I1=I1_0,A1=A1_0,Y1=Y1_0,L1=L1_0,R1=R1_0,H1=H1_0,IR1=IR1_0,S2=S2_0,V2=V2_0,E2=E2_0,I2=I2_0,A2=A2_0,Y2=Y2_0,L2=L2_0,R2=R2_0,H2=H2_0,IR2=IR2_0,S3=S3_0,V3=V3_0,E3=E3_0,I3=I3_0,A3=A3_0,Y3=Y3_0,L3=L3_0,R3=R3_0,H3=H3_0,IR3=IR3_0,IRT=IRT_0)

  x<- as.data.frame(lsoda(inits, times, SVEILR.model, params))
  output[i] <- x
  }

