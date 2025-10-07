library(deSolve)


alpha_0 = 0.76
lambda = alpha_0
sigma = 0.086
theta = 0.3
eta = 0.8
u_0 = eta #initial point
t0 = 0.0
T = 1.0
Time = c(0.0,T)
dt = 0.001 ###partition


step = 0.001
m1=1000
u0 = numeric(m1+20)
u0[1]=u_0

epsi=sqrt(.Machine$double.eps) ###### h for finite difference

M=100  


ode_func <- function(t, state, parms){
  with(as.list(c(state, parms)),{
    du <- c()
    #rate of change
    du[1] <-  (p1 + p2 - 1) * u[1] - (p1 + 3 * p2) * (u[1] ^ 2) + 2 * p2 *( u[1] ^ 3)
    #Orden 1
    for(j in 2:m1+1){
      du[j] <-  (p1 + p2 - 1) * u[j] - (p1 + 3 * p2) * (u[j] ^ 2) + 2 * p2 * (u[j] ^ 3) + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * p2 * u[1] * (1- u[1])
          
      }
    ##Orden 2 mixed
    du[m1+2] <- (p1 + p2 - 1) * u[m1+2] - (p1 + 3 * p2) * (u[m1+2] ^ 2) + 2 * p2 *( u[m1+2] ^ 3) + sqrt(2/T) * p2 * (sin(((2-1) * pi * t)/ T) * u[3] * (1- u[3]) + sin(((3-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[m1+3] <- (p1 + p2 - 1) * u[m1+3] - (p1 + 3 * p2) * (u[m1+3] ^ 2 )+ 2 * p2 *( u[m1+3] ^ 3) + p2 * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[2] * (1- u[2]))
    du[m1+4] <- (p1 + p2 - 1) * u[m1+4] - (p1 + 3 * p2) * (u[m1+4] ^ 2 )+ 2 * p2 *( u[m1+4] ^ 3) + p2 * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[4] * (1- u[4]) + sin(((4-1) * pi * t)/ T) * u[3] * (1- u[3]))
    ##Orden 2 diagonal
    du[m1+5] <- (p1 + p2 - 1) * u[m1+5] - (p1 + 3 * p2) * (u[m1+5] ^ 2) + 2 * p2 * (u[m1+5] ^ 3) + p2 * sqrt(2) *sqrt(2/T) * sin(((2-1) * pi * t)/ T) * u[2] * (1- u[2])
    du[m1+6] <- (p1 + p2 - 1) * u[m1+6] - (p1 + 3 * p2) * (u[m1+6] ^ 2) + 2 * p2 * (u[m1+6] ^ 3) + p2 * sqrt(2) *sqrt(2/T) * sin(((3-1) * pi * t)/ T) * u[3] * (1- u[3])
    #Orden 3 mixed
    du[m1+7] <- (p1 + p2 - 1) * u[m1+7] - (p1 + 3 * p2) * (u[m1+7] ^ 2) + 2 * p2 * (u[m1+7] ^ 3) + p2 * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[m1+6] * (1- u[m1+6]) + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[m1+2] * (1- u[m1+2]))
    du[m1+8] <- (p1 + p2 - 1) * u[m1+8] - (p1 + 3 * p2) * (u[m1+8] ^ 2 )+ 2 * p2 * (u[m1+8] ^ 3) + p2 * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[m1+2] * (1- u[m1+2]) + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[m1+5] * (1- u[m1+5]))
    ##Orden 3 diagonal
    du[m1+9] <- (p1 + p2 - 1) * u[m1+9] - (p1 + 3 * p2) * (u[m1+9] ^ 2 )+ 2 * p2 * (u[m1+9] ^ 3) + p2 * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+5] * (1- u[m1+5]))
    du[m1+10] <- (p1 + p2 - 1) * u[m1+10] - (p1 + 3 * p2) * (u[m1+10] ^ 2) + 2 * p2 *( u[m1+10] ^ 3) + p2 * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+6] * (1- u[m1+6]))
    #Orden 4
    du[m1+11] <- (p1 + p2 - 1) * u[m1+11] - (p1 + 3 * p2) * (u[m1+11] ^ 2) + 2 * p2 * (u[m1+11] ^ 3) + p2 * sqrt(2/T) * (sqrt(4) * sin(((2-1) * pi * t)/ T) * u[m1+9] * (1- u[m1+9]))
    du[m1+12] <- (p1 + p2 - 1) * u[m1+12] - (p1 + 3 * p2) * (u[m1+12] ^ 2 )+ 2 * p2 * (u[m1+12] ^ 3) + p2 * sqrt(2/T) * (sqrt(4) * sin(((3-1) * pi * t)/ T) * u[m1+10] * (1- u[m1+10]))
    #Orden 5
    du[m1+13] <- (p1 + p2 - 1) * u[m1+12] - (p1 + 3 * p2) * (u[m1+13] ^ 2) + 2 * p2 * (u[m1+13] ^ 3) + p2 * sqrt(2/T) * (sqrt(5) * sin(((2-1) * pi * t)/ T) * u[m1+11] * (1- u[m1+11]))
    du[m1+14] <- (p1 + p2 - 1) * u[m1+14] - (p1 + 3 * p2) * (u[m1+14] ^ 2) + 2 * p2 * (u[m1+14] ^ 3) + p2 * sqrt(2/T) * (sqrt(5) * sin(((3-1) * pi * t)/ T) * u[m1+12] * (1- u[m1+12]))
    #Orden 6
    du[m1+15] <- (p1 + p2 - 1) * u[m1+15] - (p1 + 3 * p2) * (u[m1+15] ^ 2) + 2 * p2 * (u[m1+15] ^ 3) + p2 * sqrt(2/T) * (sqrt(6) * sin(((2-1) * pi * t)/ T) * u[m1+13] * (1- u[m1+13]))
    du[m1+16] <- (p1 + p2 - 1) * u[m1+16] - (p1 + 3 * p2) * (u[m1+16] ^ 2) + 2 * p2 * (u[m1+16] ^ 3) + p2 * sqrt(2/T) * (sqrt(6) * sin(((3-1) * pi * t)/ T) * u[m1+14] * (1- u[m1+14]))
    #Orden 7
    du[m1+17] <- (p1 + p2 - 1) * u[m1+17] - (p1 + 3 * p2) * (u[m1+17] ^ 2) + 2 * p2 * (u[m1+17] ^ 3) + p2 * sqrt(2/T) * (sqrt(7) * sin(((2-1) * pi * t)/ T) * u[m1+15] * (1- u[m1+15]))
    du[m1+18] <- (p1 + p2 - 1) * u[m1+18] - (p1 + 3 * p2) * (u[m1+18] ^ 2 )+ 2 * p2 * (u[m1+18] ^ 3) + p2 * sqrt(2/T) * (sqrt(7) * sin(((3-1) * pi * t)/ T) * u[m1+16] * (1- u[m1+16]))
    #Orden 8
    du[m1+19] <- (p1 + p2 - 1) * u[m1+19] - (p1 + 3 * p2) * (u[m1+19] ^ 2) + 2 * p2 * (u[m1+19] ^ 3) + p2 * sqrt(2/T) * (sqrt(8) * sin(((2-1) * pi * t)/ T) * u[m1+17] * (1- u[m1+17]))
    du[m1+20] <-  (p1 + p2 - 1) * u[m1+20] - (p1 + 3 * p2) * (u[m1+20] ^ 2) + 2 * p2 * (u[m1+20] ^ 3) + p2 * sqrt(2/T) * (sqrt(8) * sin(((3-1) * pi * t)/ T) * u[m1+18] * (1- u[m1+18]))
    
    
    #return the rate of change
    list(c(du))  }) 
}
state <- u0
u=u0
times=seq(0,1,step)

solutionp <- function(state,times,alpha,sigma)
{
  sol_ode <- ode(y= state, times=times, func = ode_func,
                 parms = list(p1=alpha,p2=sigma))
  sol_ode
}

solu1=solutionp(state,times,alpha_0,sigma)



###one dimensional Brownian motion
BM <- function(n,TiempoC){
  MBM = matrix(0, n, length(TiempoC))
  l=length(TiempoC)
  delta=TiempoC[2]-TiempoC[1]
  for(i in 1:n){
    for(j in 2:l){
      MBM[i,j] = MBM[i,j-1]+rnorm(1)*sqrt(delta) 
    }
  }
  t(MBM)
}




solution1 <- function(Data1,M){
  n3=length(Data1[,1])  ###partition
  m3=length(Data1[,2])
  
  uu1=matrix(0, n3, M)
  for (k in 1:M) # Loop for trajectories
  {
    Xi = numeric(m1+20)
    Bt = BM(m1,Data1[,1])
    
    
    w1=matrix(0, n3, m1)
    for( i in 2:n3){
      for (j in 1:m1){
        w1[i,j]=sqrt(2/T) * sin((  pi*( (i)* Data1[i,1]))/ T)*(Bt[i,j]-Bt[i-1,j])
      }
    }
    z1=colSums(w1,na.rm=TRUE)
    
    Xi[1] = 1
    #order 1
    Xi[1:m1+1] = z1
    #order 2 mixed
    Xi[m1+2]= sqrt(2)*(z1[1]) * (z1[2])
    Xi[m1+3]= sqrt(2)*(z1[1]) * (z1[3])
    Xi[m1+4]= sqrt(2)*(z1[2]) * (z1[3])
    #order 2 diagonal
    a=((z1[1:2]^2)-1) 
    Xi[m1+5]=a[1]
    Xi[m1+6]=a[2]
    #order 3 mixed 
    Xi[m1+7]= sqrt(6)*(z1[1]) * ((z1[2]) ^ 2 - 1) /sqrt(2)
    Xi[m1+8]= sqrt(6)*(((z1[1]) ^ 2 - 1) * (z1[2]) )/sqrt(2)
    #order 3 diagonal
    b=(z1[1:2])*(((z1[1:2]^2)-3) )
    Xi[m1+9]= b[1]
    Xi[m1+10]= b[2]
    #order 4 diagonal
    c=((z1[1:2]^4)-(6*(z1[1:2]^2))+3) 
    Xi[m1+11]= c[1]
    Xi[m1+12]= c[2]
    #order 5 diagonal
    d=(z1[1:2]*(-(10*(z1[1:2]^2))+(z1[1:2]^4)+15))
    Xi[m1+13]= d[1]
    Xi[m1+14]= d[2]
    #order 6 diagonal
    e=(45*(z1[1:2]^2)-15*(z1[1:2]^4)+(z1[1:2]^6)-15)
    Xi[m1+15]= e[1]
    Xi[m1+16]= e[2]
    #order 7 diagonal
    f=(-105*z1[1:2]+105*(z1[1:2]^3)-21*(z1[1:2]^5)+(z1[1:2]^7))
    Xi[m1+17]= f[1]
    Xi[m1+18]= f[2] 
    #order 8 diagonal
    g=(-(420*(z1[1:2]^2))+(210*(z1[1:2]^4))-28*(z1[1:2]^6)+(z1[1:2]^8)+105)
    Xi[m1+19]= g[1]
    Xi[m1+20]= g[2] 
    
    for (i in 1:n3){
      #solution is X_0+ sum X_m int_0^T e_j dB
      uu1[i,k]=sum(Data1[i,2:length(Data1[1,])]*Xi,na.rm=TRUE) 
    }
  }
  return(uu1)
}

obs1=solution1(solu1,100)

variation_tox<-function(obs_1,dt){
  n=length(obs_1[,1])  ##rows
  m=length(obs_1[1,]) ##column
  z=matrix(0,n,m)
  z1=matrix(0,n,m)
  for (i in 2:n){
    z[i,]=(obs_1[i-1,]-obs_1[i,])^2
    z1[i,]=((obs_1[i,]-(obs_1[i,]^2))^2)*(dt) 
  }
  a=apply(z,1,cumsum)/apply(z1,1,cumsum)
  a=t(a)
  a[1,]=numeric(m)
  return(a)
  
}


a_var=variation_tox(obs1,0.001)


sigma_est=mean(sqrt(a_var[1001,]))

gradientcalcu<-function(obs,M,state,times,alpha, sigma){
  data_2=solutionp(state,times,alpha,sigma)
  ###we obtain the solution of the propgator of the alpha and sigma of the parameters
  
  ###we can calculate their observations but we don't use them in firther calculations
  #obs2=observations(M,data_2)
  
  ###calculate diffference between X_0^k- x^k    ,i.e., solution of the propagator of the new parameters versus the observations of the goal parameters
  #dif=reshape(repeat(data_2[:,2],M),length(data_2[:,2]),M)-obs
  me_obs=rowSums(obs, na.rm=TRUE)/M
  dif=data_2[,2]-me_obs
  
  #####generating vectors for the calculations
  data_21=solutionp(state,times,alpha*(1+epsi),sigma)
  data_22=solutionp(state,times,alpha,sigma*(1+epsi))
  derivative_w_alpha=(data_21-data_2)/(alpha*epsi)
  derivative_w_sigma=(data_22-data_2)/(sigma*epsi)
  #w=zeros(Float64,M)
  #z=zeros(Float64,M)
  
  #####-tX_0   partial derivative of X_0 with respecto alpha
  a=derivative_w_alpha[,2]
  
  
  a1=derivative_w_sigma[,2]
  #####-tX_m+sigma e^{-at}int_0^t s e^{as} sin(js/T) ds  partial derivative of X_m with respecto alpha
  #b=-data_2[:,1].*data_2[:,Not(1,2)]+y1
  b=derivative_w_alpha[,3:length(derivative_w_alpha[1,])]
  
  c=derivative_w_sigma[,3:length(derivative_w_sigma[1,])]
  #c=data_2[:,Not(1,2)]/sigma
  ######interval for the ingral calculation of int_0^T
  dt1=data_2[2,1]-data_2[1,1]
  
  
  vari=variation_tox(obs,dt1)
  mevari=rowSums(vari, na.rm=TRUE)/M
  #apply(vari,2,sum)/M
  
  
  w=sum((a*dif)*(dt1))+sum((rowSums(data_2[,3:length(data_2[1,])]^2,na.rm = TRUE)-mevari)*rowSums(data_2[,3:length(data_2[1,])]*b,na.rm = TRUE)*(dt1))
  z=sum((a1*dif)*(dt1))+sum((rowSums(data_2[,3:length(data_2[1,])]^2,na.rm = TRUE)-mevari)*rowSums(data_2[,3:length(data_2[1,])]*c,na.rm = TRUE)*(dt1))
  z1=sum((dif*dif)*(dt1))+sum(((rowSums(data_2[,3:length(data_2[1,])]^2,na.rm = TRUE)-mevari)^2)*(dt1))
  
  #end
  ###we include them in a vector that would be the gradient
  z_f=c(w,z,z1/2)*(2)
  return(z_f)
}
gradientcalcu(obs1,100,state,times,0.3,sigma_est)


gradientdescentv1 <- function(theta_1,n,obs,M,u0,Time)
{
  thetas=matrix(0,6,n)
  
  ####set the initial thetas for the gradient descent
  thetas[1:2,1]=theta_1[1:2]
  
  i=1
  gamma=1/2
  
  while(i<n){
    thetas[4:6,i]=gradientcalcu(obs,M,u0,Time,thetas[1,i], thetas[2,i])
    
    ###if loop to set the gamma,
    ##initial gamma is 0.5 and then according to expresion
    if(i==1)
    {gamma=1/2}  ###ver si con 1 o 1/2
    else {
      ###difference of gradient
      a_1=thetas[4:5,i]-thetas[4:5,i-1]
      ###difference of thetas
      b_1=thetas[1:2,i]-thetas[1:2,i-1]
      gamma=abs(sum(a_1*b_1)/sum(a_1*a_1))
    }
    
    ###we set gamma
    thetas[3,i]=gamma
    #####calculate the new thetas based on the gamma and gradient already calculated
    thetas[1:2,i+1]=thetas[1:2,i]-c(thetas[3,i],thetas[3,i])*thetas[4:5,i]
    #thetas[2,i+1]=theta_1[2]
    i=i+1
  } 
  #&& gamma>0.0001   possible additional condition to stop the while loop
  
  ###setting the gradient based on the observations and the initial thetas
  thetas[4:6,n]=gradientcalcu(obs,M,u0,Time,thetas[1,n], thetas[2,n])
  ###returns a matrix of 5xn n is the steps for the gradient descent
  return(thetas)  
}



resultado1=gradientdescentv1(c(0.4,sigma_est),6,obs1,100,state,times)