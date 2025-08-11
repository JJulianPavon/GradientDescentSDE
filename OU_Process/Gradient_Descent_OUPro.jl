#################Packages

using DifferentialEquations
using Plots
using DifferentialEquations.EnsembleAnalysis
using CSV, IterableTables, DataTables
using DataFrames
using Random, Distributions

####################### Parameters
alpha_0 = 1.7  #drift
sigma = 0.15 ###sigma
u_0 = 0.6  #####initial point
t0 = 0.0
T = 1.0
Time = (0.0,T)
dt = 0.001 ###partition


step = 0.001
t_s = range(0.0,T, step=step)
n=1000
u0 = zeros(n+1)
u0[1]=u_0

epsi=sqrt(eps(Float32)) ###### h for finite difference

M=100   #######Number of observations
############functions to generate  observations

#Propagator order 1
function Prop_order1_BM10(du,u,p,t)
    du[1] =-(p[1] ) * u[1] 
    #Orden 1
    for j in 2:n+1
        #####propagator for |m|=1
        du[j] = -(p[1])*u[j]  + sqrt(2/T) * sin(((j-1) * pi * t)/ T) * p[2]
    end

end


###Function to calculate solution based on the propagator
function solutionp(u0,Time,lambda,sigma)
    p2=[lambda,sigma]
    problem = ODEProblem(Prop_order1_MB1000,u0,Time,p2)
        Solution = solve(problem,AutoTsit5(Rosenbrock23()),saveat = step)
        df = DataFrame(Solution)
        Data = Matrix(df)
    return(Data)
end

###one dimensional Brownian motion
function BM(n,TiempoC)
    MBM = zeros(n,length(TiempoC))
    l=length(TiempoC)
    delta=TiempoC[2]-TiempoC[1]
    for i in 1:n
          for j in 2:l
              MBM[i,j] = MBM[i,j-1]+rand(Normal(0,1))*sqrt(delta)
          end
      end
    return(transpose(MBM))
end


###Wiener chaos expansion
function solution1(Data1)
    ###brownian motion
    My_BM=BM(1,Data1[:,1])
    n1=size(Data1)[2]
    m1=size(Data1)[1]
    w1=zeros(Float64,(m1,n1-2))
    #w will be the sci of wiener chaos expansion collection, but given the sci are stochastic integrals
    ###we calculate  w[j,i]=e_j (t_i) (B_i-B_{i-1})
    ####therefore, sum(w[j,:])=int_0^T e_j dB approxime by sum_i e_j (t_i) (B_i-B_{i-1})
    
    for i in 1:n1-2
        for j in 2:m1
            w1[j,i]=sqrt(2/T) * sin(((i) * pi * Data1[j,1])/ T)*(My_BM[j]-My_BM[j-1])
        end
    end
    z1=sum(eachrow(w1))
    u1=zeros(m1)
    for i in 2:m1
        #solution is X_0+ sum X_m int_0^T e_j dB
        u1[i]=Data1[i,2].+sum(Data1[i,Not(1,2)].*z1)
    end
    u1[1]=Data1[1,2]
    return(u1)
end

function observations(M,Data)
    w=zeros(Float64,(size(Data)[1],M))
    for i in 1:M
        w[:,i]=solution1(Data)
    end
    return(w)
end
###



###Quadratic variation calculation
function variation_to(obs_1,dt)
    n=length(obs_1[:,1])  ##rows
    m=length(obs_1[1,:]) ##column
    z=zeros((n,m))
    z1=zeros((n,m))
    for i in 1:n-1
        z[i+1,:]=(obs_1[i,:]-obs_1[i+1,:]).^2
        z1[i+1,:]=repeat([dt],m)
    end
    a=cumsum(z,dims=1)./cumsum(z1,dims=1)
    a[1,:]=zeros(m)
    return(a)
    
end
#####Gradient calculator and loss function on given point
function gradientcalcu(data,obs,M,u0,Time,alpha, sigma)
    data_2=solutionp(u0,Time,alpha,sigma)
    ###we obtain the solution of the propagator of the alpha and sigma of the parameters

    ###we can calculate their observations but we don't use them in further calculations
    #obs2=observations(M,data_2)

    ###calculate diffference between X_0^k- avg( x^k)    ,i.e., solution of the propagator of the new parameters versus the average of the observations of the goal parameters
    me_obs=sum(eachcol(obs))/M
    dif=data_2[:,2]-me_obs


    #####generating vectors for the calculations
    ###finite difference to calculate derivatives with respect both parameters
    data_21=solutionp(u0,Time,alpha*(1+epsi),sigma)
    derivative_w_alpha=(data_21-data_2)/(alpha*epsi)

    
    a=derivative_w_alpha[:,2]

    b=derivative_w_alpha[:,Not(1,2)]

    ######interval for the ingral calculation of int_0^T
    dt1=data_2[2,1]-data_2[1,1]

    ##gradient calculate on each paramaeter and loss function
  
    w=sum((a.*dif)*(dt1))
    
    fun=sum((dif.*dif)*(dt1))
  ###we include them in a vector that would be the gradient
    z_f=[w,fun/2]*(2)
    return(z_f)
end

##gradient descent function
function gradientdescentv1(theta_1,n,data,obs,M,u0,Time)
    thetas=zeros(Float64,(5,n))
    
    ####set the initial thetas for the gradient descent
    thetas[1:2,1]=theta_1

    i=1
    gamma=1/2
    
    while i<n 
        ###setting the gradient based on the observations and the initial thetas
        thetas[4:5,i]=gradientcalcu(data,obs,M,u0,Time,thetas[1,i], thetas[2,i])

        ###if loop to set the gamma,
        ##initial gamma is 0.5 and then according to expresion
        if i==1 
            gamma=1/2  
        else 
            ###difference of gradient
            a_1=thetas[4,i]-thetas[4,i-1]
            ###difference of thetas
            b_1=thetas[1,i]-thetas[1,i-1]
            gamma=(a_1*b_1)/(a_1*a_1)
        end

        ###we set gamma
        thetas[3,i]=gamma
        #####calculate the new thetas based on the gamma and gradient already calculated
        thetas[1,i+1]=thetas[1,i]-thetas[3,i]*thetas[4,i]
        thetas[2,i+1]=thetas[2,i]
        i=i+1
    end
    thetas[4:5,n]=gradientcalcu(data,obs,M,u0,Time,thetas[1,n], thetas[2,n])
    ###returns a matrix of 5xn n is the steps for the gradient descent
    return(thetas)
end


##############################Example

thetas_obj=[alpha_0,sigma] ##setting observation parameters
###Calculate propagator of observations
data_1=solutionp(u0,Time,thetas_obj[1],thetas_obj[2])
###observations
obs1=solution1(data_1,M)


###estimated sigma based on quadratic variation
a=variation_to(obs1,dt)
sigma_est=mean(sqrt.(a[1001,:]))

####Gradient descent example
gradientdescentv1([1.01,sigma_est],6,data_1,obs_11,M,u0,Time)

