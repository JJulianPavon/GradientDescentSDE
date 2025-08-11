#################Packages

using DifferentialEquations
using Plots
using DifferentialEquations.EnsembleAnalysis
using CSV, IterableTables, DataTables
using DataFrames
using Random, Distributions

####################### Parameters
alpha_0 = 0.9 #drift
lambda = alpha_0
sigma = 0.115 ##sigma
theta = 1.3
eta = 0.8 
u_0 = eta #initial point
t0 = 0.0
T = 1.0
Time = (0.0,T)
dt = 0.001  ###partition



step = 0.001
t_s = range(0.0,T, step=step)

m1=1000

u0 = zeros(m1+20)
u0[1]=u_0

epsi=sqrt(eps(Float32)) ###### h for finite difference

M=100   #######Number of observations
############functions to generate  observations

#Propagator order 8 
function Prop_order8_BM1000(du,u,p,t)
    du[1] = p[1] * u[1]
    #Orden 1
    for j in 2:m1+1
        du[j] = p[1] * u[j] +  sqrt(2/T) * sin(((j-1) * pi * t)/ T) * p[2] * u[1]
    end
    ##Orden 2 mixed
    du[m1+2] = p[1] * u[m1+2] + p[2] * sqrt(2/T) * (sin(((3-1) * pi * t)/ T) * u[2] + sin(((2-1) * pi * t)/ T) * u[3])
    du[m1+3] = p[1] * u[m1+3] + p[2] * sqrt(2/T) * (sin(((4-1) * pi * t)/ T) * u[2]+ sin(((2-1) * pi * t)/ T) * u[4])
    du[m1+4] = p[1] * u[m1+4] + p[2] * sqrt(2/T) * (sin(((4-1) * pi * t)/ T) * u[3] + sin(((3-1) * pi * t)/ T) * u[4])
    ##Orden 2 diagonal
    for j in 2:4
        du[m1+3+j] = p[1] * u[m1+3+j] + p[2] * sqrt(2) *sqrt(2/T) * sin(((j-1) * pi * t)/ T) * u[j]
    end
    #Orden 3 mixed
    du[m1+8] = p[1] * u[m1+8] + p[2] * sqrt(2/T) * (sin(((2-1) * pi * t)/ T) * u[m1+6] + sqrt(2)* sin(((3-1) * pi * t)/ T) * u[m1+2])
    du[m1+9] = p[1] * u[m1+9] + p[2] * sqrt(2/T) * (sqrt(2) * sin(((2-1) * pi * t)/ T) * u[m1+2] + sqrt(1)* sin(((3-1) * pi * t)/ T) * u[m1+5])
    ##Orden 3 diagonal
    du[m1+10] = p[1] * u[m1+10] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+5])
    du[m1+11] = p[1] * u[m1+11] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+6])
    #Orden 4
    du[m1+12] = p[1] * u[m1+12] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+10])
    du[m1+13] = p[1] * u[m1+13] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+11])
    #Orden 5
    du[m1+14] = p[1] * u[m1+14] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+12])
    du[m1+15] = p[1] * u[m1+15] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+13])
    #Orden 6
    du[m1+16] = p[1] * u[m1+16] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+14])
    du[m1+17] = p[1] * u[m1+17] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+15])
    #Orden 7
    du[m1+18] = p[1] * u[m1+18] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+16])
    du[m1+19] = p[1] * u[m1+19] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+17])
    #Orden 8
    du[m1+20] = p[1] * u[m1+20] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((2-1) * pi * t)/ T) * u[m1+18])
    du[m1+21] = p[1] * u[m1+21] + p[2] * sqrt(2/T) * (sqrt(3) * sin(((3-1) * pi * t)/ T) * u[m1+19])
end


###Function to calculate solution based on the propagator
function solutionp(u0,Time,lambda,sigma)
    p2=[lambda,sigma]
    problem = ODEProblem(Prop_order8_MB1000,u0,Time,p2)
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
function solution1(Data1,M)
    n3=size(Data1)[1]  ###partition
    m3=size(Data1)[2]
    
    u1=zeros(Float64,(n3,M))
    for k in 1:M # Loop for trajectories
        Xi = zeros(m1+21)
        Bt = BM(m1,Data1[:,1])
        

        w1=zeros(Float64,(n3,m1))
        for i in 2:n3
            for j in 1:m1
                w1[i,j]=sqrt(2/T) * sin((  pi*( (i)* Data1[i,1]))/ T)*(Bt[i,j]-Bt[i-1,j])
            end
        end
        z1=sum(eachrow(w1))

        Xi[1] = 1
        #order 1
        Xi[2:m1+1] = z1
        #order 2 mixed
        Xi[m1+2]= sqrt(2)*(z1[1]) * (z1[2])
        Xi[m1+3]= sqrt(2)*(z1[1]) * (z1[3])
        Xi[m1+4]= sqrt(2)*(z1[2]) * (z1[3])
        #order 2 diagonal 
        Xi[m1+5:m1+7]= ((z1[1:3].^2).-1) 
        #order 3 mixed 
        Xi[m1+8]= sqrt(6)*(z1[1]) * ((z1[2]) ^ 2 - 1) /sqrt(2)
        Xi[m1+9]= sqrt(6)*(((z1[1]) ^ 2 - 1) * (z1[2]) )/sqrt(2)
        #order 3 diagonal
        Xi[m1+10:m1+11]= (z1[1:2]).*(((z1[1:2].^2).-3) )
        #order 4 diagonal
        Xi[m1+12:m1+13]=  ((z1[1:2].^4).-(6*(z1[1:2].^2)).+3) 
        #order 5 diagonal
        Xi[m1+14:m1+15]= (z1[1:2].*(-(10*(z1[1:2].^2)).+(z1[1:2].^4).+15))
        #order 6 diagonal
        Xi[m1+16:m1+17]= (45*(z1[1:2].^2).-15*(z1[1:2].^4).+(z1[1:2].^6).-15)
        #order 7 diagonal
        Xi[m1+18:m1+19]= (-105*z1[1:2].+105*(z1[1:2].^3).-21*(z1[1:2].^5).+(z1[1:2].^7))
        #order 8 diagonal
        Xi[m1+20:m1+21]= (-(420*(z1[1:2].^2)).+(210*(z1[1:2].^4)).-28*(z1[1:2].^6).+(z1[1:2].^8).+105)

        

        for i in 1:n3
            #solution is X_0+ sum X_m int_0^T e_j dB
            u1[i,k]=sum(Data1[i,Not(1)].*Xi)
        end
        
            
    end
    return(u1)
end


###Quadratic variation calculation
function variation_to(obs_1,dt)
    n=length(obs_1[:,1])  ##rows
    m=length(obs_1[1,:]) ##column
    z=zeros((n,m))
    z1=zeros((n,m))
    for i in 1:n-1
        z[i+1,:]=(obs_1[i,:]-obs_1[i+1,:]).^2
        z1[i+1,:]=((obs_1[i+1,:].-(obs_1[i+1,:].^2)).^2)*(dt)
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
    data_22=solutionp(u0,Time,alpha,sigma*(1+epsi))
    derivative_w_alpha=(data_21-data_2)/(alpha*epsi)
    derivative_w_sigma=(data_22-data_2)/(sigma*epsi)

    
    a=derivative_w_alpha[:,2]


    a1=derivative_w_sigma[:,2]
    b=derivative_w_alpha[:,Not(1,2)]

    c=derivative_w_sigma[:,Not(1,2)]
    
    ######interval for the ingral calculation of int_0^T
    dt1=data_2[2,1]-data_2[1,1]

    ####quadratic variation

    vari=variation_to(obs,dt1)
    mevari=sum(eachcol(vari))/M
    ##gradient calculate on each paramaeter and loss function
  
    w=sum((a.*dif)*(dt1))+sum(sum(eachcol(((sum(eachcol(data_2[:,Not(1,2)].^2))-mevari).*(sum(eachcol(data_2[:,Not(1,2)].*b)))))*(dt1)))
    z=sum((a1.*dif)*(dt1))+sum(sum(eachcol(((sum(eachcol(data_2[:,Not(1,2)].^2))-mevari).*(sum(eachcol(data_2[:,Not(1,2)].*c)))))*(dt1)))
    fun=sum((dif.^2)*(dt1))*sum(sum(eachcol(((sum(eachcol(data_2[:,Not(1,2)].^2))-mevari).^2)))*(dt1))
    ###we include them in a vector that would be the gradient
    z_f=[w,z,fun/2]*(2)
    return(z_f)
end

##gradient descent function
function gradientdescentv1(theta_1,n,data,obs,M,u0,Time)
    thetas=zeros(Float64,(6,n))
    
    ####set the initial thetas for the gradient descent
    thetas[1:2,1]=theta_1

    i=1
    gamma=1/2
    
    while i<n 
        ###setting the gradient based on the observations and the initial thetas
        thetas[4:6,i]=gradientcalcu(data,obs,M,u0,Time,thetas[1,i], thetas[2,i])

        ###if loop to set the gamma,
        ##initial gamma is 0.5 and then according to expresion
        if i==1 
            gamma=1/2  
        else 
            ###difference of gradient
            a_1=thetas[4:5,i]-thetas[4:5,i-1]
            ###difference of thetas
            b_1=thetas[1:2,i]-thetas[1:2,i-1]
            gamma=abs((a_1'b_1)/(a_1'a_1))
        end

        ###we set gamma
        thetas[3,i]=gamma
        #####calculate the new thetas based on the gamma and gradient already calculated
        thetas[1:2,i+1]=thetas[1:2,i]-[thetas[3,i],thetas[3,i]].*thetas[4:5,i]
        i=i+1
    end
    thetas[4:6,n]=gradientcalcu(data,obs,M,u0,Time,thetas[1,n], thetas[2,n])
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
gradientdescentv1([0.4,sigma_est],6,data_1,obs1,M,u0,Time)
