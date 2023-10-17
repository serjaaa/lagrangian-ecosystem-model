

#
#
#
from __future__ import division
import numpy as np
import sys
#
#
#


## ###################### GENERAL PARAMETERS ######################

simulation_name = ''
output_dir = ''


## ###################### TEMPORAL PARAMETERS ######################

## integration time step duration dt [days]
dt = 0

## initial time [days]
ini_t = 0.0

## final time [days]
fin_t = 0


## ###################### INITIAL CONDITIONS ######################

## initial position of the patch center [km]
X0 = np.array([0.0, 0.0])

## initial patch length, width, thickness [km]
W0 = 0
L0 = 0
T0=0

## initial values of strain, horizontal diffusion (given as internal parameters or command line arguments)
gamma_ini = float(sys.argv[1])
kappa_h_ini = float(sys.argv[2])

## Fe:C and N:C ratio [moles/moles]
fe_to_c_ratio = 0
n_to_c_ratio = 0

## Chl:C ratio [mg/mg]
c_to_chl_ratio = 0

## initial iron distributions mean [micro-moles / m^3] in iron currency
ini_mean1 = 0

## initial phyto distributions mean [micro-moles / m^3] in iron currency 
ini_mean2 = 0

## initial distribution variances [micro-moles / m^3]*[micro-moles / m^3]
ini_var1_h = 1e-8
ini_var2_h = 1e-8

## initial distribution covariance [micro-moles / m^3]*[micro-moles / m^3]
ini_cov_h = 1e-8


## ###################### EXTERNAL FIELDS ######################


#### PHYSICS #### 

## calculating constants for strain and horizontal diffusion at time 0 following expected scaling relationships
alpha = gamma_ini * ( (W0+L0)**(2/3) )
beta = kappa_h_ini / (W0+L0)

## velocity at the patch center [km * days^-1]
def vel(time,position,size_horizontal):
    center_vel = np.array([0.0, 0.0])  
    return(center_vel)

## horizontal strain rate [day^-1]
def gamma(time,position,size_horizontal):
    strain = alpha * (S_h**(-2/3))
    return(strain)

## horizontal diffusion  [km^2 * day^-1]
def kappa_h(time,position,size_horizontal):
    diffusion = beta * S_h
    return(diffusion)

## internal horizontal homogeneization rate [day^-1]
def lambd_h(time,position,size_horizontal):
    homog_rate =  kappa_h(t,X,S_h) / (S_h**2)
    return(homog_rate)



#### SURROUNDING DISTRIBUTIONS #### 

## surrounding distributions means [micro-moles / m^3]
def suh_mean1(time,position,size_horizontal):
    sur_mean = 0
    return(sur_mean)

def suh_mean2(time,position,size_horizontal):
    sur_mean = 0
    return(sur_mean)


## surrounding distributions variances  [micro-moles / m^3]*[micro-moles / m^3]
def suh_var1(time,position,size_horizontal):
    sur_var =  1e-8
    return(sur_var)

def suh_var2(time,position,size_horizontal):
    sur_var =  1e-8
    return(sur_var)

## surrounding distributions covariance [micro-moles / m^3]*[micro-moles / m^3]
def suh_cov(time,position,size_horizontal,size_vertical):
    sur_covariance = 1e-8
    return(sur_covariance)


#### BIO-GEOCHEMICAL DYNAMICS ####

## external source/sink of resource (i.e. type 1) [micro-moles / (m^3*day)]
def e(time,position,size_horizontal):
    res_source = 0.0
    return(res_source)

## maximum growth rate for the consumer (i.e. type 2) [day^-1]
def nu(time,position,size_horizontal):
    max_grow = 0
    return(max_grow)

## mortality rate for the consumer (i.e. type 2) [day^-1] (note that if quadratic it changes dimension to [(m^3/micro-moles)/day])
def m(time,position,size_horizontal):
    mort = 0
    return(mort)

## half saturation constant for resource and consumer (i.e. type 1 and 2) [micro-moles / m^3] 
def k(time,position,size_horizontal):
    half_sat = 0
    return(half_sat)





#######
#######
#######
#######




## ###################### VARIABLES INITIALIZATION ######################

## time [day]
t = ini_t

## patch center [km]
X = X0

## length width and thickness [km]
W = W0
L = L0
T = T0

## sizes and volume [km]
S_h = W0 + L0
S_v = 2*T0
V = (4 * np.pi * W0 * L0 * T0)/3

## distributions means [micro-moles / m^3]
mean1 = ini_mean1
mean2 = ini_mean2

## distributions variances [micro-moles / m^3]*[micro-moles / m^3]
var1_h = ini_var1_h 
var2_h = ini_var2_h 
var1 = var1_h 
var2 = var2_h

## distribution covariance [micro-moles / m^3]*[micro-moles / m^3]
cov_h = ini_cov_h
cov =  cov_h 


## ###################### OPENING OUTPUT FILE ######################

## output path
output_path = str(output_dir + simulation_name + '_'
                  + 'gamma' + '{:.1e}'.format(gamma_ini) + '_'
                  + 'kappah' + '{:.1e}'.format(kappa_h_ini) + '_'
                  + '.dat')

print(output_path)

## open file
file = open(output_path, "w")

## writing file header
file.write("t" + " " )
file.write("W" + " " )
file.write("L" + " " )
file.write("V" + " " )
file.write("S_h" + " " )
file.write("mean1" + " " )
file.write("mean2" + " " )
file.write("var1_h" + " " )
file.write("var2_h" + " " )
file.write("cov" + " " )
file.write('\n')




## ###################### EVOLUTION LOOP ######################


## while temporal loop
while t <= fin_t:


    #### PHYSICS #### 
    
    ## patch center equation
    X[0] = X[0] + vel(t,X,S_h,S_v)[0] * dt
    X[1] = X[1] + vel(t,X,S_h,S_v)[1] * dt
    
    ## width equation
    W = W + ( kappa_h(t,X,S_h)/W - gamma(t,X,S_h)*W ) * dt

    ## lenght equation
    L = L + ( kappa_h(t,X,S_h)/L + gamma(t,X,S_h)*L ) * dt

    ## horizontal volume derivative
    V_dev_h = ((4 * np.pi)/3) * kappa_h(t,X,S_h) * T * ( (W*W + L*L) / (W*L) )
    V_dev = V_dev_h

    ## volume equation [km^3]
    V = V + V_dev * dt

    ## sizes equation [km]
    S_h = W + L
    S_v = 2*T

    
    #### TRACERS ####

    ## equation for mean type 1
    mean1_sur_term = (1/V) * ( V_dev_h*suh_mean1(t,X,S_h) - V_dev*mean1 )
    
    mean1_bio_term = ( - nu(t,X,S_h) * ( (mean1*mean2) / (mean1+k(t,X,S_h)) )
                  + nu(t,X,S_h)*k(t,X,S_h) * ( (mean2*var1)/((mean1+k(t,X,S_h))**3) )
                  - nu(t,X,S_h)*k(t,X,S_h) * ( (cov)/((mean1+k(t,X,S_h))**2) )
                  + e(t,X,S_h) )
    
    mean1 = mean1 + (mean1_sur_term + mean1_bio_term)*dt

    
    ## equation for mean type 2
    mean2_sur_term = (1/V) * ( V_dev_h*suh_mean2(t,X,S_h,S_v) - V_dev*mean2 )

    mean2_bio_term = ( nu(t,X,S_h) * ( (mean1*mean2) / (mean1+k(t,X,S_h)) )
                  - nu(t,X,S_h)*k(t,X,S_h) * ( (mean2*var1)/((mean1+k(t,X,S_h))**3) )
                  + nu(t,X,S_h)*k(t,X,S_h) * ( (cov)/((mean1+k(t,X,S_h))**2) )
                  - m(t,X,S_h)*mean2 )

    mean2 = mean2 + (mean2_sur_term + mean2_bio_term)*dt


    ## positiveness conditions    
    if mean1 < 0 :
        mean1 = 1e-8
    if mean2 < 0 :
        mean2 = 1e-8
    if var1_h < 0 :
        var1_h = 1e-8
    if var2_h < 0 :
        var2_h = 1e-8
    
    ## equation for variance type 1
    var1_sur_h_term = (1/V) * V_dev_h * ( suh_var1(t,X,S_h) - var1 + (suh_mean1(t,X,S_h) - mean1)**2 )
        
    var1_h = var1_h + (var1_sur_h_term)*dt
    
    var1_homo_h_term = - lambd_h(t,X,S_h) * var1_h
        
    var1_h = var1_h + (var1_homo_h_term)*dt
    
    var1 = var1_h 

    var1_bio_term = ( - 2*nu(t,X,S_h)*k(t,X,S_h) * ( (mean2*var1) / ((mean1+k(t,X,S_h))**2) )
                      - 2*nu(t,X,S_h) * ( (mean1*cov) / (mean1+k(t,X,S_h)) ) )

    var1_h = var1_h * ( (var1 + var1_bio_term*dt)/var1 )
                                            
    var1 = var1_h                                        

    ## equation for variance type 2
    var2_sur_h_term = (1/V) * V_dev_h * ( suh_var2(t,X,S_h) - var2 + (suh_mean2(t,X,S_h) - mean2)**2 )
    
    var2_h = var2_h + (var2_sur_h_term)*dt
    
    var2_homo_h_term = - lambd_h(t,X,S_h) * var2_h
        
    var2_h = var2_h + (var2_homo_h_term)*dt
    
    var2 = var2_h 

    var2_bio_term = ( 2*nu(t,X,S_h) * ( (mean1*var2) / (mean1+k(t,X,S_h)) )
                      + 2*nu(t,X,S_h)*k(t,X,S_h) * ( (mean2*cov) / ((mean1+k(t,X,S_h))**2) )
                      - 2*m(t,X,S_h)*var2 )

    var2_h = var2_h * ( (var2 + var2_bio_term*dt)/var2 )
                                             
    var2 = var2_h 
    
    ## equation for covariance
    cov_sur_h_term = ( (1/V) * V_dev_h * ( suh_cov(t,X,S_h) - cov + (suh_mean1(t,X,S_h) - mean1)*(suh_mean2(t,X,S_h) - mean2) ) )
    
    cov_h = cov_h + (cov_sur_h_term)*dt
    
    cov_homo_h_term = - lambd_h(t,X,S_h) * cov_h
            
    cov_h = cov_h + (cov_homo_h_term)*dt
    
    cov = cov_h
    
    cov_bio_term = ( nu(t,X,S_h) *  ( (mean1) / (mean1+k(t,X,S_h)) ) * ( cov - var2 )
                + nu(t,X,S_h)*k(t,X,S_h) * ( (mean2) / ((mean1+k(t,X,S_h))**2) )  * ( var1 - cov)
                - m(t,X,S_h)*cov  )
        
    cov_h = cov_h * ( (cov + cov_bio_term*dt)/cov )

    cov = cov_h 
    
    ## positiveness conditions    
    if mean1 < 0 :
        mean1 = 1e-8
    if mean2 < 0 :
        mean2 = 1e-8
    if var1_h < 0 :
        var1_h = 1e-8
    if var2_h < 0 :
        var2_h = 1e-8



    #### WRITING OUTPUT #### 
    
    ## writing variables for the current time
    file.write(str(t) + " " )
    file.write(str(W) + " " )
    file.write(str(L) + " " )
    file.write(str(V) + " " )
    file.write(str(S_h) + " " )
    file.write(str(mean1) + " " )
    file.write(str(mean2) + " " )
    file.write(str(var1_h) + " " )
    file.write(str(var2_h) + " " )
    file.write(str(cov) + " " )
    file.write('\n')


    
    #### STEPPING TIME #### 
    
    ## time increment
    t = t + dt

    
## closing output file
file.close()


