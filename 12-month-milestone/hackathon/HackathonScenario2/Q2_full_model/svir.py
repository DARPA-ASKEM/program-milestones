import numpy as np
from get_beta import *
from get_vaccine_rate import *


def svir(y,t,beta,beta_v1,beta_v2,beta_R,gamma,nu_v1,nu_v2,nu_R,mu,mu_I,mu_IV): 
    
    # retrieve current populations
    n_var = 3;  # number of variants simulated
    n_vax = 2;  # number of vaccines
    ind = 0
    
    # original susceptible
    S = y[ind].reshape((1,1))
    ind = ind + 1

    # lost immunity after vaccination or recovery
    SVR = y[ind].reshape((1,1))
    ind = ind + 1

    # one-dose vaccination
    V1 = y[ind:ind+n_vax].reshape((n_vax,1))
    ind = ind + n_vax

    # fully vaccinated
    V2 = y[ind:ind+n_vax].reshape((n_vax,1))
    ind = ind + n_vax

    # infected
    I = y[ind:ind+n_var].reshape((n_var,1))
    ind = ind + n_var

    # infected even with vaccination
    IV = y[ind:ind+n_var].reshape((n_var,1))
    ind = ind + n_var

    # infected again after recovery from a different variant
    IR = y[ind:ind+n_var].reshape((n_var,1))
    ind = ind + n_var

    # recovered
    R = y[ind:ind+n_var].reshape((n_var,1))
    ind = ind + n_var

    # recovered after getting both variants
    R2 = y[ind].reshape((1,1))
    ind = ind + 1

    # get time-dependent parameters
    (vr1,vr2) = get_vaccine_rate(t)
    beta_scale = get_beta(t)
    beta = beta * beta_scale
    beta_v1 = beta_v1 * beta_scale
    beta_v2 = beta_v2 * beta_scale
    beta_R = beta_R * beta_scale

    # compute time derivatives
    dSdt =  -np.sum(beta*S*(I+IV+IR)) - np.sum(vr1*S) + mu*(1-S)
    dSVRdt = nu_v1*np.sum(V1) + nu_v2*np.sum(V2) + np.sum(nu_R*R) + nu_R*R2 - np.sum(beta*SVR*(I+IV+IR)) - mu*SVR
    dV1dt  = vr1*S - vr2*V1 - nu_v1*V1 - np.sum(beta_v1*(V1@(I+IV+IR).T),1, keepdims = True) - mu*V1
    dV2dt  = vr2*V1 - nu_v2*V2 - np.sum(beta_v2*(V2@(I+IV+IR).T),1, keepdims = True) - mu*V2
    dIdt   = beta*S*(I+IV+IR) + beta*SVR*(I+IV+IR) - gamma*I - mu_I*I
    dIVdt  = np.sum(beta_v1*(V1@(I+IV+IR).T),0, keepdims = True).T + np.sum(beta_v2*(V2@(I+IV+IR).T),0,keepdims = True).T - gamma*IV - mu_IV*IV
    dIRdt  = beta_R*np.flip(R)*(I+IV+IR) - gamma*IR - mu*IR
    dRdt   = gamma*(I+IV+IR) - nu_R*R - beta_R*np.flip(R)*(I+IV+IR) - mu*R
    dR2dt  = np.sum(gamma*IR) - nu_R*R2 - mu*R2

    # simulate mutation of WT into variants and importation
    if (t>300) and (t<360): # alpha appears
        dWTtoA = 0.0 * dIdt[0] # don't assume mutation
        dIdt[0] = dIdt[0] - dWTtoA
        dIdt[1] = dIdt[1] + dWTtoA + ((t - 250) / 60) * 10 / 14570000
    elif (t>385) and (t<445):
        dIdt[2] = ((t - 385) / 60) * 40 / 14570000

    yp = np.vstack((dSdt,dSVRdt, dV1dt, dV2dt, dIdt, dIVdt, dIRdt, dRdt, dR2dt))
    yp = yp.reshape(len(yp),)
    
    return yp


