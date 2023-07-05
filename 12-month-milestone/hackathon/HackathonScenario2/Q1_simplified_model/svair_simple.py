import numpy as np
import numpy.matlib
from get_beta import *

def svair_simple(y, t, beta,beta_v1,beta_v2,beta_R,ai_beta_ratio,gamma,nu_v1,nu_v2,nu_R,ai,ai_V,ai_R,mu,mu_I,mu_IV): 
    # retrieve current populations
    n_var = 1;  # number of variants simulated
    n_vax = 1;  # number of vaccines
    ind = 0
    
    # original susceptible
    S = y[ind] 
    ind = ind + 1
    
    # lost immunity after vaccination or recovery
    SVR = y[ind]
    ind = ind + 1
    
    # one-dose vaccination
    V1 = y[ind:ind+n_vax]
    ind = ind + n_vax
    
    # fully vaccinated
    V2 = y[ind:ind+n_vax]
    ind = ind + n_vax
    
    # infected
    I = y[ind:ind+n_var]
    ind = ind + n_var
    
    # infected even with vaccination
    IV = y[ind:ind+n_var]
    ind = ind + n_var
    
    # infected again after recovery from a different variant
    IR = y[ind:ind+n_var]
    ind = ind + n_var
    
    # asymptomatic infections
    A = y[ind:ind+n_var]
    ind = ind + n_var
    
    # asymptomatic infections after recovery from a different variant
    AR = y[ind:ind+n_var]
    ind = ind + n_var
    
    # recovered
    R = y[ind:ind+n_var]
    ind = ind + n_var
    
    # recovered after getting both variants
    R2 = y[ind]
    ind = ind + 1
    
    # get time-dependent parameters
    # no vaccine for the first year
    # vr1,vr2 = get_vaccine_rate(t)
    vr1 = 0
    vr2 = 0

    beta_scale = get_beta(t)
    beta = beta * beta_scale
    beta_v1 = beta_v1 * beta_scale
    beta_v2 = beta_v2 * beta_scale
    beta_R = beta_R * beta_scale
    
    # total infectious population
    I_total = I + IV + IR + ai_beta_ratio*(A + AR)
    # need the following to compute infection of recovered from another variant
    # mm = np.ones((nv + 1,nv + 1)) - diag(np.ones((1,nv + 1)))
    # mv = np.multiply(mm,np.matlib.repmat(R,1,nv + 1))
    Rv = 0 #np.transpose(sum(mv)); in simple case, don't consider infection from a different variant
    
    # compute time derivatives
    dSdt = - np.sum(beta*S*I_total) - np.sum(vr1*S) + mu * (1 - S)
    dSVRdt = + nu_v1*np.sum(V1) + nu_v2*np.sum(V2) + np.sum(nu_R*R) + nu_R*R2 - np.sum(beta*SVR*I_total) - mu*SVR;
    dV1dt  = + vr1*(S+np.sum(A)) - vr2*V1 - nu_v1*V1 - np.sum(beta_v1*V1*I_total.T) - mu*V1;
    dV2dt  = + vr2*V1 - nu_v2*V2 - np.sum(beta_v2*V2*I_total.T) - mu*V2;
    dIdt   = (1-ai)*(beta*S*I_total + beta*SVR*(I_total)) - gamma*I - mu_I*I;
    dIVdt  = (1-ai_V)*(np.sum(beta_v1*V1*I_total.T) + np.sum(beta_v2*(V2*I_total.T))) - gamma*IV - mu_IV*IV;
    dIRdt  = (1-ai_R)*(+beta_R*Rv*(I_total)) - gamma*IR - mu*IR;
    dAdt   = ai*(beta*S*I_total + beta*SVR*I_total) + ai_V*(np.sum(beta_v1*V1*I_total.T) + np.sum(beta_v2*V2*I_total.T)) - vr1*A - gamma*A - mu*A;
    dARdt  = ai_R*(beta_R*Rv*I_total) - gamma*AR - mu*AR;
    dRdt   = + gamma*I_total - nu_R*R - beta_R*Rv*I_total - mu*R;
    dR2dt  = + np.sum(gamma*(IR+AR)) - nu_R*R2 - mu*R2;
    
    yp = np.vstack((dSdt,dSVRdt, dV1dt, dV2dt, dIdt, dIVdt, dIRdt, dAdt, dARdt, dRdt, dR2dt))
    yp = yp.reshape(len(yp),)
    
    # yp = np.array([[dSdt],[dSVRdt],[dV1dt],[dV2dt],[dIdt],[dIVdt],[dIRdt],[dAdt],[dARdt],[dRdt],[dR2dt]])
    return yp