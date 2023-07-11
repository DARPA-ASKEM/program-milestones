import numpy as np
from datetime import date
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('svg')
from get_beta import *
from get_vaccine_rate import *
from svair import *

# define integration interval
dt = 1
tstart = 0
tfinal = 2.0 * 365
t = np.arange(tstart, tfinal, dt) 

# Ontario population, 14.57M 2019
N = 14570000

# define parameters
# these are infection rates without lockdown
beta = np.array([[3.3e-09],[5.5e-09],[7.6e-09]]) * N

# beta_v1 = [0.3;0.5].*beta; # infection rate, first dose
# beta_v2 = 0.05*beta; # infection rate, both doses
# https://www.gov.uk/government/news/vaccines-highly-effective-against-b-1-617-2-variant-after-2-doses
beta_v1 = np.array([[0.2,0.5,0.67],[0.2,0.5,0.67]])*np.vstack((beta.T,beta.T))
beta_v2 = np.array([[0.05,0.07,0.12],[0.05,0.34,0.4]])*np.vstack((beta.T,beta.T))
beta_R = 0.05 * beta

ai_beta_ratio = np.array([[3],[3],[3]])

# vaccination rates now given as function of time in get_vaccine_rates.m
# vr1  = 1e-3; # vaccination rates (per day)
# vr2  = 1/21; # 21 days delay
gamma = 1 / 28
nu_v1 = 2 * 0.25 / 182
nu_v2 = 2 * 0.125 / 365
nu_R = 2 * 0.125 / 365
ai = np.array([[0.5],[0.5],[0.5]])
ai_V = np.array([[0.85],[0.85],[0.85]])
ai_R = np.array([[0.85],[0.85],[0.85]])
mu = 109019 / N / 365
mu_I = 1.75 * np.array([[9255 / 555927 * gamma],[1.6 * 9255 / 555927 * gamma],[1.8 * 9255 / 555927 * gamma]])
mu_IV = 0.15 * mu_I

# parameters for new killer variant, will replace wild-type after fall 2021
new_beta = beta[2]
new_beta_v1 = np.array([[0.5],[0.5]]) * new_beta
new_beta_v2 = np.array([[0.2],[0.2]]) * new_beta
new_beta_R = 0.05 * new_beta
new_ai = 0.8

day0 = date(2020, 1, 1) #'1-jan-2020'
t_new_voc = (date(2022, 9, 1)-day0).days #daysact('1-jan-2020','1-sep-2022')

# define initial population fractions
I0 = np.array([[1e-06],[0],[0]]) # infected
A0 = np.array([[0],[0],[0]]) 
S0 = 1 - sum(I0+A0) # original susceptible
SVR0 = 0  # lost immunity after vaccination or recovery
V10 = np.array([[0],[0]]) # one-dose vaccination
V20 = np.array([[0],[0]]) # fully vaccinated
IV0 = np.array([[0],[0],[0]]) # infected even with vaccination
IR0 = np.array([[0],[0],[0]]) # infected again after recovery from a different variant
AR0 = np.array([[0],[0],[0]]) # asymptomatic infection after recovery from a different variant
R0 = np.array([[0],[0],[0]]) # recovered
R20 = 0 # recovered after getting both variants

y0 = np.vstack((S0,SVR0,V10,V20,I0,IV0,IR0,A0, AR0, R0,R20))
y0 = y0.reshape(len(y0),)


########## SIMULATION
tol = 1e-8
sim = odeint(svair, y0, t, args = (beta,beta_v1,beta_v2,beta_R,ai_beta_ratio,gamma,nu_v1,nu_v2,nu_R,ai,ai_V,ai_R,mu,mu_I,mu_IV,new_beta,new_beta_v1,new_beta_v2,new_beta_R,new_ai,t_new_voc), rtol = tol, atol = tol)
S   = sim[:,0]
SVR = sim[:,1]
V1PF= sim[:,2] 
V1AZ= sim[:,3]
V2PF= sim[:,4]
V2AZ= sim[:,5]
IP  = sim[:,6]  
IA  = sim[:,7]  
ID  = sim[:,8]
IPV = sim[:,9] 
IAV = sim[:,10] 
IDV = sim[:,11]
IPR = sim[:,12] 
IAR = sim[:,13] 
IDR = sim[:,14]
AP  = sim[:,15] 
AA  = sim[:,16] 
AD  = sim[:,17]
APR = sim[:,18] 
AAR = sim[:,19] 
ADR = sim[:,20]
RP  = sim[:,21] 
RA  = sim[:,22] 
RD  = sim[:,23]
R2  = sim[:,24]

######### PLOT RESULTS
fig, axs = plt.subplots(3, 3)
axs[0, 0].plot(t, S, t, SVR)
axs[0, 0].set_title('Susceptible Population')
axs[0, 0].legend('Susceptible','SVR')
# axs[0, 0].xticks([90, 181, 273, 365, 455, 546, 638], ...
#   ['31 Mar 2020','30 Jun 2020','30 Sep 2020','31 Dec 2020','31 Mar 2021','30 Jun 2021','30 Sep 2021'])
axs[0, 0].set_xlabel('Time');
axs[0, 0].set_ylabel('Fraction of Population')

axs[0, 1].plot(t, V1PF, t, V1AZ, t, V2PF, t, V2AZ)
axs[0, 1].set_title('Vaccinated Population')
axs[0, 1].legend(['V1PF','V1AZ','V2PF','V2AZ'])
axs[0, 1].set_xlabel('Time');
axs[0, 1].set_ylabel('Fraction of Population')

axs[0, 2].plot(t, IP, t, IA, t, ID)
axs[0, 2].set_title('I Population')
axs[0, 2].legend(['$I^P$','$I^A$','$I^D$'])
axs[0, 2].set_xlabel('Time');
axs[0, 2].set_ylabel('Fraction of Population')

axs[1, 0].plot(t, IPV, t, IDV, t, IAV)
axs[1, 0].set_title('IV Population')
axs[1, 0].legend(['$I^P_V$','$I^A_V$','$I^D_V$'])
axs[1, 0].set_xlabel('Time')
axs[1, 0].set_ylabel('Fraction of Population')

axs[1, 1].plot(t, IPR, t, IDR, t, IAR)
axs[1, 1].set_title('IR Population')
axs[1, 1].legend(['$I^P_R$','$I^A_R$','$I^D_R$'])
axs[1, 1].set_xlabel('Time')
axs[1, 1].set_ylabel('Fraction of Population')

axs[1, 2].plot(t, AP, t, AD, t, AA)
axs[1, 2].set_title('Asymptomatic Population')
axs[1, 2].legend(['$A^P$','$A^A$','$A^D$'])
axs[1, 2].set_xlabel('Time')
axs[1, 2].set_ylabel('Fraction of Population')

axs[2, 0].plot(t, APR, t, ADR, t, AAR)
axs[2, 0].set_title('AR Population')
axs[2, 0].legend(['$A^P_R$','$A^A_R$','$A^D_R$'])
axs[2, 0].set_xlabel('Time')
axs[2, 0].set_ylabel('Fraction of Population')

axs[2, 1].plot(t, RP, t, RD, t, RA, t, R2)
axs[2, 1].set_title('Recovered Population')
axs[2, 1].legend(['$R^P$','$R^A$','$R^D$','R2'])
axs[2, 1].set_xlabel('Time')
axs[2, 1].set_ylabel('Fraction of Population')
# plt.show()
plt.savefig("voc2driverplots.png")