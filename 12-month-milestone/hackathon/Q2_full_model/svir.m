function yp = svir(t,y,beta,beta_v1,beta_v2,beta_R,gamma,nu_v1,nu_v2,nu_R,mu,mu_I,mu_IV)

% retrieve current populations
n_var = 3;  % number of viruses simulated
n_vax = 2; % number of vaccines simulated
ind  = 1;
S    = y(ind); ind=ind+1;  % original susceptible
SVR  = y(ind); ind=ind+1;  % lost immunity after vaccination or recovery
V1   = y(ind:ind+n_vax-1); ind=ind+n_vax;  % one-dose vaccination
V2   = y(ind:ind+n_vax-1); ind=ind+n_vax;  % fully vaccinated
I    = y(ind:ind+n_var-1); ind=ind+n_var;  % infected
IV   = y(ind:ind+n_var-1); ind=ind+n_var;  % infected even with vaccination
IR   = y(ind:ind+n_var-1); ind=ind+n_var;  % infected again after recovery from a different variant
R    = y(ind:ind+n_var-1); ind=ind+n_var;  % recovered
R2   = y(ind); ind=ind+1;  % recovered after getting both variants

% get time-dependent parameters
[vr1, vr2] = get_vaccine_rate (t);
beta_scale = get_beta (t);
beta       = beta*beta_scale;
beta_v1    = beta_v1*beta_scale;
beta_v2    = beta_v2*beta_scale;
beta_R     = beta_R*beta_scale;

% compute time derivatives
dSdt   = - sum(beta.*S.*(I+IV+IR)) - sum(vr1.*S) + mu*(1-S);
dSVRdt = + nu_v1.*sum(V1) + nu_v2.*sum(V2) + sum(nu_R.*R) + nu_R.*R2 - sum(beta.*SVR.*(I+IV+IR)) - mu*SVR;
dV1dt  = + vr1.*S - vr2.*V1 - nu_v1.*V1 - sum(beta_v1.*(V1*(I+IV+IR)'),2) - mu*V1;
dV2dt  = + vr2.*V1 - nu_v2.*V2 - sum(beta_v2.*(V2.*(I+IV+IR)'),2) - mu*V2;
dIdt   = + beta.*S.*(I+IV+IR) + beta.*SVR.*(I+IV+IR) - gamma.*I - mu_I.*I;
dIVdt  = + sum(beta_v1.*(V1*(I+IV+IR)'))' + sum(beta_v2.*(V2.*(I+IV+IR)'))' - gamma.*IV - mu_IV.*IV;
dIRdt  = + beta_R.*flip(R).*(I+IV+IR) - gamma.*IR - mu*IR;
dRdt   = + gamma.*(I+IV+IR) - nu_R.*R - beta_R.*flip(R).*(I+IV+IR) - mu*R;
dR2dt  = + sum(gamma.*IR) - nu_R.*R2 - mu*R2;

% simulate mutation of WT into variants and importation
if (t>300 && t<360)  % alpha appears
    dWTtoA  = 0e-5*dIdt(1);  % don't assume mutation
    dIdt(1) = dIdt(1)-dWTtoA;
    dIdt(2) = dIdt(2)+dWTtoA + ((t-250)/60)*10/14570000;  % people enter Ontario with alpha

elseif (t>385 && t<445)
    dIdt(3) = ((t-385)/60)*40/14570000;  % delta was born
end


yp = [dSdt;dSVRdt;dV1dt;dV2dt;dIdt;dIVdt;dIRdt;dRdt;dR2dt];

