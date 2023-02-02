# Evaluation Scenario 3

```@example evalscenario3
using EasyModelAnalysis, LinearAlgebra, CSV
using Catlab, AlgebraicPetri
using Catlab.CategoricalAlgebra
```

## Question 1

## Setup Model

>  1. Begin with a basic SIR model without vital dynamics. Calibrate the model parameters using data on cases during the ‘training period’. Then evaluate the model during the out-of-sample ‘test period’.

To get started with the code, we implemented the basic SIR without vital dynamics directly in ModelingToolkit.jl. This is a version
that was written by an epidemiologist at Microsoft Pandemic, Simon Frost, who has become a fan of the TA3 automated simulation tools
and wrote an entire repository of tutorials for this software. It is found at https://github.com/epirecipes/sir-julia.

In there is an SIR without vital dynamics which we took in full.

```@example evalscenario3
sir = read_json_acset(LabelledPetriNet, "sir.json")
sys = ODESystem(sir)
sys = complete(sys)
@unpack S, I, R, inf, rec = sys
@parameters N = 1
param_sub = [
    inf => inf / N,
]
sys = substitute(sys, param_sub)
defs = ModelingToolkit.defaults(sys)
defs[S] = 990
defs[I] = 10
defs[R] = 0.0
defs[N] = sum(x -> defs[x], (S, I, R))
defs[inf] = 0.5
defs[rec] = 0.25
tspan = (0.0, 40.0)
prob = ODEProblem(sys, [], tspan);
sol = solve(prob);
```

```@example evalscenario3
plot(sol)
```

### Perform Model Calibration

#### Model Calibration Unit Test

As a unit test of the model calibration tools, we generated data at the default parameters, then ran the global optimization,
to see how well the parameters were recovered.

```@example evalscenario3
dataset = solve(prob, saveat = 0.1)
t_train = dataset.t[1:201]
t_test = dataset.t[202:end]
data_train = [S => dataset[S][1:201], I => dataset[I][1:201], R => dataset[R][1:201]]
data_test = [S => dataset[S][202:end], I => dataset[I][202:end], R => dataset[R][202:end]]
```

```@example evalscenario3
fitparams = global_datafit(prob, [inf => [0.2, 2.0], rec => [0.05, 0.5]],
                           t_train, data_train)
```

This then gives the forecasts in the test data:

```@example evalscenario3
_prob = remake(prob, p = fitparams)
sol = solve(_prob, saveat = t_test);
plot(sol, idxs = S)
plot!(t_test, data_test[1][2])
```

```@example evalscenario3
plot(sol, idxs = I)
plot!(t_test, data_test[2][2])
```

```@example evalscenario3
plot(sol, idxs = R)
plot!(t_test, data_test[3][2])
```

This looks very good and matches the original data, confirming that the inverse problem functionality is functional.

Now we train on data from June 1 2021 to September 30 2021.

#### Application to Real Data from TA1

```@example evalscenario3
using CSV, DataFrames, Downloads

# Infectious/Recovered day by day:
url = "https://raw.githubusercontent.com/DARPA-ASKEM/program-milestones/data-h-d-breakdown/6-month-milestone/evaluation/scenario_3/ta_4/usa-IRDVHN_age_HD_breakdown.csv"
file = CSV.File(Downloads.download(url))
df_raw = DataFrame(file)

start_train = 171
stop_train = 171 + 121
start_test = 171 + 122
stop_test = 171 + 122 + 92

df_train = df_raw[start_train:stop_train, :]
df_test = df_raw[start_test:stop_test, :]

t_train = collect(0:(size(df_train, 1) - 1))
t_test = collect(0:(size(df_test, 1) - 1))

N_total = 334998398 # assumed to be constant from (https://github.com/DARPA-ASKEM/program-milestones/blob/main/6-month-milestone/evaluation/scenario_3/ta_1/usa-2021-population-age-stratified.csv)
#S = N_total - R - I
data_train = [S => N_total .- df_train.I .- df_train.R, I => df_train.I, R => df_train.R]
data_test = [S => N_total .- df_test.I .- df_test.R, I => df_test.I, R => df_test.R]

u0s = [S => N_total - df_train.I[1] - df_train.R[1], I => df_train.I[1], R => df_train.R[1]]
_prob = remake(prob, u0 = u0s, tspan = (t_train[1], t_train[end]), p = [N => N_total])

fitparams = global_datafit(_prob, [inf => [0, 1.0], rec => [0.0, 1.0]], t_train, data_train,
                           maxiters = 1_000_000)
```

```@example evalscenario3
# Plot training fit
_prob_train = remake(_prob, p = fitparams)
sol = solve(_prob_train, saveat = t_train);
```

```@example evalscenario3
plot(map(data_train) do (var, num)
         plot(sol, idxs = var)
         plot!(t_train, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("train_fit_S3_Q1.png")
```

#### Why is that the best fit?

At first glance it may look like the system was incorrect, i.e. that it did not find the global optima for this problem.
However, upon further inspection we can show that this is truly the global optima. To see this, we have to inspect
against the "intuitive" solution. The intuitive solution would be to simply place the peak of the infections at the
right spot.

```@example evalscenario3
_prob_train = remake(_prob, p = [inf => 0.363, rec => 0.29])
sol = solve(_prob_train, saveat = t_train);
plot(sol, idxs = I)
plot!(t_train, data_train[2][2], lab = "I_train")
```

However, if we check the L2 loss of this fit we will see it's a lot higher.

```@example evalscenario3
pkeys = [inf, rec]
EasyModelAnalysis.l2loss([0.363, 0.29], (_prob, pkeys, t_train, data_train))
```

```@example evalscenario3
EasyModelAnalysis.l2loss([fitparams[1][2], fitparams[2][2]],
                         (_prob, pkeys, t_train, data_train))
```

The reason for this is because the fits of "making the infected maximum have the correct peak" forces the susceptible
population to be far off. This then introduces a much larger total error.

```@example evalscenario3
_prob_train = remake(_prob, p = [inf => 0.363, rec => 0.29])
sol = solve(_prob_train, saveat = t_train);
plot(sol, idxs = S)
plot!(t_train, data_train[1][2], lab = "S_train")
```

```@example evalscenario3
_prob_train = remake(_prob, p = [inf => 0.363, rec => 0.29])
sol = solve(_prob_train, saveat = t_train);
plot(sol, idxs = R)
plot!(t_train, data_train[2][2], lab = "R_train")
```

The reason that it is off is because the onset of this pandemic has a delay, i.e. it is flat for a bit before taking off.
That cannot be the case for the SIR model. If `inf > rec`, then the onset of the pandemic is at time zero.

This motivates fitting in terms of not the L2 norm but the relative L2 norm, i.e. with values weighted in terms of the
relative size of S. This would make the much larger values of S not dominate the overall loss due to the relative
difference in units.

```@example evalscenario3
fitparams = global_datafit(_prob, [inf => [0, 1.0], rec => [0.0, 1.0]], t_train, data_train,
                           maxiters = 1_000_000, loss = EasyModelAnalysis.relative_l2loss)
```

which makes no substantial difference to the result. This is because the "intuitive solution" is also worse in terms
of relative loss:

```@example evalscenario3
EasyModelAnalysis.relative_l2loss([0.363, 0.29], (_prob, pkeys, t_train, data_train))
```

```@example evalscenario3
EasyModelAnalysis.relative_l2loss([fitparams[1][2], fitparams[2][2]],
                                  (_prob, pkeys, t_train, data_train))
```

In other words, while one may wish to fit the infected spike, doing so would cause the susceptible and recovered values
to be so far off that it leads to more error than a bad fit of the infected. The SIR model is simply not a good fit
to this data.

Another way to see this result is to notice that both the number of susceptible individuals and recovered individuals
are both dropping exponentially at a growing rate at the end of the time after the peak of the infection, which is
incompatible with the SIR model's assumptions that the rate of S -> I and I -> R would both drop after the infection's
peak.

### SIR Forecasting Plots

Demonstrated are the forecasts with the best fitting SIR parameters

```@example evalscenario3
# Plot testing fit
u0s = [S => N_total - df_test.I[1] - df_test.R[1], I => df_test.I[1], R => df_test.R[1]]
_prob_test = remake(_prob, p = fitparams, u0 = u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob_test, saveat = t_test);

plot(map(data_test) do (var, num)
         plot(sol, idxs = var)
         plot!(t_test, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("test_fit_S3_Q1.png")
```

## Question 2: Add Hospitalizations and Deaths

This expands the original SIR model to explore a model space comprising SIRD, SIRH, and SIRHD.

```@example evalscenario3
sird = read_json_acset(LabelledPetriNet, "sird.json")
sirh = read_json_acset(LabelledPetriNet, "sirh.json")
sirhd = read_json_acset(LabelledPetriNet, "sirhd.json")
sirhd_sys = ODESystem(sirhd)
sirhd_sys = complete(sirhd_sys)
@unpack S, I, R, H, D, inf, rec, ideath, death, hosp, hrec = sirhd_sys
@parameters N = 1
param_sub = [
    inf => inf / N,
]
sirhd_sys = substitute(sirhd_sys, param_sub)
defs = ModelingToolkit.defaults(sirhd_sys)
defs[S] = N_total - 10
defs[I] = 10
defs[H] = 0
defs[D] = 0
defs[R] = 0.0
defs[N] = N_total
defs[inf] = 0.5
defs[rec] = 0.25
defs[ideath] = 0.25
defs[death] = 0.25
defs[hosp] = 0.25
defs[hrec] = 0.25
tspan = (0.0, 40.0)
sirhd_prob = ODEProblem(sirhd_sys, [], tspan)
sirhd_sol = solve(sirhd_prob)
plot(sirhd_sol)
```

Question 2 involves doing the same analysis as question one but on the SIR model with hospitalizations and deaths included.

#### Model Calibration

```@example evalscenario3
data_train = [
    S => N_total .- df_train.I .- df_train.R .- df_train.D .- df_train.H,
    I => df_train.I, R => df_train.R, H => df_train.H, D => df_train.D,
]
data_test = [
    S => N_total .- df_test.I .- df_test.R .- df_test.D .- df_test.H,
    I => df_test.I, R => df_test.R, H => df_test.H, D => df_test.D,
]

u0s = [
    S => N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1],
    I => df_train.I[1], R => df_train.R[1], H => df_train.H[1], D => df_train.D[1],
]
_prob2 = remake(sirhd_prob, u0 = u0s, tspan = (t_train[1], t_train[end]),
                p = [N => N_total])

param_bounds = [inf => [0.0, 1.0]
                rec => [0.0, 1.0]
                death => [0.0, 5.0]
                ideath => [0.0, 5.0]
                hosp => [0.0, 5.0]
                hrec => [0.0, 5.0]]
fitparams2 = global_datafit(_prob2, param_bounds, t_train, data_train,
                            maxiters = 500_000, loss = EasyModelAnalysis.relative_l2loss)
```

```@example evalscenario3
# Plot training fit
_prob2_train = remake(_prob2, p = fitparams2)
sol = solve(_prob2_train, saveat = t_train);
plot(map(data_train) do (var, num)
         plot(sol, idxs = var)
         plot!(t_train, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("train_fit_S3_Q2.png")
```

```@example evalscenario3
u0s = [
    S => N_total - df_test.I[1] - df_test.R[1] - df_test.H[1] - df_test.D[1],
    I => df_test.I[1], R => df_test.R[1], H => df_test.H[1], D => df_test.D[1],
]
_prob2_test = remake(_prob2, p = fitparams2, u0 = u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob2_test, saveat = t_test);
plot(map(data_test) do (var, num)
         plot(sol, idxs = var)
         plot!(t_test, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("test_fit_S3_Q2.png")
```

#### Explanation of the Fit

Once again it's clear that the model is unable to fit the data well. The same issues apply to the new data as well.
The S, R, and D data all increase the rate of growth after the infection's peak, which is impossible to occur in
the SIRHD model and thus suggests that the infection peak might be an anomoly of the data. In either case, it's
impossible to both fit the non-decreasing derivatives of these 3 data series while also fitting the peak of the
infection, with this model. Additionally. for no parameters does the SIRHD model emit oscillatory solutions as seen
in the data, which suggests a model deficiency.

### Evaluate Model Forecasts

In order to evaluate the model forecasts, we developed a functional which does the forecasting part with multiple models
and puts a score on the forecast result. This score is calculated using the L2 norm.

```@example evalscenario3
norm(solve(_prob, saveat = t_test)[S] - data_test[1][2]) +
norm(solve(_prob, saveat = t_test)[I] - data_test[2][2]) +
norm(solve(_prob, saveat = t_test)[R] - data_test[3][2])
```

```@example evalscenario3
norm(solve(_prob2, saveat = t_test)[S] - data_test[1][2]) +
norm(solve(_prob2, saveat = t_test)[I] - data_test[2][2]) +
norm(solve(_prob2, saveat = t_test)[R] - data_test[3][2]) +
norm(solve(_prob2, saveat = t_test)[H] - data_test[4][2]) +
norm(solve(_prob2, saveat = t_test)[D] - data_test[5][2])
```

Overall the SIRHD model gives a better forecast. Though for no values of the model can an SIR or SIRHD model predict
a second wave, all of these models can only have a singular peak in the infections, and thus we would say that neither
of the models are adequate predictors of this data.

## Question 3: Add Vaccinations

This expands the previous SIRHD model to add vaccination.

```@example evalscenario3
using ModelingToolkit.Symbolics: variable
using ModelingToolkit: toparam
sirhd_vax = read_json_acset(LabelledPetriNet, "sirhd_vax.json")
sirhd_vax_sys = structural_simplify(ODESystem(sirhd_vax))
sirhd_vax_sys = complete(sirhd_vax_sys)
names = string.(ModelingToolkit.getname.(states(sirhd_vax_sys)))
sts_names = Symbol.(getindex.(names, 6), :_, getindex.(names, 11))
@variables t
@parameters N
sts = map(n -> variable(n, T = SymbolicUtils.FnType)(t), sts_names)
names = split.(string.(parameters(sirhd_vax_sys)), "\"")
ps_names = Symbol.(getindex.(split.(getindex.(names, 3), "\\"), 1), :_,
                   getindex.(split.(getindex.(names, 5), "\\"), 1))
ps = map(n -> toparam(variable(n)), ps_names)
nps = findall(x -> occursin("inf", string(x)), ps)
ups = findall(x -> !occursin("inf", string(x)), ps)
subs = [parameters(sirhd_vax_sys)[nps] .=> ps[nps] ./ N
        parameters(sirhd_vax_sys)[ups] .=> ps[ups]
        states(sirhd_vax_sys) .=> sts]
sirhd_vax_sys = substitute(sirhd_vax_sys, subs)
@unpack S_U, I_U, R_U, H_U, D_U, S_V, I_V, R_V, H_V, D_V = sirhd_vax_sys
S, I, R, H, D = S_U, I_U, R_U, H_U, D_U
Sv, Iv, Rv, Hv, Dv = S_V, I_V, R_V, H_V, D_V
@unpack id_vax, inf_infuu, inf_infuv, hosp_id, ideath_id, rec_id, hrec_id, death_id, inf_infvu, inf_infvv = sirhd_vax_sys
sys3 = sirhd_vax_sys
prob3 = ODEProblem(sys3, nothing, (0, 1.0))
```

Question 3 is the same analysis as questions 1 and 2 done on a model with vaccination added. In order to build unit tests for
the analysis and functionality, we started by building the model with vaccine by hand, awaiting a swap to the version from
TA2.

#### Model Calibration

```@example evalscenario3
data_train = [(S + Sv) => N_total .- df_train.I .- df_train.R .- df_train.H .- df_train.D,
    (I + Iv) => df_train.I, (R + Rv) => df_train.R, H => df_train.H_unvac,
    Hv => df_train.H_vac, D => df_train.D_unvac, Dv => df_train.D_vac]
data_test = [(S + Sv) => N_total .- df_test.I .- df_test.R .- df_test.H .- df_test.D,
    Sv => 0,
    (I + Iv) => df_test.I, (R + Rv) => df_test.R, H => df_test.H_unvac, Hv => df_test.H_vac,
    D => df_test.D_unvac, Dv => df_test.D_vac]

vac_rate = df_train.H_vac[1] / (df_train.H_vac[1] + df_train.H_unvac[1])
# 52% of hospitalizations are vaccinated, we do not have data for vaccination rates for other compartments,
# so we assume that the vaccination rate is the same for all compartments.
```

```@example evalscenario3
u0s = [
    S => (1 - vac_rate) *
         (N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1]),
    I => (1 - vac_rate) * df_train.I[1],
    R => (1 - vac_rate) * df_train.R[1],
    H => df_train.H_unvac[1],
    D => df_train.D_unvac[1],
    Sv => vac_rate *
          (N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1]),
    Iv => vac_rate * df_train.I[1],
    Rv => vac_rate * df_train.R[1],
    Hv => df_train.H_vac[1],
    Dv => df_train.D_vac[1],
]
_prob3 = remake(prob3, u0 = u0s, tspan = (t_train[1], t_train[end]),
                p = [ps .=> 0.1; N => N_total])

fitparams3 = global_datafit(_prob3, ps .=> ([0, 5.0],), t_train, data_train,
                            maxiters = 200_000, loss = EasyModelAnalysis.relative_l2loss)
```

```@example evalscenario3
u0s = [
    S => (1 - vac_rate) *
         (N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1]),
    I => (1 - vac_rate) * df_train.I[1],
    R => (1 - vac_rate) * df_train.R[1],
    H => df_train.H_unvac[1],
    D => df_train.D_unvac[1],
    Sv => vac_rate *
          (N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1]),
    Iv => vac_rate * df_train.I[1],
    Rv => vac_rate * df_train.R[1],
    Hv => df_train.H_vac[1],
    Dv => df_train.D_vac[1],
]
_prob3 = remake(_prob3, p = fitparams3, u0 = u0s, tspan = (t_train[1], t_train[end]))
sol = solve(_prob3, saveat = t_train);
plot(map(data_train) do (var, num)
         plot(sol, idxs = [var])
         plot!(t_train, num)
     end...)
```

```@example evalscenario3
u0s = [
    S => (1 - vac_rate) *
         (N_total - df_test.I[1] - df_test.R[1] - df_test.H[1] - df_test.D[1]),
    I => (1 - vac_rate) * df_test.I[1],
    R => (1 - vac_rate) * df_test.R[1],
    H => df_test.H_unvac[1],
    D => df_test.D_unvac[1],
    Sv => vac_rate * (N_total - df_test.I[1] - df_test.R[1] - df_test.H[1] - df_test.D[1]),
    Iv => vac_rate * df_test.I[1],
    Rv => vac_rate * df_test.R[1],
    Hv => df_test.H_vac[1],
    Dv => df_test.D_vac[1],
]
_prob3 = remake(_prob3, p = fitparams3, u0 = u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob3, saveat = t_test);
plot(map(data_test) do (var, num)
         plot(sol, idxs = [var])
         plot!(t_test, num)
     end...)
```

## Question 4: Age-Stratified Model

## Question 5: Add Reinfection

Question 5 is the same analysis as questions 1, 2, 3, and 4 on a model with
reinfection added. We will use SIRHD with reinfection.

```@example evalscenario3
sirhd = read_json_acset(LabelledPetriNet, "sirhd_renew.json")
sirhd_sys = ODESystem(sirhd)
sirhd_sys = complete(sirhd_sys)
@unpack S, I, R, H, D, inf, rec, ideath, death, hosp, hrec, renew = sirhd_sys
@parameters N = 1
param_sub = [
    inf => inf / N,
]
sirhd_sys = substitute(sirhd_sys, param_sub)
defs = ModelingToolkit.defaults(sirhd_sys)
defs[S] = 100 - 10
defs[I] = 10
defs[H] = 0
defs[D] = 0
defs[R] = 0.0
defs[N] = 100
defs[inf] = 0.5
defs[rec] = 0.25
defs[ideath] = 0.25
defs[death] = 0.25
defs[hosp] = 0.25
defs[hrec] = 0.25
defs[renew] = 0.25
tspan = (0.0, 40.0)
sirhd_prob = ODEProblem(sirhd_sys, [], tspan)
sirhd_sol = solve(sirhd_prob)
plot(sirhd_sol)
```

#### Model Calibration

```@example evalscenario3
start_train = 171
stop_train = 171 + 213
start_test = 171 + 214
stop_test = 171 + 214 + 151

df_train = df_raw[start_train:stop_train, :]
df_test = df_raw[start_test:stop_test, :]

t_train = collect(0:(size(df_train, 1) - 1))
t_test = collect(0:(size(df_test, 1) - 1))

data_train = [
    S => N_total .- df_train.I .- df_train.R .- df_train.D .- df_train.H,
    I => df_train.I, R => df_train.R, H => df_train.H, D => df_train.D,
]
data_test = [
    S => N_total .- df_test.I .- df_test.R .- df_test.D .- df_test.H,
    I => df_test.I, R => df_test.R, H => df_test.H, D => df_test.D,
]

u0s = [
    S => N_total - df_train.I[1] - df_train.R[1] - df_train.H[1] - df_train.D[1],
    I => df_train.I[1], R => df_train.R[1], H => df_train.H[1], D => df_train.D[1],
]
_prob5 = remake(sirhd_prob, u0 = u0s, tspan = (t_train[1], t_train[end]),
                p = [N => N_total])

param_bounds = [inf => [0.0, 1.0]
                rec => [0.0, 1.0]
                death => [0.0, 5.0]
                ideath => [0.0, 5.0]
                hosp => [0.0, 5.0]
                hrec => [0.0, 5.0]
                renew => [0.0, 1.0]]
fitparams5 = global_datafit(_prob5, param_bounds, t_train, data_train,
                            maxiters = 500_000, loss = EasyModelAnalysis.relative_l2loss)
```

```@example evalscenario3
# Plot training fit
_prob5_train = remake(_prob5, p = fitparams5)
sol = solve(_prob5_train, saveat = t_train);
plot(map(data_train) do (var, num)
         plot(sol, idxs = var)
         plot!(t_train, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("train_fit_S3_Q5.png")
```

```@example evalscenario3
u0s = [
    S => N_total - df_test.I[1] - df_test.R[1] - df_test.H[1] - df_test.D[1],
    I => df_test.I[1], R => df_test.R[1], H => df_test.H[1], D => df_test.D[1],
]
_prob5_test = remake(_prob5, p = fitparams5, u0 = u0s, tspan = (t_test[1], t_test[end]))
sol = solve(_prob5_test, saveat = t_test);
plot(map(data_test) do (var, num)
         plot(sol, idxs = var)
         plot!(t_test, num)
     end..., dpi = 300)
```

```@example evalscenario3
savefig("test_fit_S3_Q5.png")
```

## Question 6: New Data

Question 6 is currently awaiting data from TA3

## Question 7: Analysis

For each model, summarize your conclusions about the following:

 1. Do parameters fit from data seem reasonable and fall within typical ranges you might see in the broader literature? Provide references to support your conclusions.
 2. Describe how well the fitted model compares against historical data, both for the ‘training’ and ‘test’ periods.

### Answer

#### Question 2

The fits against the historical data, both training and test periods, are unremarkable but for clear reasons as described
in the text. The SIR/SIRHD/etc. models all have clear limitations which can easily be proven from their analytical form
this includes:

 1. Only a single infection peak or monotonic results in I.
 2. Rate decreases required after the infection peak for S, R, and D.
 3. H also can only have simple peaking results, or monotonic results.

These 3 points are directly incompatible with the real-world data. This is because the real-world data has the properties
that:

 1. The rate of decrease in S, R, and D do not decrease after the peak of the infection.
 2. The H data is non-monotonic.

Since these facts are directly incompatible with the SIR and SIRHD models, there is no set of parameters that is able to
fit all of the qualitative features of the model. Instead, the model fit tends to fit S and R, or S, R, and D with disreguard
to the other features as a way to minimize the L2 loss, which leads to a non-qualitative fit of the infection peak and the
hospitalization non-monotonicity.

These facts are compounded in the forecasting, which attempts to forecast a new wave of the infection. The quadratic models
are unable to model equations with dual infection peaks, and thus this is impossible to predict from any of the models except
the final which includes reinfection.

#### Question 1

The fits from the data seem reasonable according to literature on the analytical results of the SIR and SIRHD type models.
https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-022-07403-5 describes how interventions are required in the
model in order to be able to model second wave events, which makes it clear that a second wave cannot be modeled in these
models and thus the forecast predictions being incorrect is a known defficiency of the model.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7997702/ provides a full derivation for the infection peak and shows that the
calculation has a unique zero derivative point, which directly proves the assertion that there cannot be a second wave
in the SIR and SIRHD models.
