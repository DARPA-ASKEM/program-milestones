# Evalution Scenario 3

```@example evalscenario3
using EasyModelAnalysis, LinearAlgebra
```

## Question 1

## Setup Model

>1.	Begin with a basic SIR model without vital dynamics. Calibrate the model parameters using data on cases during the ‘training period’. Then evaluate the model during the out-of-sample ‘test period’. 

#### Sample Model, Swap in for TA1/TA2 model

```@example evalscenario3
@parameters t β=0.05 c=10.0 γ=0.25
@variables S(t)=990.0 I(t)=10.0 R(t)=0.0
∂ = Differential(t)
N=S+I+R # This is recognized as a derived variable
eqs = [∂(S) ~ -β*c*I/N*S,
       ∂(I) ~ β*c*I/N*S-γ*I,
       ∂(R) ~ γ*I];

@named sys = ODESystem(eqs);
```

```@example evalscenario3
tspan = (0.0,40.0)
prob = ODEProblem(sys,[],tspan);
sol = solve(prob);
```

### Sample Model End

```@example evalscenario3
plot(sol)
```

### Perform Model Calibration 

#### Make a fake dataset, swap in for TA1/TA2 code

```@example evalscenario3
dataset = solve(prob, saveat = 0.1)
t_train = dataset.t[1:201]
t_test = dataset.t[202:end]
data_train = [S => dataset[S][1:201],I => dataset[I][1:201],R => dataset[R][1:201]]
data_test = [S => dataset[S][202:end],I => dataset[I][202:end],R => dataset[R][202:end]]
```

#### Fake data code end

```@example evalscenario3
fitparams = global_datafit(prob, [β=>[0.03,0.15], c=>[9.0,13.0], γ=>[0.05,0.5]], t_train, data_train)
```

```@example evalscenario3
_prob = remake(prob, p=fitparams)
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

## Question 2: Add Hospitalizations and Deaths

```@example evalscenario3
@parameters t β=0.1 c=10.0 γ=0.25 ρ=0.1 h=0.1 d=0.1 r=0.1
@variables S(t)=990.0 I(t)=10.0 R(t)=0.0 H(t)=0.0 D(t)=0.0
∂ = Differential(t)
N=S+I+R+H+D # This is recognized as a derived variable
eqs = [∂(S) ~ -β*c*I/N*S,
       ∂(I) ~ β*c*I/N*S-γ*I-h*I-ρ*I,
       ∂(R) ~ γ*I + r*H,
       ∂(H) ~ h*I - r*H - d*H, 
       ∂(D) ~ ρ*I + d*H];

@named sys2 = ODESystem(eqs);
```

```@example evalscenario3
prob2 = ODEProblem(sys2,[],tspan);
sol = solve(prob2);
```

```@example evalscenario3
plot(sol)
```

```@example evalscenario3
fitparams2 = global_datafit(prob2, [β=>[0.03,0.15], c=>[9.0,13.0], γ=>[0.05,0.5]], t_train, data_train)
```

```@example evalscenario3
_prob2 = remake(prob2, p=fitparams2)
sol = solve(_prob2, saveat = t_test);
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

### Evaluate Model Forecasts

```@example evalscenario3
norm(solve(_prob, saveat = t_test)[S] - data_test[1][2]) + norm(solve(_prob, saveat = t_test)[I] - data_test[2][2]) + norm(solve(_prob, saveat = t_test)[R] - data_test[3][2])
```
```@example evalscenario3
norm(solve(_prob2, saveat = t_test)[S] - data_test[1][2]) + norm(solve(_prob2, saveat = t_test)[I] - data_test[2][2]) + norm(solve(_prob2, saveat = t_test)[R] - data_test[3][2])
```

## Question 3: Add Vaccinations

```@example evalscenario3
@parameters t β=0.1 c=10.0 γ=0.25 ρ=0.1 h=0.1 d=0.1 r=0.1 v=0.1
@parameters t β2=0.1 c2=10.0 ρ2=0.1 h2=0.1 d2=0.1 r2=0.1
@variables S(t)=990.0 I(t)=10.0 R(t)=0.0 H(t)=0.0 D(t)=0.0
@variables Sv(t)=990.0 Iv(t)=10.0 Rv(t)=0.0 Hv(t)=0.0 Dv(t)=0.0
@variables I_total(t)

∂ = Differential(t)
N=S+I+R+H+D+Sv+Iv+Iv+Hv+Dv # This is recognized as a derived variable
eqs = [∂(S) ~ -β*c*I_total/N*S - v*Sv,
       ∂(I) ~ β*c*I_total/N*S-γ*I-h*I-ρ*I,
       ∂(R) ~ γ*I + r*H,
       ∂(H) ~ h*I - r*H - d*H, 
       ∂(D) ~ ρ*I + d*H,
       
       ∂(Sv) ~ -β2*c2*I_total/N*Sv + v*Sv,
       ∂(Iv) ~ β2*c2*I_total/N*Sv-γ*I-h2*I-ρ2*I,
       ∂(Rv) ~ γ*I + r2*H,
       ∂(Hv) ~ h2*I - r2*H - d2*H, 
       ∂(Dv) ~ ρ2*I + d2*H,
       
       I_total ~ I + Iv
       ];

@named sys3 = ODESystem(eqs)
sys3 = structural_simplify(sys3)
```

```@example evalscenario3

```

## Question 4: Age-Stratified Model

#### Same as the other questions, just on a different model

## Question 5: Add Reinfection

#### Same as the other questions, just on a different model

## Question 6: New Data

#### Same as the other questions, just on new data

## Question 7: Analysis

For each model, summarize your conclusions about the following:

1.	Do parameters fit from data seem reasonable and fall within typical ranges you might see in the broader literature? Provide references to support your conclusions.
2.	Describe how well the fitted model compares against historical data, both for the ‘training’ and ‘test’ periods.

### Answer