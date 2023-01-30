# Evaluation Scenario 2

## Question 1

### Ingest SIDARTHE

```@example scenario2
cd(@__DIR__)
using OrdinaryDiffEq, ModelingToolkit, EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Plots

function sub_cont_ev(ev, rules)
    ModelingToolkit.SymbolicContinuousCallback(substitute(ev.eqs, rules),
                                               substitute(ev.affect, rules))
end
fn = "Giordano2020.xml"

myread(fn) = readSBML(fn, doc -> begin
                          set_level_and_version(3, 2)(doc)
                          convert_promotelocals_expandfuns(doc)
                      end)

m = myread(fn)
rn = ReactionSystem(m)
sys = convert(ODESystem, rn)
eqs = equations(sys)
defs_ = ModelingToolkit.defaults(sys)
defs = deepcopy(defs_)
evs = ModelingToolkit.get_continuous_events(sys)

# these are the constant=false params 
params_to_sub = unique(ModelingToolkit.lhss(vcat(map(x -> x.affect, evs)...)))

@unpack alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta, theta = sys
ps = [alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta]

@unpack Infected, Healed, Extinct, Diagnosed, Ailing, Recognized, Susceptible, Threatened = sys2
subd = Dict(params_to_sub .=> ps)
newevs = map(x -> sub_cont_ev(x, subd), evs)

sys2 = ODESystem(eqs, ModelingToolkit.get_iv(sys), states(sys), parameters(sys);
                 continuous_events = newevs, defaults = defs, name = nameof(sys))
sys2 = structural_simplify(sys2)
```

```@example scenario2
prob = ODEProblem(sys2, [], (0.0, 1000.0))
sol = solve(prob, Tsit5())
plot(sol, idxs = Infected)
```

### Unit Tests

#### Unit Test 1

> Set the initial values and parameters, as described in the Supplementary Methods section of the publication (pg. 9 of the pdf): 

> Initial Values: I = 200/60e6, D = 20/60e6, A = 1/60e6, R = 2/60e6, T = 0, H = 0, E = 0; S = 1 – I – D – A – R – T – H – E. Let total population = 60e6.

> Parameters: α = 0.570, β = δ = 0.011, γ = 0.456, ε = 0.171, θ = 0.371, ζ = η = 0.125, μ = 0.017, ν = 0.027, τ = 0.01, λ = ρ = 0.034 and κ = ξ = σ = 0.017.

> Simulate for 100 days, and determine the day and level of peak total infections (sum over all the infected states I, D, A, R, T). Expected output: The peak should occur around day 47, when ~60% of the population is infected.

#### Unit Test 2

> Now update the parameters to reflect various interventions that Italy implemented during the first wave, as described in detail on pg. 9.  Simulate for 100 days, reproduce the trajectories in Fig. 2B, and determine the day and level of peak total infections (sum over all the infected states I, D, A, R, T). Expected output: Trajectories in Fig. 2B, peak occurs around day 50, with ~0.2% of the total population infected.

### Sensitivity Analysis

> The difference between 1.b.i and 1.b.ii are changes in some parameter values over time. Describe the difference in outcomes between b.i and b.ii. Perform a sensitivity analysis to understand the sensitivity of the model to parameter variations and determine which parameter(s) were most responsible for the change in outcomes. 

```@example scenario2
pbounds = [param => [0.5*ModelingToolkit.defaults(sys)[param],2*ModelingToolkit.defaults(sys)[param]] for param in parameters(sys2)]
sensres = get_sensitivity(prob, 100.0, Infected, pbounds; samples = 200)
```

```@example scenario2
create_sensitivity_plot(prob, 100.0, Infected, pbounds; samples = 200)
```

### Mininmum Parameter Threshold

> Now return to the situation in b.i (constant parameters that don’t change over time). Let’s say we want to increase testing, diagnostics, and contact tracing efforts (implemented by increasing the detection parameters ε and θ). Assume that θ >= 2* ε, because a symptomatic person is more likely to be tested. What minimum constant values do these parameters need to be over the course of a 100-day simulation, to ensure that the total infected population (sum over all the infected states I, D, A, R, T) never rises above 1/3 of the total population?

```@example scenario2
threshold_observable = (Infected + Diagnosed + Ailing + Recognized + Threatened) / sum(states(sys))
cost = -(eta + theta)
EasyModelAnalysis.optimal_parameter_intervention_for_threshold(prob, threshold_observable, 0.33, 
                                             cost, [eta,theta], [0.0,0.0], 
                                             3 .* [ModelingToolkit.defaults(sys)[eta], ModelingToolkit.defaults(sys)[theta]]; 
                                             maxtime=60)
```

## Question 2

### Ingest SIDARTHE-V

```@example scenario2

```

### Setup the Parameters

> Set the same initial values and parameter settings in 1.b.i. Let V(t=0) = 0, τ (in SIDARTHE) = τ2 (in SIDDARTHE-V), and τ1 = (1/3)\*τ2 (reflecting the fact that the mortality rate for critical conditions (state T), will always be larger than for other infected states). Assume that the vaccination rate psi is 0 to start with. The SIDARTHE-V model allows for three main types of interventions: (1) Those that impact the transmission parameters (α, β, γ and δ) – social distancing, masking, lockdown; (2) Those that impact the detection parameters (ε, θ) – testing and contact tracing; (3) Those that impact the vaccination rate psi – vaccination campaigns. Assume previously stated constraints: θ >= 2* ε, and τ1 = (1/3)*τ2.

```@example scenario2

```

### b.i

> Let’s say our goal is to ensure that the total infected population (sum over all the infected states I, D, A, R, T) never rises above 1/3 of the total population, over the course of the next 100 days. If you could choose only a single intervention (affecting only one parameter), which intervention would let us meet our goal, with minimal change to the intervention parameter? Assume that the intervention will be implemented after one month (t = day 30), and will stay constant after that, over the remaining time period (i.e. the following 70 days). What are equivalent interventions of the other two intervention types, that would have the same impact on total infections?

```@example scenario2
intervention_parameters = [theta] # Need to figure out what these should be
[p => EasyModelAnalysis.optimal_parameter_intervention_for_threshold(prob, threshold_observable, 0.33, 
                                             p - ModelingToolkit.defaults(sys)[p], [p], [0.0], 
                                             3 .* [ModelingToolkit.defaults(sys)[p]],
                                             (30.0,100.0); 
                                             maxtime=60) for p in intervention_parameters]
```

### b.ii

> Let’s say our goal is to get the reproduction number R0 below 1.0, at some point within the next 100 days. Are there interventions that will allow us to meet our goal? If there are multiple options, which single intervention would have the greatest impact on R0 and let us meet our goal with minimal change to the intervention parameter? Assume that the intervention will be implemented after one month (t = day 30), and will stay constant after that, over the remaining time period (i.e. the following 70 days).

```@example scenario2
R0 = Infected # how is R0 defined from the states?
```

```@example scenario2
intervention_parameters = [theta] # Need to figure out what these should be
[p => EasyModelAnalysis.optimal_parameter_intervention_for_threshold(prob, R0, 1.0, 
                                             p - ModelingToolkit.defaults(sys)[p], [p], [0.0], 
                                             3 .* [ModelingToolkit.defaults(sys)[p]],
                                             (30.0,100.0); 
                                             maxtime=60) for p in intervention_parameters]
```