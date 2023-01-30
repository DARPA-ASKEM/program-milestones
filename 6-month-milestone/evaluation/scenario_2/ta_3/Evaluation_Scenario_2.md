# Evaluation Scenario 2

## Question 1

### Ingest SIDARTHE

```@example scenario2
cd(@__DIR__)
using OrdinaryDiffEq, ModelingToolkit, EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Test
import Plots

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

@unpack Infected, Healed, Extinct, Diagnosed, Ailing, Recognized, Susceptible, Threatened = sys
@unpack alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta, theta, tau = sys
ps = [alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta]

@parameters t
D = Differential(t)
eqs2 = deepcopy(eqs)
append!(eqs2, D.(ps) .~ 0)

sys2 = ODESystem(eqs2, ModelingToolkit.get_iv(sys), states(sys), parameters(sys);
                 continuous_events = evs, defaults = defs, name = nameof(sys))
ssys = structural_simplify(sys2)
```

### Unit Tests

#### Unit Test 1

> Set the initial values and parameters, as described in the Supplementary Methods section of the publication (pg. 9 of the pdf): 

> Initial Values: I = 200/60e6, D = 20/60e6, A = 1/60e6, R = 2/60e6, T = 0, H = 0, E = 0; S = 1 – I – D – A – R – T – H – E. Let total population = 60e6.

> Parameters: α = 0.570, β = δ = 0.011, γ = 0.456, ε = 0.171, θ = 0.371, ζ = η = 0.125, μ = 0.017, ν = 0.027, τ = 0.01, λ = ρ = 0.034 and κ = ξ = σ = 0.017.

> Simulate for 100 days, and determine the day and level of peak total infections (sum over all the infected states I, D, A, R, T). Expected output: The peak should occur around day 47, when ~60% of the population is infected.

```@example

```

#### Unit Test 2

> Now update the parameters to reflect various interventions that Italy implemented during the first wave, as described in detail on pg. 9.  Simulate for 100 days, reproduce the trajectories in Fig. 2B, and determine the day and level of peak total infections (sum over all the infected states I, D, A, R, T). Expected output: Trajectories in Fig. 2B, peak occurs around day 50, with ~0.2% of the total population infected.

```@example
ITALY_POPULATION = 60e6
u0s = [
    Infected => 200 / ITALY_POPULATION,
    Diagnosed => 20 / ITALY_POPULATION,
    Ailing => 1 / ITALY_POPULATION,
    Recognized => 2 / ITALY_POPULATION,
    Threatened => 0,
    Healed => 0,
    Extinct => 0,
]
push!(u0s, Susceptible => 1 - sum(last.(u0s)))

# The resulting basic reproduction number is R0 = 2.38.
pars = [alpha => 0.570, beta => 0.011, delta => 0.011, gamma => 0.456, epsilon => 0.171,
    theta => 0.371,
    zeta => 0.125, eta => 0.125, mu => 0.017, nu => 0.027, tau => 0.01,
    lambda => 0.034, rho => 0.034, kappa => 0.017, xi => 0.017, sigma => 0.017]
prob = ODEProblem(ssys, [], (0, 100))
prob_test1 = remake(prob, u0 = u0s, p = pars)
solt1 = solve(prob_test1, Tsit5(); saveat = 0:100)
og_states = states(sys)[1:8]
idart = [Infected, Diagnosed, Ailing, Recognized, Threatened]
plot(solt1; idxs = Infected)
plot(solt1; idxs = Diagnosed)
plot(solt1; idxs = idart)
@test solt1[Infected + Healed] == solt1[Infected] + solt1[Healed]
plot(solt1.t, solt1[sum(idart)] * ITALY_POPULATION; label = "IDART absolute")
plot(solt1.t, solt1[sum(idart)]; label = "IDART percent")

xmax, xmaxval = get_max_t(prob_test1, sum(idart))

@test isapprox(xmax, 47; atol = 4)
@test isapprox(xmaxval, 0.002;atol=0.01)

```
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
sysv = eval(quote
                var"##iv#608" = (@variables(t))[1]
                var"##sts#609" = (collect)(@variables(Infected(t), Healed(t), Extinct(t),
                                                      Diagnosed(t), Ailing(t),
                                                      Recognized(t), Susceptible(t),
                                                      Threatened(t), Vaccinated(t),
                                                      alpha(t), epsilon(t), gamma(t),
                                                      beta(t), delta(t), mu(t), nu(t),
                                                      lambda(t), rho(t), kappa(t), xi(t),
                                                      sigma(t), zeta(t), eta(t)))
                var"##ps#610" = (collect)(@parameters(ModelValue_21, epsilon_modifier, tau,
                                                      theta, ModelValue_19, ModelValue_20,
                                                      Event_trigger_Fig4b,
                                                      Event_trigger_Fig3d, ModelValue_16,
                                                      Event_trigger_Fig3b, alpha_modifier,
                                                      Event_trigger_Fig4d, ModelValue_17,
                                                      ModelValue_18, Italy, tau1, phi))
                var"##eqs#611" = [(~)((Differential(t))(Infected),
                                      (+)((*)((+)((/)((*)(Infected, alpha), Italy),
                                                  (/)((*)(Diagnosed, beta), Italy),
                                                  (/)((*)(Ailing, gamma), Italy),
                                                  (/)((*)(Recognized, delta), Italy)),
                                              Susceptible), (*)(-1 // 1, Infected, epsilon),
                                          (*)(-1 // 1, Infected, lambda),
                                          (*)(-1 // 1, Infected, zeta)))
                                  (~)((Differential(t))(Healed),
                                      (+)((*)(Ailing, kappa), (*)(Diagnosed, rho),
                                          (*)(Infected, lambda), (*)(Recognized, xi),
                                          (*)(Threatened, sigma)))
                                  (~)((Differential(t))(Extinct),
                                      (*)(tau, Threatened) + tau1 * Recognized)
                                  (~)((Differential(t))(Diagnosed),
                                      (+)((*)(Infected, epsilon),
                                          (*)(-1 // 1, Diagnosed, eta),
                                          (*)(-1 // 1, Diagnosed, rho)))
                                  (~)((Differential(t))(Ailing),
                                      (+)((*)(Infected, zeta), (*)(-1 // 1, theta, Ailing),
                                          (*)(-1 // 1, Ailing, kappa),
                                          (*)(-1 // 1, Ailing, mu)))
                                  (~)((Differential(t))(Recognized),
                                      (+)((*)(theta, Ailing), (*)(Diagnosed, eta),
                                          (*)(-1 // 1, Recognized, nu),
                                          (*)(-1 // 1, Recognized, xi)))
                                  (~)((Differential(t))(Susceptible),
                                      (*)((+)((/)((*)(-1, Infected, alpha), Italy),
                                              (/)((*)(-1, Diagnosed, beta), Italy),
                                              (/)((*)(-1, Ailing, gamma), Italy),
                                              (/)((*)(-1, Recognized, delta), Italy)) -
                                          phi * Susceptible,
                                          Susceptible))
                                  (~)((Differential(t))(Threatened),
                                      (+)((*)(Ailing, mu), (*)(Recognized, nu),
                                          (*)(-1 // 1, tau, Threatened),
                                          (*)(-1 // 1, Threatened, sigma)))
                                  Differential(t)(Vaccinated) ~ phi * Susceptible
                                  (~)((Differential(t))(alpha), -0.0);
                                  (~)((Differential(t))(epsilon), -0.0);
                                  (~)((Differential(t))(gamma), -0.0);
                                  (~)((Differential(t))(beta), -0.0);
                                  (~)((Differential(t))(delta), -0.0);
                                  (~)((Differential(t))(mu), -0.0);
                                  (~)((Differential(t))(nu), -0.0);
                                  (~)((Differential(t))(lambda), -0.0);
                                  (~)((Differential(t))(rho), -0.0);
                                  (~)((Differential(t))(kappa), -0.0);
                                  (~)((Differential(t))(xi), -0.0);
                                  (~)((Differential(t))(sigma), -0.0);
                                  (~)((Differential(t))(zeta), -0.0);
                                  (~)((Differential(t))(eta), -0.0)]
                var"##defs#612" = (Dict)((Pair)(delta, 0.011), (Pair)(xi, 0.017),
                                         (Pair)(Diagnosed, 3.33333333e-7),
                                         (Pair)(Event_trigger_Fig3b, 0.0),
                                         (Pair)(Extinct, 0.0), (Pair)(kappa, 0.017),
                                         (Pair)(zeta, 0.125), (Pair)(eta, 0.125),
                                         (Pair)(nu, 0.027), (Pair)(Healed, 0.0),
                                         (Pair)(Infected, 3.33333333e-6),
                                         (Pair)(ModelValue_16, 0.0),
                                         (Pair)(alpha_modifier, 1.0), (Pair)(Italy, 1.0),
                                         (Pair)(Event_trigger_Fig3d, 0.0),
                                         (Pair)(ModelValue_20, 1.0), (Pair)(sigma, 0.017),
                                         (Pair)(Threatened, 0.0), (Pair)(lambda, 0.034),
                                         (Pair)(alpha, 0.57),
                                         (Pair)(Event_trigger_Fig4b, 0.0),
                                         (Pair)(ModelValue_17, 0.0),
                                         (Pair)(Event_trigger_Fig4d, 0.0),
                                         (Pair)(Susceptible, 0.9999963),
                                         (Pair)(beta, 0.011),
                                         (Pair)(Recognized, 3.33333333e-8),
                                         (Pair)(rho, 0.034), (Pair)(mu, 0.017),
                                         (Pair)(epsilon, 0.171),
                                         (Pair)(Ailing, 1.66666666e-8),
                                         (Pair)(gamma, 0.456), (Pair)(ModelValue_19, 0.0),
                                         (Pair)(ModelValue_21, 1.0), (Pair)(theta, 0.371),
                                         (Pair)(epsilon_modifier, 1.0), (Pair)(tau, 0.01),
                                         (Pair)(ModelValue_18, 0.0),
                                         Vaccinated => 0,
                                         tau1 => 0.0200,
                                         phi => 0.0025)
                var"##iv#613" = (@variables(t))[1]
                (ODESystem)(var"##eqs#611", var"##iv#613", var"##sts#609", var"##ps#610";
                            defaults = var"##defs#612", name = Symbol("##SBML#530"),
                            checks = false)
            end)
# todo set the event flags
# todo validate the new params 
probv = ODEProblem(sysv, [], (0, 100))
solv = solve(probv, Tsit5())
plot(solv)
plot(solv, idxs = [og_states; Vaccinated])
plot(solt1; idxs = sum(idart))

xmax, xmaxval = get_max_t(probv, sum(idart) * ITALY_POPULATION)
xmax, xmaxval = get_max_t(probv, sum(idart))

@test isapprox(xmax, 47; atol = 5)
@test isapprox(xmaxval, 0.6; atol = 0.1)
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