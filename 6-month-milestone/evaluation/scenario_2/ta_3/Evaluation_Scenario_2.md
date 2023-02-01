# Evaluation Scenario 2

## Question 1

### Ingest SIDARTHE

This is the version directly from the SBML file. This will be replaced with a version from TA1/TA2 when available. This used our
[SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) library which reads SBML into ModelingToolkit and generates TeX'd versions
of the equations so we could read the resulting model and confirm it is correct against the paper description.

```@example scenario2
using EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Test

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

The SBML model that was given had a few oddities. First, it made use of `delay` blocks. These are usually used to describe a
[delay differential equation](https://en.wikipedia.org/wiki/Delay_differential_equation). While our
[simulator does have the capability to solve delay differential equations](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dde_example/)
and [their inverse problems](https://docs.sciml.ai/SciMLSensitivity/dev/examples/dde/delay_diffeq/), it turns out that this was an issue
with the SBML writing of the model as all of the delay values were zero. Thus as a simplification, we manually deleted the delay blocks
to give a standard ODE representation (since a delay of 0 on all states is mathematically equivalent).

The paper and SBML model also described time-dependent parameters. These are parameters that would discretely change at pre-specified
time points. However, we believe that the evaluators and/or paper must have used an SBML reading system that incorrectly handled these
time-dependent parameters. This is because dropping the time-dependency and treating the parameters as constant gives the requested
results of the unit test 1. A demonstration of this is as follows:

```@example scenario2
sysne = ODESystem(eqs2, ModelingToolkit.get_iv(sys), states(sys), parameters(sys);
                  defaults = defs, name = nameof(sys))
ssysne = structural_simplify(sysne)
probne = ODEProblem(ssysne, [], (0.0, 100.0))
solne = solve(probne, Tsit5())
plot(solne)
```

```@example scenario2
ITALY_POPULATION = 60e6
idart = [Infected, Diagnosed, Ailing, Recognized, Threatened]
xmax, xmaxval = get_max_t(probne, sum(idart))
@test isapprox(xmax, 47; atol = 0.5)
```

```@example scenario2
@test isapprox(xmaxval, 0.6, atol = 0.01)
```

#### Full Analysis of the Effect of Events

The SBML model already contains the changes of parameters requested for Unit Test #2 in the form of events. In the SBML these are enclosed by the element `listOfEvents`. The `id`s correspond to the days of the introduction of government intervention as outlined in the paper. "On day 4, R0 = 1.66 as a result of the introduction of basic social distancing, awareness of the epidemic, hygiene and behavioral recommendations, and early measures by the Italian government (for example, closing schools). At day 12, ... ". When we interpret the instructions "Unit Test #1: Set the initial values and parameters, as desribed in the Supplementary Methods..." as "Set the initial values and set and maintain the paramters (i.e. simulate without government restrictions)..." we have to remove the events from the model. If we do so, the unit tests pass.

```@example scenario2
solne = solve(probne, Tsit5())
p = plot(solne, vars = idart)
```

```@example scenario2
p = plot(solne.t, solne[sum(idart)])
```

#### Unit Test 2

> Now update the parameters to reflect various interventions that Italy implemented during the first wave, as described in detail on pg. 9.  Simulate for 100 days, reproduce the trajectories in Fig. 2B, and determine the day and level of peak total infections (sum over all the infected states I, D, A, R, T). Expected output: Trajectories in Fig. 2B, peak occurs around day 50, with ~0.2% of the total population infected.

This unit test was a straightforward implementation of the scenario, requiring a named reparameterization. It makes use of the new
ModelingToolkit feature designed for ASKEM, `remake(prob, u0 = u0s, p = pars)`, which allows for a new ODE to be generated from the
old ODE simply by mentioning which parameters need to be changed (all others are kept constant). The approximation tests on the
bottom demonstrate that the results in Fig 2B are obtained.

```@example scenario2
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
```

```@example scenario2
plot(solt1; idxs = Diagnosed)
```

```@example scenario2
plot(solt1; idxs = idart)
```

```@example scenario2
@test solt1[Infected + Healed] == solt1[Infected] + solt1[Healed]
```

```@example scenario2
plot(solt1.t, solt1[sum(idart)] * ITALY_POPULATION; label = "IDART absolute")
```

```@example scenario2
plot(solt1.t, solt1[sum(idart)]; label = "IDART percent")
```

```@example scenario2
xmax, xmaxval = get_max_t(prob_test1, sum(idart))

@test isapprox(xmax, 47; atol = 4)
```

```@example scenario2
@test isapprox(xmaxval, 0.002; atol = 0.01)
```

### Sensitivity Analysis

> The difference between 1.b.i and 1.b.ii are changes in some parameter values
> over time. Describe the difference in outcomes between b.i and b.ii. Perform a
> sensitivity analysis to understand the sensitivity of the model to parameter
> variations and determine which parameter(s) were most responsible for the
> change in outcomes.

This analysis was a straightforward application of the `get_sensitivity` function in EasyModelAnalysis. The only issue was the creation
of the bounds for the parameters, which was not given by the metadata from TA1/TA2. Thus we made a modeling choice that the viable
parameter set is 50% below and 100% above the starting parameter choice. Future iterations of the modeling platform should preserve
parameter bound data which would make this a one line analysis.

A utility was added (https://github.com/SciML/EasyModelAnalysis.jl/pull/134) to make it so the sensitivity values did not need to
be recreated for the plotting process. This was just a minor performance and "niceity" improvement. Polish.

The sensitivity analysis needed 1000 samples, we reduced it to 200 due to memory limitations of our documentation building
compute server.

```@example scenario2
pbounds = [param => [
               0.5 * ModelingToolkit.defaults(sys)[param],
               2 * ModelingToolkit.defaults(sys)[param],
           ] for param in parameters(sys2)]
sensres = get_sensitivity(probne, 100.0, Infected, pbounds; samples = 200)
sensres_vec = collect(sensres)
sort(filter(x -> endswith(string(x[1]), "_first_order"), sensres_vec), by = x -> abs(x[2]),
     rev = true)
```

```@example scenario2
sort(filter(x -> endswith(string(x[1]), "_second_order"), sensres_vec), by = x -> abs(x[2]),
     rev = true)
```

```@example scenario2
sort(filter(x -> endswith(string(x[1]), "_total_order"), sensres_vec), by = x -> abs(x[2]),
     rev = true)
```

```@example scenario2
create_sensitivity_plot(sensres, pbounds)
```

### Minimum Parameter Threshold

> Now return to the situation in b.i (constant parameters that don’t change over
> time). Let’s say we want to increase testing, diagnostics, and contact tracing
> efforts (implemented by increasing the detection parameters ε and θ). Assume
> that θ >= 2* ε, because a symptomatic person is more likely to be tested. What
> minimum constant values do these parameters need to be over the course of a
> 100-day simulation, to ensure that the total infected population (sum over all
> the infected states I, D, A, R, T) never rises above 1/3 of the total
> population?

This scenario demonstrates the
[lazily defined observables](https://sciml.github.io/EasyModelAnalysis.jl/dev/getting_started/#Lazily-Defining-Observables)
functionality that persists throughout our simulation and analysis libraries. When one solves an equation with ModelingToolkit
symbolic values, `sol[x]` gives the solution with respect to `x` by name. While that improves code legibility, `sol[x+y]` is
also allowed, and will automatically generate the solution of `x(t) + y(t)` on demand. Since this functionality is directly
handled by the solution representation, this means that all functions built on the solution have this functionality. Thus
without having to make any other changes, we can change our minimization to the complex form
`(Infected + Diagnosed + Ailing + Recognized + Threatened) / sum(states(sys))` required by the scenario.

However, this scenario also required making a modeling choice. In order to perform this minimization we needed, we needed
to define the comparative cost between the different intervention parameters, `epsilon` and `theta`. We have made the assumption
that the cost of interventions on these two parameters are the same, and have made requests to TA1/TA2 about the interpretation
of these parameters for further information.

```@example scenario2
threshold_observable = (Infected + Diagnosed + Ailing + Recognized + Threatened) /
                       sum(states(sys))
cost = -(epsilon + theta)
ineq_cons = [2 * epsilon - theta]
opt_p, sol_opt_p, ret = optimal_parameter_threshold(probne, threshold_observable,
                                                    0.33,
                                                    cost, [epsilon, theta],
                                                    [0.0, 0.0],
                                                    3 .* [
                                                        ModelingToolkit.defaults(sys)[epsilon],
                                                        ModelingToolkit.defaults(sys)[theta],
                                                    ];
                                                    maxtime = 60,
                                                    ineq_cons);
opt_p
```

```@example scenario2
plot(sol_opt_p, idxs = [threshold_observable], lab = "total infected", leg = :topright)
```

## Question 2

### Ingest SIDARTHE-V

This is a handwritten verison of the SIDARTHE-V model, built from the exported SIDARTHE SBML and then manually handcorrected to be
the SIDARTHE-V model. This should swap to the TA1/TA2 model form when available.

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
sysv = complete(sysv)
```

```@example scenario2
probv = ODEProblem(sysv, [], (0, 100))
solv = solve(probv, Tsit5())
plot(solv)
```

```@example scenario2
plot(solv, idxs = [og_states; Vaccinated])
```

```@example scenario2
plot(solt1; idxs = sum(idart))
```

```@example scenario2
xmax, xmaxval = get_max_t(probv, sum(idart) * ITALY_POPULATION)
xmax, xmaxval = get_max_t(probv, sum(idart))

@test isapprox(xmax, 47; atol = 5)
```

```@example scenario2
@test isapprox(xmaxval, 0.6; atol = 0.1)
```

### Setup the Parameters

> Set the same initial values and parameter settings in 1.b.i. Let V(t=0) = 0, τ
> (in SIDARTHE) = τ2 (in SIDDARTHE-V), and τ1 = (1/3)\*τ2 (reflecting the fact
> that the mortality rate for critical conditions (state T), will always be
> larger than for other infected states). Assume that the vaccination rate psi
> is 0 to start with. The SIDARTHE-V model allows for three main types of
> interventions: (1) Those that impact the transmission parameters (α, β, γ and
> δ) – social distancing, masking, lockdown; (2) Those that impact the detection
> parameters (ε, θ) – testing and contact tracing; (3) Those that impact the
> vaccination rate psi – vaccination campaigns. Assume previously stated
> constraints: θ >= 2* ε, and τ1 = (1/3)*τ2.

```@example scenario2
# TODO: double check. I am assuming our `tau` is `tau1` and `tau1` is `tau2`.
defs_v2 = deepcopy(ModelingToolkit.defaults(sysv))
defs_v2[sysv.tau1] = defs[tau]
defs_v2[sysv.tau] = defs[tau] / 3
defs_v2[sysv.phi] = 0
probv2 = remake(probv; p = defs_v2)
solv2 = solve(probv2)
plot(solv2)
```

### b.i

> Let’s say our goal is to ensure that the total infected population (sum over
> all the infected states I, D, A, R, T) never rises above 1/3 of the total
> population, over the course of the next 100 days. If you could choose only a
> single intervention (affecting only one parameter), which intervention would
> let us meet our goal, with minimal change to the intervention parameter?
> Assume that the intervention will be implemented after one month (t = day 30),
> and will stay constant after that, over the remaining time period (i.e. the
> following 70 days). What are equivalent interventions of the other two
> intervention types, that would have the same impact on total infections?

This is a straightforward usage of the `EasyModelAnalysis.optimal_parameter_intervention_for_threshold` function designed during
the ASKEM hackathon. It was able to be used without modification. However, a modeling decision had to be made to define
what the "intervention parameters" are. A data request back to TA1/TA2 has been made to define which parameters should be in this
set.

This example revealed a typo in our function (https://github.com/SciML/EasyModelAnalysis.jl/pull/135) which had to be fixed.

```@example scenario2
threshold_observable = (Infected + Diagnosed + Ailing + Recognized + Threatened) /
                       sum(states(sysv))
plot(solv2, idxs = [threshold_observable], lab = "total infected")
hline!([1 / 3], lab = "limit")
```

```@example scenario2
intervention_p = phi # Need to figure out what these should be
cost = intervention_p - defs_v2[intervention_p]
opt_p, solv2_s, ret = optimal_parameter_intervention_for_threshold(probv2,
                                                                   threshold_observable,
                                                                   0.33,
                                                                   cost,
                                                                   [intervention_p], [0.0],
                                                                   [1.0],
                                                                   (30.0, 100.0);
                                                                   maxtime = 10);
opt_p
```

Note that the optimization solution is trivial, i.e. there's no intervention at
all. This is expected because the model without any intervention would already
have less than 1/3 of the population infected.

### b.ii

> Let’s say our goal is to get the reproduction number R0 below 1.0, at some
> point within the next 100 days. Are there interventions that will allow us to
> meet our goal? If there are multiple options, which single intervention would
> have the greatest impact on R0 and let us meet our goal with minimal change to
> the intervention parameter? Assume that the intervention will be implemented
> after one month (t = day 30), and will stay constant after that, over the
> remaining time period (i.e. the following 70 days).

In order to do this scenario a modeling decision for how to represent R0 in
terms of the states was required. This needed expert information, which we
called out for and documented the results in https:
//github.com/ChrisRackauckas/ASKEM_Evaluation_Staging/issues/20. This led us to
a definition of the instantanious R0 as defined in https:
//www.ncbi.nlm.nih.gov/pmc/articles/PMC7325187/ equation 1. However, the R
computation requires the mean duration of infectiousness, we will use `D=20` for
now to have a non-trivial optimization. Thus using this definition of R0 and our
intervention functionality designed to find parameters to keep a value below a
threshold, we were able to solve for the intervention.

Another modeling decision required here was the definition of intervention parameters, which we decided to use the same parameters
as b.i.

```@example scenario2
D = 20
R0 = sysv.alpha * sysv.Susceptible * D # double check
plot(solv2, idxs = [R0])
```

```@example scenario2
intervention_parameters = [sysv.theta => (2 * defs_v2[sysv.eta], 1) # 𝜃 >= 2 * 𝜀
                           sysv.eta => (0, defs_v2[sysv.theta] / 2)
                           sysv.phi => (0, 1)]
opt_results = map(intervention_parameters) do (intervention_p, bounds)
    cost = intervention_p - defs_v2[intervention_p]
    optimal_parameter_intervention_for_reach(probv2,
                                             R0,
                                             1.0,
                                             cost,
                                             [intervention_p], [bounds[1]], [bounds[2]],
                                             (30.0, 100.0);
                                             maxtime = 10)
end;
map(first, opt_results)
```

```@example scenario2
plts = map(opt_results) do opt_result
    title = only(collect(opt_result[1]))
    title = title[1] => round(title[2], sigdigits = 3)
    plot(opt_result[2][2]; idxs = [R0], lab = "R0", title)
    hline!([1], lab = "limit")
end
plot(plts...)
```
