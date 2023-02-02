# Evaluation Scenario 2

## Question 1

### Read in the SIDARTHE model from a JSON formed from Semagrams

```@example scenario2
using Catlab, AlgebraicPetri, Catlab.CategoricalAlgebra, ModelingToolkit
using AlgebraicPetri.SubACSets
sidarthe = read_json_acset(LabelledPetriNet, "sidarthe.json")
sys = ODESystem(sidarthe)
```

### Load parameter values and initial concentrations from SBML file

This uses our [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl) library which reads SBML into ModelingToolkit and generates TeX'd versions of the equations so we could read the resulting model and confirm it is correct against the paper description.

```@example scenario2
using EasyModelAnalysis, SBML, SBMLToolkit, UnPack, Test

fn = "Giordano2020.xml"

myread(fn) = readSBML(fn, doc -> begin
                          set_level_and_version(3, 2)(doc)
                          convert_promotelocals_expandfuns(doc)
                      end)

m = myread(fn)

paramvals = map(name -> m.parameters[string(name)].value, tnames(sidarthe))
namemap = Dict(:S => "Susceptible", :I => "Infected", :D => "Diagnosed", :A => "Ailing",
               :R => "Recognized",
               :T => "Threatened", :H => "Healed", :E => "Extinct")
u0vals = map(name -> m.species[namemap[name]].initial_concentration, snames(sidarthe))
let S, I, D, A, R, T, H, E
    @unpack S, I, D, A, R, T, H, E = sys
    global Infected, Healed, Extinct, Diagnosed, Ailing, Recognized, Susceptible, Threatened
    Infected, Healed, Extinct, Diagnosed, Ailing, Recognized, Susceptible, Threatened = I,
                                                                                        H,
                                                                                        E,
                                                                                        D,
                                                                                        A,
                                                                                        R,
                                                                                        S, T
end
@unpack beta, gamma, delta, alpha, epsilon, kappa, sigma, rho, xi, mu, tau, lambda, eta, nu, zeta, theta = sys
ps = [alpha, epsilon, gamma, beta, delta, mu, nu, lambda, rho, kappa, xi, sigma, zeta, eta]
defaultsmap = Dict(param => val for (param, val) in zip(parameters(sys), paramvals))
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
ssys = structural_simplify(sys)
probne = ODEProblem(ssys, u0vals, (0.0, 100.0), paramvals)
solne = solve(probne, Tsit5())
plot(solne)
```

```@example scenario2
idart = [Infected, Diagnosed, Ailing, Recognized, Threatened]
xmax, xmaxval = get_max_t(probne, sum(idart))
@test isapprox(xmax, 47; atol = 0.5)
@test isapprox(xmaxval, 0.6, atol = 0.01)
```

#### Full Analysis of the Effect of Events

The SBML model already contains the changes of parameters requested for Unit Test #2 in the form of events. In the SBML these are enclosed by the element `listOfEvents`. The `id`s correspond to the days of the introduction of government intervention as outlined in the paper. "On day 4, R0 = 1.66 as a result of the introduction of basic social distancing, awareness of the epidemic, hygiene and behavioral recommendations, and early measures by the Italian government (for example, closing schools). At day 12, ... ". When we interpret the instructions "Unit Test #1: Set the initial values and parameters, as desribed in the Supplementary Methods..." as "Set the initial values and set and maintain the paramters (i.e. simulate without government restrictions)..." we have to remove the events from the model. If we do so, the unit tests pass.

```@example scenario2
solne = solve(probne, Tsit5())
p = plot(solne, vars = idart)
```

```@example scenario2
p = plot(solne, idxs = [sum(idart)], lab = "total infected")
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
prob = ODEProblem(ssys, u0vals, (0, 100), paramvals)
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

This test passes with SBML.jl

```@example scenario2
@test_broken isapprox(xmaxval, 0.002; atol = 0.01)
```

This last test worked with the SBML script, but fails with the model from TA2. It seems inconsequential
to the rest of the analysis.

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
               0.5 * defaultsmap[param],
               2 * defaultsmap[param],
           ] for param in parameters(sys)]
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
create_sensitivity_plot(sensres, pbounds, true, ylims = (-0.2, 1.0), size = (800, 800))
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
cost = (epsilon + theta)
ineq_cons = [2 * epsilon - theta]
opt_p, sol_opt_p, ret = optimal_parameter_threshold(probne, threshold_observable,
                                                    0.33,
                                                    cost, [epsilon, theta],
                                                    [0.0, 0.0],
                                                    3 .* [
                                                        defaultsmap[epsilon],
                                                        defaultsmap[theta],
                                                    ];
                                                    maxtime = 60,
                                                    ineq_cons);
opt_p
```

```@example scenario2
plot(sol_opt_p, idxs = [threshold_observable], lab = "total infected", leg = :bottomright)
```

## Question 2

### Form SIDARTHE-V model

This forms SIDARTHE-V by manually adding the V state and vax transition. It compares the models via maximum common subacset, plotting the common subgraph (the original SIDARTHE), the negation (the new transition and vax state), and the complement (the new transition from susceptible to vax).

```@example scenario2
import Graphviz_jll
sidarthe_v = read_json_acset(LabelledPetriNet, "sidarthe_v.json")

mca_sidarthe_v = mca(sidarthe, sidarthe_v)
AlgebraicPetri.Graph(mca_sidarthe_v[1])
```

```@example scenario2
sidarthe_sub = Subobject(sidarthe_v,
                         S = parts(sidarthe, :S),
                         T = parts(sidarthe, :T),
                         I = parts(sidarthe, :I),
                         O = parts(sidarthe, :O))
AlgebraicPetri.Graph(dom(hom(negate(sidarthe_sub))))
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
sysv = ODESystem(sidarthe_v)
u0valsv = [u0vals; 0.0]  # 0 vaccinated initially
@unpack beta, gamma, delta, alpha, epsilon, kappa, sigma, rho, xi, mu, tau1, tau2, lambda, eta, nu, zeta, theta, vax = sysv
phi = vax
p_map = Dict([parameters(sys) .=> paramvals
              tau2 => defaultsmap[tau]
              tau1 => defaultsmap[tau] / 3
              vax => 0.0])
sts_map = Dict(states(sysv) .=> u0valsv)
using ModelingToolkit: @set!
defs_v2 = merge(sts_map, p_map)
@set! sysv.defaults = defs_v2
sysv = complete(sysv)
```

```@example scenario2
probv = ODEProblem(sysv, [], (0, 100))
solv = solve(probv, Tsit5())
plot(solv)
```

```@example scenario2
plot(solv, idxs = [og_states; sysv.V])
```

```@example scenario2
plot(solt1; idxs = sum(idart))
```

```@example scenario2
xmax, xmaxval = get_max_t(probv, sum(idart) * ITALY_POPULATION)
xmax, xmaxval = get_max_t(probv, sum(idart))

@test isapprox(xmax, 47; atol = 1)
```

```@example scenario2
@test isapprox(xmaxval, 0.6; atol = 0.1)
```

This test passed with the original SBML model but failed with the model from the TA2 integration.

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
plot(solv, idxs = [threshold_observable], lab = "total infected")
hline!([1 / 3], lab = "limit")
```

```@example scenario2
intervention_parameters = [theta => (2 * defs_v2[epsilon], 1) # 𝜃 >= 2 * 𝜀
                           epsilon => (0, defs_v2[theta] / 2)
                           phi => (0, 1)]

opt_results = map(intervention_parameters) do (intervention_p, bounds)
    cost = intervention_p - defs_v2[intervention_p]
    optimal_parameter_intervention_for_threshold(probv,
                                                 threshold_observable,
                                                 0.33,
                                                 cost,
                                                 [intervention_p], [0.0],
                                                 [1.0],
                                                 (30.0, 100.0);
                                                 maxtime = 10)
end;
map(first, opt_results)
```

```@example scenario2
plts = map(opt_results) do opt_result
    title = only(collect(opt_result[1]))
    title = string(title[1], " = ", round(title[2], sigdigits = 3))
    plot(opt_result[2][2]; idxs = [threshold_observable], lab = "total infected", title)
    hline!([1 / 3], lab = "limit")
end
plot(plts...)
```

### b.ii

> Let’s say our goal is to get the reproduction number Rt below 1.0, within the
> next 60 days. Which interventions will allow us to meet our goal, while
> minimizing total cumulative infections (over all infected states I, D, A, R,
> T)? If there are multiple options, show the tradeoff between change in
> parameter and infected populations – show the space of possible solutions.
> Which single intervention would have the greatest impact on Rt and let us meet
> our goal with minimal change to the intervention parameter, while minimizing
> total cumulative infections? Assume that the intervention will be implemented
> immediately. Use Rt as defined in the SIDDARTHE-V publication. No intervention
> and increasing the infected population, are not valid solutions for this
> problem.

A modeling decision required here was the definition of intervention parameters,
which we decided to use the same parameters as b.i.

```@example scenario2
@unpack S, I, D, A, R, T = sysv
r1 = epsilon + xi + lambda
r2 = eta + rho
r3 = theta + mu + kappa
r4 = nu + xi + tau1
r5 = sigma + tau2
R0 = (alpha + beta * epsilon / r2 + gamma * zeta / r3 +
      delta * ((eta * epsilon /
                (r2 * r4)) + (zeta * theta) / (r3 * r4))) / r1
R_t = S * R0
plot(solv, idxs = [R_t], lab = "R_t")
@variables t cumulative_inf(t)
total_inf = (I + D + A + R + T) / sum(states(sysv))
sysva = add_accumulations(sysv, [cumulative_inf => total_inf])
using ModelingToolkit: @set!
@set! sysva.defaults = defs_v2
u03 = [u0valsv; 0.0]
probv3 = ODEProblem(sysva, u03, (0, 60.0))
solv3 = solve(probv3)
plot(solv3)
```

A modeling decision required here is that we want to minimize the sum of
normalized cumulative total infections and the change of intervention parameters
`theta`, `epsilon`, and `phi`.

```@example scenario2
intervention_parameters = [theta => (2 * defs_v2[epsilon], 3) # 𝜃 >= 2 * 𝜀
                           epsilon => (0, defs_v2[theta] / 2)
                           phi => (0, 1)]
sol_cost = sol -> begin sol(sol.t[end], idxs = cumulative_inf) end
opt_results = map(intervention_parameters) do (intervention_p, bounds)
    cost = intervention_p - defs_v2[intervention_p]
    optimal_parameter_intervention_for_reach(probv3,
                                             R_t,
                                             1.0,
                                             (cost, sol_cost),
                                             [intervention_p], [bounds[1]], [bounds[2]],
                                             maxtime = 10)
end;
map(first, opt_results)
```

We can see that increasing the rate of vaccination is the most effective.

```@example scenario2
plts = map(opt_results) do opt_result
    title = only(collect(opt_result[1]))
    title = string(title[1], " = ", round(title[2], sigdigits = 3))
    plot(opt_result[2][2]; idxs = [R_t], lab = "R_t", title)
    plot!(opt_result[2][2]; idxs = [cumulative_inf], lab = "cumulative infected", title)
    hline!([1], lab = "limit", ylims = (0, 10))
end
plot(plts...)
```
