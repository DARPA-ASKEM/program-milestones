# Evaluation Scenario 1

This document is the solution set for Scenario 1 of the Jan 2023 DARPA
ASKEM program evaluation. Quoted text is the scenario specfification provided
by DARPA/MITRE. Unquoted text, code and output are intended to address
the question posed in the scenario document, while providing use case
illustration for the libraries in question.

In this scenario, we investigate the effects of different age population distributions on the effect of viral epidemics using a simple SIR model.

First, we load various packages for model exploration and data loading,
which we will use later.

```@example scenario1
using EasyModelAnalysis, LinearAlgebra
using EasyModelAnalysis.ModelingToolkit: toparam
using EasyModelAnalysis.ModelingToolkit.Symbolics: FnType, variables
using XLSX, CSV, DataFrames, Plots
using Catlab, Catlab.CategoricalAlgebra, Catlab.Programs, AlgebraicPetri, AlgebraicPetri.TypedPetri
using Base: splat
```

## Stratified SIR

> Scenario Ask: In order to consider more nuanced interventions, we would like for models to account for different age groups and their contact dynamics. Start with a basic SIR model without vital dynamics, and stratify it according to the following questions.

We begin by creating a basic function that manually creates a stratified SIR model given a list of population buckets.

```@example scenario1
tf = 600
const k = 1000
const γ = 1/14
const R₀ = 5
const β = R₀ * γ

"""
    make_stratified_model(pops, contactmat; infectedfrac = nothing)

Given a list of population buckets `pops` (of length `n`), create a stratified SIR model with interaction
incidence given by `contactmat`. The SIR parameters γ and R0 (and thus β) are inherited from global scope.
Returns the system of differential equations that can be used for futher simulations.

The optional keyword argument `infectedfrac` specifies the fraction of each group that should be
considered infected initially.
"""
function make_stratified_model(pops, contactmat; infectedfrac = nothing)
    types = LabelledPetriNet([:Pop],
                             :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
                             :disease=>(:Pop=>:Pop),
                             :strata=>(:Pop=>:Pop))
    sir_uwd = @relation () where (S::Pop, I::Pop, R::Pop) begin
        infect(S,I,I,I)
        disease(I,R)
    end
    sir_typed = oapply_typed(types, sir_uwd, [:inf, :rec])
    totalpop = sum(pops)
    n = length(pops)
    I₀ = something(infectedfrac, n / totalpop)
    sir_paramd = add_params(sir_typed, Dict{Symbol,Float64}(:S => 1 - I₀, :I => I₀, :R => 0),
                            Dict(:inf => β, :rec => γ))

    ages = pairwise_id_typed_petri(types, :Pop, :infect, [Symbol("A$i") for i = 1:n],
                                   pops, contactmat ./ pops,
                                   codom_net = codom(sir_paramd))
    ages = add_reflexives(ages, repeat([[:disease]], n), types)
    return typed_product(sir_paramd, ages)
end

"""
    scenario1(pops, contactmat; infectedfrac, numinfected)

Run a scenario 1 simulation with the given population buckets and contact matrix.
`infectedfrac` has the same meaning as in `make_stratified_model`.
`numinfected` allows specifying a fixed number of individuals in each age group
considered infected initially.
"""
function scenario1(pops, mat; infectedfrac = nothing, numinfected = nothing)
    sir_strat = flatten_labels(make_stratified_model(pops, mat; infectedfrac).dom)
    sys = ODESystem(sir_strat)
    U₀ = map(splat(*), sir_strat[:, :concentration])
    if !isnothing(numinfected)
        U₀[2:2:end] .= numinfected
    end
    P = map(splat(*), sir_strat[:, :rate])
    prob = ODEProblem(sys, U₀, (0, tf), P)
    solve(prob)
end
```

## Question 1

> Start with a simple stratification with three age groups: young, middle-aged, and old.

### Sub-question 1.a.
> Begin with a situation where the population size across each age group is uniform: N_young = 2k, N_middle = 2k, N_old = 2k. Assume only one person in each age group is infectious at the beginning of the simulation. Let gamma = 1/14 days, and let R0 = 5. Assume gamma, beta, and R0 are the same for all age groups.

> i. Simulate this model for the case where the 3x3 contact matrix is uniform (all values in matrix are 0.33)

N.B.: Uniform `1/n_strata` is the default in our model creation function above.

```@example scenario1
sol = scenario1([2k, 2k, 2k], fill(1/3, 3, 3), numinfected = 1)
plt_a1 = plot(sol, leg = :topright)
```

> ii. Simulate this model for the case where there is significant in-group contact
> preference – you may choose the numbers in the matrix to represent this in-
> group preference.

We pick a contact matrix with significant (0.4) in-group interaction
and somewhat weak (but differening, to make the plots more interesting)
off-diagonal interactions.

```@example scenario1
contact_matrix = [0.4  0.05  0.1
                  0.05 0.4   0.15          
                  0.1  0.15  0.4]
```

We now use this updated contact matrix to re-run the simulation.

```@example scenario1
sol = scenario1([2k, 2k, 2k], contact_matrix, numinfected = 1)
plt_a2 = plot(sol, leg = :topright)
```

> iii. Simulate this model for the case where there is no contact between age groups.
> You may choose the numbers in the matrix, but ensure it meets the requirement
> of no contact between age groups.

```@example scenario1
sol = scenario1([2k, 2k, 2k], Diagonal(contact_matrix), numinfected = 1)
plt_a3 = plot(sol, leg = :topright)
```

> Simulate social distancing by scaling down the uniform contact matrix by a
> factor (e.g. multiply by 0.5)

```@example scenario1
uniform_matrix = fill(0.33, (3, 3))
sol = scenario1([2k, 2k, 2k], 0.5 * uniform_matrix, numinfected = 1)
plt_a4 = plot(sol, leg = :topright)
```

> Repeat 1.a.iv for the scenario where the young population has poor compliance
> with social distancing policies, but the old population is very compliant.

```@example scenario1
scaling = Diagonal([0.9, 0.8, 0.4])
sol = scenario1([2k, 2k, 2k], scaling * uniform_matrix, numinfected = 1)
plt_a5 = plot(sol, leg = :topright)
```

Now we combine all the plots into one to allow easy comparison.

```@example scenario1
plot(plt_a1, plt_a2, plt_a3, plt_a4, plt_a5, size = (1000, 500))
```

> Repeat 1.a for a younger-skewing population: `N_young = 3k, N_middle = 2k, N_old = 1k`

```@example scenario1
sol = scenario1([3k, 2k, 1k], fill(1/3, 3, 3), numinfected = 1)
plt_b1 = plot(sol, leg = :topright, title = "i")

sol = scenario1([3k, 2k, 1k], contact_matrix, numinfected = 1)
plt_b2 = plot(sol, leg = :topright, title = "ii")

sol = scenario1([3k, 2k, 1k], Diagonal(contact_matrix), numinfected = 1)
plt_b3 = plot(sol, leg = :topright, title = "iii")

sol = scenario1([3k, 2k, 1k], 0.5 * uniform_matrix, numinfected = 1)
plt_b4 = plot(sol, leg = :topright, title = "iv")

sol = scenario1([3k, 2k, 1k], scaling * uniform_matrix, numinfected = 1)
plt_b5 = plot(sol, leg = :topright, title = "v")

plot(plt_b1, plt_b2, plt_b3, plt_b4, plt_b5, size = (1000, 500))
```

> Repeat 1.a for an older-skewing population: `N_young = 1k, N_middle = 2k, N_old = 3k`

```@example scenario1
sol = scenario1([1k, 2k, 3k], fill(1/3, 3, 3), numinfected = 1)
plt_c1 = plot(sol, leg = :topright, title = "i")

sol = scenario1([1k, 2k, 3k], contact_matrix, numinfected = 1)
plt_c2 = plot(sol, leg = :topright, title = "ii")

sol = scenario1([1k, 2k, 3k], Diagonal(contact_matrix), numinfected = 1)
plt_c3 = plot(sol, leg = :topright, title = "iii")

sol = scenario1([1k, 2k, 3k], 0.5 * uniform_matrix, numinfected = 1)
plt_c4 = plot(sol, leg = :topright, title = "iv")

sol = scenario1([1k, 2k, 3k], scaling * uniform_matrix, numinfected = 1)
plt_c5 = plot(sol, leg = :topright, title = "v")

plot(plt_c1, plt_c2, plt_c3, plt_c4, plt_c5, size = (1000, 500))
```

> d. Compare simulation outputs from 1a-c, and describe any takeaways/conclusions.

The most difference between age demographics can be seen in the case
where there are social-distancing compliance differences. We complare
these plots here:

```@example scenario1
plot(plt_a5, plt_b5, plt_c5)
```

## Question 2

> Now find real contact matrix data and stratify the basic SIR model with the appropriate number of age groups to match the data found. To simulate the model with realistic initial values, find data on population distribution by age group. As in question 1, let gamma = 1/14 days, and let R0 = 5. Assume gamma, beta, and R0 are the same for all age groups.

TA1 provided the data from ["Projecting social contact matrices in 152 countries using contact surveys and demographic data"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697) by Prem, et al.
This paper comes with [10 Excel files](https://doi.org/10.1371/journal.pcbi.1005697.s002) that provide contact matrices
for 152 countries at work, home, school and other locations (plus
a data set of the sum of these locations). The data values in these
excel files are raw, averaged survey results.

!!! note
    We make no attempt to normalize or otherwise adjust the contact matrices. The interpretation of the contact matrices needs to be
    consistent with the SIR parameter `β`, which is fixed in our the
    given scenario. In a real world case, care would need to be taken
    to match the units of the contact matrix to the units of `β`.

```@example scenario1
# Load contact matrices
xf_all_locations1 = XLSX.readxlsx("data/MUestimates_all_locations_1.xlsx")
xf_all_locations2 = XLSX.readxlsx("data/MUestimates_all_locations_2.xlsx")
xf_work1 = XLSX.readxlsx("data/MUestimates_work_1.xlsx")
xf_work2 = XLSX.readxlsx("data/MUestimates_work_2.xlsx")
xf_school1 = XLSX.readxlsx("data/MUestimates_school_1.xlsx")
xf_school2 = XLSX.readxlsx("data/MUestimates_school_2.xlsx")
xf_home1 = XLSX.readxlsx("data/MUestimates_home_1.xlsx")
xf_home2 = XLSX.readxlsx("data/MUestimates_home_2.xlsx")
xf_other1 = XLSX.readxlsx("data/MUestimates_other_locations_1.xlsx")
xf_other2 = XLSX.readxlsx("data/MUestimates_other_locations_2.xlsx")

xfs1 = (; all = xf_all_locations1, work = xf_work1, school = xf_school1,
        home = xf_home1, other = xf_other1)
xfs2 = (; all = xf_all_locations2, work = xf_work2, school = xf_school2,
        home = xf_home2, other = xf_other2)

to_cm(sheet) = Float64[sheet[i, j] for i in 2:17, j in 1:16]
```

We begin by loading up the relevant data for Belgium
and quickly visualizing the contract matrix.

```@example scenario1
# Load Belgium contact matrix
cm_belg = to_cm(xf_all_locations1["Belgium"])
heatmap(cm_belg, yflip=true)
```

Next we load the population distribution data.

```@example scenario1
pop_belg = collect(values(CSV.read("data/2022_ Belgium_population_by_age.csv", DataFrame,
                           header = 3)[1, 2:17]))
bar(1:length(pop_belg), collect(pop_belg), permute=(:x, :y), xlabel="Age (5 year buckets)", ylabel="Total # of people", leg=:none)
```

```@example scenario1
# Set up model
# Per MITRE: Assume that the same fixed fraction of the population in each stratum is initially infected. Here: 0.01%
sol = scenario1(pop_belg, cm_belg, infectedfrac = 0.0001)
plot(sol, leg = :topright)
```

> If the data you’ve found supports this, compare the situation for a country with significant multi-generational contact beyond two generations (as indicated by multiple contact matrix diagonal bandings), and for a country without.

TA1 advises that India has significant multi-generational contact, while Belgium does not. We repeat the exercise for India.

```@example scenario1
# Load India contact matrix
cm_india = to_cm(xf_all_locations1["India"])
hm = heatmap(cm_india, yflip=true)

# Load India population distribution
pop_india = collect(values(CSV.read("data/2016_india_population_by_age.csv", DataFrame)[1, 3:18]))
bar_india = bar(1:length(pop_india), collect(pop_india), permute=(:x, :y), xlabel="Age (5 year buckets)", ylabel="Total # of people", leg=:none)
plot(hm, bar_india)
```

```@example scenario1
sol = scenario1(pop_india, cm_india, infectedfrac = 0.0001)
plot(sol, leg = :topright)
```

> If the data supports this, try implementing interventions like: (1) School closures (2) Social distancing at work and other locations, but not at home.

> (1) School closures

> Prem et al Supplementary info, page 20

```@example scenario1
function cm_school(xfs, country)
    to_cm(xfs[:home][country]) + to_cm(xfs[:work][country]) + to_cm(xfs[:other][country])
end # no school

cm_belgium_school_closure = cm_school(xfs1, "Belgium")
sol = scenario1(pop_belg, cm_belgium_school_closure, infectedfrac = 0.0001)
plot(sol, leg = :topright)
```

> (2) Social distancing at work and other locations, but not at home.

> Prem et al sets social distancing to reduce contacts by half

```@example scenario1
function cm_social_dist(xfs, country)
    to_cm(xfs[:home][country]) + 0.5 * to_cm(xfs[:work][country]) +
    0.5 * to_cm(xfs[:school][country]) + 0.5 * to_cm(xfs[:other][country])
end

cm_belgium_social_dist = cm_social_dist(xfs1, "Belgium")
sol = scenario1(pop_belg, cm_belgium_social_dist, infectedfrac = 0.0001)
plot(sol, leg = :topright)
```
