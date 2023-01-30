# Evaluation Scenario 1

## Stratified SIR

```@example scenario1
using EasyModelAnalysis, LinearAlgebra
using EasyModelAnalysis.ModelingToolkit: toparam
using EasyModelAnalysis.ModelingToolkit.Symbolics: FnType, variables
tf = 600
@parameters γ=1 / 14 R0=5
β = R0 * γ
k = 1000

function make_statified_model(pops; pop_assumption = (stratum, pop)->1)
    @variables t S(t) I(t) R(t)
    D = Differential(t)
    n_stratify = length(pops)
    Ns = map(toparam, variables(:N, 1:n_stratify))
    Ss = map(v -> v(t), variables(:S, 1:n_stratify, T = FnType))
    Is = map(v -> v(t), variables(:I, 1:n_stratify, T = FnType))
    Rs = map(v -> v(t), variables(:R, 1:n_stratify, T = FnType))
    C = map(toparam, variables(:C, 1:n_stratify, 1:n_stratify))
    uniform_contact_matrix = fill(1/n_stratify, (n_stratify, n_stratify))
    defs = Dict()

    for (i, nn) in enumerate(pops)
        defs[Ns[i]] = nn
        # ASSUMPTION: Assume only one person in each age group is infectious at the beginning of the simulation
        Ii = pop_assumption(i, nn)
        defs[Ss[i]] = nn - Ii
        defs[Is[i]] = Ii
        defs[Rs[i]] = 0
    end
    for i in eachindex(C)
        defs[C[i]] = uniform_contact_matrix[i]
    end
    eqs = [D.(Ss) .~ -β ./ Ns .* Ss .* (C * Is)
        D.(Is) .~ β ./ Ns .* Ss .* (C * Is) .- γ .* Is
        @. D(Rs) ~ γ * Is
        S ~ sum(Ss)
        I ~ sum(Is)
        R ~ sum(Rs)]
    @named model = ODESystem(eqs; defaults = defs)
    sys = structural_simplify(model)
    (C, sys)
end

(C, sys) = make_statified_model((2k, 2k, 2k))
prob = ODEProblem(sys, [], (0, tf))
sol = solve(prob)
plot(sol, leg = :topright)
```

## Question 1

> Simulate this model for the case where there is significant in-group contact
> preference – you may choose the numbers in the matrix to represent this in-
> group preference.

```@example scenario1
contact_matrix = fill(0.1, (3, 3)) + 0.2I
prob = ODEProblem(sys, [], (0, tf), vec(C .=> contact_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Simulate this model for the case where there is no contact between age groups.
> You may choose the numbers in the matrix, but ensure it meets the requirement
> of no contact between age groups.

```@example scenario1
prob = ODEProblem(sys, [], (0, tf), vec(C .=> Diagonal(contact_matrix)))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Simulate social distancing by scaling down the uniform contact matrix by a
> factor (e.g. multiply by 0.5)

```@example scenario1
prob = ODEProblem(sys, [], (0, tf), vec(C .=> 0.5 * contact_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Repeat 1.a.iv for the scenario where the young population has poor compliance
> with social distancing policies, but the old population is very compliant.

```@example scenario1
scaling = Diagonal([0.9, 0.8, 0.4])
prob = ODEProblem(sys, [], (0, tf), vec(C .=> scaling * contact_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Repeat 1.a for a younger-skewing population: `N_young = 3k, N_middle = 2k, N_old = 1k`

```@example scenario1
(C, sys) = make_statified_model((3k, 2k, 1k))
prob = ODEProblem(sys, [], (0, tf), ps)
sol = solve(prob)
plt1 = plot(sol, leg = :topright, title = "i")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> contact_matrix)])
sol = solve(prob)
plt2 = plot(sol, leg = :topright, title = "ii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> Diagonal(contact_matrix))])
sol = solve(prob)
plt3 = plot(sol, leg = :topright, title = "iii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> 0.5 * contact_matrix)])
sol = solve(prob)
plt4 = plot(sol, leg = :topright, title = "iv")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> scaling * contact_matrix)])
sol = solve(prob)
plt5 = plot(sol, leg = :topright, title = "v")
plot(plt1, plt2, plt3, plt4, plt5, size = (1000, 500))
```

> Repeat 1.a for an older-skewing population: `N_young = 1k, N_middle = 2k, N_old = 3k`

```@example scenario1
ps = Pair[]
pop = (1k, 2k, 3k)
for i in 1:n_stratify
    nn = pop[i]
    push!(ps, Ns[i] => nn)
    push!(ps, Ss[i] => nn - 1)
    push!(ps, Is[i] => 1)
    push!(ps, Rs[i] => 0)
end
prob = ODEProblem(sys, [], (0, tf), ps)
sol = solve(prob)
plt1 = plot(sol, leg = :topright, title = "i")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> contact_matrix)])
sol = solve(prob)
plt2 = plot(sol, leg = :topright, title = "ii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> Diagonal(contact_matrix))])
sol = solve(prob)
plt3 = plot(sol, leg = :topright, title = "iii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> 0.5 * contact_matrix)])
sol = solve(prob)
plt4 = plot(sol, leg = :topright, title = "iv")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> scaling * contact_matrix)])
sol = solve(prob)
plt5 = plot(sol, leg = :topright, title = "v")
plot(plt1, plt2, plt3, plt4, plt5, size = (1000, 500))
```

## Question 2

> Now find real contact matrix data and stratify the basic SIR model with the appropriate number of age groups to match the data found. To simulate the model with realistic initial values, find data on population distribution by age group. As in question 1, let gamma = 1/14 days, and let R0 = 5. Assume gamma, beta, and R0 are the same for all age groups.

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

xfs1 = (:all => xf_all_locations1, :work => xf_work1, :school => xf_school1, :home => xf_home1, :other => xf_other1)
xfs2 = (:all => xf_all_locations2, :work => xf_work2, :school => xf_school2, :home => xf_home2, :other => xf_other2)

to_cm(sheet) = Float64[sheet[i,j] for i = 2:17, j = 1:16]

# Load Belgium contact matrix
cm_belg = to_cm(xf_all_locations1["Belgium"])

# Load Belgium population distribution
pop_belg = values(CSV.read("data/2022_ Belgium_population_by_age.csv", DataFrame, header=3)[1, 2:17])

# Set up model
# Per MITRE: Assume that the same fixed fraction of the population in each stratum is initially infected. Here: 0.01%
pop_assumption(_, nn) = nn*0.0001
(C, sys_belg) = make_statified_model(pop_belg; pop_assumption)

prob = ODEProblem(sys_belg, [], (0, tf), vec(C .=> cm_belg))
sol = solve(prob)
plot(sol, leg=:topright)
```

> If the data you’ve found supports this, compare the situation for a country with significant multi-generational contact beyond two generations (as indicated by multiple contact matrix diagonal bandings), and for a country without. 

```@example scenario1
# Load India contact matrix
cm_india = to_cm(xf_all_locations1["India"])

# Load India population distribution
pop_india = values(CSV.read("data/2016_india_population_by_age.csv", DataFrame)[1, 3:18])

# Set up model
(C, sys_india) = make_statified_model(pop_india; pop_assumption)
prob = ODEProblem(sys_india, [], (0, tf), vec(C .=> cm_india))
sol = solve(prob)
plot(sol, leg=:topright)
```

> If the data supports this, try implementing interventions like: (1) School closures (2) Social distancing at work and other locations, but not at home.

> (1) School closures

> Prem et al Supplementary info, page 20

```@example scenario1
cm_school(xfs) = to_cm(xfs[:home][country]) + to_cm(xfs[:work][country]) + to_cm(xfs[:other][country]) # no school

cm_belgium_school_closure = cm_school(xfs1, "Belgium")
prob = ODEProblem(sys, [], (0, tf), vec(C .=> cm_belgium_school_closure))
sol = solve(prob)
plot(sol, leg=:topright)
```

> (2) Social distancing at work and other locations, but not at home.

> Prem et al sets social distancing to reduce contacts by half

```@example scenario1
cm_social_dist(xfs) = to_cm(xfs[:home][country]) + 0.5*to_cm(xfs[:work][country]) + 0.5*to_cm(xfs[:school][country]) + 0.5*to_cm(xfs[:other][country]) 

cm_belgium_social_dist = cm_social_dist(xfs1, "Belgium")
prob = ODEProblem(sys, [], (0, tf), vec(C .=> cm_belgium_social_dist))
sol = solve(prob)
plot(sol, leg=:topright)
```
