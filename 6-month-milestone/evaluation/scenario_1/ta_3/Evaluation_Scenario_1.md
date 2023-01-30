# Evaluation Scenario 1

## Stratified SIR

```@example scenario1
using EasyModelAnalysis, LinearAlgebra
using EasyModelAnalysis.ModelingToolkit: toparam
using EasyModelAnalysis.ModelingToolkit.Symbolics: FnType, variables
tf = 600
@parameters γ=1 / 14 R0=5
β = R0 * γ
@variables t S(t) I(t) R(t)
D = Differential(t)
n_stratify = 3
Ns = map(toparam, variables(:N, 1:n_stratify))
Ss = map(v -> v(t), variables(:S, 1:n_stratify, T = FnType))
Is = map(v -> v(t), variables(:I, 1:n_stratify, T = FnType))
Rs = map(v -> v(t), variables(:R, 1:n_stratify, T = FnType))
C = map(toparam, variables(:C, 1:n_stratify, 1:n_stratify))
k = 1000
uniform_contect_matrix = fill(0.33, (3, 3))
defs = Dict()

for i in 1:n_stratify
    nn = 2k
    defs[Ns[i]] = nn
    defs[Ss[i]] = nn - 1
    defs[Is[i]] = 1
    defs[Rs[i]] = 0
end
for i in eachindex(C)
    defs[C[i]] = uniform_contect_matrix[i]
end
eqs = [D.(Ss) .~ -β ./ Ns .* Ss .* (C * Is)
       D.(Is) .~ β ./ Ns .* Ss .* (C * Is) .- γ .* Is
       @. D(Rs) ~ γ * Is
       S ~ sum(Ss)
       I ~ sum(Is)
       R ~ sum(Rs)]
@named model = ODESystem(eqs; defaults = defs)
sys = structural_simplify(model)
prob = ODEProblem(sys, [], (0, tf))
sol = solve(prob)
plot(sol, leg = :topright)
```

## Question 1

> Simulate this model for the case where there is significant in-group contact
> preference – you may choose the numbers in the matrix to represent this in-
> group preference.

```@example scenario1
using Random
Random.seed!(123)
contect_matrix = 0.33 * rand(3, 3)
prob = ODEProblem(sys, [], (0, tf), vec(C .=> contect_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Simulate this model for the case where there is no contact between age groups.
> You may choose the numbers in the matrix, but ensure it meets the requirement
> of no contact between age groups.

```@example scenario1
prob = ODEProblem(sys, [], (0, tf), vec(C .=> Diagonal(contect_matrix)))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Simulate social distancing by scaling down the uniform contact matrix by a
> factor (e.g. multiply by 0.5)

```@example scenario1
prob = ODEProblem(sys, [], (0, tf), vec(C .=> 0.5 * contect_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Repeat 1.a.iv for the scenario where the young population has poor compliance
> with social distancing policies, but the old population is very compliant.

```@example scenario1
scaling = Diagonal([0.9, 0.8, 0.4])
prob = ODEProblem(sys, [], (0, tf), vec(C .=> scaling * contect_matrix))
sol = solve(prob)
plot(sol, leg = :topright)
```

> Repeat 1.a for a younger-skewing population: `N_young = 3k, N_middle = 2k, N_old = 1k`

```@example scenario1
ps = Pair[]
pop = (3k, 2k, 1k)
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

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> contect_matrix)])
sol = solve(prob)
plt2 = plot(sol, leg = :topright, title = "ii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> Diagonal(contect_matrix))])
sol = solve(prob)
plt3 = plot(sol, leg = :topright, title = "iii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> 0.5 * contect_matrix)])
sol = solve(prob)
plt4 = plot(sol, leg = :topright, title = "iv")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> scaling * contect_matrix)])
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

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> contect_matrix)])
sol = solve(prob)
plt2 = plot(sol, leg = :topright, title = "ii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> Diagonal(contect_matrix))])
sol = solve(prob)
plt3 = plot(sol, leg = :topright, title = "iii")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> 0.5 * contect_matrix)])
sol = solve(prob)
plt4 = plot(sol, leg = :topright, title = "iv")

prob = ODEProblem(sys, [], (0, tf), [ps; vec(C .=> scaling * contect_matrix)])
sol = solve(prob)
plt5 = plot(sol, leg = :topright, title = "v")
plot(plt1, plt2, plt3, plt4, plt5, size = (1000, 500))
```

## Question 2

> Now find real contact matrix data and stratify the basic SIR model with the appropriate number of age groups to match the data found. To simulate the model with realistic initial values, find data on population distribution by age group. As in question 1, let gamma = 1/14 days, and let R0 = 5. Assume gamma, beta, and R0 are the same for all age groups.

> If the data you’ve found supports this, compare the situation for a country with significant multi-generational contact beyond two generations (as indicated by multiple contact matrix diagonal bandings), and for a country without. 

Do simulations with India vs Belgium

> If the data supports this, try implementing interventions like: (1) School closures (2) Social distancing at work and other locations, but not at home.

Missing data right now