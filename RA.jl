### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 70f48930-71d3-11ef-0553-13b51e04efa2
begin
  using Pkg
  Pkg.activate(@__DIR__)
  Pkg.instantiate()
end

# ╔═╡ 6abe7abd-c8af-4477-b661-6f8f2195487f
begin

using ReachabilityAnalysis, Plots

# initial-value problem specification
p = @ivp(x' = -x, x(0) ∈ Interval(1, 2))

# flowpipe computation
sol = solve(p, T=5)

# post-processing or plotting
plot(sol, vars=(0, 1), xlab="t", ylab="x(t)")
	
end

# ╔═╡ 9bbebfba-0f20-48e7-855e-548122908eaa
md"""## Caso decaimiento exponencial"""

# ╔═╡ 780a2dca-a646-4980-8501-c2f9b53cdcca
md"""## Caso oscilador armonico"""

# ╔═╡ 996c9055-5deb-4eba-b5be-683efd7511e0
begin

# initial-value problem specification
A = [0 1; -6 0]
X₀ = concretize(Singleton([1.0, 0.0]) ⊕ BallInf(zeros(2), 0.1))
problem = @ivp(x' = A * x, x(0) ∈ X₀)

# flowpipe computation
solp = solve(problem, T=3)

# post-processing or plotting
plot(solp, vars=(1, 2), xlab="x(t)", ylab="y(t)", ratio=.2)
	
end

# ╔═╡ Cell order:
# ╠═70f48930-71d3-11ef-0553-13b51e04efa2
# ╟─9bbebfba-0f20-48e7-855e-548122908eaa
# ╠═6abe7abd-c8af-4477-b661-6f8f2195487f
# ╟─780a2dca-a646-4980-8501-c2f9b53cdcca
# ╠═996c9055-5deb-4eba-b5be-683efd7511e0
