### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 0377899e-70a1-11ef-1823-7d1e830cb5b8
begin
  using Pkg
  Pkg.activate(@__DIR__)
  Pkg.instantiate()
end

# ╔═╡ 3ef48295-e259-473f-aa0d-166e71983a00
using LazySets, Plots, PlutoUI

# ╔═╡ a497e55b-db26-4ab4-b26b-fbb9f246101d
md"""## Convex hull"""

# ╔═╡ e5d68b73-7ac3-48de-81a4-1df4d4cddab2
ycoord = @bind ycoord Slider(-5:0.1:5, default=3)

# ╔═╡ 0614fbf8-5e6d-4d54-a623-089b94da97f1
begin
  X = Ball2(zeros(2), 1.0)
  Y = BallInf([3, ycoord], 1.0)

  plot(ConvexHull(X, Y), alpha=0.5, lw=5.0, ls=:dash, lc=:black, ratio=1.)
  plot!(X)
  plot!(Y)
  xlims!(-8, 8)
  ylims!(-6, 6)
end

# ╔═╡ d58c1e2c-9fef-4a9d-b0db-beb6325e383e
md"## Reduccion de zonotopes"

# ╔═╡ 938b0838-c01a-49ca-b620-21d3f30a15bd
Ngens_red = @bind Ngens_red Slider(2:1:10, default=10)

# ╔═╡ 19f10d31-ec54-4789-b408-f6c2a1cfa422
begin
	Z1 = rand(Zonotope, num_generators=10, seed=1234)
	Z1_red = reduce_order(Z1, Ngens_red/2)
	plot(Z1_red, ls=:dash, lw=2.0)
	plot!(Z1)
end

# ╔═╡ 65db68e0-c8a0-4335-bcf6-5ff89cc71337
md"""## Zonotope overapproximation"""

# ╔═╡ 1b55cbf6-f814-4882-9335-46780c8b4219
begin
 Hrand = rand(HPolygon, seed=1234)
 Zoa = overapproximate(Hrand, Zonotope, OctDirections(2))

 # Custom directions example: 
 # dirs = CustomDirections([c.a for c in constraints_list(Hrand)])
 #Zoa_aligned = overapproximate(Hrand, Zonotope, dirs)
 #plot(Zoa_aligned, ls=:dash)
 
 plot(Zoa)
 plot!(Hrand)
end

# ╔═╡ e4b19660-1b82-4fd6-a9ce-800a197562b1
md"""## Hausdorff distance"""

# ╔═╡ e187cc2e-6d6a-408f-b62f-92bcd96157e6
I2_inf = @bind I2_inf Slider(-1:0.01:6, default=3)

# ╔═╡ 488a6748-1ff5-4a12-8eef-51cbdc32d14b
begin
	I1 = Interval(1, 3)
	I2 = Interval(I2_inf, 6)
	plot(I1, lw=5.0)
	plot!(I2, lw=5.0)
	ylims!(-1, 1)
end

# ╔═╡ 6c6ff784-2e99-4d10-988b-69481e776937
@show hausdorff_distance(I1, I2)

# ╔═╡ 0d1a6480-949e-43a7-a127-801e592bfb19
begin
  xvals = -1:0.01:6
  dist = [hausdorff_distance(Interval(1, 3), Interval(xx, 6)) for xx in -1:0.01:6]
  plot(xvals, dist, lab="")
end

# ╔═╡ 8c9023ca-47b5-4177-b4c9-6fcbe4019ce9
md"""## Polyomial zonotopes"""

# ╔═╡ 54033924-4850-4a83-8e81-67eadcda0962
PZ_nsdiv = @bind PZ_nsdiv Slider(1:1:30, default=5)

# ╔═╡ 90f459c3-2f6e-44f2-b15a-e235d8789f66
begin
	c = zeros(2)
	G = [4 2 1 2; 4. 0 2 2]
	GI = hcat([1; 0.])
	E = [0 1 0 3; 0 0 1 1]
    PZ = SparsePolynomialZonotope(c, G, GI, E)
	plot(PZ, nsdiv=PZ_nsdiv, lw=0.0)
end

# ╔═╡ Cell order:
# ╠═0377899e-70a1-11ef-1823-7d1e830cb5b8
# ╠═3ef48295-e259-473f-aa0d-166e71983a00
# ╟─a497e55b-db26-4ab4-b26b-fbb9f246101d
# ╠═0614fbf8-5e6d-4d54-a623-089b94da97f1
# ╠═e5d68b73-7ac3-48de-81a4-1df4d4cddab2
# ╟─d58c1e2c-9fef-4a9d-b0db-beb6325e383e
# ╠═19f10d31-ec54-4789-b408-f6c2a1cfa422
# ╠═938b0838-c01a-49ca-b620-21d3f30a15bd
# ╟─65db68e0-c8a0-4335-bcf6-5ff89cc71337
# ╠═1b55cbf6-f814-4882-9335-46780c8b4219
# ╟─e4b19660-1b82-4fd6-a9ce-800a197562b1
# ╠═488a6748-1ff5-4a12-8eef-51cbdc32d14b
# ╠═e187cc2e-6d6a-408f-b62f-92bcd96157e6
# ╠═6c6ff784-2e99-4d10-988b-69481e776937
# ╠═0d1a6480-949e-43a7-a127-801e592bfb19
# ╟─8c9023ca-47b5-4177-b4c9-6fcbe4019ce9
# ╠═90f459c3-2f6e-44f2-b15a-e235d8789f66
# ╠═54033924-4850-4a83-8e81-67eadcda0962
