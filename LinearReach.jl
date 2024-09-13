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

# ╔═╡ da84839f-d379-486a-bac3-01c09b8c0062
begin
  using Pkg
  Pkg.activate(@__DIR__)
  Pkg.instantiate()
end

# ╔═╡ 59a8d99a-9da7-44d7-9458-18cdf9a81a5a
using LazySets, Plots, LinearAlgebra, PlutoUI

# ╔═╡ c684fa0c-3d84-4302-90d1-b7147ac2fd9f
md"""
## Modelo

Consideremos el oscilador armonico simple, $u''(t) + \omega^2 u(t) = 0$ con $\omega = 4\pi$. Con el cambio de variables $x(t) = [u(t), u'(t)]$ obtenemos:

$x'(t) = Ax(t)$ con $A = \begin{pmatrix} 0 & 1 \\ -(4\pi)^2 & 0 \end{pmatrix}$
"""

# ╔═╡ 3013fa43-1876-4a2d-98e8-5e1c7f6d6839
md"""
## Caso 1: condición inicial puntual

Sea $u(0) = 1$, $v(0) = 0$.
"""

# ╔═╡ 622f5552-8e8c-4dd2-b099-3a9d03e6fe60
τ = @bind τ Slider(0.001:0.005:0.1, default=0.02)

# ╔═╡ 19979b38-715f-11ef-38e4-fdc5840341bb
begin
	# modelo
	A = [0 1; -6 0]
	X₀ = Singleton([1.0, 0.0])

	# discretizacion, ver GLG09
	Φ = exp(A * τ)
	W = ZeroSet(2)
	Anorm = opnorm(A, Inf)
	RW = norm(W, Inf)
	RX₀ = norm(X₀, Inf)
    α = (exp(τ * Anorm) - 1 - τ * Anorm) * (RX₀ + RW / Anorm)
	B = BallInf(zeros(2), 1.0)
	
	# aproximacion inicial
	S₀ = CH(X₀, Φ * X₀ ⊕ τ * W ⊕ α * B)
end

# ╔═╡ 985010de-201f-4c15-b4ce-305f267a829b
begin
	plot(Hyperrectangle(low=[0.9, -0.2], high=[1.1, 0.2]), color=:white, lc=:white)
	plot!(S₀, lab="S₀", c=:lightblue)
	plot!(X₀, lab="X₀", legend=:topright, alpha=1.)
	plot!(Φ * X₀, lab="Φ * X₀", alpha=1.)
end

# ╔═╡ 33f13241-3e25-4255-a5bf-f50b08f52265
N = @bind N Slider(10:10:50, default=20)

# ╔═╡ f25b2b58-e166-4235-9338-fefaf0afa6a1
begin
   β = (exp(τ * Anorm) - 1 - τ * Anorm) * RW / Anorm
   successor(S) = Φ * S ⊕ τ * W ⊕ β * B

  function compute_successors(S₀, N)
	   result = LazySet[S₀]
	   for k in 1:(N-1)
		   Sk = successor(last(result))
		   push!(result, Sk)
	   end
	   result
  end
  @time result = compute_successors(S₀, N)
end

# ╔═╡ a21cd24a-66d9-4a74-b46e-554259dd9ea6
last(result)

# ╔═╡ 7ed4eefb-3b6a-483f-8286-5886a50c35eb
N

# ╔═╡ d648f8e6-8714-4e0a-bf91-d8be64aa1de9
begin
	plt = plot()
	[plot!(plt, res) for res in result]
	plt
end

# ╔═╡ 8fbe74db-b773-4782-9140-b96e6e8ce844
md"""## Reduciendo en cada paso temporal"""

# ╔═╡ 0acbeee7-9025-44b0-9016-f1f341252aad
begin
  function compute_successors_red(S₀, N)
	   result = LazySet[S₀]
	   for k in 1:(N-1)
		   Sk = concretize(successor(last(result)))
		   Sk = overapproximate(Sk, Zonotope, OctDirections(2))
		   push!(result, Sk)
	   end
	   result
  end
  @time result_red = compute_successors_red(S₀, N)
end

# ╔═╡ ac10c391-6554-45a3-aa7c-3363b87ba354
begin
	plt_red = plot()
	[plot!(plt_red, res) for res in result_red]
	plt_red
end

# ╔═╡ ba9617a9-51e3-4fa3-826b-96795077bcaa
md"""## Esquema sin efecto de wrapping"""

# ╔═╡ d5d5f30b-9a48-4ec3-aa13-c383e632e2b3
begin
  function compute_successors_nowrap(S₀, N)	   
	   Z0 = CH(X₀, Φ * X₀ ⊕ τ * W ⊕ α * B)
	   V0 = CH(X₀, Φ * X₀ ⊕ τ * W ⊕ α * B)
	   Y0 = τ * W ⊕ β * B

	   Zk_sequence = LazySet[Z0]
	   Vk_sequence = LazySet[V0]
	   Yk_sequence = LazySet[Y0]

	  result = LazySet[]
      push!(result, Z0 ⊕ V0)
      dirs = OctDirections(2)
	  
	  for k in 1:(N-1)
		   Zk = Φ * last(Zk_sequence)
		   Vk = Φ * last(Vk_sequence)
		   Yk = last(Yk_sequence) ⊕ last(Vk_sequence)

		   push!(Zk_sequence, Zk)
		   push!(Vk_sequence, Vk)
		   push!(Yk_sequence, Yk)
		   Sk = Zk ⊕ Vk
		   Sk = overapproximate(Sk, Zonotope, dirs)
		   push!(result, overapproximate(Zk ⊕ Vk, Zonotope, dirs))
	   end
	   result
  end
  @time result_nowrap = compute_successors_nowrap(S₀, N)
end

# ╔═╡ 7dd87faa-1d43-4aac-a6e9-6b1d056ca69a
begin
	plt_nowrap = plot()
	[plot!(plt_nowrap, res) for res in result_nowrap]
	plt_nowrap
end

# ╔═╡ a2c0d00e-0d9b-4135-ae16-c428ef5fc410
result_nowrap[end]

# ╔═╡ 8f811510-77e1-4238-9ce3-e1c56fd6c1e4
md"""
## Caso 2: condición inicial no puntual

Consideramos como condicion inicial una caja centrada en [1, 0].

Alcanza con volver a la celda superior y definir

```
X₀ = Singleton([1.0, 0.0]) ⊕ BallInf(zeros(2), 1e-1)
```
Variar el radio para ver el compartamiento segun el diametro de la inicial.
"""

# ╔═╡ Cell order:
# ╠═da84839f-d379-486a-bac3-01c09b8c0062
# ╠═59a8d99a-9da7-44d7-9458-18cdf9a81a5a
# ╟─c684fa0c-3d84-4302-90d1-b7147ac2fd9f
# ╟─3013fa43-1876-4a2d-98e8-5e1c7f6d6839
# ╠═19979b38-715f-11ef-38e4-fdc5840341bb
# ╠═622f5552-8e8c-4dd2-b099-3a9d03e6fe60
# ╠═985010de-201f-4c15-b4ce-305f267a829b
# ╠═f25b2b58-e166-4235-9338-fefaf0afa6a1
# ╠═a21cd24a-66d9-4a74-b46e-554259dd9ea6
# ╠═7ed4eefb-3b6a-483f-8286-5886a50c35eb
# ╠═33f13241-3e25-4255-a5bf-f50b08f52265
# ╠═d648f8e6-8714-4e0a-bf91-d8be64aa1de9
# ╟─8fbe74db-b773-4782-9140-b96e6e8ce844
# ╠═0acbeee7-9025-44b0-9016-f1f341252aad
# ╠═ac10c391-6554-45a3-aa7c-3363b87ba354
# ╟─ba9617a9-51e3-4fa3-826b-96795077bcaa
# ╠═d5d5f30b-9a48-4ec3-aa13-c383e632e2b3
# ╠═7dd87faa-1d43-4aac-a6e9-6b1d056ca69a
# ╠═a2c0d00e-0d9b-4135-ae16-c428ef5fc410
# ╟─8f811510-77e1-4238-9ce3-e1c56fd6c1e4
