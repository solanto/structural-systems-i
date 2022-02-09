### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 21236c1f-c66c-46ee-80d0-fac24e625716
using
	Chain,
	Unitful,
	SparseArrays,
	Rotations,
	PaddedViews,
	CoordinateTransformations,
	LinearAlgebra

# ╔═╡ 955d1dd4-1e3e-4c3d-a845-26ad335a9cfd
md"# ⚙️ notebook setup"

# ╔═╡ 29708c92-615c-4e8e-a582-cb1e2329ccfe
# units
begin
	° = u"°"
	m = u"m"
	kN = u"kN"
	inch = u"inch"
	ksi = 1000u"psi"
	lbf = u"lbf"
	ft = u"ft"
	mm = u"mm"
	cm = u"cm"
end;

# ╔═╡ 7b699be4-73a5-4203-92cb-869f8164b692
md"""

# 🖼️ setting the scene

In this analysis, we'll consider the diagonal members of the bottom of the truss spine from Homework 1. Instead of inserting a horizontal member between the bototm two nodes, we fix the bottom two nodes. With nodes $\text{i}$, $\text{j}$, and $\text{k}$ as well wind load with magnitude $P$, the setup is as follows:

```
 P →k
   / \
  i   j
  ^   ^
```

The wind load selected is the 100-year wind load calculated in Homework 1, amounting. As well, the members' provided cross-sectional area and standard steel modulus of elasticity are recorded.

"""

# ╔═╡ 166c57dd-f02e-4dbc-b074-da6db827e235
begin
	P = 15.30kN * 7
	E = 29000ksi
	A = 66.9cm^2
end;

# ╔═╡ 15b33dd2-d93b-4fb8-baa0-42e087b87e5b
md"Additionally, we will use the building's measurements as necessary."

# ╔═╡ 41f91e06-5261-49a3-b92e-f26220f5f625
begin
	buildingwidth = 18m
	buildingheight = 22.2m
	trusswidth = 5m
end;

# ╔═╡ f184e7dc-ee4c-4368-969b-6c111027f459
md"""

# 🧮 devising an algorithm

## representing and scaling degrees of freedom

Member $\text{ij}$'s untransformed degrees of freedom $\bf{k}$, as an example, are:

| | $\delta_{\text{i},x}$ | $\delta_{\text{i},y}$ | $\delta_{\text{j},x}$ | $\delta_{\text{j},y}$ | 
|:-:|:-:|:-:|:-:|:-:|
| $F_{\text{i},x}$ | $\frac{EA}{L}$ | $\cdot$ | $-\frac{EA}{L}$ | $\cdot$ |
| $F_{\text{i},y}$ | $\cdot$ | $\cdot$ | $\cdot$ | $\cdot$ |
| $F_{\text{j},x}$ | $-\frac{EA}{L}$ | $\cdot$ | $\frac{EA}{L}$ | $\cdot$ |
| $F_{\text{j},y}$ | $\cdot$ | $\cdot$ | $\cdot$ | $\cdot$ |

With member length $L$. Factoring out $\frac{EA}{L}$, each of the matrix's four quadrants denotes a single degree of freedom at which a displacement and force in the axial direction are the only related quantities.

$$\bf{d}=\begin{bmatrix} 1 & 0 \\ 0 & 0 \end{bmatrix}$$

"""

# ╔═╡ 5e2c1a32-5ecb-4be3-b864-1b1b50679176
node_dofs = [1 0
			 0 0]

# ╔═╡ ee8c9021-7b7b-4305-89cc-a6286a8eb451
md"""

Thus:

$${\bf{k}} = \frac{EA}{L} \begin{bmatrix} \bf{d} & -\bf{d} \\ -\bf{d} & \bf{d} \end{bmatrix}$$

Such degrees of freedom can be transformed into a global coordinate system (so we can consider many angled members in a single coordinate system) by rotation. Here, we use a rotation matrix $\bf{T}$ on each axial degree of freedom.

$${\bf{TdT}}^{-1} = {\bf{TdT}}^\intercal$$

"""

# ╔═╡ adf29f84-c825-450e-b101-c2236ba662dd
rotate(system, angle) =
	let T = RotMatrix(angle)
		T * system * T'
	end

# ╔═╡ 121852ab-a833-445e-b4f2-7a6b912f48a4
md"""

## representing the structure as a graph

A system of members can be represented as a graph of nodes and links, with the additional information of the members' secion properties for practical analyses.

Here, we define a vector of nodes as well as a vector of their links. Links are defined as a three-way `Pair`. For instance, a link connecting nodes $\text{i}$ and $\text{j}$ separated by 1 unit of distance and at an angle 0° from the $x$-axis might be represented by either of the following:

```julia

:i => [1; 0] => :j

:i => Polar(1, 0) => :j

```

Finally, we store information about the members' cross sections.

"""

# ╔═╡ 38b6985e-ff8b-44f2-be97-6c3baaebc3ac
struct Section
	E::Number
	A::Number
end

# ╔═╡ fe05aa25-966b-4e61-8bdf-2bbc2a594ef7
struct Geometry{T, U <: Number, V <: Union{Vector{U},Polar}}
	nodes::Vector{T}
	section::Section
	links::Vector{Pair{T,Pair{V,T}}}
end

# ╔═╡ db05d34b-5ac2-423f-b10b-b4ece9e348ca
md"""

## constructing system matrices from graphs

First, we consider taking one directed link of the graph. We scale and transform each of the four degree of freedom quadrants, and we arrange (offset) each one to coincide with the appropriate degree of freedom in the final system matrix.

"""

# ╔═╡ 70e835a5-344a-4834-9361-b8e62212fdd1
direct(
	nodes::Vector{T},
	section::Section,
	vec::U,
	a::T,
	b::T
) where {T, U <: Union{Vector, Polar}} =
	let (; r, θ) = convert(Polar, vec)
		(; E, A) = section

		# scale DOFs
		scalednodaldofs = node_dofs * (E * A / r)

		#rotate DOFs
		dofs = rotate(scalednodaldofs, θ)

		# allow each to be offset based on node order
		size = 2 * length(nodes)
		findindex(node) = 2 * findfirst(nodes .== node) - 1
		indexA = findindex(a)
		indexB = findindex(b)
		dofsat(loc) = PaddedView(0 * unit(dofs[1]), dofs, (size, size), loc)
		sumdofsat(locs...) = locs .|> dofsat |> sum

		# offset and overlay terms, negating antidiagonal terms
		sumdofsat(
			(indexA, indexA),
			(indexB, indexB)
		) - sumdofsat(
			(indexA, indexB),
			(indexB, indexA)
		)
	end

# ╔═╡ 9ec4e489-0cbe-4243-8ae2-5bd156dc0f2c
md"We perform this transformation on all links in a structural system."

# ╔═╡ 2b971d1e-6f2b-42b8-8437-150d66db7820
connectiondofs(geometry::Geometry) =
	[direct(geometry.nodes, geometry.section, vec, origin, target)
	 for (origin, (vec, target))
	 in geometry.links]

# ╔═╡ 4174ba98-da7d-4c31-a3c8-97298c6f79a9
md"The system's final degree of freedom matrix is the sum of all member matrices."

# ╔═╡ bf16d56c-26f0-42c9-ac8b-9c28068b0b8f
dofs = sum ∘ connectiondofs

# ╔═╡ fda1cae4-cbe4-470e-a70c-f6b991a9e2e7
md"We constrain nodes in the $x$ and $y$ directions by setting their effects on the system to zero. As this is generally necessary to creating a solvable system, the system's geometric (graph) definitions are not transformed to their matrix representation until this step."

# ╔═╡ 9e8af36f-b232-4830-bdb3-df904598b7a7
constrain(
	geometry::Geometry{T,U},
	xcs=[]::Vector{T},
	ycs=[]::Vector{T}
) where {T, U} =
	let dofs = dofs(geometry) # transform graph to matrix

		# account for units
		dofunits = unit(dofs[1])

		# locate constrained rows and columns
		findindex(node) = findfirst(geometry.nodes .== node)
		xlocs = 2 .* findindex.(xcs) .- 1
		ylocs = 2 .* findindex.(ycs)
		locs = [xlocs..., ylocs...]

		# create a mask matrix to block out constrained terms
		s = size(dofs)
		mask = [a in locs || b in locs ? 0 : 1
				for (a, b)
				in getfield.(CartesianIndices(dofs), :I)]
		invertedmask = mask .|> Bool .|> !

		# apply constraints
		(Matrix(I, s...) .* invertedmask * dofunits) + (mask .* dofs)
	end

# ╔═╡ 4924c4e5-d837-4b12-b8e9-f9b0e710e22f
md"# ✏️ solving the problem"

# ╔═╡ 71e996d7-0b61-46bb-9fc1-6e6a7a16fef7
section = Section(E, A)

# ╔═╡ fc94f286-1a36-49e7-9f21-9bd5e07e54b1
system = @chain [:i, :j, :k] begin
	Geometry(section, [
		:i => [trusswidth / 2; buildingheight / 7] => :k,
		:j => [-trusswidth / 2; buildingheight / 7] => :k
	])
	constrain([:i, :j], [:i, :j])
end

# ╔═╡ a1ec62ae-71b9-4a0b-85c2-e755602271ab
loading = P * sparsevec([5], [1], 6)

# ╔═╡ c0d02e11-f2cb-49ce-a0df-da1456c7a1bc
system^-1 * loading .|> mm

# ╔═╡ 2445f198-baaa-47a3-94e9-553903844520
md"""

Node $\text{k}$ is displaced 0.422mm horizontally subject to loading due to Aruba's 100-year wind.

According to [STRUCTURE Magazine](https://www.structuremag.org/?p=15358#:~:text=Typical%20wind%20drift%20limits%20in%20common%20usage%20vary%20from%20H/100%20to%20H/600%20for%20total%20building%20drift%20and%20h/200%20to%20h/600%20for%20interstory%20drift), this is well below the typical range of allowable maximum interstory wind drift values—usually between 1/600 to 1/200 of the building's height.

"""

# ╔═╡ a30df75c-4c26-47c3-8cfe-4a9fd5cade2d
buildingheight / 7 ./ [600, 200] .|> mm

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PaddedViews = "5432bcbf-9aad-5242-b902-cca2824c8663"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Chain = "~0.4.10"
CoordinateTransformations = "~0.6.2"
PaddedViews = "~0.5.11"
Rotations = "~1.2.0"
Unitful = "~1.10.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra"]
git-tree-sha1 = "adf644ef95a5e26c8774890a509a55b7791a139f"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "405148000e80f70b31e7732ea93288aecb1793fa"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.2.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e6bf188613555c78062842777b116905a9f9dd49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "a635a9333989a094bddc9f940c04c549cd66afcf"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─955d1dd4-1e3e-4c3d-a845-26ad335a9cfd
# ╠═21236c1f-c66c-46ee-80d0-fac24e625716
# ╠═29708c92-615c-4e8e-a582-cb1e2329ccfe
# ╟─7b699be4-73a5-4203-92cb-869f8164b692
# ╠═166c57dd-f02e-4dbc-b074-da6db827e235
# ╟─15b33dd2-d93b-4fb8-baa0-42e087b87e5b
# ╠═41f91e06-5261-49a3-b92e-f26220f5f625
# ╟─f184e7dc-ee4c-4368-969b-6c111027f459
# ╠═5e2c1a32-5ecb-4be3-b864-1b1b50679176
# ╟─ee8c9021-7b7b-4305-89cc-a6286a8eb451
# ╠═adf29f84-c825-450e-b101-c2236ba662dd
# ╟─121852ab-a833-445e-b4f2-7a6b912f48a4
# ╠═fe05aa25-966b-4e61-8bdf-2bbc2a594ef7
# ╠═38b6985e-ff8b-44f2-be97-6c3baaebc3ac
# ╟─db05d34b-5ac2-423f-b10b-b4ece9e348ca
# ╠═70e835a5-344a-4834-9361-b8e62212fdd1
# ╟─9ec4e489-0cbe-4243-8ae2-5bd156dc0f2c
# ╠═2b971d1e-6f2b-42b8-8437-150d66db7820
# ╟─4174ba98-da7d-4c31-a3c8-97298c6f79a9
# ╠═bf16d56c-26f0-42c9-ac8b-9c28068b0b8f
# ╟─fda1cae4-cbe4-470e-a70c-f6b991a9e2e7
# ╠═9e8af36f-b232-4830-bdb3-df904598b7a7
# ╟─4924c4e5-d837-4b12-b8e9-f9b0e710e22f
# ╠═71e996d7-0b61-46bb-9fc1-6e6a7a16fef7
# ╠═fc94f286-1a36-49e7-9f21-9bd5e07e54b1
# ╠═a1ec62ae-71b9-4a0b-85c2-e755602271ab
# ╠═c0d02e11-f2cb-49ce-a0df-da1456c7a1bc
# ╟─2445f198-baaa-47a3-94e9-553903844520
# ╠═a30df75c-4c26-47c3-8cfe-4a9fd5cade2d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
