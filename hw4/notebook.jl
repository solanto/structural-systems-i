### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 28965851-aa66-479d-8828-4b901271df45
using
	Unitful,
	QuadGK,
	Calculus,
	LinearAlgebra,
	Rotations,
	Chain,
	SparseArrays,
	Rotations,
	PaddedViews,
	CoordinateTransformations,
	LinearAlgebra

# ╔═╡ d55eacdf-fb9f-417d-99ff-f168165b11ee
md"""
# ⚙️ notebook setup
"""

# ╔═╡ 5c751d6b-076e-4b68-ad35-4fa2214717e5
begin
	m = u"m"
	° = u"°"
	kN = u"kN"
	cm = u"cm"
	ksi = 1000u"psi"
	mm = u"mm"
end;

# ╔═╡ 153fc91e-711b-40e1-9efd-0f71d540b0fd
integral = first ∘ quadgk

# ╔═╡ fb63659d-e8fe-4c27-b029-ba28bd31f6cd
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(
				 include(mapexpr::Function, x) = 
				 	$(Expr(:top, :include))(mapexpr, $name, x)
			 ),
             :(include($path))
		)
	)
	
	return m
end

# ╔═╡ dd751d0d-1263-40d1-bcdc-ab3cf83011bf
(;
	Geometry,
	constrain,
	Section
) = ingredients("../hw2/notebook.jl")

# ╔═╡ 73b4dcf9-6a8b-4380-9545-d57c9692eb6f
md"""
# 📐 defining geometry

The bottom chord of the structure can be modeled with the following equation of horizontal distance in meters.
"""

# ╔═╡ 056d80e0-9f08-11ec-3780-4d0e6312b10a
 bottomchord(x) = (21 - 4) * (1 - (2x / 73 - 1)^2)

# ╔═╡ 8edc697c-8669-40c3-aa00-fe834efaa790
md"""
Elements join at equally spaced nodes along each arc's length.
"""

# ╔═╡ 830b3b37-ee00-4903-9678-b32324e59e3d
function arclength(f, a, b)
	dfdx = derivative(f)

	integral(x -> √(1 + dfdx(x)^2), a, b)
end

# ╔═╡ 195442e5-45f4-4971-a56e-33d4c7a15c66
bottomlength(x) = arclength(bottomchord, 0, x)

# ╔═╡ 611b1e1d-13dd-404c-a1c4-02b2cbfa3d0a
elementlength = bottomlength(73) / 36

# ╔═╡ 29282c7d-74c9-4bf2-b0c6-397acab72c48
md"""
We iterate through the segments that make up the chord and represent each element as a vector.
"""

# ╔═╡ bb394427-c819-4f10-b946-1a24a998c752
struct Segments
	curve
	dcurve
	length

	Segments(curve, length) = new(curve, derivative(curve), length)
end

# ╔═╡ 4b309137-8869-4414-834a-11df0970a77c
function Base.iterate(
	(; curve, dcurve, length)::Segments,
	state = 0
)
	segment = normalize([1; dcurve(state)]) * length

	(
		segment,
		state + segment[1]
	)
end

# ╔═╡ f4c773da-8311-41c0-8500-5ee49e932bd4
Base.length(T::Iterators.Take{Segments}) = T.n

# ╔═╡ 44765efa-edb7-482d-9cfc-107424796aee
bottomelementvecs =
	Iterators.take(Segments(bottomchord, elementlength), 36) |> collect

# ╔═╡ 9d79907f-96f0-43b3-a812-d9bb789d997c
md"""
To use these in our structural analysis toolkit, we embed them in connection definitions, numbering each node along the bottom of the chord in ascending order left to right.
"""

# ╔═╡ 92da2895-aea5-4211-b06a-3d4d6e4278a7
bottomelements = [a => vec => b for (a, vec, b) in zip(1:36, bottomelementvecs, 2:37)]

# ╔═╡ 67b19c0c-5f76-4e3c-81c8-cb1dfcc382b0
md"""
Diagonals extending upward from these elements form equilaterial triangles on each element. These are numbered in ascending order left to right just as the nodes to the left of them, but with 100 added to their node IDs.
"""

# ╔═╡ 20c8ead1-6aa4-48b7-b499-36c62a72d7fb
leftdiags = [a => RotMatrix(60°) * vec => a + 100
			 for (a, (vec, b)) in bottomelements]

# ╔═╡ 3bb84c92-9733-48c9-9f45-8bff2c77b44c
rightdiags = [a + 100 => RotMatrix(-60°) * vec => b
			  for (a, (vec, b)) in bottomelements]

# ╔═╡ e488c9aa-a248-46f6-8c8e-b257f1a5a165
md"""
Finally, the elements connecting the tops of the triangles are fully constrained by the known geometry.
"""

# ╔═╡ cababe6e-dce1-4ae2-bb3b-ab752f86c1b2
topelements = let
	topnodes = first.(rightdiags)
	rightdiagvecs = [vec for (a, (vec, b)) in rightdiags]
	info = zip(
		topnodes[1:end-1],
		topnodes[2:end],
		bottomelementvecs[1:end-1],
		rightdiagvecs[1:end-2],
		rightdiagvecs[2:end-1]
	)

	[a => bottom + leftb - lefta => b
	 for (a, b, bottom, lefta, leftb) in info]
end

# ╔═╡ 4cf14548-a5a7-4dff-9366-9ed3014ee1f9
md"""
This yields the full set of elements.
"""

# ╔═╡ 15efa6c4-7ecb-4606-879f-dbf3a496b321
elements =
	[a => Vector(vec * 1m) => b
	 for (a, (vec, b)) in [bottomelements; leftdiags; rightdiags; topelements]]

# ╔═╡ 630e3b46-deaa-4922-b7fb-ba77046761ec
md"""
And the corresponding set of nodes.
"""

# ╔═╡ 31037e05-334a-464e-8dab-5afce71ade59
nodes = Iterators.flatten(
	(a, b) for (a, (vec, b)) in elements
) |> Set |> collect

# ╔═╡ 9d6f13f4-68cf-4cd3-87b8-3e124b3fd2aa
md"""
# 🔨 defining physical properties

From the provided paper, the force densities (weight per unit volume) of concrete and steel are as follows:
"""

# ╔═╡ 21048f0d-bde7-491d-b067-4d74826881e5
begin
	f_c = 25kN/m^3
	f_s = 78.5kN/m^3
end;

# ╔═╡ cfbd1674-39be-4e78-a7bf-452ed193e43c
md"""
The using the numbers of each type of element and the paper's provided section properties, we calculate the average area of the elements used in the structure.
"""

# ╔═╡ dbfa860c-086c-497b-a3ae-4e3a2c9d3b38
bararea = (
	(
		length(topelements) * 0.47m * 0.3m +
		length(bottomelements) * 0.35m * 0.35m +
		(length(rightdiags) + length(leftdiags)) * 0.1m * 0.1m
	) / (
		length(topelements) +
		length(bottomelements) +
		length(rightdiags) + 
		length(leftdiags)
	)
)

# ╔═╡ 2d48f4cc-e2ed-41a6-beb8-c579bb1deb43
md"""
From this, we define the section properties of the structure.
"""

# ╔═╡ a72a472c-b351-41e8-85cc-7cda07d65087
section = Section(29000ksi, bararea)

# ╔═╡ b4553f1a-f0e3-4c5d-9170-7b718ea222a5
md"""
Using the number of vaults in the structure, the structure's dimensions, and the concrete shell's thickness: the self-weight of each vault is calculated.
"""

# ╔═╡ d39b4a96-42ee-4826-9e9f-070744afaf1f
numvaults = 28

# ╔═╡ 23f88358-709a-4e2a-a763-842d2e51f2aa
vaultweight = kN(
	f_c * 5cm * bottomlength(73) * 1m * 230m +
	f_s * length(elements) * elementlength * 1m * bararea
) / numvaults

# ╔═╡ b19c76b8-8322-4134-944e-9108ec00a640
md"""
# 🧮 solving the problem
"""

# ╔═╡ 2019757e-30e5-4e08-9a0e-ea9f57d0c4da
system = @chain nodes begin
	Geometry(section, elements)
	constrain([1, 36], [1, 36])
end

# ╔═╡ f394e161-ae14-4167-b53c-40258c618c24
loading = (
	[[0kN; -vaultweight / length(nodes)] for _ in 1:length(nodes)]
	|> Iterators.flatten
	|> collect
)

# ╔═╡ ae12c0b7-b0c2-48d3-8c83-b92b954ba9ed
nodes .=> Iterators.partition(system^-1 * loading .|> mm, 2) .|> collect

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Calculus = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PaddedViews = "5432bcbf-9aad-5242-b902-cca2824c8663"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Calculus = "~0.5.1"
Chain = "~0.4.10"
CoordinateTransformations = "~0.6.2"
PaddedViews = "~0.5.11"
QuadGK = "~2.4.2"
Rotations = "~1.3.0"
Unitful = "~1.11.0"
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
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

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

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

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
git-tree-sha1 = "90b158083179a6ccbce2c7eb1446d5bf9d7ae571"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.7"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

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
git-tree-sha1 = "3f7cb7157ef860c637f3f4929c8ed5d9716933c6"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.7"

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
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "522770af103809e8346aefa4b25c31fbec377ccf"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "a167638e2cbd8ac41f9cd57282cab9b042fa26e6"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.0"

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
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

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
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

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
# ╟─d55eacdf-fb9f-417d-99ff-f168165b11ee
# ╠═28965851-aa66-479d-8828-4b901271df45
# ╠═5c751d6b-076e-4b68-ad35-4fa2214717e5
# ╠═153fc91e-711b-40e1-9efd-0f71d540b0fd
# ╟─fb63659d-e8fe-4c27-b029-ba28bd31f6cd
# ╠═dd751d0d-1263-40d1-bcdc-ab3cf83011bf
# ╟─73b4dcf9-6a8b-4380-9545-d57c9692eb6f
# ╠═056d80e0-9f08-11ec-3780-4d0e6312b10a
# ╟─8edc697c-8669-40c3-aa00-fe834efaa790
# ╠═830b3b37-ee00-4903-9678-b32324e59e3d
# ╠═195442e5-45f4-4971-a56e-33d4c7a15c66
# ╠═611b1e1d-13dd-404c-a1c4-02b2cbfa3d0a
# ╟─29282c7d-74c9-4bf2-b0c6-397acab72c48
# ╠═bb394427-c819-4f10-b946-1a24a998c752
# ╠═4b309137-8869-4414-834a-11df0970a77c
# ╠═f4c773da-8311-41c0-8500-5ee49e932bd4
# ╠═44765efa-edb7-482d-9cfc-107424796aee
# ╟─9d79907f-96f0-43b3-a812-d9bb789d997c
# ╠═92da2895-aea5-4211-b06a-3d4d6e4278a7
# ╟─67b19c0c-5f76-4e3c-81c8-cb1dfcc382b0
# ╠═20c8ead1-6aa4-48b7-b499-36c62a72d7fb
# ╠═3bb84c92-9733-48c9-9f45-8bff2c77b44c
# ╟─e488c9aa-a248-46f6-8c8e-b257f1a5a165
# ╠═cababe6e-dce1-4ae2-bb3b-ab752f86c1b2
# ╟─4cf14548-a5a7-4dff-9366-9ed3014ee1f9
# ╠═15efa6c4-7ecb-4606-879f-dbf3a496b321
# ╟─630e3b46-deaa-4922-b7fb-ba77046761ec
# ╠═31037e05-334a-464e-8dab-5afce71ade59
# ╟─9d6f13f4-68cf-4cd3-87b8-3e124b3fd2aa
# ╠═21048f0d-bde7-491d-b067-4d74826881e5
# ╟─cfbd1674-39be-4e78-a7bf-452ed193e43c
# ╠═dbfa860c-086c-497b-a3ae-4e3a2c9d3b38
# ╟─2d48f4cc-e2ed-41a6-beb8-c579bb1deb43
# ╠═a72a472c-b351-41e8-85cc-7cda07d65087
# ╟─b4553f1a-f0e3-4c5d-9170-7b718ea222a5
# ╠═d39b4a96-42ee-4826-9e9f-070744afaf1f
# ╠═23f88358-709a-4e2a-a763-842d2e51f2aa
# ╟─b19c76b8-8322-4134-944e-9108ec00a640
# ╠═2019757e-30e5-4e08-9a0e-ea9f57d0c4da
# ╠═f394e161-ae14-4167-b53c-40258c618c24
# ╠═ae12c0b7-b0c2-48d3-8c83-b92b954ba9ed
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
