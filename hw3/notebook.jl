### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 03bfe070-8dbb-11ec-22ce-f9e9c204e33e
using
	Chain,
	Unitful,
	SparseArrays,
	Rotations,
	PaddedViews,
	CoordinateTransformations,
	LinearAlgebra

# ╔═╡ 1632f500-9126-4b29-8407-ace14db0148e
md"""

# ⚙️ notebook setup

We import packages used in homework assignment 2 as well, as we will be importing definitions from the previous notebook.

"""

# ╔═╡ 5748a0cc-8ac9-47cd-a1f4-e1f2dbdc717a
md"

We import the definitions in homework assignment two's notebook for reuse. Since `import` statements don't work well with Pluto, we use the `ingredients` function defined in a relevant [GitHub issue](https://github.com/fonsp/Pluto.jl/issues/115#issuecomment-661722426).

"

# ╔═╡ 16d58eff-e6db-43fb-96c2-a859b511bd5c
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

# ╔═╡ ed02bdf6-6a99-4797-b337-34a38d8c0222
(;
	Geometry,
	constrain,
	section,
	trusswidth,
	buildingheight,
	A
) = ingredients("../hw2/notebook.jl")

# ╔═╡ 2f994004-70b1-4820-9f35-d80791c346af
md"As well, we define units."

# ╔═╡ 3ef2b81b-111d-47b0-832d-cba65220457e
begin
	mm = u"mm"
	kN = u"kN"
	MPa = u"MPa"
end;

# ╔═╡ 567dae34-2e24-4832-99aa-7d652d76e5a9
md"""

# 🖼️ setting the scene

We have the following structure, as from homework assignment two:

```
 P →k
   / \
  i   j
  ^   ^
```

We'll use the displacement field we found after applying the 100-year wind load.

"""

# ╔═╡ 03ea3a33-6d30-4f67-b11d-0f5868846025
displacement = 0.422mm * sparsevec([5], [1], 6)

# ╔═╡ 47b13e20-93fa-431b-80f9-9e9b5dd1a2a6
md"""

# 🧮 determining internal forces

Since $\text{ij}$ and $\text{jk}$ are fully constrained at their bottoms, their internal forces are the opposites (negatives) of their reactions. While the system of external forces is defined asserting that those forces have zero effect, the system of internal forces does not have those restrictions. Thus, our system graph is the same as in notebook two with its original constraints removed and its effects negated.

"""

# ╔═╡ 9ebe7b44-35b5-486e-85fd-5bd41f4949b1
system = @chain [:i, :j, :k] begin
	Geometry(section, [
		:i => [trusswidth / 2; buildingheight / 7] => :k,
		:j => [-trusswidth / 2; buildingheight / 7] => :k
	])
	constrain
	_ * -1
end

# ╔═╡ dc665c51-1122-4ff4-a1de-f3ee8570d24e
md"The internal forces are the result of the system subject to the displacement field."

# ╔═╡ 88569649-01a6-4846-a98a-7106b2356ef0
nodals = system * displacement .|> kN

# ╔═╡ 2ae24376-88ef-48cd-9a9f-213e072bd998
md"# ✅ ensuring tensile diagonal safety"

# ╔═╡ ac6bfbbf-7688-42b3-a908-956c91387172
md"""

The internal forces in $\text{ij}$ and $\text{jk}$ are opposite, since the structure is symmetric and subject to a symmetric load.

Since forces are represented by $x$ and $y$ component vectors, we calculate the members' (identical) absolute magnitudes of internal force. This is $P_r$: required load.

"""

# ╔═╡ a4cd6723-206a-420f-85ac-fcd9e3d682ac
P_r = (sqrt ∘ sum)(nodals[1:2].^2)

# ╔═╡ 49f519d2-08be-4de3-b2cc-0202cf67ccd0
md"""

This required load value matches the load calculated in homework one.

According to the load resistance factor design guidelines described in [AISC 360-16 (16.1-83)](https://www.aisc.org/globalassets/aisc/publications/standards/a360-16w-rev-june-2019.pdf)—with critical tensile load $P_c$, tensile resistance factor $\phi_t$, nominal tensile strength $P_n$, yield stress $F_y$, and gross cross-sectional area $A_g$:

$$P_c=\phi_t P_n=\phi_t F_y A_g$$

$$\phi_t=0.90$$

According to [MatWeb](http://www.matweb.com/search/datasheet.aspx?matguid=9ced5dc901c54bd1aef19403d0385d7f&ckck=1):

$$F_y=345\text{MPa}$$

"""

# ╔═╡ 22e29638-2d05-41cf-9347-040fb28b3533
P_c = 0.90 * 345MPa * A |> kN

# ╔═╡ b22c5817-9ece-4133-965b-e27c9ebf5c7d
md"""

According to H3-6 in the AISC booklet—with no moment, shear, or torsion—$P_r$ is allowable given $P_r$ for:

$$\frac{P_r}{P_c} \le 1$$

In other words, the required load can not exceen 100% of the critical load.

"""

# ╔═╡ 8bcbe91c-594d-443e-ae5f-e3eab25cb254
P_r / P_c

# ╔═╡ bb3bf6f0-f15c-4cad-928b-d0c35b13b3e2
md"""

The required tensile load is only 4.17% of the critical tensile load, indicating that the load is allowable.

"""

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
git-tree-sha1 = "8d0c8e3d0ff211d9ff4a0c2307d876c99d10bdf1"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "95c6a5d0e8c69555842fc4a927fc485040ccc31c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.5"

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
# ╟─1632f500-9126-4b29-8407-ace14db0148e
# ╠═03bfe070-8dbb-11ec-22ce-f9e9c204e33e
# ╟─5748a0cc-8ac9-47cd-a1f4-e1f2dbdc717a
# ╟─16d58eff-e6db-43fb-96c2-a859b511bd5c
# ╠═ed02bdf6-6a99-4797-b337-34a38d8c0222
# ╟─2f994004-70b1-4820-9f35-d80791c346af
# ╠═3ef2b81b-111d-47b0-832d-cba65220457e
# ╟─567dae34-2e24-4832-99aa-7d652d76e5a9
# ╠═03ea3a33-6d30-4f67-b11d-0f5868846025
# ╟─47b13e20-93fa-431b-80f9-9e9b5dd1a2a6
# ╠═9ebe7b44-35b5-486e-85fd-5bd41f4949b1
# ╟─dc665c51-1122-4ff4-a1de-f3ee8570d24e
# ╠═88569649-01a6-4846-a98a-7106b2356ef0
# ╟─2ae24376-88ef-48cd-9a9f-213e072bd998
# ╟─ac6bfbbf-7688-42b3-a908-956c91387172
# ╠═a4cd6723-206a-420f-85ac-fcd9e3d682ac
# ╟─49f519d2-08be-4de3-b2cc-0202cf67ccd0
# ╠═22e29638-2d05-41cf-9347-040fb28b3533
# ╟─b22c5817-9ece-4133-965b-e27c9ebf5c7d
# ╠═8bcbe91c-594d-443e-ae5f-e3eab25cb254
# ╟─bb3bf6f0-f15c-4cad-928b-d0c35b13b3e2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
