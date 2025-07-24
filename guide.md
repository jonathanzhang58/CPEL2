<!-- 2 from Eleanor 

finite diff 2.7.0

NLSolversBase 7.7.1

StatsFuns@0.9.5

add MultipleTesting@0.4.1

add Distributions@0.23.12


(Profile absent in Eleanor's file and mine, in CPel)



############################################################### -->
<!-- Correct solution:  -->
Install Julia 1.3.1 by going to: 
```https://julialang.org/downloads/oldreleases/```

and downloading the appropriate `tar.gz` file and moving that to your home directory in the JHPCE using your method choice.



Replace text of file at 

```
/yourhome/directory/.julia/environments/v1.3/Mainfest.toml
```

With magical file: 
<details>
<summary>Click to expand Manifest.toml</summary>


```
# This file is machine-generated - editing it directly is not advised

[[ArrayInterface]]
deps = ["LinearAlgebra", "Requires", "SparseArrays"]
git-tree-sha1 = "920136b6b8ae5bd28a3c45d68e55421dec156d7f"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "2.13.6"

[[Artifacts]]
deps = ["Pkg"]
git-tree-sha1 = "c30985d8821e0cd73870b17b0ed0ce6dc44cb744"
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.3.0"

[[Automa]]
deps = ["DataStructures", "Printf", "Random", "Test", "TranscodingStreams"]
git-tree-sha1 = "c81526bf5f6fb4616b4e22a3cd62ac20e255fd3c"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.0"

[[BED]]
deps = ["Automa", "BGZFStreams", "BioCore", "BufferedStreams", "ColorTypes", "FixedPointNumbers", "GenomicFeatures", "Indexes"]
git-tree-sha1 = "a80beaf68b53e09da91be3f0ec32a285cb7c3ebe"
uuid = "8e4a8c10-cb6b-11e8-08d2-83478d609d67"
version = "0.1.0"

[[BGZFStreams]]
deps = ["CodecZlib"]
git-tree-sha1 = "3aca54d25f8c30056577aa37ea68184da68df685"
uuid = "28d598bf-9b8f-59f1-b38c-5a06b4a0f5e6"
version = "0.3.2"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "9e62e66db34540a0c919d72172cc2f642ac71260"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "0.5.0"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[BioAlignments]]
deps = ["BioGenerics", "BioSequences", "BioSymbols", "IntervalTrees", "LinearAlgebra"]
git-tree-sha1 = "f610a3a965f187890edb0b1fdef4f30d77852edd"
uuid = "00701ae9-d1dc-5365-b64a-a3a3ebf5695e"
version = "2.0.0"

[[BioCore]]
deps = ["Automa", "BufferedStreams", "YAML"]
git-tree-sha1 = "476edbf4ef94594fff430a84ca96f86cb2327a71"
uuid = "37cfa864-2cd6-5c12-ad9e-b6597d696c81"
version = "2.0.5"

[[BioGenerics]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "017562e86afcd2a6a2a9220606a40b54604887c9"
uuid = "47718e42-2ac5-11e9-14af-e5595289c2ea"
version = "0.1.5"

[[BioSequences]]
deps = ["BioGenerics", "BioSymbols", "Combinatorics", "IndexableBitVectors", "Printf", "Random", "StableRNGs", "Twiddle"]
git-tree-sha1 = "255fbe7cd0f5fa346d549ba70f8dd255499e27e2"
uuid = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
version = "2.0.6"

[[BioSymbols]]
deps = ["Automa"]
git-tree-sha1 = "ec77888ac3e78f9d372c2b533bdb52668f9e2b09"
uuid = "3c28c6f8-a34d-59c4-9654-267d177fcfa9"
version = "4.0.4"

[[BufferedStreams]]
deps = ["Compat", "Test"]
git-tree-sha1 = "5d55b9486590fdda5905c275bb21ce1f0754020f"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.0.0"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "32ad4ece064a61855a35bdc34e3da0b495e01169"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.2"

[[ChangesOfVariables]]
deps = ["InverseFunctions", "LinearAlgebra", "Test"]
git-tree-sha1 = "3aa4bf1532aa2e14e0374c4fd72bed9a9d0d0f6c"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.10"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b9de8dc6106e09c79f3f776c27c62360d30e5eb8"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.9.1"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "d476eaeddfcdf0de15a67a948331c69a585495fa"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.47.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8e695f735fca77e9708e795eda62afdb869cbb70"
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.3.4+0"

[[CpelTdm]]
deps = ["BED", "BenchmarkTools", "Calculus", "Combinatorics", "Dates", "DelimitedFiles", "Distributed", "Distributions", "FASTX", "GenomicFeatures", "LinearAlgebra", "MultipleTesting", "Optim", "Random", "StatsBase", "Test", "XAM"]
git-tree-sha1 = "3aad17dca4e7b56a0171643312f25abacd4fa254"
repo-rev = "master"
repo-url = "https://github.com/jordiabante/CpelTdm.jl.git"
uuid = "63cc84f1-2b0d-464a-9e1f-381981f7b89b"
version = "0.1.0"

[[DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[DataStructures]]
deps = ["InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "88d48e133e6d3dd68183309877eac74393daa7eb"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.17.20"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "eb0c34204c8410888844ada5359ac8b96292cfd1"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.0.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "StaticArrays", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "501c11d708917ca09ce357bed163dbaf0f30229f"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.23.12"

[[DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[FASTX]]
deps = ["Automa", "BioGenerics", "BioSequences", "BioSymbols", "TranscodingStreams"]
git-tree-sha1 = "a980d6ac14c84c3ed17d0d07a0963a1c6d074b34"
uuid = "c2308a5c-f048-11e8-3e8a-31650f418d12"
version = "1.1.3"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "502b3de6039d5b78c76118423858d981349f3823"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.9.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "51c8f36c81badaa0e9ec405dcbabaf345ed18c84"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.1"

[[FixedPointNumbers]]
git-tree-sha1 = "4aaea64dd0c30ad79037084f8ca2b94348e65eaa"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.7.1"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "e2af66012e08966366a43251e1fd421522908be6"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.18"

[[GenomicFeatures]]
deps = ["BioGenerics", "DataStructures", "IntervalTrees"]
git-tree-sha1 = "51b0906aab4a9ae1a4095c27994361a582c5ebe1"
uuid = "899a7d2d-5c61-547b-bef9-6698a8d05446"
version = "2.1.0"

[[IndexableBitVectors]]
deps = ["Random", "Test"]
git-tree-sha1 = "b7f5e42dc867b8a8654a5f899064632dac05bc82"
uuid = "1cb3b9ac-1ffd-5777-9e6b-a3d42300664d"
version = "1.0.0"

[[Indexes]]
deps = ["BGZFStreams", "BioGenerics", "GenomicFeatures", "TranscodingStreams"]
git-tree-sha1 = "275bce824b40fd2e70358e0a652ba1b34172f240"
uuid = "4ffb77ac-cb80-11e8-1b35-4b78cc642f6d"
version = "0.1.3"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalTrees]]
git-tree-sha1 = "dc3b97bb5c9cb7c437f74027309f2c2f09a82aaf"
uuid = "524e6230-43b7-53ae-be76-1e9e4d08d11b"
version = "1.1.0"

[[InverseFunctions]]
deps = ["Dates", "Test"]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[JLLWrappers]]
git-tree-sha1 = "7cec881362e5b4e367ff0279dd99a06526d51a55"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.1.2"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[LibGit2]]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "b211c553c199c111d998ecdaf7623d1b89b69f93"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.12"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MultipleTesting]]
deps = ["Distributions", "Random", "SpecialFunctions", "StatsBase", "Test"]
git-tree-sha1 = "d6c79f17fbcdf492ca969ab1d3c3c402fd25d323"
uuid = "f8716d33-7c4a-5097-896f-ce0ecbd3ef6b"
version = "0.4.1"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "39d6bc45e99c96e6995cbddac02877f9b61a1dd1"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.7.1"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9db77584158d0ab52307f8c04f8e7c08ca76b5b3"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.3+4"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "3286df38aba45acf7445f3acd87b7b57b7c7feb7"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.2.4"

[[OrderedCollections]]
git-tree-sha1 = "d78db6df34313deaca15c5c0b9ff562c704fe1ab"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.5.0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "95a4038d1011dfdbde7cecd2ad0ac411e53ab1bc"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.10.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[Pkg]]
deps = ["Dates", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Test", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "e237232771fdafbae3db5c31275303e056afaa9f"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.10.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "28faf1c963ca1dc3ec87f166d92982e3c4a1f66d"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.0"

[[Rmath]]
deps = ["BinaryProvider", "Libdl", "Random", "Statistics"]
git-tree-sha1 = "2bbddcb984a1d08612d0c4abb5b4774883f6fa98"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.6.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "f6cb12bae7c2ecff6c4986f28defff8741747a9b"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.2"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "d8d8b8a9f4119829410ecd706da4cc8594a1e020"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.10.3"

[[StableRNGs]]
deps = ["Random", "Test"]
git-tree-sha1 = "b57c4216b6c163a3a9d674f6b9f7b99cdccdb959"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "0.1.2"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "b0e70de949b55b1f8d2b9b69dcb378a02e5e9e00"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.12.6"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics"]
git-tree-sha1 = "7bab7d4eb46b225b35179632852b595a3162cb61"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.2"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[Twiddle]]
git-tree-sha1 = "29509c4862bfb5da9e76eb6937125ab93986270a"
uuid = "7200193e-83a8-5a55-b20d-5d36d44a0795"
version = "1.1.2"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[XAM]]
deps = ["Automa", "BGZFStreams", "BioAlignments", "BioGenerics", "BioSequences", "GenomicFeatures", "Indexes", "Printf", "TranscodingStreams"]
git-tree-sha1 = "1be615462b8614e0012442dee0fe148b85732809"
uuid = "d759349c-bcba-11e9-07c2-5b90f8f05f7c"
version = "0.2.8"

[[YAML]]
deps = ["Base64", "Dates", "Printf"]
git-tree-sha1 = "209c033ada051007a934f7ab4738a4776bc041c3"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.2"

[[Zlib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "320228915c8debb12cb434c59057290f0834dbf6"
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.11+18"
```
</details>

\
Then run: 

```
Pkg.resolve()
Pkg.instantiate()
```

Hopefully it should work lmao. 
