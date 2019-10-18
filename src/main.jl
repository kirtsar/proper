module ProperFamilies
    using Combinatorics
    using Base.Iterators   # for product
    using LightGraphs
    using LightGraphs.Experimental
    using DataStructures  # DefaultDict




    include("zhegalkin.jl")
    include("family.jl")
    include("proper.jl")
    include("pairing.jl")
    include("latsquares.jl")
    include("orthogonality.jl")
    include("minimal_repr.jl")
    include("plot_family.jl")
    include("assoc_triples.jl")

    export Monom, ZhegFun, Family, LatinSquare, Pairing
    export orb, stab, permute, is_proper
    export latin, latin!, assoc, triples
    export is_ortho, has_ortho, generate_ortho_family
    export plot

end