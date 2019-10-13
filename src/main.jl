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

end