# check whether given family of functions is proper or not
function is_proper(fam :: Family)
    size = length(fam.fs)
    maxNo = 2^size - 1
    for x in 0 : maxNo
        for y in 0 : maxNo
            if xor(x, y) != 0
                fx = fam(x)
                fy = fam(y)
                res = xor(x, y) & ~(xor(fx, fy))
                if res == 0
                    return false
                end
            end
        end
    end
    return true
end


# try to construct counterexample in order
# to show that the family is not proper
function counterexample(fam :: Family)
    size = length(fam.fs)
    maxNo = 2^size - 1
    for x in 0 : maxNo
        for y in 0 : maxNo
            if xor(x, y) != 0
                fx = fam(x)
                fy = fam(y)
                res = xor(x, y) & ~(xor(fx, fy))
                if res == 0
                    println("x   : ", bitstring(x)[end - size + 1 : end])
                    println("y   : ", bitstring(y)[end - size + 1 : end])
                    println("xor : ", bitstring(xor(x, y))[end - size + 1 : end])
                    println("fx  : ", bitstring(fx)[end - size + 1 : end])
                    println("fy  : ", bitstring(fy)[end - size + 1 : end])
                    return false
                end
            end
        end
    end
    return true
end


# generate all proper families of given size
# brute force approach
function generate_all_proper(numVars)
    # generate all possible coefficients
    # for i-th function
    # i-th function must NOT depend on x_i
    possibleCoefs = [Int[] for i in 1 : numVars]
    # we can exclude 0
    # because f_i -> f_i + a_i
    # is also a proper family
    for i in 1 : (2^numVars - 1)
        bits = reverse(bitstring(i))
        for k in 1 : numVars
            if bits[k] == '0'
                append!(possibleCoefs[k], i)
            end
        end
    end
    # generate all possible subsets of coefficients
    funSets = []
    for i in 1 : numVars
        push!(funSets, powerset(possibleCoefs[i]))
    end
    numProper = 0
    proper_families = Array{Family}([])

    for funArray in product(funSets...)
        fam = Family(ZhegFun.(funArray)...)
        if is_proper(fam)
            numProper += 1
            push!(proper_families, fam)
        end
    end
    return (2^numVars * numProper, proper_families)
end



# permute proper family of functions
# fi(x1, ... xn) -> f_pi(i)(x_pi(1), ... x_pi(n))
# this family is also proper
# matrix graph is the same
function permute(f :: ZhegFun, perm)
    g = ZhegFun(Monom.(zeros(Int, length(f.ANF))))
    for (ind, mon) in enumerate(f.ANF)
        g.ANF[ind] = permute(mon, perm)
    end
    sort!(g.ANF, by = x -> x.val)
    return g
end


# given a family of functions
# apply all possible permutations to this family
# permutation of the type I
# fi(x1, ... xn) -> f_pi(i)(x_pi(1), ... x_pi(n))
# obtain orbit of the family
function orb(fam :: Family{T}) where T
    n = length(fam.fs)
    perms = permutations(collect(1 : n))
    resOrbit = Set{Family{T}}([])
    for perm in perms
        funArray = Vector{ZhegFun{T}}([])
        for i in 1 : n
            fi = permute(fam.fs[i], perm)
            push!(funArray, fi)
        end
        push!(resOrbit, Family(funArray[perm]...))
    end
    return resOrbit
end


function stab(fam :: Family{T}) where T
    n = length(fam)
    perms = permutations(collect(1 : n))
    stabilizers = Vector{Vector{Int}}([])
    for perm in perms
        funArray = Vector{ZhegFun{T}}([])
        for i in 1 : n
            fi = permute(fam.fs[i], perm)
            push!(funArray, fi)
        end
        fam_new = Family(funArray[perm]...)
        if fam_new == fam 
            push!(stabilizers, perm)
        end
    end
    return stabilizers
end


function as_vec!(funVec, numb :: T) where T <: Integer
    n = length(funVec)
    for i in n : -1 : 1
        funVec[i] = numb & 0x1
        numb >>= 1
    end
    return nothing
end



function shift_orb(fam :: Family{T}) where T
    n = length(fam)
    shiftOrb = Vector{Family{T}}()
    currFam = deepcopy(fam)
    shiftVec = zeros(Int, n)
    for shift in 0 : 2^n - 1
        as_vec!(shiftVec, shift)
        for j in 1 : n
            if (shiftVec[j] == 1) && !(Monom(0) in currFam[j].ANF)
                push!(currFam[j].ANF, Monom(0))
                sort!(currFam[j].ANF, by = (x -> x.val))
            elseif (shiftVec[j] == 0) && (Monom(0) in currFam[j].ANF)
                filter!(x -> x != Monom(0), currFam[j].ANF)
                sort!(currFam[j].ANF, by = (x -> x.val))
            end
        end
        push!(shiftOrb, deepcopy(currFam))
    end
    return shiftOrb
end





include("pairing.jl")
# derive permutation x -> x + y + fam o pi(x, y)
# as array: perm[i] = pi(i)
# fixed point : perm[i] == i
function as_perm(fam :: Family, pi :: Pairing, y)
    n = size(fam)
    all_elems = Set(collect(0 : 2^n - 1))
    permutation = Vector{Array{Int}}([])
    while !isempty(all_elems)
        elstart = minimum(all_elems)
        push!(permutation, [elstart])
        el = xor(xor(elstart, y), fam(pi(elstart, y)))
        while el != elstart
            push!(permutation[end], el)
            el = xor(xor(el, y), fam(pi(el, y)))
        end
        all_elems = setdiff(all_elems, permutation[end])
    end
    return permutation
end


# returns mapping (i-1 -> fam(i-1))
function cyclic(fam :: Family)
    n = length(fam)
    mapping = zeros(Int, 2^n)
    for i in 1 : 2^n
        mapping[i] = fam(i - 1)
    end
    return mapping
end


function visual(map)
    n = length(map)
    for i in 1 : n
        println(i - 1, " -> ", map[i])
    end
end


