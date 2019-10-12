include("show_graphs.jl")

function bits(val; len = 0)
    if len == 0
        len = Int(floor(log2(val))) + 1  
    end  
    res = zeros(Int, len)
    for i in 1 : len
        #res[len - i + 1] = val & 0x1
        res[i] = val & 0x1
        val >>= 1
    end
    return res
end

function from_bits(bitarr :: Vector{T}) where T
    n = length(bitarr)
    res = 0
    for i in 1 : n
        res <<= 1
        # res |= bitarr[i]
        res |= bitarr[n - i + 1]
    end
    return res
end


# asynchronous state graph
function ASG(fam :: Family{T}) where T
    n = length(fam)
    verts = collect(0 : 2^n - 1)
    adj_table = zeros(Int, (2^n, 2^n))
    for v in verts
        bitrepr = bits(v, len = n)
        for i in 1 : n
            if bitrepr[i] != fam[i](v)
                adj_repr = bits(v, len = n)
                adj_repr[i] = fam[i](v)
                adj_v = from_bits(adj_repr)
                adj_table[v + 1, adj_v + 1] = 1
            end
        end
    end
    return adj_table
end



# state graph
function SG(fam :: Family{T}) where T
    n = length(fam)
    verts = collect(0 : 2^n - 1)
    adj_table = zeros(Int, (2^n, 2^n))
    for v in verts
        bitrepr = bits(v; len = n)
        for i in 1 : n
            adj_repr = bits(v, len = n)
            adj_repr[i] = 1 - bitrepr[i] #fam[i](v)
            adj_v = from_bits(adj_repr)
            if bitrepr[i] == fam[i](v)
                adj_table[v + 1, adj_v + 1] = 1
            else #if bitrepr[i] > fam[i](v)
                adj_table[adj_v + 1, v + 1] = 1
            end
        end
    end
    return adj_table
end


function fx_graph(fam :: Family{T}) where T
    n = length(fam)
    verts = collect(0 : 2^n - 1)
    adj_table = zeros(Int, (2^n, 2^n))
    for v in verts
        vnext = fam(v)
        adj_table[v + 1, vnext + 1] = 1
    end
    return adj_table
end


function gshow(fam :: Family{T}) where T
    gshow(fx_graph(fam))
end


function fx_preimages(fam :: Family{T}) where T
    adjmat = fx_graph(fam)
    indegs = sum(adjmat, dims = 1)
    res = 0
    for v in indegs
        res += binomial(v, 2)
    end
    return res
end


function xfx_preimages(fam :: Family{T}) where T
    adjmat = xfx_graph(fam)
    indegs = sum(adjmat, dims = 1)
    res = 0
    for v in indegs
        res += binomial(v, 2)
    end
    return res
end


function xfx_graph(fam :: Family{T}) where T
    n = length(fam)
    verts = collect(0 : 2^n - 1)
    adj_table = zeros(Int, (2^n, 2^n))
    for v in verts
        vnext = xor(fam(v), v)
        adj_table[v + 1, vnext + 1] = 1
    end
    return adj_table
end


# given a family of functions
# return the matrix of essentials
function deps_graph(fam :: Family{T}) where T
    n = length(fam)
    graphMatrix = zeros(Int8, n, n)
    for i in 1 : n
        f = fam[i]
        for numBit in get_essentials(f)
            graphMatrix[numBit, i] = 1
        end
    end
    return graphMatrix
end


function deps_graph(pr :: Pairing{T}) where T
    return deps_graph(pr.funs)
end


function is_full(fam :: Family{T}) where T
    n = length(fam)
    fullmat = ones(T, (n, n))
    for i in 1 : n
        fullmat[i, i] = 0
    end
    fam_mat = deps_graph(fam)
    return fam_mat == fullmat
end

