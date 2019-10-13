struct LatinSquare
    sq :: Matrix{Int}
end


function latin(fam :: Family, pair :: Pairing)
    n = length(fam)
    m = 2^n - 1
    table = zeros(Int, (m + 1, m + 1))
    ls = LatinSquare(table)
    latin!(ls, fam, pair)
    return LatinSquare(table)
end


# in-place constructor of latin square
function latin!(LS :: LatinSquare, fam :: Family, pair :: Pairing)
    n = length(fam)
    m = 2^n - 1
    # table = zeros(Int, (m + 1, m + 1))
    for i in 0 : m
        for j in 0 : m
            # (i, j) -> i + j + f o pi (i, j)
            LS.sq[i + 1, j + 1] = xor(xor(i, j), fam(pair(i, j)))
        end
    end
    # return LatinSquare(table)
end


function (ls :: LatinSquare)(i :: Int, j :: Int)
    return ls.sq[i + 1, j + 1]
end


import Base.length
function length(ls :: LatinSquare)
    return size(ls.sq)[1]
end


function assoc(ls :: LatinSquare)
    m = length(ls)
    assocNum = 0
    for i in 0 : m - 1
        for j in 0 : m - 1
            for k in 0 : m - 1
                # (ij)k == i(jk) ?
                if ls(ls(i, j), k) == ls(i, ls(j, k))
                    assocNum += 1
                end
            end
        end
    end
    return assocNum
end


function triples(ls :: LatinSquare)
    m = length(ls)
    assoc_triples = Vector{Tuple{Int, Int, Int}}([])
    for i in 0 : m - 1
        for j in 0 : m - 1
            for k in 0 : m - 1
                # (ij)k == i(jk) ?
                if ls(ls(i, j), k) == ls(i, ls(j, k))
                    #assocNum += 1
                    push!(assoc_triples, (i, j, k))
                end
            end
        end
    end
    return assoc_triples
end




function test_LS()
    zeroZh = Int[]
    fam = Family(ZhegFun.([zeroZh, zeroZh, zeroZh]))
    # does not depend on pairing
    pair = Pairing(Family(ZhegFun.([[1, 2], [1, 2], [1, 2]])))
    LS = latin(fam, pair)
    return LS.sq
end



    