

function test2(fam :: Family)
    tmp = []
    orbFam = orb(fam)
    res = zeros(length(orbFam))
    for (i, ff) in enumerate(orbFam)
        total_assoc = assoc(ff)
        res[i] = sum(total_assoc) / length(total_assoc)
        push!(tmp, countmap(total_assoc))
    end 
    return res, tmp
end

function test3(fams :: Vector{Family{T}}) where T
    min_assoc = zeros(Int, length(fams))
    mean_assoc = zeros(length(fams))
    for (i, fam) in enumerate(fams)
        total_assoc = assoc(fam)
        min_assoc[i] = minimum(total_assoc)
        mean_assoc[i] = sum(total_assoc) / length(total_assoc)
    end
    return min_assoc, mean_assoc
end   


function test5(assoc4)
    res = [[] for i in 1 : 56]
    for i in 1 : 56
        res[i] = min_pairings(assoc4[1][i])
    end
    return res
end

function test6(assoc3)
    res = [[] for i in 1 : 2]
    for i in 1 : 2
        res[i] = min_pairings(assoc3[1][i])
    end
    return res
end


n = length(prs3)
for i in 1 : n
    prs = prs3[i]
    for pr in prs
        print(switch(pr) in prs)
        print()
    end
    println()
end

Family(
[[2,3], [2,4], [2,5], [3,4], [3,5], [4,5]],
[[1], [1,3], [1,4], [1,5], [3,4], [3,5], [4,5]],
[[1], [2], [1,2], [1,4], [1,5], [2,4], [2,5], [4,5]],
[[1], [2], [3], [1,2], [1,3], [1,5], [2,3], [2,5], [3,5]],
[[1], [2], [3], [4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
)


function generate_ortho_family3()
    fun1 = ZhegFun([6, 4])
    fun2 = ZhegFun([5, 1])
    fun3 = ZhegFun([3, 2])
    return Family([fun1, fun2, fun3])
end


function generate_ortho_family4()
    fun1 = ZhegFun([12, 14])
    fun2 = ZhegFun([ 9, 13])
    fun3 = ZhegFun([ 3, 11])
    fun4 = ZhegFun([ 6,  7])
    return Family([fun1, fun2, fun3, fun4])
end


function monom(arr)
    res = 0
    for x in arr
        res += 2^(x-1)
    end
    return Monom(res)
end


function generate_ortho_family5()
    mon5 = [2,3,4,5,1,2,3,4,5]
    funs = Vector{ZhegFun}([])
    for i in 1 : 5
        m1 = monom(mon5[i:i+3])
        m2 = monom(mon5[i+1:i+3])
        push!(funs, ZhegFun([m1, m2]))
    end
    return Family(funs...)
end