function shift_vec!(shiftVec, shift)
    n = length(shiftVec)
    for i in n : -1 : 1
        shiftVec[i] = shift & 0x1
        shift >>= 1
    end
    return nothing
end


function as_vec!(funVec, f :: ZhegFun{T}) where T
    n = length(funVec)
    for i in 0 : (n - 1)
        funVec[i + 1] = f(i)
    end
    return nothing
end

function not!(funVec)
    for i in 1 : length(funVec)
        funVec[i] = 1 - funVec[i]
    end
    return nothing
end


function is_ortho(fam :: Family{T}) where T
    shiftVec = zeros(Int, length(fam))
    return is_ortho(fam, shiftVec)
end

function is_ortho(fam :: Family{T}, shiftVec) where T
    n = length(fam)
    fVec = zeros(T, 2^n)
    gVec = zeros(T, 2^n)
    res = true
    for i in 1 : n
        as_vec!(fVec, fam[i])
        if shiftVec[i] == 1
            not!(fVec)
        end
        for j in (i + 1) : n
            as_vec!(gVec, fam[j])
            if shiftVec[j] == 1
                not!(gVec)
            end
            res &= is_ortho(fVec, gVec)
        end
    end
    return res
end



function is_ortho(fVec, gVec)
    res = 0
    for i in 1 : length(fVec)
        res += fVec[i] * gVec[i]
    end
    return (res == 0)
end


function has_ortho(fam :: Family{T}) where T
    n = length(fam)
    shiftVec = zeros(Int, n)
    for i in 0 : 2^n - 1
        shift_vec!(shiftVec, i)
        if is_ortho(fam, shiftVec)
            return true
        end
    end
    return false
end

function monom(arr)
    res = 0
    for x in arr
        res += 2^(x-1)
    end
    return Monom(res)
end


function generate_ortho_family(n)
    mon = vcat(collect(2 : n), collect(1 : n))
    funs = Vector{ZhegFun}([])
    for i in 1 : n
        m1 = monom(mon[i : i+(n-2)])
        m2 = monom(mon[i+1 : i+(n-2)])
        push!(funs, ZhegFun([m1, m2]))
    end
    return Family(funs...)
end