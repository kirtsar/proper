BitType = Union{UInt8, UInt16, UInt32, Int}
using AutoHashEquals  # checking for equality

################################
# structures
# monoms in zhegalkin polys
struct Monom{T}
    val :: T
end


# zhegalkin poly 
@auto_hash_equals mutable struct ZhegFun{T}
    ANF :: Vector{Monom{T}}
end


# zhegalkin poly from array
function ZhegFun(arr :: Vector{T}) where T <: BitType
    return ZhegFun(Monom.(arr))
end


function ZhegFun(monArr :: Vector{Vector{Int}})
    mons = zeros(Int, length(monArr))
    for (i, monRepr) in enumerate(monArr)
        val = 0
        for bit in monRepr
            val += 2^(bit - 1)
        end
        mons[i] = val
    end
    sort!(mons)
    return ZhegFun(mons)
end



################################
# APPLICATION OF DIFFERENT TYPES
# apply monom
function (m :: Monom)(x)
    return (m.val & x == m.val)
end


# apply zhegalkin poly

function (f :: ZhegFun)(x :: T) where T <: BitType
    res = 0
    for m in f.ANF
        res = xor(res, m(x))
    end
    return res
end





################################
# nice printing
# helper 
# find non-zero bits in number
function nzBits(number)
    res = []
    i = 1
    while number > 0
        bit = number & 0x1
        if bit != 0
            append!(res, i)
        end
        number >>= 1
        i += 1
    end
    return res
end


# nice printing of monoms
function Base.show(io::IO, m::Monom)
    if m.val == 0
        print(io, "1")
    else
        for u in nzBits(m.val)
            print(io, "x$u")
        end
    end
end


# nice printing of zhegalkin polys
function Base.show(io::IO, f::ZhegFun)
    if isempty(f.ANF)
        print(io, "0")
    else
        for i in 1 : length(f.ANF) - 1
            mon = f.ANF[i]
            print(io, mon)
            print(io, " + ")
        end
        print(io, f.ANF[end])
    end
end


# find all essential variables of Zhegalkin Poly
function get_essentials(f :: ZhegFun)
    essentials = 0
    for coef in f.ANF
        essentials |= coef.val
    end
    return nzBits(essentials)
end



function to_bits!(m :: Monom, bits)
    curr = m.val
    i = 1
    for i in 1 : length(bits)
        bits[i] = curr & 0x1
        curr >>= 1
    end
    return nothing
end

function from_b(bits)
    res = 0
    for i in length(bits) : -1 : 1
        res <<= 1
        res |= bits[i]
    end
    return Monom(res)
end


function permute(m :: Monom, perm)
    n = length(perm)
    bits = zeros(UInt8, n)
    to_bits!(m, bits)
    return from_b(bits[perm])
end






