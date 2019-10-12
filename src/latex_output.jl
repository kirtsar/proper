function latexify(m :: Monom)
    if m.val == 0
        res = "1"
    else
        res = ""
        inds = nzBits(m.val)
        for i in 1 : length(inds) - 1
            res *= "x_$(inds[i]) \\cdot "
        end
        res *= "x_$(inds[end])"
    end
    return res
end


function latexify(f :: ZhegFun)
    if !isempty(f.ANF)
        res = reduce(*, map(x -> x * " \\oplus ", latexify.(f.ANF)[1:end-1]))
        res *= latexify(f.ANF[end])
    else
        res = "0"
    end
    return res
end


function latexify(fam :: Family)
    return latexify.(fam.fs)
end