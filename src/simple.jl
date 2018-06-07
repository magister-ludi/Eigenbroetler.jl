
for name in (:real, :imag, :angle, :abs, :conj)
    @eval begin
        import Base.$name
        $name(eb::Eigenbrot) = Eigenbrot(Array{Complex128}($name(eb.vals)))
    end
end

phase(eb::Eigenbrot) = angle(eb)

import Base: (+), (-), (*), (/), (^)

(+)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .+ eb2.vals)
(+)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals + n)
(+)(n::Number, eb1::Eigenbrot) = eb1 + n
(-)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .- eb2.vals)
(-)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals .- n)
(-)(n::Number, eb1::Eigenbrot) = Eigenbrot(n .- eb1.vals)

(*)(eb::Eigenbrot, x::Number) = Eigenbrot(x * eb)
(*)(x::Number, eb::Eigenbrot) = Eigenbrot(x .* eb.vals)
(/)(eb::Eigenbrot, x::Number) = Eigenbrot(eb.vals ./ x)

Base.broadcast(::typeof(*), eb1::Eigenbrot, eb2::Eigenbrot) =
    Eigenbrot(eb1.vals .* eb2.vals)

Base.broadcast(::typeof(/), eb1::Eigenbrot, eb2::Eigenbrot) =
    Eigenbrot(eb1.vals ./ eb2.vals)

#=
Utility method used by pad() methods, not exported.
=#
function _pad(v::Matrix{Complex128}, value::Complex,
               left::Integer, top::Integer, right::Integer, bottom::Integer)
    rows, cols = size(v)
    cols_out = cols + left + right
    rows_out = rows + top + bottom
    (cols_out ≤ 0 || rows_out ≤ 0) &&
        error("Invalid resize: $(cols)x$rows -> $(cols_out)x$rows_out")

    col_dest1 = left > 0 ? 1 + left : 1
    col_src1 = left > 0 ? 1 : 1 - left
    col_src2 = right > 0 ? cols : cols + right
    col_dest2 = col_dest1 + col_src2 - col_src1

    row_dest1 = bottom > 0 ? 1 + bottom : 1
    row_src1 = bottom > 0 ? 1 : 1 - bottom
    row_src2 = top > 0 ? rows : rows + top
    row_dest2 = row_dest1 + row_src2 - row_src1

    v2 = fill(complex(value), rows_out, cols_out)
    v2[row_dest1:row_dest2, col_dest1:col_dest2] = v[row_src1:row_src2, col_src1:col_src2]
    return v2
end

"""
    pad(eb::Eigenbrot, margin::Integer)
Pad edges of `eb` by `margin` pixels. If `margin` is negative,
extend by that amount, filling the new pixels with zero.

    pad(eb::Eigenbrot, value, margin::Integer)
Pad edges of `eb` by `margin` pixels. If `margin` is negative,
extend by that amount, filling the new pixels with `value`.

    pad(eb::Eigenbrot, value = 0; left::Integer = 0, top::Integer = 0, right::Integer = 0, bottom::Integer = 0)
Pad the named edges of `eb` by the number of pixels given. Negative
edge values will extend by the relevant amount, filling the new pixels
with `value`.
"""
pad(eb::Eigenbrot, value::Number = 0.0im; left::Integer = 0, top::Integer = 0,
        right::Integer = 0, bottom::Integer = 0) =
    Eigenbrot(_pad(eb.vals, Complex128(value), left, top, right, bottom), eb.fft)

pad(eb::Eigenbrot, value::Real, margin::Integer) =
    Eigenbrot(_pad(eb.vals, Complex128(value), margin, margin, margin, margin), eb.fft)

pad(eb::Eigenbrot, value::Complex, margin::Integer) =
    Eigenbrot(_pad(eb.vals, value, margin, margin, margin, margin), eb.fft)

pad(eb::Eigenbrot, margin::Integer) =
   pad(eb, 0.0im, margin)

"""
    pow2pad(eb::Eigenbrot)
Create a new `Eigenbrot` which has a copy of the data
padded by zeros so that its dimensions are powers of 2
(see also `pow2pad!`).
"""
function pow2pad(eb::Eigenbrot)
    rows, cols = size(eb.vals)
    rows2 = nextpow2(rows)
    cols2 = nextpow2(cols)
    eb2 = if rows == rows2 && cols == cols2
        copy(eb)
    else
        pad(eb, top = rows2 - rows, right = cols2 - cols)
    end
    return eb2
end

"""
    pow2pad!(eb::Eigenbrot)
Replace the data of `eb` with the original data
padded by zeros so that its dimensions are powers of 2
(see also `pow2pad`).
"""
function pow2pad!(eb::Eigenbrot)
    eb2 = pow2pad(eb)
    eb.vals = eb2.vals
    return eb
end

function flipver(eb::Eigenbrot)
    flip = copy(eb)
    h = size(flip, 1)
    for r in 1:h
        reverse!(@view flip.vals[r, :])
    end
    return flip
end

function fliphor(eb::Eigenbrot)
    flip = copy(eb)
    w = size(flip, 2)
    for c in 1:w
        reverse!(@view flip.vals[:, c])
    end
    return flip
end

function swapxy(eb::Eigenbrot, xflip = false)
    h, w = size(eb)
    swap = Eigenbrot(w, h)
    if xflip
        for c in 1:h
            swap.vals[:, c] = @view eb.vals[c, :]
        end
    else
        for c in 1:h
            swap.vals[w:-1:1, h-c+1] = @view eb.vals[c, :]
        end
    end
    return swap
end
