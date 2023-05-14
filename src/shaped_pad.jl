
#=
Utility method used by extend() methods, not exported.
=#
function _extend(
    v::Matrix{ComplexF64},
    Λ,
    left::Integer,
    top::Integer,
    right::Integer,
    bottom::Integer,
)
    @assert abs(Λ) ≤ 1
    @assert 0 ≤ left
    @assert 0 ≤ right
    @assert 0 ≤ top
    @assert 0 ≤ bottom
    rows, cols = size(v)
    cols_out = cols + left + right
    rows_out = rows + top + bottom
    v2 = Matrix{ComplexF64}(undef, rows_out, cols_out)
    v2[bottom:(rows + bottom - 1), left:(cols + left - 1)] .= v
    for c = left:(cols + left - 1)
        for r = (bottom - 1):-1:1
            v2[r, c] = Λ * v2[r + 1, c]
        end
        for r = (rows + bottom):rows_out
            v2[r, c] = v2[r - 1, c]
        end
    end
    for c = (left - 1):-1:1
        @views v2[:, c] = Λ * v2[:, c + 1]
    end
    for c = (cols + left):cols_out
        @views v2[:, c] = Λ * v2[:, c - 1]
    end
    return v2
end

"""
    extend(eb::Eigenbrot, Λ::Number, margin::Integer)

Extend edges of `eb` by `margin` pixels. Extension
multiplies original edge pixels successively by the value `Λ` which
must satisfy `abs(Λ) ≤ 1`.
"""
extend(eb::Eigenbrot, Λ::Number, margin::Integer)

"""
    extend(eb::Eigenbrot, Λ::Number; left::Integer = 0, top::Integer = 0, right::Integer = 0, bottom::Integer = 0)

Extend given edges of `eb` by the number of pixels indicated. Extension
multiplies original edge pixels successively by the value `Λ` which
must satisfy `abs(Λ) ≤ 1`.
"""
extend(
    eb::Eigenbrot,
    Λ::Number;
    left::Integer = 0,
    top::Integer = 0,
    right::Integer = 0,
    bottom::Integer = 0,
)

extend(
    eb::Eigenbrot,
    Λ::Number;
    left::Integer = 0,
    top::Integer = 0,
    right::Integer = 0,
    bottom::Integer = 0,
) = Eigenbrot(_extend(eb.vals, Λ, left, top, right, bottom), eb.fft)

extend(eb::Eigenbrot, Λ::Number, margin::Integer) =
    Eigenbrot(_extend(eb.vals, Λ, margin, margin, margin, margin), eb.fft)
