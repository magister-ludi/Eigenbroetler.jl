
"""
    alt_xy(r::Integer, c::Integer)
    alt_xy{T <: Integer, U <: Integer}(rc::Tuple{T, U})
Construct a matrix with `r` rows and `c` columns, with elements
alternating between `+1.0` and `-1.0`.
The first element is +1.
"""
function alt_xy(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], outer = div(r + 1, 2))
    isodd(r) && (a = a[1:end - 1])
    b = repeat(vcat(a, -a), outer = div(c + 1, 2))
    isodd(c) && (b = b[1:r * c])
    reshape(b, r, c)
end

alt_xy{T <: Integer, U <: Integer}(rc::Tuple{T, U}) = alt_xy(rc[1], rc[2])

"""
    alt_x(r::Integer, c::Integer)
    alt_x{T <: Integer, U <: Integer}(rc::Tuple{T, U})
Construct a matrix with `r` rows and `c` columns, with columns
alternating between `+1.0` and `-1.0`.
The first column is +1.
"""
function alt_x(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], inner = r)
    b = repeat(a, outer = div(c + 1, 2))
    isodd(c) && (b = b[1:r * c])
    reshape(b, r, c)
end

alt_x{T <: Integer, U <: Integer}(rc::Tuple{T, U}) = alt_x(rc[1], rc[2])

"""
    alt_y(r::Integer, c::Integer)
    alt_y{T <: Integer, U <: Integer}(rc::Tuple{T, U})
Construct a matrix with `r` rows and `c` columns, with rows
alternating between `+1.0` and `-1.0`.
The first row is +1.
"""
function alt_y(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], outer = div(r + 1, 2))
    isodd(r) && (a = a[1:end-1])
    b = repeat(a, outer = c)
    reshape(b, r, c)
end

alt_y{T <: Integer, U <: Integer}(rc::Tuple{T, U}) = alt_y(rc[1], rc[2])

"""
    fft(eb[, recentre = true])
Returns the two-dimensional DFT of Eigenbrot `eb`. If
`recentre` is `true` the DFT is stored so that the value
at the coordinate origin is at the central location of
the data.
"""
function fft(eb::Eigenbrot, recentre::Bool = true)
    w = width(eb)
    h = height(eb)
    scale = 1.0 / sqrt(h * w)
    alt = alt_xy(h, w)
    data = recentre ? alt .* eb.vals : copy(eb.vals)
    eb.fft ? fft!(data) : bfft!(data)
    trf = Eigenbrot(scale * (recentre ? alt .* data : data), !eb.fft)
    return trf
end

"""
    fftx(eb[, recentre = true])
Returns the one-dimensional DFT of Eigenbrot `eb` in
the `x`-direction. If `recentre` is `true` the DFT is
stored so that the value at the coordinate origin is
at the central location of the data.
"""
function fftx(eb::Eigenbrot, recentre::Bool = true)
    w = width(eb)
    h = height(eb)
    scale = 1.0 / sqrt(w)
    alt = alt_x(h, w)
    data = recentre ? alt .* eb.vals : copy(eb.vals)
    eb.fft ? fft!(data, 2) : bfft!(data, 2)
    return Eigenbrot(scale * (recentre ? alt .* data : data), !eb.fft)
end

"""
    ffty(eb[, recentre = true])
Returns the one-dimensional DFT of Eigenbrot `eb` in
the `y`-direction. If `recentre` is `true` the DFT is
stored so that the value at the coordinate origin is
at the central location of the data.
"""
function ffty(eb::Eigenbrot, recentre::Bool = true)
    w = width(eb)
    h = height(eb)
    scale = 1.0 / sqrt(h)
    alt = alt_y(h, w)
    data = recentre ? alt .* eb.vals : copy(eb.vals)
    eb.fft ? fft!(data, 1) : bfft!(data, 1)
    return Eigenbrot(scale * (recentre ? alt .* data : data), !eb.fft)
end

"""
    removeDC(eb)
Return an Eigenbrot that is the Eigenbrot `eb`
with the zero-frequency component removed.
"""
function removeDC(eb::Eigenbrot)
  r = fft(eb, false)
  r[1, 1] = 0
  return fft(r, false)
end
