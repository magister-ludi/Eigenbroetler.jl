
"""
    alt_xy(r::Integer, c::Integer)
    alt_xy(rc::Tuple{T, U}) where {T <: Integer, U <: Integer}

Construct a matrix with `r` rows and `c` columns, with elements
alternating between `+1.0` and `-1.0`.
The first element is +1.
"""
function alt_xy(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], outer = div(r + 1, 2))
    isodd(r) && (a = a[1:(end - 1)])
    b = repeat(vcat(a, -a), outer = div(c + 1, 2))
    isodd(c) && (b = b[1:(r * c)])
    reshape(b, r, c)
end

alt_xy(rc::Tuple{T, U}) where {T <: Integer, U <: Integer} = alt_xy(rc[1], rc[2])

"""
    alt_x(r::Integer, c::Integer)
    alt_x(rc::Tuple{T, U}) where {T <: Integer, U <: Integer}

Construct a matrix with `r` rows and `c` columns, with columns
alternating between `+1.0` and `-1.0`.
The first column is +1.
"""
function alt_x(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], inner = r)
    b = repeat(a, outer = div(c + 1, 2))
    isodd(c) && (b = b[1:(r * c)])
    reshape(b, r, c)
end

alt_x(rc::Tuple{T, U}) where {T <: Integer, U <: Integer} = alt_x(rc[1], rc[2])

"""
    alt_y(r::Integer, c::Integer)
    alt_y(rc::Tuple{T, U}) where {T <: Integer, U <: Integer}

Construct a matrix with `r` rows and `c` columns, with rows
alternating between `+1.0` and `-1.0`.
The first row is +1.
"""
function alt_y(r::Integer, c::Integer)
    a = repeat([1.0, -1.0], outer = div(r + 1, 2))
    isodd(r) && (a = a[1:(end - 1)])
    b = repeat(a, outer = c)
    reshape(b, r, c)
end

alt_y(rc::Tuple{T, U}) where {T <: Integer, U <: Integer} = alt_y(rc[1], rc[2])

"""
    fft(eb[, recentre = true])

Returns the two-dimensional DFT of Eigenbrot `eb`. If
`recentre` is `true` the DFT is stored so that the value
at the coordinate origin is at the central location of
the data.
"""
function FFTW.fft(eb::Eigenbrot, recentre::Bool = true)
    w = width(eb)
    h = height(eb)
    scale = 1.0 / sqrt(h * w)
    if recentre
        alt = alt_xy(h, w)
        data = alt .* eb.vals
        eb.fft ? fft!(data) : bfft!(data)
        Eigenbrot(scale * (alt .* data), !eb.fft)
    else
        data = copy(eb.vals)
        eb.fft ? fft!(data) : bfft!(data)
        Eigenbrot(scale * data, !eb.fft)
    end
end

"""
    fftx(eb[, recentre = true])

Returns the one-dimensional DFT of Eigenbrot `eb` in
the `x`-direction. If `recentre` is `true` the DFT is
stored so that the value at the coordinate origin is
at the central location of the data.
"""
function fftx!(eb::Eigenbrot, recentre::Bool = true)
    if recentre
        w = width(eb)
        h = height(eb)
        scale = 1.0 / sqrt(h)
        alt = alt_y(h, w)
        eb.vals .*= alt
        eb.fft ? fft!(eb.vals, 2) : bfft!(eb.vals, 2)
        eb.vals .= scale * (alt .* eb.vals)
    else
        eb.fft ? fft!(eb.vals, 2) : bfft!(eb.vals, 2)
    end
    eb.fft = !eb.fft
    return eb
end

fftx(eb::Eigenbrot, recentre::Bool = true) = fftx!(copy(eb), recentre)

"""
    ffty(eb[, recentre = true])

Returns the one-dimensional DFT of Eigenbrot `eb` in
the `y`-direction. If `recentre` is `true` the DFT is
stored so that the value at the coordinate origin is
at the central location of the data.
"""
function ffty!(eb::Eigenbrot, recentre::Bool = true)
    if recentre
        w = width(eb)
        h = height(eb)
        scale = 1.0 / sqrt(h)
        alt = alt_y(h, w)
        eb.vals .*= alt
        eb.fft ? fft!(eb.vals, 1) : bfft!(eb.vals, 1)
        eb.vals .= scale * (alt .* eb.vals)
    else
        eb.fft ? fft!(eb.vals, 1) : bfft!(eb.vals, 1)
    end
    eb.fft = !eb.fft
    return eb
end

ffty(eb::Eigenbrot, recentre::Bool = true) = ffty!(copy(eb), recentre)

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
