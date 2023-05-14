
"""
    hilbert_x(eb::Eigenbrot)

Return the Hilbert x-transform of `eb`.
"""
function hilbert_x(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c = 1:w
        for r = 1:h
            sign = if c == centre_c
                r > centre_r ? im : -im
            elseif r == centre_r
                c > centre_c ? im : -im
            else
                c > centre_c ? im : -im
            end
            nd[r, c] *= sign
        end
    end
    return nd
end

"""
    hilbert_y(eb::Eigenbrot)

Return the Hilbert x-transform of `eb`.
"""
function hilbert_y(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c = 1:w
        for r = 1:h
            sign = if c == centre_c
                r > centre_r ? im : -im
            elseif r == centre_r
                c > centre_c ? im : -im
            else
                r > centre_r ? im : -im
            end
            nd[r, c] *= sign
        end
    end
    return nd
end

"""
    fft_xshear(eb::Eigenbrot, xShear::Real, expand::Bool = true)

Shear `eb` by a factor 'xShear' in the x-direction using the Fourier
shift theorem, and return the result. If `expand` is true the input
is padded so that the shear doesn't wrap the data as a result of
aliasing in the Fourier domain.
"""
function fft_xshear(eb::Eigenbrot, xShear::Real, expand::Bool = true)
    xform = if expand
        extra = ceil(Int, abs(xShear) * height(eb) / 2)
        fftx(pad(eb, left = extra, right = extra), false)
    else
        fftx(eb, false)
    end
    h, w = size(xform)
    halfHeight = h / 2
    lim = iseven(w) ? w >> 1 : (w + 1) >> 1
    for row = 1:h
        δ = xShear * (halfHeight - row + 1)
        for i = 1:lim
            θ = (2.0 * π * δ * i) / w
            cs = cos(θ) + im * sin(θ)
            xform[row, i + 1] *= cs
            xform[row, w - i + 1] *= cs'
        end
        if iseven(w)
            i2 = trunc(Int, abs(δ))
            δ < 0 && (i2 = -i2)
            θ = (2.0 * π * δ * i2) / w
            xform[row, lim] *= cos(θ) + im * sin(θ)
        end
    end

    return fftx(xform, false)
end

"""
    fft_yshear(eb::Eigenbrot, xShear::Real, expand::Bool = true)

Shear `eb` by a factor 'yShear' in the y-direction using the Fourier
shift theorem, and return the result. If `expand` is true the input
is padded so that the shear doesn't wrap the data as a result of
aliasing in the Fourier domain.
"""
function fft_yshear(eb::Eigenbrot, yShear::Real, expand::Bool = true)
    xform = if expand
        extra = ceil(Int, abs(yShear * width(eb) / 2))
        ffty(pad(eb, top = extra, bottom = extra), false)
    else
        ffty(eb, false)
    end
    h, w = size(xform)
    halfWidth = w / 2
    lim = iseven(h) ? h >> 1 : (h + 1) >> 1
    for col = 1:w
        δ = yShear * (halfWidth - col + 1)
        for i = 1:lim
            θ = (2.0 * π * δ * i) / h
            cs = cos(θ) + im * sin(θ)
            xform[i + 1, col] *= cs
            xform[h - i + 1, col] *= cs'
        end
        if iseven(h)
            i2 = trunc(Int, abs(δ))
            δ < 0.0 && (i2 = -i2)
            θ = (2.0 * π * i2 * lim) / h
            xform[lim + 1, col] *= cos(θ) + im * sin(θ)
        end
    end

    return ffty(xform, false)
end

"""
    rotate90(eb::Eigenbrot, n::Integer)

Return the result of rotating `eb` through an angle of
`n` * π /2 (i.e. `n` right angles).
"""
function rotate90(eb::Eigenbrot, n::Integer)
    if n == 0
        return copy(eb)
    elseif n == 1
        h, w = size(eb)
        e2 = Eigenbrot(w, h, eb.fft)
        for r = 1:w
            e2.vals[r, 1:h] = eb.vals[h:-1:1, r]
        end
        return e2
    elseif n == 2
        e2 = Eigenbrot(height(eb), width(eb), eb.fft)
        l = length(eb.vals)
        e2.vals[1:l] = eb.vals[l:-1:1]
        return e2
    elseif n == 3
        h, w = size(eb)
        e2 = Eigenbrot(w, h, eb.fft)
        for r = 1:w
            e2.vals[w - r + 1, h:-1:1] = eb.vals[h:-1:1, r]
        end
        return e2
    else
        return rotate90(eb, mod(n, 4))
    end
end

"""
    fourierRotation(eb::Eigenbrot, θ::Real, expand::Bool = true)

Return the result of rotating `eb` through an angle `θ`,
using the algorithm described by Larkin, Oldfield and Klemm
(Optics Communications 139 (1997) 99-106). If `expand` is true,
the size of the output image is increased to contain the rotated image,
otherwise the aliasing described in Larkin et al. occurs.
"""
function fourierRotation(eb::Eigenbrot, θ::Real, expand::Bool = true)
    h, w = size(eb)
    θ = -θ  # reverse rotation direction

    # Deconstruct θ as the sum of an integral number, n90,
    # of right angles, and a remainder, ϕ, with absolute
    # value less than π/4...
    n90, ϕ = divrem(θ, π / 2)
    if abs(ϕ) > π / 4
        s = sign(ϕ)
        n90 += s
        ϕ = ϕ - s * π / 2
    end
    n90 = trunc(Int, mod(n90, 4))
    result = rotate90(eb, n90)
    ϕ ≈ 0 && return result

    xshear = tan(ϕ / 2)
    yshear = -sin(ϕ)
    if expand
        hh, ww = size(result)
        hm, wm = 0, 0
        msx = [1.0 xshear; 0.0 1.0]
        msy = [1.0 0.0; yshear 1.0]
        cnrs = [0 ww 0 ww; 0 0 hh hh]
        for mm in [msx, msy, msx]
            cnrs = mm * cnrs
            wm = max(wm, maximum(cnrs[1, :]) - minimum(cnrs[1, :]))
            hm = max(hm, maximum(cnrs[2, :]) - minimum(cnrs[2, :]))
        end
        lr = ceil(Int, (wm - ww) / 2)
        tb = ceil(Int, (hm - hh) / 2)
        result = pad(result, left = lr, right = lr, top = tb, bottom = tb)
    end

    result = fft_xshear(result, xshear, false)
    result = fft_yshear(result, yshear, false)
    result = fft_xshear(result, xshear, false)
    return result
end
