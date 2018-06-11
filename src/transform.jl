
"""
    hilbertX(eb::Eigenbrot)
Return the Hilbert x-transform of `eb`.
"""
function hilbertX(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c in 1:w
        for r in 1:h
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
    hilbertY(eb::Eigenbrot)
Return the Hilbert x-transform of `eb`.
"""
function hilbertY(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c in 1:w
        for r in 1:h
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
    for row in 1:h
        δ = xShear * (halfHeight - row + 1)
        for i in 1:lim
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
    for col in 1:w
        δ = yShear * (halfWidth - col + 1)
        for i in 1:lim
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
