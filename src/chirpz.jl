struct Placeholder{T <: Number}
    r::T
    c::T
end

struct scaleDataType
    newDim::Placeholder{Int}
    nTiles::Placeholder{Int}
    paddedDim::Placeholder{Int}
    offset::Placeholder{Int}
    scale::Float64
    alpha::Placeholder{Float64}
    isFFT::Bool
    transform::Function
end

chirpZScale(eb::Eigenbrot, xfactor::Real, yfactor::Real) =
    chirpZScaleY(chirpZScaleX(eb,xfactor), yfactor)

function chirpZScaleX(eb::Eigenbrot, factor::Real)
    factor == 1 && return copy(eb)
    h, w = size(eb)
    newDim = Placeholder(h, round(Int, w * factor))
    sc = scaleDataType(
        newDim,                                  # newDim
        Placeholder(1, ceil(Int, newDim.c / w)), # nTiles
        Placeholder(h, 2 * w),                   # paddedDim
        Placeholder(0, w >> 1),                  # offset
        sqrt(2.0),                               # scale
        Placeholder(0.0, π / (w * factor)),      # alpha
        eb.fft,                                  # isFFT
        fftx                                     # transform
    )
    scale_1d(eb, sc)
end

function chirpZScaleY(eb::Eigenbrot, factor::Real)
    factor == 1 && return copy(eb)
    h, w = size(eb)
    newDim = Placeholder(round(Int, h * factor), w)
    sc = scaleDataType(
        newDim,                                  # newDim
        Placeholder(ceil(Int, newDim.r / h), 1), # nTiles
        Placeholder(2 * h, w),                   # paddedDim
        Placeholder(h >> 1, 0),                  # offset
        sqrt(2.0),                               # scale
        Placeholder(π / (h * factor), 0.0),      # alpha
        eb.fft,                                  # isFFT
        ffty                                     # transform
    )
    scale_1d(eb, sc)
end

function scale_1d(input::Eigenbrot, sc::scaleDataType)
    h, w = size(input)
    result = Eigenbrot(sc.newDim.r, sc.newDim.c, sc.isFFT)
    quadProd1 = Eigenbrot(sc.paddedDim.r, sc.paddedDim.c, !sc.isFFT)
    fill!(quadProd1.vals, 0.0)
    inputFFT = sc.transform(input, true)
    inputFFT.fft = sc.isFFT

    for cc in 1:w
        r2x = (div(w, 2) - cc + 1) ^ 2 * sc.alpha.c
        for rr in 1:h
            r2y = (div(h, 2) - rr + 1) ^ 2 * sc.alpha.r
            quadProd1[rr + sc.offset.r, cc + sc.offset.c] =
                inputFFT[rr, cc] * (cos(r2x + r2y) - im * sin(r2x + r2y))
        end
    end
    quadFFT1 = sc.transform(quadProd1, true)
    quadProd2 = Eigenbrot(sc.paddedDim.r, sc.paddedDim.c, !sc.isFFT)
    fftProd = Eigenbrot(sc.paddedDim.r, sc.paddedDim.c, sc.isFFT)
    for deltaHor in 1:sc.nTiles.c
        firstXpixel = (deltaHor == sc.nTiles.c && mod(sc.newDim.c, w) != 0) ?
            w - mod(sc.newDim.c, w) + 1 : 1
        for deltaVer in 1:sc.nTiles.r
            firstYpixel = (deltaVer == sc.nTiles.r && mod(sc.newDim.r, h) != 0) ?
                h - mod(sc.newDim.r, h) + 1 : 1

            for cc in 1:sc.paddedDim.c
                r2x = (div(sc.paddedDim.c, 2) - cc + 1 + (deltaHor - 1) * w) ^ 2 * sc.alpha.c
                for rr in 1:sc.paddedDim.r
                    r2y = (div(sc.paddedDim.r, 2) - rr + 1 + (deltaVer - 1) * h) ^ 2 * sc.alpha.r
                    quadProd2[rr, cc] = cos(r2x + r2y) + im * sin(r2x + r2y)
                end
            end
            quadFFT2 = sc.transform(quadProd2, true)
            fftProd = quadFFT2 * quadFFT1
            conv = sc.transform(fftProd, true)
            for cc in firstXpixel:w
                r2x = (div(w, 2) - cc + 1 + (deltaHor - 1) * w) ^ 2 * sc.alpha.c
                col = mod1(cc + div(sc.newDim.c - w, 2) - (deltaHor - 1) * w + sc.newDim.c,
                           sc.newDim.c)
                for rr in firstYpixel:h
                    r2y = (div(h, 2) - rr + 1 + (deltaVer - 1) * h) ^ 2 * sc.alpha.r
                    row = mod1(rr + div(sc.newDim.r - h, 2) - (deltaVer - 1) * h + sc.newDim.r,
                               sc.newDim.r)
                    newval = sc.scale * conv[rr + sc.offset.r, cc + sc.offset.c] *
                        (cos(r2x + r2y) - im * sin(r2x + r2y))
                    result[row,col] =
                        sc.scale * conv[rr + sc.offset.r, cc + sc.offset.c] *
                        (cos(r2x + r2y) - im * sin(r2x + r2y))
                end
            end
        end
    end
    return result
end
