
struct ScaleData{F <: Function}
    newDim::NTuple{2, Int}
    nTiles::NTuple{2, Int}
    paddedDim::NTuple{2, Int}
    offset::NTuple{2, Int}
    scale::Float64
    alpha::NTuple{2, Float64}
    isFFT::Bool
    transform::F
end

scale_chirpz(eb::Eigenbrot, xfactor::Real, yfactor::Real) =
    scale_y_chirpz(scale_x_chirpz(eb, xfactor), yfactor)

function scale_x_chirpz(eb::Eigenbrot, factor::Real)
    factor == 1 && return copy(eb)
    h, w = size(eb)
    newDim = (h, round(Int, w * factor))
    sc = ScaleData(
        newDim,                                  # newDim
        (1, ceil(Int, newDim[2] / w)),           # nTiles
        (h, 2 * w),                              # paddedDim
        (0, w >> 1),                             # offset
        sqrt(2.0),                               # scale
        (0.0, π / (w * factor)),                 # alpha
        eb.fft,                                  # isFFT
        fftx!,                                   # transform
    )
    scale_1d(eb, sc)
end

function scale_y_chirpz(eb::Eigenbrot, factor::Real)
    factor == 1 && return copy(eb)
    h, w = size(eb)
    newDim = (round(Int, h * factor), w)
    sc = ScaleData(
        newDim,                                  # newDim
        (ceil(Int, newDim[1] / h), 1),           # nTiles
        (2 * h, w),                              # paddedDim
        (h >> 1, 0),                             # offset
        sqrt(2.0),                               # scale
        (π / (h * factor), 0.0),                 # alpha
        eb.fft,                                  # isFFT
        ffty!,                                   # transform
    )
    scale_1d(eb, sc)
end

function scale_1d(input::Eigenbrot, sc::ScaleData)
    h, w = size(input)
    inputFFT = sc.transform(copy(input), true)
    inputFFT.fft = sc.isFFT

    result = Eigenbrot(sc.newDim[1], sc.newDim[2], sc.isFFT)
    quad1 = Eigenbrot(fill(zero(ComplexF64), sc.paddedDim[1], sc.paddedDim[2]), !sc.isFFT)

    for cc = 1:w
        r2x = (div(w, 2) - cc + 1)^2 * sc.alpha[2]
        for rr = 1:h
            r2y = (div(h, 2) - rr + 1)^2 * sc.alpha[1]
            quad1[rr + sc.offset[1], cc + sc.offset[2]] = inputFFT[rr, cc] * cis(-r2x - r2y)
        end
    end
    sc.transform(quad1, true)
    quad2 = Eigenbrot(sc.paddedDim[1], sc.paddedDim[2])
    fftProd = Eigenbrot(sc.paddedDim[1], sc.paddedDim[2])

    for deltaHor = 1:(sc.nTiles[2])
        firstXpixel =
            (deltaHor == sc.nTiles[2] && mod(sc.newDim[2], w) != 0) ? w - mod(sc.newDim[2], w) + 1 :
            1
        for deltaVer = 1:(sc.nTiles[1])
            firstYpixel =
                (deltaVer == sc.nTiles[1] && mod(sc.newDim[1], h) != 0) ?
                h - mod(sc.newDim[1], h) + 1 : 1

            for cc = 1:(sc.paddedDim[2])
                r2x = (div(sc.paddedDim[2], 2) - cc + 1 + (deltaHor - 1) * w)^2 * sc.alpha[2]
                for rr = 1:(sc.paddedDim[1])
                    r2y = (div(sc.paddedDim[1], 2) - rr + 1 + (deltaVer - 1) * h)^2 * sc.alpha[1]
                    quad2[rr, cc] = cis(r2x + r2y)
                end
            end
            quad2.fft = !sc.isFFT
            fftProd.vals .= sc.transform(sc.transform(quad2, true) * quad1, true).vals

            for cc = firstXpixel:w
                r2x = (div(w, 2) - cc + 1 + (deltaHor - 1) * w)^2 * sc.alpha[2]
                col = mod1(
                    cc + div(sc.newDim[2] - w, 2) - (deltaHor - 1) * w + sc.newDim[2],
                    sc.newDim[2],
                )
                for rr = firstYpixel:h
                    r2y = (div(h, 2) - rr + 1 + (deltaVer - 1) * h)^2 * sc.alpha[1]
                    row = mod1(
                        rr + div(sc.newDim[1] - h, 2) - (deltaVer - 1) * h + sc.newDim[1],
                        sc.newDim[1],
                    )
                    result[row, col] =
                        sc.scale * fftProd[rr + sc.offset[1], cc + sc.offset[2]] * cis(-r2x - r2y)
                end
            end
        end
    end
    return result
end
