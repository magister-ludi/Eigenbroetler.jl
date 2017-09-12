
"""
    Eigenbrot(filename)
Construct an Eigenbrot from a file, which should be a FITS
file (including a gzipped FITS file), or an image type that
can be loaded by ImageMagick.

    Eigenbrot(data::Matrix{Complex128})
Construct an Eigenbrot from a matrix of complex values.

    Eigenbrot(rows::Integer, cols::Integer)
Construct an unitialised Eigenbrot of size `rows×cols`.

    Eigenbrot(rows::Integer, cols::Integer, f::Function)
Construct an Eigenbrot of size `rows×cols`, using `f(x, y)` to
calculate the data values.
"""
mutable struct Eigenbrot
    vals::Matrix{Complex128}
    errorString::AbstractString
    fft::Bool
    have_min_max::Bool
    maxCmp::Float64
    minCmp::Float64
    maxMag::Float64
    Eigenbrot(data::Matrix{Complex128}, fft = false) = new(data, "", fft, false)
    Eigenbrot(rows::Integer, cols::Integer, fft = false) =
        new(Matrix{Complex128}(rows, cols), "", fft, false)
end

@enum Scale Linear Log Root
@enum Component RealPart ImagPart Magn Phase

struct ImageSetting
    c::Component
    s::Scale
    p::Int
    ImageSetting(component::Component, scale::Scale, power::Integer = 2) =
        new(component, scale, power)
end

function Eigenbrot(rows::Integer, cols::Integer, f::Function)
    xMax = div(cols, 2) - 1
    xMin = xMax - cols + 1
    yMax = div(rows, 2) - 1
    yMin = yMax - rows + 1
    return Eigenbrot([Complex128(f(x, y)) for x in xMax:-1:xMin, y in yMin:yMax])
end

function read_fits(file::AbstractString)
    fits = FITS(file)
    hdu = nothing
    for u in fits
        if isa(u, ImageHDU)
            hdu = u
            break
        end
    end
    if hdu == nothing
        error("No image found in FITS file '$file'")
    end
    data = read(hdu)

    if ndims(data) > 2
        w, h = size(data)
        cdata = Matrix{Complex128}(h, w)
        for c in 1:w
            for r in 1:h
                cdata[h - r + 1, c] = Complex128(data[c, r, 1], data[c, r, 2])
            end
        end
    elseif ndims(data) == 2
        cdata = map(x -> Complex128(x), data)
        reshape(cdata, size(data, 1), size(data, 2))
    else
        cdata = map(x -> Complex128(x), data)
    end
    return cdata
end

function read_image(file::AbstractString)
    img = nothing
    try
        img = load(file)
    catch
        error("Can't read file $file.")
    end
    #=
    Images.jl greyscale conversion uses the
    values from CCIR 601, as does the C++-version of
    Eigenbroetler. See
    https://en.wikipedia.org/wiki/Grayscale
    =#
    img = Gray.(img)
    cdata = Vector{Complex128}(length(img))
    for i in 1:length(img)
        cdata[i] = Complex128(round(UInt8, 255 * img[i].val))
    end
    return reshape(cdata, size(img))
end

const fits_magic = b"SIMPLE "
const gzip_magic = [0x1f, 0x8b]
const max_magic_len = max(length(fits_magic), length(gzip_magic))

function Eigenbrot(filename::AbstractString)
    t = open(filename)
    by = read(t, max_magic_len)
    close(t)
    if by[1:length(fits_magic)] == fits_magic || by[1:length(gzip_magic)] == gzip_magic
        data = read_fits(filename)
    else
        data = read_image(filename)
    end
    return Eigenbrot(data)
end

Images.width(eb::Eigenbrot) = size(eb.vals, 2)
Images.height(eb::Eigenbrot) = size(eb.vals, 1)
Base.size(eb::Eigenbrot) = size(eb.vals)
isValid(eb::Eigenbrot) = isempty(eb.errorString)
isFFT(eb::Eigenbrot) = eb.fft

getindex(eb::Eigenbrot, r::Integer, c::Integer) =
    getindex(eb.vals, r, c)
getindex(eb::Eigenbrot, i::Integer) =
    getindex(eb.vals, i)
getindex{T<:Real, U<:Real}(eb::Eigenbrot, rc::Tuple{T, U}) =
    getindex(eb.vals, rc[1], rc[2])
setindex!(eb::Eigenbrot, v::Number, r::Integer, c::Integer) =
    setindex!(eb.vals, r, c)
setindex!(eb::Eigenbrot, v::Number, i::Integer) =
    setindex!(eb.vals, i)
setindex!{T<:Real, U<:Real}(eb::Eigenbrot, rc::Tuple{T, U}) =
    setindex!(eb.vals, rc[1], rc[2])
reset!(eb::Eigenbrot) = (eb.have_min_max = false)

"""
    fill!(eb, x)
Fill Eigenbrot `eb` with the number `x`.
"""
fill!(eb::Eigenbrot, x::Number) = fill!(eb.vals, x)

"""
    save(filename, eb)
Write Eigenbrot `eb` to a file in FITS format

    save(filename, eb, cmp::ImageSetting; colours::Symbol)
Save a bitmap representation of a component of Eigenbrot `eb` to
`filename`. The file format will be determined by the file name.
The component of the data and the scaling method are determined by `cmp`.

    save(filename, eb, cmp::ImageSettingb, cmp2::ImageSetting; colours::Symbol)
Save a bitmap representation of two components of Eigenbrot `eb` to
`filename`. The file format will be determined by the file name.
The components of the data and the scaling method are determined by `cmp`
and `cmp2`.
"""
function save(filename::AbstractString, eb::Eigenbrot)
    keys = ["COMMENT", "ISFFT"]
    values = [nothing, eb.fft]
    comments = ["CREATED BY AN EIGENBROETLER",
                eb.fft ? "image in Fourier space?" : "image in real space?"]
    hdr = FITSHeader(keys, values, comments)
    fits = FITS(filename, "w")
    w = width(eb)
    h = height(eb)
    hlen = w * h
    lpData = Vector{Float64}(hlen * 2)
    real = 1
    imag = 1 + hlen
    for r in h:-1:1
        for c in 1:w
            v = eb.vals[r, c]
            lpData[real] = v.re
            lpData[imag] = v.im
            real += 1
            imag += 1
        end
    end
    lpData = reshape(lpData, w, h, 2)
    write(fits, lpData, header=hdr, name=nothing, ver=nothing)
    close(fits)
end

"""
    image(filename, eb, cmp::ImageSetting; colours::Symbol)
Construct a bitmap representation of a component of Eigenbrot `eb`.
The component of the data and the scaling method are determined by `cmp`.

    image(filename, eb, cmp::ImageSettingb, cmp2::ImageSetting; colours::Symbol)
Construct a bitmap representation of two components of Eigenbrot `eb`.
The components of the data and the scaling method are determined by `cmp`
and `cmp2`.
"""
function image(eb::Eigenbrot,
               cmp::ImageSetting,
               cmp2::Union{Void, ImageSetting} = nothing;
               colours::Union{Symbol, AbstractString} = :grey)
    w = w2 = width(eb)
    h = height(eb)
    hlen = w * h
    calculate_ranges!(eb)
    (cmp2 != nothing) && (w2 *= 2)
    pal = palette(colours)
    img = Matrix{RGB{N0f8}}(h, w2)
    calculate_pixels!(img, eb, cmp, pal, 0)
    (cmp2 != nothing) && calculate_pixels!(img, eb, cmp2, pal, h * w)
    return img
end

save(filename::AbstractString, eb::Eigenbrot,
     cmp::ImageSetting,
     cmp2::Union{Void, ImageSetting} = nothing;
     colours::Union{Symbol, AbstractString} = :grey) =
         save(filename, image(eb, cmp, cmp2, colours = colours))

scaleLogarithmic(x::Float64, ::Int) = log10(x)
scaleLinear(x::Float64, ::Int) = x
scaleRoot(x::Float64, p::Int) = x ^ (1.0 / p)

function calculate_pixels!(img::Matrix{RGB{N0f8}}, eb::Eigenbrot,
                           cmp::ImageSetting, palette::Palette, idx::Integer)
    if cmp.c == Phase
        # Phase part
        minValue = -π
        maxValue = +π
        scaleFactor = (NUM_COLOURS - 1) / (maxValue - minValue)
        for cval in eb.vals
            idx += 1
            k = trunc(Int, scaleFactor * (angle(cval) - minValue))
            img[idx] = palette[1 + k]
        end
    elseif cmp.c == Magn
        # Magnitude part
        scaler = cmp.s == Linear ? scaleLinear : (cmp.s == Log ? scaleLogarithmic : scaleRoot)
        minValue = 0.0
        maxValue = eb.maxMag
        offset = (cmp.s == Log) ? 1.0 : 0.0
        if isapprox(maxValue, minValue)
            scaleFactor = (NUM_COLOURS >> 1) - 1.0
            minValue -= offset
            maxValue -= offset
        else
            offset = minValue - offset
            scaleFactor = (NUM_COLOURS - 1) /
                (scaler(maxValue - offset, cmp.p) - scaler(minValue - offset, cmp.p))
        end
        for cval in eb.vals
            idx += 1
            k = trunc(Int, scaler(abs(cval) - offset, cmp.p) * scaleFactor)
            img[idx] = palette[1 + k]
        end
    else
        scaler = scaleLinear
        minValue = eb.minCmp
        maxValue = eb.maxCmp
        if isapprox(maxValue, minValue)
            scaleFactor = (NUM_COLOURS >> 1) - 1.0
            offset = 0.0
        else
            offset = minValue
            scaleFactor = (NUM_COLOURS - 1) /
                (scaler(maxValue - offset, cmp.p) - scaler(minValue - offset, cmp.p))
        end
        select = (cmp.c == RealPart) ? real : imag
        for cval in eb.vals
            idx += 1
            k = trunc(Int, scaler(select(cval) - offset, cmp.p) * scaleFactor)
            img[idx] = palette[1 + k]
        end
    end
end

function show(io::IO, eb::Eigenbrot)
    print(io, "Eigenbrot(", height(eb), "x", width(eb), ", ",
          isFFT(eb) ? "Fourier" : "Real", " space)")
 end

"""
    calculate_ranges!(eb::Eigenbrot)
Calculate minimum and maximum values in an `Eigenbrot`.
Used for calculating pixel values in an image.
"""
function calculate_ranges!(eb::Eigenbrot)
    if eb.have_min_max
        return
    end
    eb.maxCmp = typemin(Float64)
    eb.minCmp = typemax(Float64)
    eb.maxMag = 0.0
    for x in eb.vals
        eb.maxCmp = max(eb.maxCmp, real(x), imag(x))
        eb.minCmp = min(eb.minCmp, real(x), imag(x))
        eb.maxMag = max(eb.maxMag, abs(x))
    end
    eb.have_min_max = true
    return
end

similar(eb::Eigenbrot) = Eigenbrot(height(eb), width(eb))

copy(eb::Eigenbrot) = Eigenbrot(copy(eb.vals), eb.fft)
