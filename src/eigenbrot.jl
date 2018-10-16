
"""
    Eigenbrot(filename)
Construct an Eigenbrot from a file, which should be a FITS
file (including a gzipped FITS file), or an image type that
can be loaded by ImageMagick.

    Eigenbrot(data::Matrix{ComplexF64}, fft = false)
Construct an Eigenbrot from a matrix of complex values.

    Eigenbrot(rows::Integer, cols::Integer, fft = false)
Construct an unitialised Eigenbrot of size `rows×cols`.

    Eigenbrot(f::Function, rows::Integer, cols::Integer, fft = false)
Construct an Eigenbrot of size `rows×cols`, using `f(x, y)` to
calculate the data values.
"""
mutable struct Eigenbrot
    vals::Matrix{ComplexF64}
    fft::Bool
    have_min_max::Bool
    maxCmp::Float64
    minCmp::Float64
    maxMag::Float64
    Eigenbrot(data::Matrix{T}, fft = false) where {T <: Number} =
        new(ComplexF64.(data), fft, false)
end

Eigenbrot(rows::Integer, cols::Integer, fft = false) =
        Eigenbrot(Matrix{ComplexF64}(undef, rows, cols), fft)

@enum Scale Linear Log Root
@enum Component RealPart ImagPart Magn Phase

struct ImageSetting
    c::Component
    s::Scale
    p::Int
    ImageSetting(component::Component, scale::Scale, power::Integer = 2) =
        new(component, scale, power)
end

function Eigenbrot(f::Function, rows::Integer, cols::Integer, fft = false)
    xMax = div(cols, 2) - 1
    xMin = xMax - cols + 1
    yMax = div(rows, 2) - 1
    yMin = yMax - rows + 1
    return Eigenbrot([ComplexF64(f(x, y)) for y in yMin:yMax, x in xMin:xMax], fft)
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
    isfft, comment = read_key(hdu, "ISFFT")
    data = read(hdu)

    if ndims(data) > 2
        w, h = size(data)
        cdata = Matrix{ComplexF64}(undef, h, w)
        for c in 1:w
            for r in 1:h
                cdata[r, c] = ComplexF64(data[c, r, 1], data[c, r, 2])
            end
        end
    elseif ndims(data) == 2
        cdata = map(x -> ComplexF64(x), data)
        reshape(cdata, size(data, 1), size(data, 2))
    else
        cdata = map(x -> ComplexF64(x), data)
    end
    return Eigenbrot(cdata, isfft)
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
    w, h = size(img)
    cdata = Matrix{ComplexF64}(undef, h, w)
    for c in 1:w
        for r in 1:h
            cdata[h - r + 1, c] = ComplexF64(round(UInt8, 255 * img[r, c].val))
        end
    end
    return Eigenbrot(cdata)
end

const fits_magic = b"SIMPLE "
const gzip_magic = [0x1f, 0x8b]
const max_magic_len = max(length(fits_magic), length(gzip_magic))

function Eigenbrot(filename::AbstractString)
    t = open(filename)
    by = read(t, max_magic_len)
    close(t)
    if by[1:length(fits_magic)] == fits_magic || by[1:length(gzip_magic)] == gzip_magic
        return read_fits(filename)
    else
        return read_image(filename)
    end
end

Images.width(eb::Eigenbrot) = size(eb.vals, 2)
Images.height(eb::Eigenbrot) = size(eb.vals, 1)
isFFT(eb::Eigenbrot) = eb.fft

for name in (:ndims, :length, :size)
    @eval begin
        Base.$name(eb::Eigenbrot) = $name(eb.vals)
    end
end

Base.size(eb::Eigenbrot, args...) = size(eb.vals, args...)

Base.getindex(eb::Eigenbrot, row::Integer, col::Integer) = getindex(eb.vals, row, col)
Base.getindex(eb::Eigenbrot,
              rows::AbstractRange{T},
              cols::AbstractRange{U}) where {T <: Integer, U <: Integer} =
                  getindex(eb.vals, rows, cols)
Base.getindex(eb::Eigenbrot, indices::Vector{T}) where {T <: Integer} =
    getindex(eb.vals, indices)
Base.getindex(eb::Eigenbrot, rc::Tuple{T, U}) where {T <: Integer, U <: Integer} =
                  getindex(eb.vals, rc[1], rc[2])

function Base.setindex!(eb::Eigenbrot, val::Number, row::Integer, col::Integer)
    reset!(eb)
    setindex!(eb.vals, val, row, col)
end

function Base.setindex!(eb::Eigenbrot, val::Number, indices::Vector{T}) where {T <: Integer}
    reset!(eb)
    setindex!(eb.vals, val, indices)
end

function Base.setindex!(eb::Eigenbrot, val::Number,
                        rc::Tuple{T, U}) where {T <: Integer, U <: Integer}
    reset!(eb)
    setindex!(eb.vals, val, rc[1], rc[2])
end

reset!(eb::Eigenbrot) = (eb.have_min_max = false)

"""
    fill!(eb, x)
Fill Eigenbrot `eb` with the number `x`.
"""
function Base.fill!(eb::Eigenbrot, x::Number)
    reset!(eb)
    fill!(eb.vals, x)
end

"""
    save(filename::AbstractString, eb::Eigenbrot)
Write Eigenbrot `eb` to a file in FITS format

    save(filename::AbstractString, eb::Eigenbrot, cmp::ImageSetting; colours::Symbol)
Save a bitmap representation of a component of Eigenbrot `eb` to
`filename`. The file format will be determined by the file name.
The component of the data and the scaling method are determined by `cmp`.

    save(filename, eb::Eigenbrot, cmp::ImageSetting, cmp2::ImageSetting; colours::Symbol)
Save a bitmap representation of two components of Eigenbrot `eb` to
`filename`. The file format will be determined by the file name.
The components of the data and the scaling method are determined by `cmp`
and `cmp2`.
"""
function FileIO.save(filename::AbstractString, eb::Eigenbrot)
    keys = ["COMMENT", "ISFFT"]
    values = [nothing, eb.fft]
    comments = ["CREATED BY AN EIGENBROETLER",
                eb.fft ? "image in Fourier space?" : "image in real space?"]
    hdr = FITSHeader(keys, values, comments)
    fits = FITS(filename, "w")
    w = width(eb)
    h = height(eb)
    hlen = w * h
    lpData = Vector{Float64}(undef, hlen * 2)
    real = 1
    imag = 1 + hlen
    for r in 1:h
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
    image(eb::Eigenbrot, cmp::ImageSetting; colours::Symbol)
Construct a bitmap representation of a component of Eigenbrot `eb`.
The component of the data and the scaling method are determined by `cmp`.

    image(eb::Eigenbrot, cmp::ImageSettingb, cmp2::ImageSetting; colours::Symbol)
Construct a bitmap representation of two components of Eigenbrot `eb`.
The components of the data and the scaling method are determined by `cmp`
and `cmp2`.
"""
function image(eb::Eigenbrot,
               cmp::ImageSetting,
               cmp2::Union{Nothing, ImageSetting} = nothing;
               colours::Union{Symbol, AbstractString} = :grey)
    w = w2 = width(eb)
    h = height(eb)
    hlen = w * h
    calculate_ranges!(eb)
    (cmp2 != nothing) && (w2 *= 2)
    pal = palette(colours)
    img = Matrix{RGB{N0f8}}(undef, h, w2)
    calculate_pixels!(img, eb, cmp, pal, 0)
    (cmp2 != nothing) && calculate_pixels!(img, eb, cmp2, pal, h * w)
    return img
end

FileIO.save(filename::AbstractString, eb::Eigenbrot,
     cmp::ImageSetting,
     cmp2::Union{Nothing, ImageSetting} = nothing;
     colours::Union{Symbol, AbstractString} = :grey) =
         save(filename, image(eb, cmp, cmp2, colours = colours))

scaleLogarithmic(x::Float64, ::Int) = log10(x)
scaleLinear(x::Float64, ::Int) = x
scaleRoot(x::Float64, p::Int) = x ^ (1.0 / p)

function calculate_pixels!(img::Matrix{RGB{N0f8}}, eb::Eigenbrot,
                           cmp::ImageSetting, palette::Palette, idx::Integer)
    h, w = size(eb.vals)
    if cmp.c == Phase
        # Phase part
        minValue = -π
        maxValue = +π
        scaleFactor = (NUM_COLOURS - 1) / (maxValue - minValue)
        for c in 1:w
            for r in h:-1:1
                k = trunc(Int, scaleFactor * (angle(eb.vals[r, c]) - minValue))
                idx += 1
                img[idx] = palette[1 + k]
            end
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
        for c in 1:w
            for r in h:-1:1
                k = trunc(Int, scaler(abs(eb.vals[r, c]) - offset, cmp.p) * scaleFactor)
                idx += 1
                img[idx] = palette[1 + k]
            end
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
        for c in 1:w
            for r in h:-1:1
                k = trunc(Int, scaler(select(eb.vals[r, c]) - offset, cmp.p) * scaleFactor)
                idx += 1
                img[idx] = palette[1 + k]
            end
        end
    end
end

function Base.show(io::IO, eb::Eigenbrot)
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

Base.similar(eb::Eigenbrot) = Eigenbrot(height(eb), width(eb))

Base.copy(eb::Eigenbrot) = Eigenbrot(copy(eb.vals), eb.fft)

"""
    pixel(eb::Eigenbrot, x::Real, y::Real)
Return a tuple , `(row, column)`, indexing the closest position
in `eb` to to the Cartesian location `(x, y)`. The inverse of `coords`.
    pixel(eb::Eigenbrot, xy::Tuple{Real, Real})
Return a tuple , `(row, column)`, indexing the closest position
in `eb` to to the Cartesian location given by `xy`. The inverse of `coords`.
"""
function pixel(eb::Eigenbrot, x::Real, y::Real)
    h, w = size(eb.vals)
    xMax = div(w, 2)
    xMin = xMax - w
    yMax = div(h, 2)
    yMin = yMax - h
    (round(Int, 1 + h * (y - yMin) / height(eb)),
     round(Int, 1 + w * (x - xMin) / width(eb)))
end

pixel(eb::Eigenbrot, xy::Tuple{T, U}) where {T <: Real, U <: Real} =
    pixel(eb, xy[1], xy[2])

"""
    coords(eb::Eigenbrot, r::Integer, c::Integer)
Return a tuple, `(x, y)`, of the Cartesian location
mapped by the pixel `eb[r, c]`. The inverse of `pixel`.
    coords(eb::Eigenbrot, rc::Tuple{Integer, Integer})
Return a tuple , `(x, y)`, of the Cartesian location
mapped by the pixel `eb[rc]`. The inverse of `pixel`.
"""
function coords(eb::Eigenbrot, r::Integer, c::Integer)
    h, w = size(eb.vals)
    xMin = div(w, 2) - w
    yMin = div(h, 2) - h
    (xMin + (c - 1) * width(eb) / w,
     yMin + (r - 1) * height(eb) / h)
end

coords(eb::Eigenbrot, rc::Tuple{T, U}) where {T <: Integer, U <: Integer} =
    coords(eb, rc[1], rc[2])
