
type Eigenbrot
  vals::Matrix{Complex128}
  errorString::AbstractString
  fft::Bool
  have_min_max::Bool
  maxCmp::Float64
  minCmp::Float64
  maxMag::Float64
  minMag::Float64
  Eigenbrot(data::Matrix{Complex128}) = new(data, "", false, false)
  Eigenbrot(rows::Integer, cols::Integer) = new(Matrix{Complex128}(rows, cols), "", false, false)
end

function Eigenbrot(w::Integer, h::Integer, f::Function)
    xMax = div(w, 2)
    xMin = xMax - w + 1
    yMax = div(h, 2)
    yMin = yMax - h + 1
    return Eigenbrot([Complex128(f(x, y)) for y in yMax:-1:yMin, x in xMin:xMax])
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
   #println("----\nFITS file '$file' is good")
   data = read(hdu)
   #println("    data: ", typeof(data), size(data))

   if ndims(data) > 2
     cdata = Matrix{Complex128}(size(data, 1), size(data, 2))
     for r in 1:size(data, 1)
       for c in size(data, 2)
         cdata[r, c] = Complex128(data[r, c, 1], data[r, c, 2])
       end
     end
   elseif ndims(data) == 2
     cdata = map(x -> Complex128(x), data)
     reshape(cdata, size(data, 1), size(data, 2))
   else
     cdata = map(x -> Complex128(x), data)
   end
   #println("   cdata: ", typeof(cdata), size(cdata))
   #println()
   return cdata
end

#=
The conversion from RGB to greyscale uses the
values from CCIR 601. See
https://en.wikipedia.org/wiki/Grayscale
=#
const RED_INTENSITY = 0.299
const GREEN_INTENSITY = 0.587
const BLUE_INTENSITY = 0.114

function read_image(file::AbstractString)
  img = load(file)
  data = raw(img)
  if size(data, 1) != 3
    error("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  end
  if size(data, 1) == 1
    # TODO: does this ever happen?
  elseif size(data, 1) == 3
    cdata = Matrix{Complex128}(size(data, 2), size(data, 3))
    for r in 1:size(data, 2)
      for c in 1:size(data, 3)
        cdata[r, c] = RED_INTENSITY   * data[1, r, c] +
                      GREEN_INTENSITY * data[2, r, c] +
                      BLUE_INTENSITY  * data[3, r, c]
      end
    end
  else
    error("Raw data has dimensions ", size(data), ", sorry don't know what to do.")
  end
  return cdata
end

const fits_magic = b"SIMPLE "
const gzip_magic = [0x1f, 0x8b]
const max_magic_len = max(length(fits_magic), length(gzip_magic))

# Workaround for magic number used in FileIO
if haskey(FileIO.sym2info, :PCX)
    delete!(FileIO.sym2info, :PCX)
    add_format(format"PCX", UInt8[0x0a,0x05,0x01], ".pcx", [:ImageMagick])
end

"""
    Eigenbrot(filename)
Construct an Eigenbrot from a file, which should be a FITS
file (including a gzipped FITS file), or an image type that
can be loaded by ImageMagick.
"""
function Eigenbrot(filename::AbstractString)
    t = open(filename)
    by = read(t, max_magic_len)
    close(t)
    if by[1:length(fits_magic)] == fits_magic || by[1:length(gzip_magic)] == gzip_magic
        data = read_fits(filename)
    else
        try
            data = read_image(filename)
        catch
            error("Can't read '$filename': unknown format")
        end
    end
    return Eigenbrot(data)
end

Images.width(eb::Eigenbrot) = size(eb.vals, 2)
Images.height(eb::Eigenbrot) = size(eb.vals, 1)
isValid(eb::Eigenbrot) = isempty(eb.errorString)
isFFT(eb::Eigenbrot) = eb.fft

getindex(eb::Eigenbrot, r::Integer, c::Integer) =
  getindex(eb.vals, r, c)
getindex(eb::Eigenbrot, i::Integer) =
  getindex(eb.vals, i)
getindex{T<:Real, U<:Real}(eb::Eigenbrot, rc::Tuple{T, U}) =
  getindex(eb.vals, r[1], c[2])
setindex!(eb::Eigenbrot, v::Number, r::Integer, c::Integer) =
  setindex!(eb.vals, r, c)
setindex!(eb::Eigenbrot, v::Number, i::Integer) =
  setindex!(eb.vals, i)
setindex!{T<:Real, U<:Real}(eb::Eigenbrot, rc::Tuple{T, U}) =
  setindex!(eb.vals, r[1], c[2])
reset!(eb::Eigenbrot) = (eb.have_min_max = false)
"""
    fill!(eb, x)
Fill Eigenbrot `eb` with the number `x`.
"""
fill!(eb::Eigenbrot, x::Number) = fill!(eb.vals, Complex128(x))

"""
Return a tuple containing a pair of indices corresponding to the
Cartesian coordinates (x, y). The coordinate values will be rounded
to the nearest integer.
"""
index(eb::Eigenbrot, x::Real, y::Real) =
    (1 + div(width(eb), 2) - round(Int, y), round(Int, x) + div(height(eb), 2))
"""
Return a tuple containing a pair of indices corresponding to the
Cartesian coordinates xy. The coordinate values will be rounded
to the nearest integer.
"""
index{T<:Real, U<:Real}(eb::Eigenbrot, xy::Tuple{T, U}) =
    index(eb, xy[1], xy[2])

"""
    coords(eb, row, col)
Return a tuple containing the Cartesian coordinates of the
element at the given row and column of Eigenbrot `eb`.
"""
coords(eb::Eigenbrot, row::Integer, col::Integer) =
    (col - div(height(eb), 2), 1 + div(width(eb), 2) - row)
"""
    coords(eb, rc)
Return a tuple containing the Cartesian coordinates of the
element `eb[rc]` of Eigenbrot `eb`, where `rc` is a two-tuple of integers.
"""
coords{T<:Integer, U<:Integer}(eb::Eigenbrot, rc::Tuple{T,U}) =
    coords(eb, rc[1], rc[2])
"""
    coords(eb, i)
Return a tuple containing the Cartesian coordinates of the
element `eb[i]` of Eigenbrot `eb`, where `i` is an integer.
"""
coords(eb::Eigenbrot, i::Integer) =
    coords(eb, mod1(i, height(eb)), 1 + div(i, height(eb)))

"""
    save(filename, eb)
Write Eigenbrot `eb` to a file in FITS format
"""
function save(filename::AbstractString, eb::Eigenbrot)
    keys = ["COMMENT", "ISFFT"]
    values = [nothing, true]
    comments = ["CREATED BY AN EIGENBROETLER", "image in Fourier space?"]
    hdr = FITSHeader(keys, values, comments)
    fits = FITS(filename, "w")
    w = width(eb)
    h = height(eb)
    hlen = w * h
    lpData = Vector{Float64}(hlen * 2)
    real = 1
    imag = 1 + hlen

    for i in 1:hlen
        lpData[real] = eb.vals[i].re
        lpData[imag] = eb.vals[i].im
        real += 1
        imag += 1
    end
    lpData = reshape(lpData, w, h, 2)
    write(fits, lpData, header=hdr, name=nothing, ver=nothing)
    close(fits)
    return
end

function show(io::IO, eb::Eigenbrot)
  print(io, "Eigenbrot(", width(eb), "x", height(eb), ", ", isFFT(eb) ? "Fourier" : "Real", " space)")
end
