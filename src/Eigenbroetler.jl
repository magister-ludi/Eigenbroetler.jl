
isdefined(Base, :__precompile__) && __precompile__()

module Eigenbroetler

using FITSIO
using FileIO
using Images
using FITSIO.Libcfitsio

export Eigenbrot
export width, height, size
export isValid, isFFT
export save, show
export fill!, getindex, setindex!
export index, coords
export fft, fftx, ffty, removeDC
export ImageSetting
export Linear, Log, Root
export RealPart, ImagPart, Magn, Phase
export image, pad, pow2pad, pow2pad!
export flipver, fliphor, swapxy
export pixel, coords
export chirpZScale, chirpZScaleX, chirpZScaleY

include("colourmaps.jl")
include("eigenbrot.jl")
include("simple.jl")
include("fft.jl")
include("chirpz.jl")

end
