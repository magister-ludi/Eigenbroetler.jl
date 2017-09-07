
isdefined(Base, :__precompile__) && __precompile__()

module Eigenbroetler

using FITSIO
using FileIO
using Images
using FITSIO.Libcfitsio

import Base: show,
             fft,
             fill!,
             getindex,
             setindex!,
             similar
import FileIO.save

export Eigenbrot
export width, height
export isValid, isFFT
export save, show
export fill!, getindex, setindex!
export index, coords
export fft, fftx, ffty, removeDC
export ImageSetting
export Linear, Log, Root
export RealPart, ImagPart, Magn, Phase
export image, trim, pow2pad, pow2pad!

include("colourmaps.jl")
include("eigenbrot.jl")
include("simple.jl")
include("fft.jl")

end
