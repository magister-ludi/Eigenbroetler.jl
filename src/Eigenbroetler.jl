# TODO: remove comment
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
             setindex!
import FileIO.save

export Eigenbrot
export width, height
export isValid, isFFT
export save, show, fft
export fill!, getindex, setindex!
export index, coords

include("eigenbrot.jl")
include("fft.jl")

end
