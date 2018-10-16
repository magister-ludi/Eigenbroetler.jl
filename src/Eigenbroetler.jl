
module Eigenbroetler

using FileIO
using FFTW
using FITSIO
using FITSIO.Libcfitsio
using Images

export Eigenbrot
export width, height
export size
export isFFT, save
export fft, fftx, ffty, removeDC
export ImageSetting
export Linear, Log, Root
export RealPart, ImagPart, Magn, Phase
export image, pad
export flipver, fliphor, swapxy
export pixel, coords
export chirpZScale, chirpZScaleX, chirpZScaleY
export hilbertX, hilbertY
export fft_xshear, fft_yshear, rotate90, fourierRotation

include("colourmaps.jl")
include("eigenbrot.jl")
include("simple.jl")
include("fft.jl")
include("chirpz.jl")
include("transform.jl")

end
