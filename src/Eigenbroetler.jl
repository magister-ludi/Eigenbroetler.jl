
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
export LinearScale, LogScale, RootScale
export RealPart, ImagPart, Magn, Phase
export image, pad
export flipver, fliphor, swapxy
export pixel, coords
export scale_chirpz, scale_x_chirpz, scale_y_chirpz
export hilbert_x, hilbert_y
export fft_xshear, fft_yshear, rotate90, fourierRotation

include("colourmaps.jl")
include("eigenbrot.jl")
include("simple.jl")
include("fft.jl")
include("chirpz.jl")
include("transform.jl")

end
