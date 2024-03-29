
module Eigenbroetler

using FileIO
using FFTW
using FITSIO
using FITSIO.Libcfitsio
using ImageCore
using Statistics

export Eigenbrot
export Butterworth, Gaussian, Hamming, Hann, Square
#export Parzen Welch

export width, height
export size
export isFFT, save
export fft, fftx, ffty, removeDC
export ImageSetting
export LinearScale, LogScale, RootScale
export RealPart, ImagPart, Magn, Phase
export image, pad, extend
export flipver, fliphor, swapxy
export pixel, coords
export scale_chirpz, scale_x_chirpz, scale_y_chirpz
export hilbert_x, hilbert_y
export fft_xshear, fft_yshear, rotate90, fourierRotation
export highpass, lowpass, conv, cor

include("colourmaps.jl")
include("eigenbrot.jl")
include("simple.jl")
include("shaped_pad.jl")
include("fft.jl")
include("chirpz.jl")
include("transform.jl")
include("filters.jl")
include("fourier_multiply.jl")

end
