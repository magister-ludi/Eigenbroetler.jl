#!/usr/bin/env julia

using Eigenbroetler

squarified_sinusoid(x::Number, y::Number) =
    exp(im * (x ^ 4 + y ^ 4) ^ (1 / 4))

eb = Eigenbrot(squarified_sinusoid, 512, 512)
save("squarified_sinusoid.png", eb, ImageSetting(RealPart, LinearScale))
save("squarified_sinusoid_dft.png", fft(eb), ImageSetting(Magn, LogScale))
