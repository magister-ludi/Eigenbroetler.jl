#!/usr/bin/env julia

using Eigenbroetler

squarified_sinusoid(x::Number, y::Number) =
    exp(im * (x ^ 4 + y ^ 4) ^ (1 / 4))

eb = Eigenbrot(512, 512, squarified_sinusoid)
save("squarified_sinusoid.png", eb, ImageSetting(RealPart, Linear))
save("squarified_sinusoid_dft.png", fft(eb), ImageSetting(Magn, Log))
