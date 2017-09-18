#!/usr/bin/env julia

using Eigenbroetler

function high_freq_sinusoid(x::Number, y::Number)
    # sin(100 * r)
    r = sqrt(x * x + y * y)
    return sin(100 * r)
end

eb = Eigenbrot(high_freq_sinusoid, 512, 512)
save("high_freq_sinusoid.png", eb, ImageSetting(RealPart, Linear))
save("high_freq_sinusoid_dft.png", fft(eb), ImageSetting(Magn, Log))
