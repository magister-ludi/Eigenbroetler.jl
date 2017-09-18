#!/usr/bin/env julia

using Eigenbroetler

function starry(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    phi = atan2(y, x)
    return sin(r * (3 + cos(7 * phi)) / 4)
end

eb = Eigenbrot(starry, 512, 512)
save("starry.png", eb, ImageSetting(RealPart, Linear))
save("starry_dft.png", fft(eb), ImageSetting(Magn, Log))
