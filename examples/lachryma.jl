#!/usr/bin/env julia

using Eigenbroetler

function lachryma(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    return exp(-(r / 100) ^ 2) * exp(2 * im * ((x + y) * exp(-(r / 100) ^ 2)))
end

eb = fft(Eigenbrot(lachryma, 512, 512))
save("lachryma.png", eb, ImageSetting(Magn, LogScale), colours = :rainbow)
