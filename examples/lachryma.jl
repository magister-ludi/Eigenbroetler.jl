#!/usr/bin/env julia

using Eigenbroetler

function lachryma(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    return exp(-(r / 100) ^ 2) * exp(2 * im * ((x + y) * exp(-(r / 100) ^ 2)))
end

eb = fft(Eigenbrot(512, 512, lachryma))
save("lachryma.png", eb, ImageSetting(Magn, Log), colours = :rainbow)
