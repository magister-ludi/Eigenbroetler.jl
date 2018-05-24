#!/usr/bin/env julia

using Eigenbroetler

function scaled_sinusoid(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    return sin(100 * r) / (r + 1)
end

eb = fft(Eigenbrot(scaled_sinusoid, 512, 512))
save("scaled_sinusoid.png", eb, ImageSetting(Magn, Log),
     colours = joinpath(dirname(@__FILE__), "iron.map"))
