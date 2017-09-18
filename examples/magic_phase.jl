#!/usr/bin/env julia

using Eigenbroetler

function magic_phase(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    return exp(-(r / 100) ^ 2) * exp(2 * im * x * (x ^ 2 - y ^ 2) / 200 ^ 2)
end

eb = Eigenbrot(magic_phase, 512, 512)
save("magic_phase.png", eb, ImageSetting(Phase, Linear))
save("magic_phase_dft.png", fft(eb), ImageSetting(Phase, Linear))
