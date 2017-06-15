#!/usr/bin/env julia

using Eigenbroetler

function magic_phase(x::Number, y::Number)
    r = sqrt(x * x + y * y)
    return exp(-(r / 100) ^ 2) * exp(2 * im * x * (x ^ 2 - y ^ 2) / 200 ^ 2)
end

eb = Eigenbrot(512, 512, magic_phase)
save("magic_phase.png", eb, ImageSetting(Phase, Linear))
save("magic_phase_dft.png", fft(eb), ImageSetting(Phase, Linear))
