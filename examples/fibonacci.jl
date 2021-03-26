#!/usr/bin/env julia

using Eigenbroetler

function fibonacci(x::Number, y::Number)
    # cos(30*log(r)+13*(ϕ))**2+cos(20*log(r)-11*(ϕ))**2
    if x == 0 && y == 0
        return 0.0
    end
    r = sqrt(x * x + y * y)
    ϕ = atan(y, x)
    return cos(30.0 * log(r) + 13.0 * ϕ) ^ 2 + cos(20.0 * log(r) - 11.0 * ϕ) ^ 2
end

eb = Eigenbrot(fibonacci, 512, 512)
save("fibonacci.png", eb, ImageSetting(RealPart, LinearScale))
