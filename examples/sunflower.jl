#!/usr/bin/env julia

using Eigenbroetler

function sunflower(x::Number, y::Number)
    if x == 0 && y == 0
        return 0.0
    end
    r = sqrt(x * x + y * y)
    phi = atan2(y, x)
    return cos(30 * log(r) + 13 * phi) ^ 2 + cos(20 * log(r) - 11 * phi) ^ 2
end

eb = Eigenbrot(sunflower, 512, 512)
save("sunflower.png", eb, ImageSetting(RealPart, Linear), colours = :rainbow)
