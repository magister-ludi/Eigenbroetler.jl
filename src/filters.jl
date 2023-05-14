
@enum FilterType Butterworth Gaussian Hamming Hann Parzen Square Welch

function lowpass(
    eb::Eigenbrot,
    filter::FilterType,
    r::Number;
    x0::Number = 0,
    y0::Number = 0,
    order::Int = 1,
)
    func = if filter == Butterworth
        filterLowButterworth
    elseif filter == Gaussian
        filterLowGaussian
    elseif filter == Hamming
        filterLowHamming
    elseif filter == Hann
        filterLowHann
    elseif filter == Parzen
        filterLowParzen
    elseif filter == Square
        filterLowSquare
    else  # filter == Welch
        filterLowWelch
    end
    func(eb, r, x0, y0, order)
end

function highpass(
    eb::Eigenbrot,
    filter::FilterType,
    r::Number;
    x0::Number = 0,
    y0::Number = 0,
    order::Int = 1,
)
    func = if filter == Butterworth
        filterHighButterworth
    elseif filter == Gaussian
        filterHighGaussian
    elseif filter == Hamming
        filterHighHamming
    elseif filter == Hann
        filterHighHann
    elseif filter == Parzen
        filterHighParzen
    elseif filter == Square
        filterHighSquare
    else  # filter == Welch
        filterHighWelch
    end
    func(eb, r, x0, y0, order)
end

function filterLowButterworth(eb::Eigenbrot, r::Number, x0::Number, y0::Number, order::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filterFactor = 1.0 / (1.0 + ((colSq + rowSq) / r2)^order)
            filtered[row, col] = filterFactor * eb[row, col]
        end
    end
    return filtered
end

function filterLowGaussian(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = exp(-(colSq + rowSq) / r2) * eb[row, col]
        end
    end
    return filtered
end

function filterLowHamming(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = if colSq + rowSq > r2
                0.0
            else
                (0.54 - 0.46 * cos(π * (r - sqrt(colSq + rowSq)) / r)) * eb[row, col]
            end
        end
    end
    return filtered
end

function filterLowHann(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = if colSq + rowSq > r2
                0.0
            else
                0.5 * (1.0 - cos(π * (r - sqrt(colSq + rowSq)) / r)) * eb[row, col]
            end
        end
    end
    return filtered
end

function filterLowParzen(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq > r2
                filtered[row, col] = 0
            else
                filterFactor = (r - sqrt(colSq + rowSq)) / r
                filtered[row, col] = filterFactor * eb[row, col]
            end
        end
    end
    return filtered
end

function filterLowSquare(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = copy(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq > r2
                filtered[row, col] = 0
            end
        end
    end
    return filtered
end

function filterLowWelch(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq > r2
                filtered[row, col] = 0
            else
                filterFactor = (r - sqrt(colSq + rowSq)) / r
                filtered[row, col] = filterFactor * filterFactor * eb[row, col]
            end
        end
    end
    return filtered
end

function filterHighButterworth(eb::Eigenbrot, r::Number, x0::Number, y0::Number, order::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filterFactor = 1.0 / (1.0 + (r2 / (colSq + rowSq))^order)
            filtered[row, col] = filterFactor * eb[row, col]
        end
    end
    return filtered
end

function filterHighGaussian(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = (1.0 - exp(-(colSq + rowSq) / r2)) * eb[row, col]
        end
    end
    return filtered
end

function filterHighHamming(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = if colSq + rowSq <= r2
                (0.54 - 0.46 * cos(π * sqrt(colSq + rowSq) / r)) * eb[row, col]
            else
                eb[row, col]
            end
        end
    end
    return filtered
end

function filterHighHann(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            filtered[row, col] = if colSq + rowSq <= r2
                0.5 * (1.0 - cos(π * sqrt(colSq + rowSq) / r)) * eb[row, col]
            else
                eb[row, col]
            end
        end
    end
    return filtered
end

function filterHighParzen(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = copy(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq > r2
                filtered[row, col] = eb[row, col] * sqrt(colSq + rowSq) / r
            end
        end
    end
    return filtered
end

function filterHighSquare(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = copy(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq < r2
                filtered[row, col] = 0
            end
        end
    end
    return filtered
end

function filterHighWelch(eb::Eigenbrot, r::Number, x0::Number, y0::Number, ::Int)
    filtered = similar(eb)
    r2 = r * r
    h, w = size(eb)
    centreCol = w >> 1 + x0
    centreRow = h >> 1 + y0
    for col = 1:w
        colSq = (centreCol - col)^2
        for row = 1:h
            rowSq = (centreRow - row)^2
            if colSq + rowSq <= r2
                filtered[row, col] = eb[row, col] * (colSq + rowSq) / r2
            end
        end
    end
    return filtered
end
