
function fft(eb::Eigenbrot, recentre::Bool = true)
    const w = width(eb)
    const h = height(eb)
    const scale = 1.0 / sqrt(h * w)
    data = recentre ? fftshift(eb.vals) : copy(eb.vals)
    fft!(data)
    trf = Eigenbrot(scale * (recentre ? fftshift(data) : data))
    trf.fft = !eb.fft
    return trf
end

function fftx(eb::Eigenbrot, recentre::Bool = true)
    const w = width(eb)
    const scale = 1.0 / sqrt(w)
    data = recentre ? fftshift(eb.vals, 2) : copy(eb.vals)
    fft!(data, 2)
    return Eigenbrot(scale * (recentre ? fftshift(fftdata, 2) : fftdata))
end

function ffty(eb::Eigenbrot, recentre::Bool = true)
    const h = height(eb)
    const scale = 1.0 / sqrt(h)
    data = recentre ? fftshift(eb.vals, 1) : copy(eb.vals)
    fft!(data, 1)
    return Eigenbrot(scale * (recentre ? fftshift(data, 1) : data))
end
