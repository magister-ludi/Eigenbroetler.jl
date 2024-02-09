
"""
    fft_multiply(ebL::Eigenbrot, ebR::Eigenbrot, do_cc::Bool)

Compute a pseudo-multiplication of `eb1` and `eb2`. Utility method
for convolution and correlation functions.
"""
function fft_multiply(ebL::Eigenbrot, ebR::Eigenbrot, do_cc::Bool)
    wL, wR = width(ebL), width(ebR)
    hL, hR = height(ebL), height(ebR)
    rL = lL = tL = bL = 0
    rR = lR = tR = bR = 0
    δw = wL - wR
    δh = hL - hR
    if δw > 0
        lR = δw >> 1
        rR = δw - lR
    elseif δw < 0
        lL = δw >> 1
        rL = δw - lR
    end
    if δh > 0
        tR = δh >> 1
        bR = δh - tR
    elseif δh < 0
        tL = δh >> 1
        bL = δh - tR
    end
    if any(!=(0), (rL, lL, tL, bL))
        ebL = pad(ebL, right = rL, left = lL, top = tL, bottom = bL)
    end
    if any(!=(0), (rR, lR, tR, bR))
        ebR = pad(ebR, right = rR, left = lR, top = tR, bottom = bR)
    end
    fftL = fft(ebL)
    fftR = fft(ebR)
    if do_cc
        fftR = conj(fftR)
    end
    fftprod = fftL * fftR
    return fft(fftprod)
end

"""
    cor(eb1::Eigenbrot, eb2::Eigenbrot)

Compute the correlation of `eb1` and `eb2`. Uses FFT algorithm.
"""
Statistics.cor(eb1::Eigenbrot, eb2::Eigenbrot) = fft_multiply(eb1, eb2, true)

"""
    cor(eb::Eigenbrot)

Compute the autocorrelation of `eb`. Uses FFT algorithm.
"""
Statistics.cor(eb::Eigenbrot) = cor(eb, eb)

"""
    conv(eb1::Eigenbrot, eb2::Eigenbrot)

Compute the convolution of `eb1` and `eb2`. Uses FFT algorithm.
"""
conv(eb1::Eigenbrot, eb2::Eigenbrot) = fft_multiply(eb1, eb2, false)
