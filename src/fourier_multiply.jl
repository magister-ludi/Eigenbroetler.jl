
const do_display = true

function cdisplay(a::Eigenbrot, msg = "")
    if do_display
        println(
            "***** ",
            height(a),
            "x",
            width(a),
            ", ",
            a.fft ? "Fourier" : "Real",
            " space",
            isempty(msg) ? "" : " | $msg",
        )
        println("Real part")
        for row = 1:height(a)
            for col = 1:width(a)
                if abs(real(a[row, col])) < 1e-6
                    @printf(" %+.5g", 0.0)
                else
                    @printf(" %+.5g", real(a[row, col]))
                end
            end
            println()
        end
        println("Imaginary part")
        for row = 1:height(a)
            for col = 1:width(a)
                if abs(imag(a[row, col])) < 1e-6
                    @printf(" %+.5g", 0.0)
                else
                    @printf(" %+.5g", imag(a[row, col]))
                end
            end
            println()
        end
    end
end

"""
    fft_multiply(ebL::Eigenbrot, ebR::Eigenbrot, do_cc::Bool)

Compute a pseudo-multiplication of `eb1` and `eb2`. Utility method
for convolution and correlation functions.
"""
function fft_multiply(ebL::Eigenbrot, ebR::Eigenbrot, do_cc::Bool)
    fftL = fft(ebL)
    fftR = fft(ebR)
    if do_cc
        fftR = conj(fftR)
    end
    fftprod = fftL .* fftR
    return fft(fftprod)
end

"""
    cor(eb1::Eigenbrot, eb2::Eigenbrot)

Compute the correlation of `eb1` and `eb2`. Uses FFT algorithm.
"""
Base.cor(eb1::Eigenbrot, eb2::Eigenbrot) = fft_multiply(eb1, eb2, true)

"""
    cor(eb::Eigenbrot)

Compute the autocorrelation of `eb`. Uses FFT algorithm.
"""
Base.cor(eb::Eigenbrot) = cor(eb, eb)

"""
    conv(eb1::Eigenbrot, eb2::Eigenbrot)

Compute the convolution of `eb1` and `eb2`. Uses FFT algorithm.
"""
Base.conv(eb1::Eigenbrot, eb2::Eigenbrot) = fft_multiply(eb1, eb2, false)
