
"""
    hilbertX(eb::Eigenbrot)
Return the Hilbert x-transform of `eb`.
"""
function hilbertX(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c in 1:w
        for r in 1:h
            sign = if c == centre_c
                r > centre_r ? im : -im
            elseif r == centre_r
                c > centre_c ? im : -im
            else
                c > centre_c ? im : -im
            end
            nd[r, c] *= sign
        end
    end
    return nd
end

"""
    hilbertY(eb::Eigenbrot)
Return the Hilbert x-transform of `eb`.
"""
function hilbertY(eb::Eigenbrot)
    nd = copy(eb)
    h, w = size(nd)
    centre_c = 1 + w >> 1
    centre_r = 1 + h >> 1
    for c in 1:w
        for r in 1:h
            sign = if c == centre_c
                r > centre_r ? im : -im
            elseif r == centre_r
                c > centre_c ? im : -im
            else
                r > centre_r ? im : -im
            end
            nd[r, c] *= sign
        end
    end
    return nd
end
