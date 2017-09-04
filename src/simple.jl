
for name in (:real, :imag, :angle, :abs, :conj)
    @eval begin
        import Base.$name
        $name(eb::Eigenbrot) = Eigenbrot(Array{Complex128}($name(eb.vals)))
    end
end

phase(eb::Eigenbrot) = angle(eb)

import Base: (+), (-), (*), (/), (^)

(+)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .+ eb2.vals)
(+)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals + n)
(+)(n::Number, eb1::Eigenbrot) = eb1 + n
(-)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .- eb2.vals)
(-)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals .- n)
(-)(n::Number, eb1::Eigenbrot) = Eigenbrot(n .- eb1.vals)

(*)(eb::Eigenbrot, x::Number) = Eigenbrot(x * eb)
(*)(x::Number, eb::Eigenbrot) = Eigenbrot(x .* eb.vals)
(/)(eb::Eigenbrot, x::Number) = Eigenbrot(eb.vals ./ x)

Base.broadcast(::typeof(*), eb1::Eigenbrot, eb2::Eigenbrot) =
    Eigenbrot(eb1.vals .* eb2.vals)

Base.broadcast(::typeof(/), eb1::Eigenbrot, eb2::Eigenbrot) =
    Eigenbrot(eb1.vals ./ eb2.vals)
