
for name in (:real, :imag, :angle, :abs, :conj)
    @eval begin
        import Base.$name
        $name(eb::Eigenbrot) = Eigenbrot(Array{Complex128}($name(eb.vals)))
    end
end

phase(eb::Eigenbrot) = angle(eb)

import Base: (+), (.+), (-), (.-), (*), (.*), (/), (./), (.^)

(+)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .+ eb2.vals)
(+)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals .+ n)
(+)(n::Number, eb1::Eigenbrot) = eb1 .+ n
(.+)(eb1::Eigenbrot, eb2) = eb1 + eb2
(.+)(eb1, eb2::Eigenbrot) = eb1 + eb2
(-)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .- eb2.vals)
(-)(eb1::Eigenbrot, n::Number) = Eigenbrot(eb1.vals .- n)
(-)(n::Number, eb1::Eigenbrot) = eb1 .- n
(.-)(eb1::Eigenbrot, eb2) = eb1 - eb2
(.-)(eb1::Eigenbrot, eb2::Eigenbrot) = eb1 - eb2

(.*)(eb1::Eigenbrot, eb2::Eigenbrot) = Eigenbrot(eb1.vals .* eb2.vals)
(*)(eb::Eigenbrot, x::Number) = x * eb
(*)(x::Number, eb::Eigenbrot) = Eigenbrot(x .* eb.vals)
(/)(eb::Eigenbrot, x::Number) = Eigenbrot(eb.vals ./ x)
(./)(x::Number, eb::Eigenbrot) = Eigenbrot(x ./ eb.vals)
(.^)(eb::Eigenbrot, x::Number) = Eigenbrot(eb.vals .^ x)
(.^)(x::Number, eb::Eigenbrot) = Eigenbrot(x .^ eb.vals)
