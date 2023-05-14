using Eigenbroetler
using Test

random(rows::Integer, cols::Integer) = Eigenbrot(rand(ComplexF64, rows, cols), rand(Bool))

function sample(x, y)
    θ = angle(x + y * im)
    r = hypot(x, y)
    x^2 - y^2 - im * (r * cos(θ))
end

@testset "Eigenbroetler" begin
    @testset "Constructors" begin
        eb = Eigenbrot(512, 512, rand(Bool)) do x, y
            if x == 0 && y == 0
                return 0.0
            end
            r = sqrt(x * x + y * y)
            ϕ = atan(y, x)
            return cos(30 * log(r) + 13 * ϕ)^2 + cos(20 * log(r) - 11 * ϕ)^2
        end
        @test height(eb) == 512
        @test width(eb) == 512
        @test size(eb, 1) == 512
        @test size(eb, 2) == 512
        @test length(eb) == 512 * 512
        @test ndims(eb) == 2
        @test size(eb) == (512, 512)
    end

    @testset "Files" begin
        eb = Eigenbrot(sample, 64, 64)
        testfile = tempname() * ".png"
        save(testfile, eb, ImageSetting(RealPart, LinearScale))
        pngmagic = [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
        t = open(testfile)
        by = read(t, length(pngmagic))
        close(t)
        @test by[1:length(pngmagic)] == pngmagic
    end

    @testset "Basics" begin
        eb = random(256, 512)
        @test height(eb) == 256
        @test width(eb) == 512
        @test size(eb, 1) == 256
        @test size(eb, 2) == 512
        @test length(eb) == 256 * 512
        @test ndims(eb) == 2
        @test size(eb) == (256, 512)
        testfile = tempname()
        save(testfile, eb)
        @test ispath(testfile)
        ebr = Eigenbrot(testfile)
        @test eb.vals ≈ ebr.vals
        @test eb.fft == ebr.fft
        rm(testfile, force = false)
    end

    @testset "Indexing" begin
        eb = Eigenbrot(50, 50)
        for c = 1:50
            for r = 1:50
                val = rand(ComplexF64)
                eb[r, c] = val
                @test val == eb[r, c]
            end
        end
        for c = 1:10
            for r = 1:10
                xy = x, y = coords(eb, r, c)
                @test pixel(eb, x, y) == (r, c)
                @test pixel(eb, xy) == (r, c)
                @test xy == coords(eb, (r, c))
            end
        end
    end

    @testset "Simple methods" begin
        eb = random(8, 8)
        ans = real(eb)
        @test real(eb.vals) ≈ ans.vals
        @test eb.fft == ans.fft
        ans = imag(eb)
        @test imag(eb.vals) ≈ ans.vals
        @test eb.fft == ans.fft
        ans = angle(eb)
        @test angle.(eb.vals) ≈ ans.vals
        @test eb.fft == ans.fft
        ans = abs(eb)
        @test abs.(eb.vals) ≈ ans.vals
        @test eb.fft == ans.fft
        ans = conj(eb)
        @test conj(eb.vals) ≈ ans.vals
        @test eb.fft == ans.fft
    end

    @testset "Arithmetic" begin
        a = random(8, 8)
        b = random(8, 8)
        ans = 1 + a
        @test ans.vals ≈ 1 .+ a.vals
        @test ans.fft == a.fft
        ans = a + 1
        @test ans.vals ≈ 1 .+ a.vals
        @test ans.fft == a.fft

        ans = 1 - a
        @test ans.vals ≈ 1 .- a.vals
        @test ans.fft == a.fft
        ans = a - 1
        @test ans.vals ≈ a.vals .- 1
        @test ans.fft == a.fft

        ans = 2 * a
        @test ans.vals ≈ 2 * a.vals
        @test ans.fft == a.fft
        ans = a * 2
        @test ans.vals ≈ 2 * a.vals
        @test ans.fft == a.fft
        ans = a / 2
        @test ans.vals ≈ a.vals / 2
        @test ans.fft == a.fft

        ans = a * b
        @test ans.vals ≈ a.vals .* b.vals
        @test ans.fft == a.fft
        ans = b * a
        @test ans.vals ≈ a.vals .* b.vals
        @test ans.fft == b.fft

        ans = a / b
        @test ans.vals ≈ a.vals ./ b.vals
        @test ans.fft == a.fft
        ans = b / a
        @test ans.vals ≈ b.vals ./ a.vals
        @test ans.fft == b.fft
    end

    @testset "Sample with padding" begin
        eb = Eigenbrot(sample, 64, 64)
        ebf = Eigenbrot(joinpath(@__DIR__, "sample.fit"))
        @test eb.vals ≈ ebf.vals
        @test eb.fft == ebf.fft

        l, t, r, b = 5, 8, 3, 7
        padvalue = rand(ComplexF64)
        padded = pad(eb, padvalue, left = l, top = t, right = r, bottom = b)
        @test size(padded) == size(eb) .+ (t + b, l + r)
        @test all(padded.vals[1:b, :] .== padvalue)
        @test all(padded.vals[(end - t + 1):end, :] .== padvalue)
        @test all(padded.vals[:, 1:l] .== padvalue)
        @test all(padded.vals[:, (end - r + 1):end] .== padvalue)
        @test padded.vals[(b + 1):(end - t), (l + 1):(end - r)] == eb.vals
    end

    @testset "Fourier" begin
        eb = random(50, 50)
        for ft in [fftx, ffty, fft]
            ebfft = ft(eb)
            ebfft2 = ft(ebfft)
            Δ = eb.vals .- ebfft2.vals
            @test eb.vals ≈ ebfft2.vals
            @test eb.fft == ebfft2.fft
        end
        nodc = removeDC(eb)
        nodcff = fft(nodc)
        @test isapprox(nodcff[pixel(nodcff, 0, 0)], 0, atol = 1e-15)

        xscale = 3 // 2
        yscale = 5 // 4
        eb = random(8, 8)
        ebs = scale_chirpz(eb, xscale, yscale)
        @test size(ebs) == (10, 12)
        ebs = scale_x_chirpz(eb, xscale)
        @test size(ebs) == (8, 12)
        ebs = scale_y_chirpz(eb, yscale)
        @test size(ebs) == (10, 8)
    end
end
