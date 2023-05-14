
const MAX_COLOUR = 255
const NUM_COLOURS = MAX_COLOUR + 1

struct Palette
    p::Vector{RGB{N0f8}}
end

const palettes = Dict{Symbol, Palette}()

Base.getindex(pal::Palette, idx::Integer) = pal.p[idx]

"""
    readpalette(stream::IO)
Read palette data from `stream`. The data should be
formatted as up to $NUM_COLOURS lines with three integers
each, separated by white space. The integers must fit into
a `UInt8`; they represent the RGB values of a 24-bit colour.
'#' marks the beginning of a comment which runs to the
end of a line and is ignored.
"""
function readpalette(stream::IO)
    i = 0
    pal = Vector{RGB{N0f8}}(undef, NUM_COLOURS)
    for line in readlines(stream)
        line = replace(strip(line), r"\s*#.*" => "")
        if isempty(line)
            continue
        end
        i += 1
        if i > NUM_COLOURS
            break
        end
        r, g, b = split(line)
        pal[i] = RGB{N0f8}(
            parse(Float32, r) / MAX_COLOUR,
            parse(Float32, g) / MAX_COLOUR,
            parse(Float32, b) / MAX_COLOUR,
        )
    end
    while i < NUM_COLOURS
        i += 1
        pal[i] = RGB{N0f8}(0, 0, 0)
    end
    return Palette(pal)
end

"""
    palette(filename::AbstractString)
Reads and returns a `Palette` from `filename`.
    palette(name::Symbol)
Reads and returns a palette from a user supplied file
or from a standard palette file. Defaults to a palette
of $NUM_COLOURS greys.
"""
palette(filename::AbstractString) = palette(Symbol(filename))

function palette(name::Symbol)
    paths = ["$name", joinpath(dirname(@__FILE__), "maps", "$name.map")]
    if haskey(palettes, name)
        return palettes[name]
    elseif any(ispath, paths)
        path = paths[findfirst(ispath, paths)]
        if !haskey(palettes, name)
            pal = readpalette(open(path))
            palettes[name] = pal
        end
        return palettes[name]
    else
        if !haskey(palettes, :grey)
            pal = Vector{RGB{N0f8}}(undef, NUM_COLOURS)
            for i = 1:NUM_COLOURS
                g = (i - 1) / MAX_COLOUR
                pal[i] = RGB{N0f8}(g, g, g)
            end
            palettes[:gray] = palettes[:grey] = Palette(pal)
        end
        return palettes[:grey]
    end
end
