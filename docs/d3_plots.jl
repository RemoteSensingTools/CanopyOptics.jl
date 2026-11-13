# Zero-dependency Julia -> JSON/HTML serializer for interactive D3.js v7 figures
# in CanopyOptics.jl documentation (Documenter.jl + Literate.jl).

struct D3Figure
    spec::Dict{String,Any}
end

D3Figure(; kwargs...) = D3Figure(Dict{String,Any}(String(k) => v for (k, v) in kwargs))

d3_line(x, y; kwargs...) = (; kind = "line", x = collect(Float64, x), y = collect(Float64, y), kwargs...)
d3_scatter(x, y; kwargs...) = (; kind = "scatter", x = collect(Float64, x), y = collect(Float64, y), kwargs...)
d3_bar(y; kwargs...) = (; kind = "bar", y = collect(Float64, y), kwargs...)
d3_cartesian(; series, kwargs...) = (; type = "cartesian", series, kwargs...)
d3_heatmap(x, y, z; kwargs...) = (; type = "heatmap", x = collect(Float64, x), y = collect(Float64, y), z = Matrix{Float64}(z), kwargs...)
d3_polar(vza, vaz, z; kwargs...) = (; type = "polar", vza = collect(Float64, vza), vaz = collect(Float64, vaz), z = Matrix{Float64}(z), kwargs...)

function _d3_json_escape(io::IO, s::AbstractString)
    write(io, '"')
    for c in s
        if c == '"'
            write(io, "\\\"")
        elseif c == '\\'
            write(io, "\\\\")
        elseif c == '/'
            # Prevent "</script>" from closing the JSON script tag early
            write(io, "\\/")
        elseif c == '\n'
            write(io, "\\n")
        elseif c == '\r'
            write(io, "\\r")
        elseif c == '\t'
            write(io, "\\t")
        elseif UInt32(c) < 0x20
            write(io, "\\u", lpad(string(UInt32(c), base = 16), 4, '0'))
        else
            write(io, c)
        end
    end
    write(io, '"')
end

_d3_json(io::IO, ::Nothing) = write(io, "null")
_d3_json(io::IO, b::Bool) = write(io, b ? "true" : "false")
_d3_json(io::IO, x::Integer) = print(io, x)

function _d3_json(io::IO, x::Real)
    xf = Float64(x)
    if !isfinite(xf)
        write(io, "null")
    else
        print(io, round(xf; sigdigits = 6))
    end
end

_d3_json(io::IO, s::Symbol) = _d3_json_escape(io, String(s))
_d3_json(io::IO, s::AbstractString) = _d3_json_escape(io, s)

function _d3_json(io::IO, v::Union{AbstractVector,Tuple})
    write(io, '[')
    first_item = true
    for item in v
        first_item || write(io, ',')
        first_item = false
        _d3_json(io, item)
    end
    write(io, ']')
end

function _d3_json(io::IO, m::AbstractMatrix)
    # Emit row-major nested arrays [row1, row2, ...]
    write(io, '[')
    for r in axes(m, 1)
        r == first(axes(m, 1)) || write(io, ',')
        _d3_json(io, view(m, r, :))
    end
    write(io, ']')
end

function _d3_json(io::IO, nt::NamedTuple)
    write(io, '{')
    first_kv = true
    for k in keys(nt)
        val = getfield(nt, k)
        val === nothing && continue
        first_kv || write(io, ',')
        first_kv = false
        _d3_json_escape(io, String(k))
        write(io, ':')
        _d3_json(io, val)
    end
    write(io, '}')
end

function _d3_json(io::IO, d::AbstractDict)
    write(io, '{')
    first_kv = true
    for (k, val) in d
        val === nothing && continue
        first_kv || write(io, ',')
        first_kv = false
        _d3_json_escape(io, String(k))
        write(io, ':')
        _d3_json(io, val)
    end
    write(io, '}')
end

function Base.show(io::IO, ::MIME"text/html", fig::D3Figure)
    write(io, "<div class=\"canopy-d3-plot\"><script type=\"application/json\" class=\"canopy-plot-spec\">")
    _d3_json(io, fig.spec)
    write(io, "</script></div>")
end

struct D3RawFigure
    json::String
end

function Base.show(io::IO, ::MIME"text/html", fig::D3RawFigure)
    write(io, "<div class=\"canopy-d3-plot\"><script type=\"application/json\" class=\"canopy-plot-spec\">")
    write(io, replace(fig.json, "</" => "<\\/"))
    write(io, "</script></div>")
end

d3_from_file(path::AbstractString) = D3RawFigure(read(path, String))
