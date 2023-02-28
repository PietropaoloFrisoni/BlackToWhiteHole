# dimension of spin representations
@inline dimension(x::Int64) = 2x
@inline dimension(x::HalfInteger) = twice(x)
@inline dimension(x::Float64) = round(Int64, 2x)

@inline Dimension(x::Int64) = 2x + 1
@inline Dimension(x::HalfInteger) = twice(x) + 1
@inline Dimension(x::Float64) = round(Int64, 2x + 1)

# returns number of files in a folder
function file_count(folder::String)
    files_and_dirs = readdir(folder)
    number_of_files = size(files_and_dirs)[1]
    return number_of_files
end

# allowed intertwiner range for spins (j1, j2, j3, j4) in the recoupling (j1, j2) -- (j3, j4)."
function intertwiner_range_float(j1::Float64, j2::Float64, j3::Float64, j4::Float64)

    if (2j1 + 2j2) % 2 != (2j3 + 2j4) % 2
        return (), 0
    end

    imin = max(abs(j1 - j2), abs(j3 - j4))
    imax = min(j1 + j2, j3 + j4)

    range = imax >= imin ? (imin, imax) : ()
    size = imax >= imin ? Int(imax - imin) + 1 : 0

    (range, size)

end

# retrieves from process
macro retrieve_from_process(p, obj, mod=:Main)
    quote
        remotecall_fetch($(esc(p)), $(esc(mod)), $(QuoteNode(obj))) do m, o
            Core.eval(m, o)
        end
    end
end