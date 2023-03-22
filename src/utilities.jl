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

# retrieves from process
macro RetrieveFromProcess(p, obj, mod=:Main)
    quote
        remotecall_fetch($(esc(p)), $(esc(mod)), $(QuoteNode(obj))) do m, o
            Core.eval(m, o)
        end
    end
end

# computes a Wigner 4jm symbol on-the-fly
@inline function Wigner4jm(j1::HalfInteger, j2::HalfInteger, j3::HalfInteger, j4::HalfInteger, m1::HalfInteger, m2::HalfInteger, m3::HalfInteger, m4::HalfInteger, i::HalfInteger)
    return (-1)^(i + m1 + m2) * wigner3j(j1, j2, i, m1, m2, -m1 - m2) * wigner3j(i, j3, j4, m1 + m2, m3, m4)
end

# from index to intertwiner of tuple ((i_min, i_max), i_range)
@inline function from_index_to_intertwiner(tuple, index)
    return tuple[1][1] + index - 1
  end