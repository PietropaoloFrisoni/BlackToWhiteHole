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

# computes a Wigner 4jm symbol on-the-fly efficiently for all outgoing links
@inline function Wigner4jm(j1::HalfInteger, j2::HalfInteger, j3::HalfInteger, j4::HalfInteger, m1, m2, m3, m4, i::HalfInteger)
    return (-1)^(i + m1 + m2) * wigner3j(j1, j2, i, m1, m2, -m1 - m2) * wigner3j(i, j3, j4, m1 + m2, m3, m4)
end

# computes a Wigner 4jm symbol on-the-fly with sum over m
@inline function Wigner4jm2(j1::HalfInteger, j2::HalfInteger, j3::HalfInteger, j4::HalfInteger, m1, m2, m3, m4, i::HalfInteger)
    w4j = 0.0
    for m = -i:i
        w4j += (-1)^(i - m) * wigner3j(j1, j2, i, m1, m2, m) * wigner3j(i, j3, j4, -m, m3, m4)
    end
    return w4j
end

# from index to intertwiner of tuple ((i_min, i_max), i_range)
@inline function from_index_to_intertwiner(tuple, index)
    return tuple[1][1] + index - 1
end

# logging function (flushing needed)
function log(x...)
    println("[ ", now(), " ] - ", join(x, " ")...)
    flush(stdout)
end

# release workers
function ReleaseWorkers(number_of_workers::Int64=nworkers())
    if (number_of_workers > 1)
        for i in workers()
            rmprocs(i)
        end
    end
end