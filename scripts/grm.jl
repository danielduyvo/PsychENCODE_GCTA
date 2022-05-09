using CSV, DataFrames, LinearAlgebra

function readgrmbin(prefix; alln = true, numtype = Float32)
    function sumi(i)::Integer
        return (i + 1) * i / 2
    end
    binfilename = prefix * ".grm.bin"
    nfilename = prefix * ".grm.N.bin"
    idfilename = prefix * ".grm.id"
    id = CSV.read(idfilename, DataFrame;
                  header = 0)
    n = size(id, 1)
    grm = Vector{numtype}(undef, sumi(n))
    read!(binfilename, grm)
    N = Vector{numtype}(undef, sumi(n))
    if alln
        read!(nfilename, N)
    else
        N = Vector{numtype}(undef, 1)
        read!(nfilename, N)
    end
    return (grm, N, id)
end

function writegrmbin(prefix, grm, N, id)
    function sumi(i)::Integer
        return (i + 1) * i / 2
    end
    binfilename = prefix * ".grm.bin"
    nfilename = prefix * ".grm.N.bin"
    idfilename = prefix * ".grm.id"
    id = CSV.write(idfilename, id;
                   delim = "\t", writeheader = false)
    write(binfilename, grm)
    write(nfilename, N)
    return nothing
end

function vectolowertri(vec)
    lvec = length(vec)
    n = convert(Int, 1 / 2 * (sqrt(8 * lvec + 1) - 1))
    @assert (n + 1) * n / 2 == lvec
    i = 1
    j = 1
    mat = Matrix{eltype(vec)}(undef, n, n)
    for el in vec
        mat[i, j] = el
        if j + 1 > i
            i += 1
            j = 1
        else
            j += 1
        end
    end
    return LowerTriangular(mat)
end

function lowertritovec(mat)
    @assert size(mat, 1) == size(mat, 2)
    n = size(mat, 1)
    vec = Vector{eltype(mat)}(undef, convert(Int, (n + 1) * n / 2))
    curr = 1
    for i in 1:n
        for j in 1:i
            vec[curr] = mat[i, j]
            curr += 1
        end
    end
    return vec
end

function generatesampletrimestervec(idfilename)
    tri1ids = filter(id -> id != "", split(read(idfilename, String), "\n"))
    sampletrimester = fill(2, size(id, 1))
    for i in 1:length(sampletrimester)
        if id[i, :Column2] in tri1ids
            sampletrimester[i] = 1
        end
    end
    return sampletrimester
end

function generatetri1_tri1tri2(prefix, sampletrimester::AbstractVector)
    grm, N, id = readgrmbin(prefix)
    grmmat = vectolowertri(grm)

    # Edit grm matrices
    tri1_grmmat = deepcopy(grmmat)
    for i in 1:size(tri1_grmmat, 1)
        for j in 1:i
            if sampletrimester[i] == 2 || sampletrimester[j] == 2
                tri1_grmmat[i, j] = 0
            end
        end
    end
    tri1tri2_grmmat = deepcopy(grmmat)
    for i in 1:size(tri1tri2_grmmat, 1)
        for j in 1:i
            if sampletrimester[i] != sampletrimester[j]
                tri1tri2_grmmat[i, j] = 0
            end
        end
    end
    # Convert back to grm vector and print out
    tri1_grmvec = lowertritovec(tri1_grmmat)
    tri1tri2_grmvec = lowertritovec(tri1tri2_grmmat)
    writegrmbin(prefix * "_tri1", tri1_grmvec, N, id)
    writegrmbin(prefix * "_tri1tri2", tri1tri2_grmvec, N, id)
    return nothing
end

