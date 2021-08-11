### Run SAD algorithm with SWOT data

using Sad
using DelimitedFiles
using Distributions
using NCDatasets
using JSON

const FILL = -999999999999

"""
Get reach file names.
"""
function get_reach_files(indir, reachjson)
    line = parse(Int64, ENV["AWS_BATCH_JOB_ARRAY_INDEX"]) + 1
    open(joinpath(indir, reachjson)) do jf
        data = read(jf, String)
        reachlist = JSON.parse(data)[line]
        reachlist["reach_id"], joinpath(indir, "swot", reachlist["swot"]), joinpath(indir, "sos", reachlist["sos"])
    end
end

"""
Load SWOT observations.
"""
function read_swot_obs(ncfile)
    ds = NCDataset(ncfile)
    nodes = NCDatasets.group(ds, "node")
    # S = nodes["slope2"][:]
    S = Array{Union{Float64, Missing}}(nodes["slope2"][:])
    S[S .=== missing] .= 0.0
    # S[S .≈ 9.969209968386869e36] .= missing
    # H = nodes["wse"][:]
    H = Array{Union{Float64, Missing}}(nodes["wse"][:])
    H[H .=== missing] .= 0.0
    # H[H .≈ 9.969209968386869e36] .= missing
    # W = nodes["width"][:]
    W = Array{Union{Float64, Missing}}(nodes["width"][:])
    W[W .=== missing] .= 0.0
    # W[W .≈ 9.969209968386869e36] .= missing
    close(ds)
    H, W, S
end

"""
Extract river channel chainage.
"""
function channel_chainage(H, W, S)
    # sort nodes from downstream to upstream
    # FIXME: currently uses the mean water surface elevation
    hm = [mean(skipmissing(H[:, c])) for c in 1:size(H, 2)]
    i = sortperm(hm)
    H = H[:, i]'
    W = W[:, i]'
    S = S[:, i]'
    j = [v[1] for v in findall(all(!ismissing, H, dims=2))]
    H = convert(Array{Float64, 2}, H[j, :])
    W = convert(Array{Float64, 2}, W[j, :])
    S = convert(Array{Float64, 2}, S[j, :])
    # calculate distance from downstream node
    # FIXME: CSV with geographical information does not correspond directly to the nodes in each SWOT
    # orbit pass file. We can't use the calculated slope do derive the distance between
    # cross-sections as there are rather implausible discontinuities in the test data.
    # For now, just assume that cross-sections are 100 m apart
    x = (collect(0:size(H, 1)) * 100.0)[1:end-1]
    H, W, S, x
end

"""
Write SAD output to NetCDF.
"""
function write_output(reachid, valid, outdir, A0, n, Qa, Qu)
    outfile = joinpath(outdir, "$(reachid)_sad.nc")
    out = Dataset(outfile, "c")
    out.attrib["valid"] = valid   # FIXME Determine what is considered valid in the context of a SAD run
    defDim(out, "nt", size(Qa,2))
    ridv = defVar(out, "reach_id", Int64, (), fillvalue = FILL)
    ridv[:] = reachid
    A0v = defVar(out, "A0", Float64, (), fillvalue = FILL)
    A0v[:] = A0
    nv = defVar(out, "n", Float64, (), fillvalue = FILL)
    nv[:] = n
    Qav = defVar(out, "Qa", Float64, ("nt",), fillvalue = FILL)
    Qav[:,:] = Qa
    Quv = defVar(out, "Q_u", Float64, ("nt",), fillvalue = FILL)
    Quv[:,:] = Qu
    close(out)
end

"""
Main driver routine.
"""
function main()
    indir = joinpath("/mnt", "data", "input")
    outdir = joinpath("/mnt", "data", "output")

    isempty(ARGS) ? reachfile = "reaches.json" : reachfile = ARGS[1]
    reachid, swotfile, swordfile = get_reach_files(indir, reachfile)

    H, W, S = read_swot_obs(swotfile)
    A0 = missing
    n = missing
    Qa = Array{Missing}(missing, 1, size(W,2))
    Qu = Array{Missing}(missing, 1, size(W,2))
    if all(ismissing, H) || all(ismissing, W) || all(ismissing, S)
        println("$(reachid): INVALID")
        write_output(reachid, 0, outdir, A0, n, Qa, Qu)
    else
        H, W, S, x = channel_chainage(H, W, S)
        Qₚ, nₚ, rₚ, zₚ = Sad.priors(swordfile, H, reachid)
        if ismissing(Qₚ) || ismissing(nₚ) || ismissing(rₚ) || ismissing(zₚ)
            println("$(reachid): INVALID")
            write_output(reachid, 0, outdir, A0, n, Qa, Qu)
        else
            try
                nens = 100 # default ensemble size
                A0, n, Qa, Qu = Sad.assimilate(H, W, x, maximum(W, dims=2), maximum(H, dims=2), S,
                                    Qₚ, nₚ, rₚ, zₚ, nens, [1, length(x)])
                println("$(reachid): VALID")
                write_output(reachid, 1, outdir, A0, n, Qa, Qu)
            catch e
                println("$(reachid): INVALID")
                println(e)
                write_output(reachid, 0, outdir, A0, n, Qa, Qu)
            end
        end
    end
end

main()