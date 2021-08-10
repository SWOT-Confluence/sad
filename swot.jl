### Run SAD algorithm with SWOT data

using Sad
using DelimitedFiles
using Distributions
using NCDatasets

const FILL = -999999999999

"""
Get reach file names.
"""
function get_reach_files(reachtxt)
    line = parse(Int64, ENV["AWS_BATCH_JOB_ARRAY_INDEX"]) + 1
    lines = readdlm(reachtxt, '\n', String, '\n')
    reach = lines[line]
    "$(reach)_SWOT.nc", "$(reach)_SOS.nc" 
end

"""
Load SWOT observations.
"""
function read_swot_obs(ncfile)
    ds = Dataset(ncfile)
    nodes = NCDatasets.group(ds, "node")
    S = nodes["slope2"][:]
    H = nodes["wse"][:]
    W = nodes["width"][:]
    H, W, S
end

"""
Write SAD output to NetCDF.
"""
function write_output(sosfile, outdir, A0, n, Qa, Qu)
    sos = Dataset(sosfile)
    reach_id = sos.attrib["reach_id"]
    valid = sos.attrib["valid"]    
    close(sos)
    
    outfile = joinpath(outdir, "$(reach_id)_sad.nc")
    out = Dataset(outfile, "c")
    out.attrib["reach_id"] = reach_id
    out.attrib["valid"] = valid   # FIXME Determine what is considered valid in the context of a SAD run
    defDim(out, "nt", size(Qa)[2])
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
Main driver routine.
"""
function main()
        indir = joinpath("/mnt", "data", "input")
        outdir = joinpath("/mnt", "data", "output")

        isempty(ARGS) ? reachfile = "reaches.txt" : reachfile = ARGS[1]
        swot, sos = get_reach_files(joinpath(indir, reachfile))
        swotfile = joinpath(indir, "swot", swot) # NetCDF file with SWOT observations
        swordfile = joinpath(indir, "sos", sos) # NetCDF file with SWORD database 

        H, W, S = read_swot_obs(swotfile)
        H, W, S, x = channel_chainage(H, W, S)
        Qₚ, nₚ, rₚ, zₚ = Sad.priors(swordfile, H)
        nens = 100 # default ensemble size
        A0, n, Qa, Qu = Sad.assimilate(H, W, x, maximum(W, dims=2), maximum(H, dims=2), S,
                            Qₚ, nₚ, rₚ, zₚ, nens, [1, length(x)])
        
        write_output(swordfile, outdir, A0, n, Qa, Qu)
    
end

main()
