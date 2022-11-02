using NCDatasets
using Distributions

# type alias
const SMatrix = Adjoint{Union{Missing, Float64}, Matrix{Union{Missing, Float64}}}

"""
Read prior information from SWORD from `ncfile` NetCDF file.
Requires prior information on the minimum observed water surface elevation (`minH`).
"""
function read_sword_data(ncfile::String, minH::Float64, reachid::Int64)
    ds = Dataset(ncfile)
    rg = NCDatasets.group(ds, "reaches")
    reachids = rg["reach_id"][:]
    index = findall(x -> x==reachid, reachids)
    
    mg = NCDatasets.group(ds, "model")
    Qm = mg["mean_q"][index][1] 
    if ismissing(Qm)
        close(ds)
        Qₚ, nₚ, rₚ, zₚ = missing, missing, missing, missing
    else
        Qub = ismissing(mg["max_q"][index][1]) ? Qm / 10 : mg["max_q"][index][1]
        Qlb = ismissing(mg["min_q"][index][1]) ? 10 * Qm : mg["min_q"][index][1]
        
        gg = NCDatasets.group(ds, "gbpriors").group["reach"]
        Qc = ismissing(gg["logQ_sd"][index][1]) ? rand(Uniform(0.5,2.0)) : exp(gg["logQ_sd"][index][1]) / Qm
        Qₚ = Truncated(LogNormal(log(Qm/sqrt(Qc^2)), log(1+Qc^2)), Qlb, Qub)

        n = [exp(gg[var][index][1]) for var in ["logn_hat", "logn_sd", "upperbound_logn", "lowerbound_logn"]]
        nₚ = all(ismissing, n) ? Uniform(0.01, 0.07) : Truncated(Normal(n[1], n[2]), n[4], n[3])

        r = [exp(gg[var][index][1]) for var in ["logr_hat", "logr_sd", "upperbound_logr", "lowerbound_logr"]]
        rₚ = all(ismissing, r) ? Truncated(Normal(2.5, 1.0), 0.8, 20) : Truncated(Normal(r[1], r[2]), r[4], r[3])
        
        z = [exp(gg[var][index][1]) for var in ["logDb_hat", "logDb_sd", "upperbound_logDb", "lowerbound_logDb"]]
        zₚ = all(ismissing, z) ? Uniform(minH-20, minH) : Truncated(Normal(minH-z[1], z[2]), minH-z[3], minH-z[4])

        close(ds)
    end
    Qₚ, nₚ, rₚ, zₚ
end

"""
Derive prior probability distributions from database.
"""
function priors(ncfile::String, H::SMatrix, reachid::Int64)
    minH = minimum(skipmissing(H[1, :]))
    Qₚ, nₚ, rₚ, zₚ = read_sword_data(ncfile, minH, reachid)
    Qₚ, nₚ, rₚ, zₚ
end
