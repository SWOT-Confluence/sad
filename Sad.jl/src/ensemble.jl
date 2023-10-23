import Distributions.quantile
import BlackBoxOptim.Utils.latin_hypercube_sampling

"""
Get quantiles from `p` according to weights provided by
Latin Hypercube Sampling `samples`.

"""
function get_samples(p, samples)
    out = try
        quantile.(p, samples)
    catch
        [quantile(p, s) for s in samples]
    end
    out
end

"""
Generate an ensemble of size `nens` using Latin Hypercube Sampling of the
list of distributions or collections provided as `args`.

"""
function lhs_ensemble(nens, args...)
    nvars = length(args)
    usamples = latin_hypercube_sampling([0. for _ in 1:nvars], [1. for _ in 1:nvars], nens)
    [get_samples(args[i], usamples[i, :]) for i in 1:nvars]
end
    
"""
Generate a prior ensemble of discharge, roughness coefficient, channel shape parameter, 
and downstream bed elevation from provided distributions or sample data.

"""
function prior_ensemble(x, Qₚ, nₚ, rₚ, zₚ, nens)
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    Qe, ne, re, ze
end

"""
Generate an ensemble of water height profiles from Gradually-Varied-Flow simulations 
and associated profiles of bed elevation.

- `H`: water surface elevation
- `W`: water width
- `S`: water surface slope
- `x`: channel chainage
- `hbf`: bankfull water surface elevation
- `wbf`: bankfull width
- `Qe`: ensemble discharge
- `ne`: ensemble roughness coefficient
- `re`: ensemble channel shape parameter
- `ze`: ensemble bed elevation profiles
- `perturb_slope`: treat water surface slope stochastically (default is false)

"""
function gvf_ensemble!(H, S, x, hbf, wbf, Qe, ne, re, ze)
    nens = length(Qe)
    if ndims(S) > 1
        Se = S
    else
        Se = repeat(S', outer=nens)'
    end
    for j in 2:length(x)
        ze[j, :] = ze[j-1, :] .+ Se[j, :] .* (x[j] - x[j-1]);
    end
    ybf = hbf .- ze
    he = zeros(length(x), nens)
    for i in 1:nens
        he[:, i] = gvf(Qe[i], (H[1]-ze[1, i])*re[i]/(re[i] + 1), Se[:, i], ne[i], x, wbf, ybf[:, i], [re[i] for _ in 1:length(x)])
    end
    he
end
