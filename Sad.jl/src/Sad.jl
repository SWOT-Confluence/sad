module Sad

using Distributions
using KernelDensity: kde

abstract type CrossSection end

include("crosssections.jl")
include("kalman.jl")
include("gvf.jl")
include("ensemble.jl")
include("priors.jl")

# Define some helper functions
smean(x) = mean(skipmissing(x))

"""
Estimate bed slope by assimilating observed water surface elevations.

"""
function bed_slope(S, H, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens; ϵₒ=0.01)
    S0 = [smean(S[j, :]) > 0 ? smean(S[j, :]) : smean(S[j, :][S[j,:] .> 0]) for j in 1:size(S, 1)]
    So = try
        S0[isnan.(S0)] .= minimum(S[S .> 0])
        Se = rand.(Normal(1., 0.2), length(x), nens)
        Se = [Se[i, j] .* S0[i] for i in 1:length(x), j in 1:nens]
        ze = zeros(length(x), nens)
        Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
        he = gvf_ensemble!(mean(H, dims=2), Se, x, hbf, wbf, Qe, ne, re, ze)
        i = findall(he[1, :] .> 0)
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = Se[:, i]
        XA = h[:, i]
        d = mean(H, dims=2)[:, 1]
        E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = letkf(X, d, XA, E, [collect(1:length(x))], [collect(1:length(d))], diagR=true)
        So = mean(A, dims=2)[:, 1]
        So[So .< 0] .= minimum(So[So .> 0])
        # REVIEW: Put some checks here to capture quality control issues with the observation data    
        zz = zeros(length(x))
        zz[1] = params(zₚ)[4]
        for j in 2:length(x)
            zup = zz[j-1] .+ So[j] .* (x[j] - x[j-1])
            if zup > minimum(H[j, :])
                So[j] = minimum(So)
                zz[j] = zz[j-1] .+ So[j] .* (x[j] - x[j-1])
            end
        end
        So
    catch
        S0
    end
    So
end

"""
Assimilate SWOT observations for river reach.

- `H`: water surface elevation
- `W`: water surface width
- `x`: downstream distance for each cross section
- `wbf`: bankfull width
- `hbf`: bankfull depth
- `Qₚ`: prior probability distribution for discharge
- `nₚ`: prior probability distribution for roughness coefficient
- `rₚ`: prior probability distribution for channel shape parameter
- `zₚ`: prior distribution for downstream bed elevation
- `nens`: ensemble size
- `ri`: reach definition indices
- `ϵₒ`: observation error standard deviation

"""
function assimilate(H, W, x, wbf, hbf, S, Qₚ, nₚ, rₚ, zₚ, nens, ri; ϵₒ=0.01, logQ=false)
    min_ensemble_size = 5
    Qa = zeros(length(ri)-1, size(H, 2))
    Qu = zeros(length(ri)-1, size(H, 2))
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    So = bed_slope(S, H, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens)
    for t in 1:size(H, 2)
        vxs = findall(!ismissing, H[:, t]) # valid cross sections
        if length(vxs) > 1
            # FIXME: Need a more robust way to remap the reaches within the river if there are missing data
            # rri = [ri[1:end-1]; length(x) + 1]
            rri = [1; length(x[vxs])+1]
            he = gvf_ensemble!(H[vxs, t], So[vxs], x[vxs] .- x[vxs][1], hbf[vxs], wbf[vxs], Qe, ne, re, ze[vxs, :])
            i = findall(he[1, :] .> 0)
            if length(i) < min_ensemble_size
                Sf = [S[j, t] > 0 ? S[j, t] : minimum(S[:, t][S[:, t] .> 0]) for j in vxs]
                he = [(Qe[e] .* ne[e]) ./ (W[vxs[j], t] .* Sf[j].^0.5).^(3/5) for j=1:length(vxs), e=1:nens]
            end
            i = findall(he[1, :] .> 0)
            h = ze[vxs, :] .+ he .* ((re .+ 1) ./ re)'
            X = repeat(Qe[i]', outer=length(ri)-1)
            XA = h[:, i]
            d = H[vxs, t]
            E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
            if logQ
                A = letkf(log.(X), d, XA, E, [[j] for j in 1:length(ri)-1],
                          [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
                Qa[:, t] = mean(exp.(A), dims=2)
                Qu[:, t] = std(exp.(A), dims=2)
            else
                A = letkf(X, d, XA, E, [[j] for j in 1:length(ri)-1],
                          [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
                A[A .< 0] .= 0.0
                Qa[:, t] = mean(A, dims=2)  # could also use the absolute of the mean
                Qu[:, t] = std(A, dims=2)
            end
        else
            Qa[:, t] .= NaN
            Qu[:, t] .= 0.0
        end
    end
    A0, n = estimate_Q_params(H, W, S, ri, ze, re, ne, Qa)
    A0, n, Qa, Qu
end

"""
Estimate parameters for discharge estimation, i.e. A₀ and n.

"""
function estimate_Q_params(H, W, S, ri, ze, re, ne, Qa)
    nr, nt = size(Qa)
    nens = length(ne)
    A0 = zeros(nr)
    n = zeros(nr)
    Wbf = [maximum(skipmissing(W[c, :])) for c=1:size(W, 1)]
    Hbf = [maximum(skipmissing(H[c, :])) for c=1:size(H, 1)]
    vt = findall(!ismissing, H[1, :])
    Hmin = H[:, vt[findmin(H[1, vt])[2]]]
    Wmin = W[:, vt[findmin(H[1, vt])[2]]]
    ybf = Hbf .- ze
    ymin = Hmin .- ze
    ybf[ybf .<= 0] .= minimum(ybf[ybf .> 0])
    try
        ymin[ymin .<= 0] .= minimum(ymin[ymin .> 0])
    catch
        # ymin = [0.1 for i=1:length(ybf)]
        ymin = 0.1
    end
    A0e = Wbf .* (ymin ./ ybf).^(1 ./ re') .* (re ./ (re .+ 1))' .* ymin
    # A0e = Wmin .* (re ./ (re .+ 1))' .* ymin
    A0e = [mean(A0e[ri[j]:ri[j+1], e]) for j in 1:nr, e in 1:nens]
    dA = (H .- Hmin) .* (W .+ Wmin) ./ 2
    dA = reshape([smean(dA[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    Wr = reshape([smean(W[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    Sr = reshape([smean(S[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    A = repeat(A0e, outer=nt) .+ dA
    A[A .< 0] .= 0.0
    Sr[Sr .< 0] .= minimum(Sr[Sr .> 0])
    Qe = (1 ./ ne') .* A.^(5/3) .* Wr.^(-2/3) .* Sr.^(1/2)
    X = zeros(nr*2, nens)
    X[1:2:end, :] = A0e
    X[2:2:end, :] = repeat(ne', outer=nr)
    vt = [t for t=1:size(Qe, 1) if all(.!isnan.(Qe[t, :]))]
    XA = Qe[vt, :]
    nt = length(vt)
    d = reshape(Qa[:, vt], nr*nt)
    E = rand(Normal(0.1*mean(d), 1e-6), length(d), nens) .* rand([-1, 1], length(d), nens)
    A = letkf(X, d, XA, E, [[2*j-1;2*j] for j in 1:nr], [collect(nt*j-1:nt*j) for j in 1:nr], diagR=true)
    A0 = mean(A[1:2:end, :], dims=2)
    n = mean(A[2:2:end, :], dims=2)
    A0, n
end


end
