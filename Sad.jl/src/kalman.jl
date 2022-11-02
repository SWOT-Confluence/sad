using LinearAlgebra
using Statistics

"""
Apply the Local Ensemble Transform Kalman Filter algorithm.

- `A`: state variable ensemble matrix
- `d`: observation vector
- `HA`: model-predicted observation ensemble matrix
- `E`: observation error ensemble matrix
- `xs`: indices of state vector that map to observation local patches
- `ys`: indices of observation local patches
- `ρ`: covariance inflation parameter (optional, default value of 1.05)
- `diagR`: force observation error covariance to be diagonal, i.e. no spatial correlation (optional, default is false)

"""
function letkf(A::Matrix, d::Vector, HA::Matrix, E::Matrix, xs::Vector, ys::Vector; ρ=1.05, diagR=false)
    Aa = zeros(size(A))
    ndim, nens = size(A)
    nobs = length(d)
    if size(E) == (nobs, nobs)
        R = E
    else
        if diagR
            R = diagm(0 => diag((1 / (nens - 1)) * E * E'))
        else
            R = (1 / (nens - 1)) * E * E'
        end
    end
    Y = HA .- mean(HA, dims=2)
    X = A .- mean(A, dims=2)
    for (lx, ly) in zip(xs, ys)
        Xₗ = X[lx, :]
        Yₗ = Y[ly, :]
        Rₗ = R[ly, ly]
        C = Yₗ' * pinv(Rₗ)
        P = pinv((nens - 1) / ρ * I + C * Yₗ)
        W = real(sqrt((nens - 1) * P))
        w = P * C * (d[ly] .- mean(HA[ly, :], dims=2))
        W = W .+ w
        Xₗᵃ = Xₗ * W .+ mean(A[lx, :], dims=2)
        Aa[lx, :] = Xₗᵃ
    end
    Aa
end
