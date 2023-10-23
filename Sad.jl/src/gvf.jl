using DifferentialEquations

"""
Calculate the Froude number.

"""
function froude(Q::Real, xs::CrossSection)
    g = 9.806
    Fr = Q / (depth(xs) * width(xs) * sqrt(g * depth(xs)))
end

"""
Ordinary differential equation describing the Gradually-Varied-Flow model.

"""
function dydx(y, p, x)
    i = trunc(Int, ceil(x / (p[3] / length(p[1]))))
    i = i > 0 ? i : 1
    # TODO: Change this dynamically depending on whether the flow is subcritical or supercritical
    pm = -1 # integrate upstream
    Q = p[1]
    xs = p[2][i]
    xs.Yₘ = y * (xs.r + 1) / xs.r
    S0 = xs.S0
    Sf = xs.n^2 * Q^2 / (width(xs)^2 * depth(xs)^(10/3))
    Fr = froude(Q, xs)
    pm * (S0 - Sf) / (1 - Fr^2)
end

"""
Calculate a water surface profile by solving the Gradually-Varied-Flow equation.

- `Q`: discharge
- `ybc`: downstream boundary condition for depth
- `z0`: bed elevation at downstream cross section
- `S0`: bed slope for reach
- `n`: Manning's roughness coefficient
- `x`: downstream distance for each cross section
- `b`: bottom width for each cross section
- `r`: channel geometry coefficient

"""
function gvf(Q, ybc, S0, n, x, wbf, ybf, r)
    c = [Dingman(wbf[i], ybf[i], ybc, r[i], S0[i], n) for i in 1:length(x)]
    prob = DiscreteProblem(dydx, ybc, (x[1], x[end]), (Q, c, x[end]))
    h = try
        sol = solve(prob, Tsit5(), abstol=1e-2, saveat=x)
        sol.u
    catch
        zeros(length(x)) .- 9999.
    end
    h
end

