using Sad: CrossSection

"""
Rectangular river channel cross section.

- `W⃰`: bankfull width
- `Yₘ⃰`: bankfull depth
- `Yₘ`: flow depth
- `S0`: bed slope
- `n`: Manning roughness coefficient

See also: [`Dingman`](@ref)
"""
mutable struct Rectangular <: CrossSection
    W⃰:: Real
    Yₘ⃰:: Real
    Yₘ:: Real
    S0:: Real
    n:: Real
end

"""
Dingman river channel cross section.

Uses generalized expressions for cross-section geometry based on at-a-station hydraulic geometry relations.

- `W⃰`: bankfull width
- `Yₘ⃰`: bankfull depth
- `Yₘ`: flow depth
- `r`: geometry exponent
- `S0`: bed slope
- `n`: Manning roughness coefficient

See also: [`Rectangular`](@ref)
"""
mutable struct Dingman <: CrossSection
    W⃰:: Real
    Yₘ⃰:: Real
    Yₘ:: Real
    r:: Real
    S0:: Real
    n:: Real
end

"""
Calculate water surface width of cross section.
"""
function width(xs::Dingman)
    xs.W⃰ * (xs.Yₘ / xs.Yₘ⃰)^(1 / xs.r)
end

function width(xs::Rectangular)
    xs.W⃰
end

"""
Calculate average depth of cross section.
"""
function depth(xs::Dingman)
    (xs.r / (xs.r + 1)) * xs.Yₘ
end

function depth(xs::Rectangular)
    xs.Yₘ
end

export width, depth
