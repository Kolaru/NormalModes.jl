module NormalModes

using Chain
using Distributions
using LinearAlgebra
using StatsBase
using Unitful
using UnitfulAtomic

export NormalModes
export project, normal_modes, mode_masses, frequencies, sample

# TODO add the number of releveant mode somewhere
struct NormalModes{T}
    ωs::Vector{T}
    modes::Matrix{T}
    m::Vector{T}  # Diagonal of the mass weighting matrix M
    U::Matrix{T}  # Orthonormal modes
end

function NormalModes(hessian, masses ; valid_modes = 7:3length(masses))
    m = @chain masses begin
        vcat([fill(m, 3) for m in _]...)
        1 ./ sqrt.(_)
    end

    M = Diagonal(m)

    W, U = eigen(M * hessian * M)

    # Sometimes very low frequencies are artificially negative
    ωs = sqrt.(abs.(W))
    perm = sortperm(ωs)[valid_modes]
    ωs = ωs[perm]
    U = U[:, perm]

    return NormalModes(ωs, M * U, m, U)
end

function project(nm::NormalModes, v::Array{<:Any, 3})
    projector = nm.U' * Diagonal(1 ./ nm.m)
    return projector * reshape(v, length(nm.m), :)
end

function normal_modes(nm::NormalModes)
    modes = Diagonal(nm.m) * nm.U
    μs = mode_masses(nm)
    return modes ./ reshape(μs, 1, :)
end

mode_masses(nm::NormalModes) = norm.(eachcol(nm.modes))
frequencies(nm::NormalModes) = nm.ωs

# TODO Not sure this is correct
function position_covariance(nm::NormalModes ; regularize = true)
    modes = normal_modes(nm)
    rawΣ = modes * Diagonal(nm.ωs) * modes'

    !regularize && return Symetric(rawΣ)

    evals, evecs = eigen(rawΣ ; sortby = abs)
    evals[1:6] .= evals[7] * 1e-15

    return Symmetric(evecs * Diagonal(evals) * evecs')
end

# TODO Extend to non-atomic units
function StatsBase.sample(nm::NormalModes, n_samples)
    hbar = austrip(1u"hbar")
    Δz_dist = MvNormal(Diagonal(hbar ./ 2nm.ωs))
    Δp_dist = MvNormal(Diagonal(nm.ωs / 2))
    MU = Diagonal(nm.m) * nm.U
    MU = normal_modes(nm)

    Δx = MU*rand(Δz_dist, n_samples) 
    Δp = rand(Δp_dist, n_samples)
    return Δx, Δp
end

end
