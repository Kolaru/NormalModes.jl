module NormalModes

using Chain
using Distributions
using LinearAlgebra
using PeriodicTable
using StatsBase
using Unitful
using UnitfulAtomic

export NormalDecomposition
export project, normal_modes, frequencies, wave_number, reduced_masses
export sample

# TODO add the number of releveant mode somewhere
struct NormalDecomposition{T}
    ωs::Vector{T}  # Angular frequencies
    m::Vector{T}  # Inverse square root of the masses
    U::Matrix{T}  # Orthonormal modes
end

function NormalDecomposition(hessian::Matrix, masses::Vector ; valid_modes = 1:3length(masses))
    hessian = to_atomic_units(hessian)
    masses = to_atomic_masses(masses)

    m = @chain masses begin
        vcat([fill(m, 3) for m in _]...)
        1 ./ sqrt.(_)
    end

    M = Diagonal(m)

    Ω, U = eigen(M * hessian * M)

    # Sometimes very low frequencies are artificially negative
    ωs = sqrt.(abs.(Ω))
    perm = sortperm(ωs)[valid_modes]
    ωs = ωs[perm]
    U = U[:, perm]

    return NormalDecomposition(ωs, m, U)
end

function to_atomic_units(x::VecOrMat)
    @warn "A unitless quantity was given. We assume it was atomic units."
    return x
end
to_atomic_units(x::VecOrMat{<:Quantity}) = austrip.(x)

to_atomic_masses(x) = to_atomic_units(x)
to_atomic_masses(x::VecOrMat{<:Element}) = [austrip(elem.atomic_mass) for elem in x]
to_atomic_masses(x::VecOrMat{<:Integer}) = to_atomic_masses([elements[i] for i in x])

function project(nm::NormalDecomposition, v::Array{<:Any, 3})
    projector = nm.U' * Diagonal(1 ./ nm.m)
    return projector * reshape(v, length(nm.m), :)
end

function normal_modes(nm::NormalDecomposition)
    modes = Diagonal(nm.m) * nm.U
    μs = norm.(eachcol(modes))
    return modes ./ reshape(μs, 1, :)
end

function reduced_masses(nm::NormalDecomposition)
    @chain nm.m begin
        Diagonal(_) .^ 2
        _ * (nm.U .^ 2)
        sum(_ ; dims = 1)
        vec(_)
        1 ./ _
        _ * aunit(u"kg")
        uconvert.(u"u", _)
    end
end

frequencies(nm::NormalDecomposition) = @chain nm.ωs begin
    _ ./ 2π 
    _ * aunit(u"s^-1")  # Add atomic unit of inverse time
    uconvert.(u"GHz", _)  # Standard units
end

wave_number(nm::NormalDecomposition) = @chain nm begin
    frequencies(_)
    _ ./ 1u"c"  # Divide by the speed of light
    uconvert.(u"cm^-1", _)  # Standard units
end

# TODO Not sure this is correct
function position_covariance(nm::NormalDecomposition ; regularize = true)
    modes = normal_modes(nm)
    rawΣ = modes * Diagonal(nm.ωs) * modes'

    !regularize && return Symetric(rawΣ)

    evals, evecs = eigen(rawΣ ; sortby = abs)
    evals[1:6] .= evals[7] * 1e-15

    return Symmetric(evecs * Diagonal(evals) * evecs')
end

# TODO Extend to non-atomic units
function StatsBase.sample(nm::NormalDecomposition, n_samples)
    hbar = 1
    Δz_dist = MvNormal(Diagonal(hbar ./ 2nm.ωs))
    Δp_dist = MvNormal(Diagonal(hbar * nm.ωs / 2))
    MU = Diagonal(nm.m) * nm.U
    
    Δx = MU*rand(Δz_dist, n_samples) 
    Δp = Diagonal(nm.m).^2 * MU*rand(Δp_dist, n_samples)
    return Δx, Δp
end

end
