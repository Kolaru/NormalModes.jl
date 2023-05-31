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

"""
    NormalDecomposition(hessian, elements)

Compute the normal mode decomposition of the `3n x 3n` Hessian for a system
composed of atoms of the given elements, where `n` is the number of atoms.

The elements can be given as either elements from PeriodicTable, by their
atomic number or by their mass, defaulting to atomic units when not usng a
Unitful quantity.
"""
function NormalDecomposition(hessian::Matrix, elements::Vector ; valid_modes = 1:3length(elements))
    hessian = to_atomic_units(hessian)
    masses = to_atomic_masses(elements)

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
    @warn "NormalModes: A unitless quantity was given. We assume it was in atomic units."
    return x
end
to_atomic_units(x::VecOrMat{<:Quantity}) = austrip.(x)

to_atomic_masses(x) = to_atomic_units(x)
to_atomic_masses(x::VecOrMat{<:Element}) = [austrip(elem.atomic_mass) for elem in x]
to_atomic_masses(x::VecOrMat{<:Integer}) = to_atomic_masses([elements[i] for i in x])

"""
    project(nm::NormalDecomposition, v)

Project the given geometries (3 x n_atoms x n_obs) on the normal modes.
"""
function project(nm::NormalDecomposition, v::Array{<:Any, 3})
    projector = nm.U' * Diagonal(1 ./ nm.m)
    return projector * reshape(v, length(nm.m), :)
end

"""
    normal_modes(nm::NormalDecomposition)

Real space normal modes according to a normal decomposition.
"""
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

function frequencies(nm::NormalDecomposition)
    @chain nm.ωs begin
        _ ./ 2π 
        _ * aunit(u"s^-1")  # Add atomic unit of inverse time
        uconvert.(u"GHz", _)  # Standard units
    end
end

function wave_number(nm::NormalDecomposition)
    @chain nm begin
        frequencies(_)
        _ ./ 1u"c"  # Divide by the speed of light
        uconvert.(u"cm^-1", _)  # Standard units
    end
end

"""
    sample(nm::NormalDecomposition, n_samples)

Perform Wigner sampling according to the normal decomposition.

Return the deviation from the average geometry.
"""
function StatsBase.sample(nm::NormalDecomposition, n_samples)
    hbar = 1  # Atomic units
    m = 1  # All modes have mass 1
    Δz_dist = MvNormal(Diagonal(1/2 * (hbar ./ (m * nm.ωs))))
    Δp_dist = MvNormal(Diagonal(1/2 * hbar * m * nm.ωs))
    MU = Diagonal(nm.m) * nm.U
    
    Δx = MU*rand(Δz_dist, n_samples) 
    Δp = Diagonal(nm.m).^2 * MU*rand(Δp_dist, n_samples)
    return Δx * aunit(u"m"), Δp * aunit(u"kg*m/s")
end

end
