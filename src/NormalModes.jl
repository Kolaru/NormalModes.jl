module NormalModes

using Chain
using Distributions
using LinearAlgebra
using PeriodicTable
using Printf
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

By default the first 
"""
function NormalDecomposition(hessian::AbstractMatrix, elements ; skip_modes = 6)
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
    perm = sortperm(ωs)[(skip_modes + 1):end]
    ωs = ωs[perm]
    U = U[:, perm]

    return NormalDecomposition(ωs, m, U)
end

function Base.show(io::IO, nm::NormalDecomposition)
    println(io,
        "NormalDecomposition(n_atoms = $(n_atoms(nm)), n_modes = $(size(nm.U, 2)))"
    )
    println(io,
        "[mode]   [frequency]"
    )

    for (i, k) in enumerate(wave_number(nm))
        val = @sprintf "%5d %7.0f" i ustrip(k)
        println(io, "$val $(unit(k))")
    end
end

n_atoms(nm::NormalDecomposition) = div(length(nm.m), 3)

function to_atomic_units(x)
    if !(eltype(x) <: Quantity)
        @warn "NormalModes: A unitless quantity was given. We assume it is given in atomic units."
        return x
    end
    return austrip.(x)
end

function to_atomic_masses(x)
    if eltype(x) <: Element
        return [austrip(elem.atomic_mass) for elem in x]
    elseif eltype(x) <: Integer
        return [austrip(elements[E].atomic_mass) for E in x]
    else
        return to_atomic_units(x)
    end
end

"""
    project(nm::NormalDecomposition, geometries)

Project the given geometries (3 x n_atoms x n_obs) on the normal modes.
"""
function project(nm::NormalDecomposition, geometries::AbstractArray)
    projector = nm.U' * Diagonal(1 ./ nm.m)
    return projector * reshape(geometries, length(nm.m), :)
end

"""
    normal_modes(nm::NormalDecomposition)

Real space normal modes of the normal decomposition.
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
