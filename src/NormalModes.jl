module NormalModes

using Chain
using Distributions
using LinearAlgebra
using PeriodicTable
using Printf
using Random
using StatsBase
using Unitful
using UnitfulAtomic

import PhysicalConstants.CODATA2018: k_B

export NormalDecomposition
export project_geometries, project_momenta, project_per_atom
export normal_modes, normal_mode, frequencies, wave_number, reduced_masses
export energies, omegas
export spatial_variances, momentum_variances
export spatial_fluctuation, momentum_fluctuation
export spatial_Σ, momentum_Σ
export sample, wigner_sample

const hbar = 1  # Atomic units

function mass_weight_matrix(elements)
    masses = to_atomic_masses(elements)

    m = @chain masses begin
        vcat([fill(m, 3) for m in _]...)
        1 ./ sqrt.(_)
    end

    return Diagonal(m)
end

# TODO add the number of releveant mode somewhere
struct NormalDecomposition{T}
    elements::Vector
    ωs::Vector{T}  # Angular frequencies
    M::Diagonal{T, Vector{T}}  # Inverse square root of the masses
    U::Matrix{T}  # Orthonormal modes
end

"""
    NormalDecomposition(hessian, elements)

Compute the normal mode decomposition of the `3n x 3n` Hessian for a system
composed of atoms of the given elements, where `n` is the number of atoms.

The elements can be given as either elements from PeriodicTable, by their
atomic number or by their mass, defaulting to atomic units when not usng a
Unitful quantity.

By default the 6 modes with lowest frequencies are skipped, as they are likely
to represent rotations and translations.
"""
function NormalDecomposition(hessian::AbstractMatrix, elements ; skip_modes = 6)
    hessian = to_atomic_units(hessian)
    M = mass_weight_matrix(elements)
    Ω, U = eigen(M * hessian * M)

    # Sometimes very low frequencies are artificially negative
    ωs = sqrt.(abs.(Ω))
    perm = sortperm(ωs)[(skip_modes + 1):end]
    ωs = ωs[perm]
    U = U[:, perm]

    return NormalDecomposition(elements, ωs, M, U)
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

n_atoms(nm::NormalDecomposition) = div(size(nm.M, 1), 3)

function to_atomic_units(x)
    if !(eltype(x) <: Quantity)
        @warn "NormalModes: A unitless quantity was given. We assume it was given in atomic units."
        return x
    end
    return austrip.(x)
end

function to_atomic_masses(x)
    if eltype(x) <: Element
        return [austrip(elem.atomic_mass) for elem in x]
    elseif eltype(x) <: Union{Integer, Symbol}
        return [austrip(elements[E].atomic_mass) for E in x]
    else
        return to_atomic_units(x)
    end
end

"""
    project_displacements(nm::NormalDecomposition, displacements)

Project the given displacements (3 x n_atoms x n_obs) from the equilibrium
geometry on the normal modes.
"""
function project_geometries(nm::NormalDecomposition, displacements::AbstractArray)
    projector = nm.U' * inv(nm.M)
    return projector * reshape(displacements, size(nm.M, 1), :)
end

"""
    project_momenta(nm::NormalDecomposition, geometries)

Project the given momenta (3 x n_atoms x n_obs) on the normal modes.
"""
function project_momenta(nm::NormalDecomposition, momenta::AbstractArray)
    projector = nm.U' * nm.M
    return projector * reshape(momenta, size(nm.M, 1), :)
end

function wigner_pdf(nm::NormalDecomposition, displacements, momenta, T = 0.0u"K")
    zs = project_geometries(nm, displacements)
    pzs = project_momenta(nm, momenta)

    zs_dist = MvNormal(Diagonal(spatial_variances(nm, T)))
    pz_dist = MvNormal(Diagonal(momentum_variances(nm, T)))

    return pdf(zs_dist, zs) .* pdf(pz_dist, pzs)
end

function temperature_bias(nm, displacements, momenta, T)
    w0 = wigner_pdf(nm, displacements, momenta, 0.0u"K")
    wT = wigner_pdf(nm, displacements, momenta, T)

    return wT ./ w0
end

function project_per_atom(nm::NormalDecomposition, geometries)
    n_atoms = length(nm.elements)
    xx = reshape(geometries, 3, n_atoms, :)
    uu = reshape(inv(nm.M) * nm.U, 3, n_atoms, :)

    projections = zeros(n_atoms, size(nm.U, 2), size(xx, 3))
    for (k, x) in enumerate(eachslice(xx ; dims = 3))
        projections[:, :, k] = sum(x .* uu ; dims = 1)
    end

    return projections
end

"""
    normal_modes(nm::NormalDecomposition)

Return the real space normal modes.
"""
function normal_modes(nm::NormalDecomposition)
    modes = nm.M * nm.U
    return modes
    μs = norm.(eachcol(modes))
    return modes ./ reshape(μs, 1, :)
end

"""
    normal_mode(nm::NormalDecomposition, n)

Return the n-th real space normal modes.
"""
function normal_mode(nm::NormalDecomposition, mode_number)
    return reshape(normal_modes(nm)[:, mode_number], 3, :)
end

function reduced_masses(nm::NormalDecomposition)
    @chain nm.M begin
        _ .^ 2
        _ * (nm.U .^ 2)
        sum(_ ; dims = 1)
        vec(_)
        1 ./ _
        _ * aunit(u"kg")
        uconvert.(u"u", _)
    end
end

omegas(nm::NormalDecomposition) = nm.ωs * aunit(u"1/s")

function frequencies(nm::NormalDecomposition)
    @chain nm.ωs begin
        _ ./ 2π 
        _ * aunit(u"s^-1")  # Add atomic unit of inverse time
        uconvert.(u"GHz", _)  # Standard units
    end
end

function energies(nm::NormalDecomposition)
    return hbar * nm.ωs * aunit(u"J")
end

function wave_number(nm::NormalDecomposition)
    @chain nm begin
        frequencies(_)
        _ ./ 1u"c"  # Divide by the speed of light
        uconvert.(u"cm^-1", _)  # Standard units
    end
end

"""
    thermal_widenings(nm::NormalDecomposition, T::Quantity)

Widening of the distribution for Wigner sampling in the harmonic approximation.
"""
function thermal_widenings(nm::NormalDecomposition, T::Quantity)
    T = uconvert(u"K", T)

    T < 0u"K" && throw(DomainError(T, "temperature can not be below absolute zero"))
    T == 0u"K" && return ones(length(energies(nm)))

    return 1 ./ tanh.(upreferred.(energies(nm) ./ (k_B * T)))
end

spatial_variances(nm::NormalDecomposition) = 1/2 * (hbar ./ nm.ωs)
momentum_variances(nm::NormalDecomposition) = 1/2 * hbar * nm.ωs

function spatial_variances(nm::NormalDecomposition, T::Quantity)
    return spatial_variances(nm) .* thermal_widenings(nm, T)
end

function momentum_variances(nm::NormalDecomposition, T::Quantity)
    return momentum_variances(nm) .* thermal_widenings(nm, T)
end

function spatial_fluctuation(nm::NormalDecomposition)
    n = length(nm.ωs)
    x = zeros(n, n)
    Σ = spatial_variances(nm)

    for k in 1:n
        x[k, k] = sqrt(Σ[k])
    end
    return reshape(nm.M * nm.U * x, 3, :, n) * aunit(u"m")
end

function momentum_fluctuation(nm::NormalDecomposition)
    n = length(nm.ωs)
    p = zeros(n, n)
    Σ = momentum_variances(nm)

    for k in 1:n
        p[k, k] = sqrt(Σ[k])
    end
    return reshape(inv(nm.M) * nm.U * p, 3, :, n) * aunit(u"kg * m / s")
end

function spatial_Σ(nm::NormalDecomposition)
    nm.M * nm.U * Diagonal(spatial_variances(nm)) * nm.U' * nm.M' * aunit(u"m^2")
end

function momentum_Σ(nm::NormalDecomposition)
    inv(nm.M) * nm.U * Diagonal(momentum_variances(nm)) * nm.U' * inv(nm.M)' * aunit(u"kg^2*m^2/s^2")
end


"""
    sample(nm::NormalDecomposition[, n_samples])

Perform ground state Wigner sampling according to the normal decomposition.

Return the deviation from the average geometry and the deviation from zero
momentum.
"""
function StatsBase.sample(rng::AbstractRNG, nm::NormalDecomposition, n_samples::Integer)
    Δx_dist = MvNormal(Diagonal(spatial_variances(nm)))
    Δp_dist = MvNormal(Diagonal(momentum_variances(nm)))
    
    Δx = nm.M * nm.U * rand(rng, Δx_dist, n_samples)
    Δp = inv(nm.M) * nm.U * rand(rng, Δp_dist, n_samples)

    return Δx * aunit(u"m"), Δp * aunit(u"kg*m/s")
end

function StatsBase.sample(rng::AbstractRNG, nm::NormalDecomposition)
    geometries, momenta = sample(rng, nm, 1)
    return geometries[:, 1], momenta[:, 1]
end

StatsBase.sample(nm::NormalDecomposition, n_samples) = sample(Random.GLOBAL_RNG, nm, n_samples)
StatsBase.sample(nm::NormalDecomposition) = sample(Random.GLOBAL_RNG, nm)

function wigner_sample(rng::AbstractRNG, nm::NormalDecomposition, T::Quantity, n_samples::Integer)
    Δx_dist = MvNormal(Diagonal(spatial_variances(nm, T)))
    Δp_dist = MvNormal(Diagonal(momentum_variances(nm, T)))
    
    Δx = nm.M * nm.U * rand(rng, Δx_dist, n_samples)
    Δp = inv(nm.M) * nm.U * rand(rng, Δp_dist, n_samples)

    return Δx * aunit(u"m"), Δp * aunit(u"kg*m/s")
end

function wigner_sample(rng::AbstractRNG, nm::NormalDecomposition, T::Quantity)
    geometries, momenta = wigner_sample(rng, nm, T, 1)
    return geometries[:, 1], momenta[:, 1]
end

wigner_sample(nm::NormalDecomposition, T, n_samples) = wigner_sample(Random.GLOBAL_RNG, nm, T, n_samples)
wigner_sample(nm::NormalDecomposition, T) = wigner_sample(Random.GLOBAL_RNG, nm, T)

end
