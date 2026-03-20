using DelimitedFiles
using Distributions
using LinearAlgebra
using NormalModes 
using PeriodicTable
using Statistics
using Test
using Unitful
using UnitfulAtomic

example_folder = "../example"
ħ = auconvert(1u"ħ")

begin
    systems = Dict()

    let
        folder = joinpath(example_folder, "iodopyridine")
        Z = vec(readdlm(joinpath(folder, "Z.txt"), Int))
        masses = [elements[z].atomic_mass for z in Z]
        hessian = readdlm(joinpath(folder, "hessian.dat")) * aunit(u"J*m^-2")

        nm = NormalDecomposition(hessian, Z)

        systems["iodopyridine"] = nm
    end
end

@testset "Basics" begin
    @testset "$name" for (name, nm) in systems
        xs, ps = sample(nm, 1000000)
        zs = project_geometries(nm, xs)
        pzs = project_momenta(nm, ps)
        m = 1aunit(u"kg")
        ωs = omegas(nm)

        Es = @. pzs^2 / (2m) + 1/2 * m * ωs^2 * zs^2
        Es = auconvert.(Es)
        fs = @. Es / (1/2 * ħ * ωs)

        @test all(abs.(mean(fs ; dims = 2) .- 1) .< 0.01) 
        @test all(abs.(std(fs ; dims = 2) .- 1) .< 0.01) 
    end
end

@testset "Orca" begin
    data = load_orca("../example/formic_acid/mol_opt.hess")
    nm = load_orca(NormalDecomposition, "../example/formic_acid/mol_opt.hess")

    @test data["vibrational_frequencies"][7:end] ≈ wave_numbers(nm) rtol=1e-4

    for (mode, orca_mode) in zip(eachcol(normal_modes(nm)), eachcol(data["normal_modes"][:, 7:end]))
        mode ./= norm(mode)
        mode *= sign(argmax(abs, mode))
        orca_mode *= sign(argmax(abs, orca_mode))

        @test mode ≈ orca_mode rtol=2e-2 atol=1e-3
    end
end