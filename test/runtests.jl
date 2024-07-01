using DelimitedFiles
using NormalModes 
using PeriodicTable
using Statistics
using Test
using Unitful
using UnitfulAtomic

example_folder = "../example"
ħ = auconvert(1u"ħ")

@testset "NormalModes.jl" begin
    @testset "Iodopyridine" begin
        folder = joinpath(example_folder, "iodopyridine")
        Z = vec(readdlm(joinpath(folder, "Z.txt"), Int))
        masses = [elements[z].atomic_mass for z in Z]
        hessian = readdlm(joinpath(folder, "hessian.dat")) * aunit(u"J*m^-2")

        nm = NormalDecomposition(hessian, Z)

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