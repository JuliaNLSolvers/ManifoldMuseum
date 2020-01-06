using Test
using Manifolds

@testset "Fields" begin
    @test ℝ isa Manifolds.RealNumbers
    @test Manifolds.RealNumbers() === ℝ
    @test field_dimension(ℝ) == 1

    @test ℂ isa Manifolds.ComplexNumbers
    @test Manifolds.ComplexNumbers() === ℂ
    @test field_dimension(ℂ) == 2
end
