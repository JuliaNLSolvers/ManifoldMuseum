@info "Manifolds.jl Test settings:\n\n" *
      "Testing Float32:  $(TEST_FLOAT32)\n" *
      "Testing Double64: $(TEST_DOUBLE64)\n" *
      "Testing Static:   $(TEST_STATIC_SIZED)\n\n" *
      "Test group:       $(TEST_GROUP)\n\n" *
      "These settings are stored in environment variables, see in test/header.jl"

function include_test(path)
    @info "Testing $path"
    @time include(path)  # show basic timing, (this will print a newline at end)
end

@testset "Manifolds.jl" begin
    if TEST_GROUP ∈ ["all", "test_manifolds"]
        include_test("differentiation.jl")
        # include_test("ambiguities.jl") # -> moved to Aqua.jl
        include_test("test_deprecated.jl")
        @testset "utils test" begin
            Random.seed!(42)
            @testset "usinc_from_cos" begin
                @test Manifolds.usinc_from_cos(-1) == 0
                @test Manifolds.usinc_from_cos(-1.0) == 0.0
            end
            @testset "log_safe!" begin
                n = 8
                Q = qr(randn(n, n)).Q
                A1 = Matrix(Hermitian(Q * Diagonal(rand(n)) * Q'))
                @test exp(Manifolds.log_safe!(similar(A1), A1)) ≈ A1 atol = 1e-6
                A1_fail = Matrix(Hermitian(Q * Diagonal([-1; rand(n - 1)]) * Q'))
                @test_throws DomainError Manifolds.log_safe!(similar(A1_fail), A1_fail)

                T = triu!(randn(n, n))
                T[diagind(T)] .= rand.()
                @test exp(Manifolds.log_safe!(similar(T), T)) ≈ T atol = 1e-6
                T_fail = copy(T)
                T_fail[1] = -1
                @test_throws DomainError Manifolds.log_safe!(similar(T_fail), T_fail)

                A2 = Q * T * Q'
                @test exp(Manifolds.log_safe!(similar(A2), A2)) ≈ A2 atol = 1e-6
                A2_fail = Q * T_fail * Q'
                @test_throws DomainError Manifolds.log_safe!(similar(A2_fail), A2_fail)

                A3 = exp(SizedMatrix{n,n}(randn(n, n)))
                @test A3 isa SizedMatrix
                @test exp(Manifolds.log_safe!(similar(A3), A3)) ≈ A3 atol = 1e-6
                @test exp(Manifolds.log_safe(A3)) ≈ A3 atol = 1e-6

                A3_fail = Float64[1 2; 3 1]
                @test_throws DomainError Manifolds.log_safe!(similar(A3_fail), A3_fail)

                A4 = randn(ComplexF64, n, n)
                @test exp(Manifolds.log_safe!(similar(A4), A4)) ≈ A4 atol = 1e-6
            end
            @testset "isnormal" begin
                @test !Manifolds.isnormal([1.0 2.0; 3.0 4.0])
                @test !Manifolds.isnormal(complex.(reshape(1:4, 2, 2), reshape(5:8, 2, 2)))

                # diagonal
                @test Manifolds.isnormal(diagm(randn(5)))
                @test Manifolds.isnormal(diagm(randn(ComplexF64, 5)))
                @test Manifolds.isnormal(Diagonal(randn(5)))
                @test Manifolds.isnormal(Diagonal(randn(ComplexF64, 5)))

                # symmetric/hermitian
                @test Manifolds.isnormal(Symmetric(randn(3, 3)))
                @test Manifolds.isnormal(Hermitian(randn(3, 3)))
                @test Manifolds.isnormal(Hermitian(randn(ComplexF64, 3, 3)))
                x = Matrix(Symmetric(randn(3, 3)))
                x[3, 1] += eps()
                @test !Manifolds.isnormal(x)
                @test Manifolds.isnormal(x; atol=sqrt(eps()))

                # skew-symmetric/skew-hermitian
                skew(x) = x - x'
                @test Manifolds.isnormal(skew(randn(3, 3)))
                @test Manifolds.isnormal(skew(randn(ComplexF64, 3, 3)))

                # orthogonal/unitary
                @test Manifolds.isnormal(Matrix(qr(randn(3, 3)).Q); atol=sqrt(eps()))
                @test Manifolds.isnormal(
                    Matrix(qr(randn(ComplexF64, 3, 3)).Q);
                    atol=sqrt(eps()),
                )
            end
            @testset "realify/unrealify!" begin
                # round trip real
                x = randn(3, 3)
                @test Manifolds.realify(x, ℝ) === x
                @test Manifolds.unrealify!(similar(x), x, ℝ) == x

                # round trip complex
                x2 = randn(ComplexF64, 3, 3)
                x2r = Manifolds.realify(x2, ℂ)
                @test eltype(x2r) <: Real
                @test size(x2r) == (6, 6)
                x2c = Manifolds.unrealify!(similar(x2), x2r, ℂ)
                @test x2c ≈ x2

                # matrix multiplication is preserved
                x3 = randn(ComplexF64, 3, 3)
                x3r = Manifolds.realify(x3, ℂ)
                @test x2 * x3 ≈ Manifolds.unrealify!(similar(x2), x2r * x3r, ℂ)
            end
            @testset "allocation" begin
                @test allocate([1 2; 3 4], Float64, Size(3, 3)) isa Matrix{Float64}
                @test allocate(SA[1 2; 3 4], Float64, Size(3, 3)) isa MMatrix{3,3,Float64}
                @test allocate(SA[1 2; 3 4], Size(3, 3)) isa MMatrix{3,3,Int}
                @test Manifolds.quat_promote(Float64) === Quaternions.QuaternionF64
                @test Manifolds.quat_promote(Float32) === Quaternions.QuaternionF32
                @test Manifolds.quat_promote(QuaternionF64) === Quaternions.QuaternionF64
                @test Manifolds.quat_promote(QuaternionF32) === Quaternions.QuaternionF32
            end
            @testset "eigen_safe" begin
                @test Manifolds.eigen_safe(SA[1.0 0.0; 0.0 1.0]) isa
                      Eigen{Float64,Float64,<:SizedMatrix{2,2},<:SizedVector{2}}
            end
            @testset "max_eps" begin
                x64 = randn(Float64, 2)
                x32 = randn(Float32, 2)
                z32 = randn(ComplexF32, 2)
                xi = rand(0:1, 2)
                @test Manifolds.max_eps(x64, x64) == eps()
                @test Manifolds.max_eps(x64, x32) == eps(Float32)
                @test Manifolds.max_eps(x32, x64) == eps(Float32)
                @test Manifolds.max_eps(xi, xi) == 0
                @test Manifolds.max_eps(xi, x64) == eps()
                @test Manifolds.max_eps(xi, x32) == eps(Float32)
                @test Manifolds.max_eps(xi, z32) == eps(Float32)
                @test Manifolds.max_eps(xi, x64, x32, z32) == eps(Float32)
            end
        end

        @test Manifolds.is_metric_function(flat)
        @test Manifolds.is_metric_function(sharp)

        include_test("groups/group_utils.jl")
        include_test("notation.jl")
        # starting with tests of simple manifolds
        include_test("manifolds/centered_matrices.jl")
        include_test("manifolds/circle.jl")
        include_test("manifolds/cholesky_space.jl")
        include_test("manifolds/elliptope.jl")
        include_test("manifolds/euclidean.jl")
        include_test("manifolds/fixed_rank.jl")
        include_test("manifolds/flag.jl")
        include_test("manifolds/generalized_grassmann.jl")
        include_test("manifolds/generalized_stiefel.jl")
        include_test("manifolds/grassmann.jl")
        include_test("manifolds/hamiltonian.jl")
        include_test("manifolds/hyperbolic.jl")
        include_test("manifolds/hyperrectangle.jl")
        include_test("manifolds/lorentz.jl")
        include_test("manifolds/multinomial_doubly_stochastic.jl")
        include_test("manifolds/multinomial_symmetric.jl")
        include_test("manifolds/multinomial_spd.jl")
        include_test("manifolds/positive_numbers.jl")
        include_test("manifolds/probability_simplex.jl")
        include_test("manifolds/projective_space.jl")
        include_test("manifolds/rotations.jl")
        include_test("manifolds/shape_space.jl")
        include_test("manifolds/skewhermitian.jl")
        include_test("manifolds/spectrahedron.jl")
        include_test("manifolds/sphere.jl")
        include_test("manifolds/sphere_symmetric_matrices.jl")
        include_test("manifolds/stiefel.jl")
        include_test("manifolds/symmetric.jl")
        include_test("manifolds/symmetric_positive_definite.jl")
        include_test("manifolds/spd_fixed_determinant.jl")
        include_test("manifolds/symmetric_positive_semidefinite_fixed_rank.jl")
        include_test("manifolds/symplectic.jl")
        include_test("manifolds/symplecticgrassmann.jl")
        include_test("manifolds/symplecticstiefel.jl")
        include_test("manifolds/tucker.jl")
        include_test("manifolds/unitary_matrices.jl")

        include_test("manifolds/essential_manifold.jl")
        include_test("manifolds/multinomial_matrices.jl")
        include_test("manifolds/oblique.jl")
        include_test("manifolds/torus.jl")

        #meta manifolds
        include_test("manifolds/product_manifold.jl")
        include_test("manifolds/power_manifold.jl")
        include_test("manifolds/quotient_manifold.jl")
        include_test("manifolds/fiber_bundle.jl")
        include_test("manifolds/vector_bundle.jl")
        include_test("manifolds/graph.jl")

        include_test("metric.jl")
        include_test("statistics.jl")
    end

    if TEST_GROUP ∈ ["all", "test_integration"]
        include_test("approx_inverse_retraction.jl")

        # manifolds requiring ODE solvers
        include_test("manifolds/embedded_torus.jl")
    end

    if TEST_GROUP ∈ ["test_lie_groups", "all"]
        # Lie groups and actions
        include_test("groups/groups_general.jl")
        include_test("groups/validation_group.jl")
        include_test("groups/circle_group.jl")
        include_test("groups/translation_group.jl")
        include_test("groups/general_linear.jl")
        include_test("groups/general_unitary_groups.jl")
        include_test("groups/special_linear.jl")
        include_test("groups/special_orthogonal.jl")
        include_test("groups/heisenberg.jl")
        include_test("groups/product_group.jl")
        include_test("groups/semidirect_product_group.jl")
        include_test("groups/power_group.jl")
        include_test("groups/special_euclidean.jl")
        include_test("groups/group_operation_action.jl")
        include_test("groups/rotation_action.jl")
        include_test("groups/translation_action.jl")
        include_test("groups/rotation_translation_action.jl")
        include_test("groups/connections.jl")
        include_test("groups/metric.jl")
    end
    if TEST_GROUP ∈ ["all", "test_integration"] && !Sys.isapple()
        include_test("recipes.jl")
    end
end
