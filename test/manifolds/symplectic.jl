include("../utils.jl")

@testset "Symplectic" begin
    @testset "Real" begin
        small_n = 1
        M = Symplectic(2 * small_n)
        M2 = MetricManifold(M, EuclideanMetric())
        p = [0.0 1.0/2.0; -2.0 -2.0]
        X1 = [
            -0.121212 0.121212
            0.969697 -1.0
        ]
        X2 = [
            0.0 0.0
            0.0 -1.0
        ]

        @testset "Basics" begin
            @test repr(M) == "Symplectic($(2*small_n), ℝ)"
            @test representation_size(M) == (2small_n, 2small_n)
            @test base_manifold(M) === M

            @test is_point(M, p)
            @test !is_point(M, p + I)

            @test is_vector(M, p, X1; atol=1.0e-6)
            @test is_vector(M, p, X2; atol=1.0e-12)
            @test is_vector(M, p, X1 + X2; atol=1.0e-6)
            @test !is_vector(M, p, X1 + [0.1 0.1; -0.1 0.1]; atol=1.0e-6)
        end
        @testset "Embedding and Projection" begin
            x = [0.0 1.0/2.0; -2.0 -2.0]
            y = similar(x)
            z = embed(M, x)
            @test z == x

            Y = similar(X1)
            embed!(M, Y, p, X1)
            @test Y == X1
        end
        @testset "Retractions and Exponential Mapping" begin
            q_exp = [
                -0.0203171 0.558648
                -1.6739 -3.19344
            ]
            @test isapprox(exp(M, p, X2), q_exp; atol=1.0e-5)
            @test isapprox(retract(M, p, X2, ExponentialRetraction()), q_exp; atol=1.0e-5)

            q_cay = [
                0.0 0.5
                -2.0 -3.0
            ]
            @test retract(M, p, X2) == q_cay
            @test retract(M, p, X2, CayleyRetraction()) == q_cay

            X_inv_cayley_retraction = inverse_retract(M, p, q_cay)
            X_inv_cayley_retraction_2 =
                inverse_retract(M, p, q_cay, CayleyInverseRetraction())
            @test X_inv_cayley_retraction == X_inv_cayley_retraction_2
            @test X_inv_cayley_retraction ≈ X2
        end

        @testset "Riemannian metric" begin
            X1_p_norm = 0.49259905148939337
            @test norm(M, p, X1) == X1_p_norm
            @test norm(M, p, X1) == √(inner(M, p, X1, X1))

            X2_p_norm = 1 / 2
            @test norm(M, p, X2) == X2_p_norm
            @test norm(M, p, X2) == √(inner(M, p, X2, X2))
        end
    end

    # TODO: test_manifold()
end