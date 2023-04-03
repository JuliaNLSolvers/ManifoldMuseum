include("../utils.jl")

@testset "Grassmann" begin
    @testset "Real" begin
        M = Grassmann(3, 2)
        @testset "Basics" begin
            @test repr(M) == "Grassmann(3, 2, ℝ)"
            @test representation_size(M) == (3, 2)
            @test manifold_dimension(M) == 2
            @test !is_flat(M)
            @test is_flat(Grassmann(2, 1))
            @test get_total_space(M) == Stiefel(3, 2, ℝ)
            @test get_orbit_action(M) ==
                  Manifolds.RowwiseMultiplicationAction(M, Orthogonal(2))
            @test !is_point(M, [1.0, 0.0, 0.0, 0.0])
            @test !is_vector(M, [1.0 0.0; 0.0 1.0; 0.0 0.0], [0.0, 0.0, 1.0, 0.0])
            @test_throws ManifoldDomainError is_point(M, [2.0 0.0; 0.0 1.0; 0.0 0.0], true)
            @test_throws ManifoldDomainError is_vector(
                M,
                [2.0 0.0; 0.0 1.0; 0.0 0.0],
                zeros(3, 2),
                true,
            )
            @test_throws ManifoldDomainError is_vector(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                ones(3, 2),
                true,
            )
            @test is_point(M, [1.0 0.0; 0.0 1.0; 0.0 0.0], true)
            @test_throws ManifoldDomainError is_point(
                M,
                1im * [1.0 0.0; 0.0 1.0; 0.0 0.0],
                true,
            )
            @test is_vector(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                zero_vector(M, [1.0 0.0; 0.0 1.0; 0.0 0.0]),
                true,
            )
            @test_throws ManifoldDomainError is_vector(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                1im * zero_vector(M, [1.0 0.0; 0.0 1.0; 0.0 0.0]),
                true,
            )
            @test injectivity_radius(M) == π / 2
            @test injectivity_radius(M, ExponentialRetraction()) == π / 2
            @test injectivity_radius(M, [1.0 0.0; 0.0 1.0; 0.0 0.0]) == π / 2
            @test injectivity_radius(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                ExponentialRetraction(),
            ) == π / 2
        end
        types = [Matrix{Float64}]
        TEST_STATIC_SIZED && push!(types, MMatrix{3,2,Float64,6})

        TEST_FLOAT32 && push!(types, Matrix{Float32})
        basis_types = (ProjectedOrthonormalBasis(:gram_schmidt),)
        @testset "Type $T" for T in types
            p1 = [1.0 0.0; 0.0 1.0; 0.0 0.0]
            X = [0.0 0.0; 0.0 0.0; 0.0 1.0]
            p2 = exp(M, p1, X)
            Y = [0.0 1.0; -1.0 0.0; 1.0 0.0]
            p3 = exp(M, p1, Y)
            pts = convert.(T, [p1, p2, p3])
            test_manifold(
                M,
                pts,
                test_exp_log=true,
                test_injectivity_radius=false,
                test_project_tangent=true,
                test_project_point=true,
                test_default_vector_transport=false,
                point_distributions=[Manifolds.uniform_distribution(M, pts[1])],
                test_vee_hat=false,
                test_rand_point=true,
                test_rand_tvector=true,
                retraction_methods=[PolarRetraction(), QRRetraction()],
                inverse_retraction_methods=[
                    PolarInverseRetraction(),
                    QRInverseRetraction(),
                ],
                #basis_types_vecs = basis_types,
                # investigate why this is so large on dev
                exp_log_atol_multiplier=10.0 * (VERSION >= v"1.6-DEV" ? 10.0^8 : 1.0),
                is_tangent_atol_multiplier=20.0,
                projection_atol_multiplier=10.0,
                retraction_atol_multiplier=10.0,
            )

            @testset "inner/norm" begin
                X1 = inverse_retract(M, pts[1], pts[2], PolarInverseRetraction())
                X2 = inverse_retract(M, pts[1], pts[3], PolarInverseRetraction())

                @test real(inner(M, pts[1], X1, X2)) ≈ real(inner(M, pts[1], X2, X1))
                @test imag(inner(M, pts[1], X1, X2)) ≈ -imag(inner(M, pts[1], X2, X1))
                @test imag(inner(M, pts[1], X1, X1)) ≈ 0

                @test norm(M, pts[1], X1) isa Real
                @test norm(M, pts[1], X1) ≈ sqrt(inner(M, pts[1], X1, X1))
            end
            @test riemann_tensor(M, p1, X, Y, 2 * X + Y) ≈ [0 -2; 0 1; 2 0]
            @testset "gradient and metric conversion" begin
                Y = change_metric(M, EuclideanMetric(), p1, X)
                @test Y == X
                Z = change_representer(M, EuclideanMetric(), p1, X)
                @test Z == X
            end
        end

        @testset "Distribution tests" begin
            ugd_mmatrix = Manifolds.uniform_distribution(M, @MMatrix [
                1.0 0.0
                0.0 1.0
                0.0 0.0
            ])
            @test isa(rand(ugd_mmatrix), MMatrix)
        end

        @testset "vector transport" begin
            p1 = [1.0 0.0; 0.0 1.0; 0.0 0.0]
            X = [0.0 0.0; 0.0 0.0; 0.0 1.0]
            p2 = exp(M, p1, X)
            @test vector_transport_to(M, p1, X, p2, ProjectionTransport()) ==
                  project(M, p2, X)
            @test is_vector(
                M,
                p2,
                vector_transport_to(M, p1, X, p2, ProjectionTransport()),
                true;
                atol=10^-15,
            )
        end

        @testset "default_* functions" begin
            p = [1.0 0.0; 0.0 1.0; 0.0 0.0]
            pS = StiefelPoint(p)
            @test default_vector_transport_method(M, typeof(p)) == ProjectionTransport()
            @test default_vector_transport_method(M, typeof(pS)) == ProjectionTransport()
        end
    end

    @testset "Complex" begin
        M = Grassmann(3, 2, ℂ)
        @testset "Basics" begin
            @test repr(M) == "Grassmann(3, 2, ℂ)"
            @test representation_size(M) == (3, 2)
            @test manifold_dimension(M) == 4
            @test !is_flat(M)
            @test !is_point(M, [1.0, 0.0, 0.0, 0.0])
            @test !is_vector(M, [1.0 0.0; 0.0 1.0; 0.0 0.0], [0.0, 0.0, 1.0, 0.0])
            @test Manifolds.allocation_promotion_function(M, exp!, (1,)) == complex
            @test_throws ManifoldDomainError is_point(M, [2.0 0.0; 0.0 1.0; 0.0 0.0], true)
            @test_throws ManifoldDomainError is_vector(
                M,
                [2.0 0.0; 0.0 1.0; 0.0 0.0],
                zeros(3, 2),
                true,
            )
            @test_throws ManifoldDomainError is_vector(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                ones(3, 2),
                true,
            )
            @test is_vector(
                M,
                [1.0 0.0; 0.0 1.0; 0.0 0.0],
                1im * zero_vector(M, [1.0 0.0; 0.0 1.0; 0.0 0.0]),
            )
            @test is_point(M, [1.0 0.0; 0.0 1.0; 0.0 0.0])
            @test injectivity_radius(M) == π / 2
        end
        types = [Matrix{ComplexF64}]
        @testset "Type $T" for T in types
            p1 = [0.5+0.5im 0.5+0.5im; 0.5+0.5im -0.5-0.5im; 0.0 0.0]
            X = [0.0 0.0; 0.0 0.0; 0.0 1.0]
            p2 = exp(M, p1, X)
            Y = [0.0 1.0; -1.0 0.0; 1.0 0.0]
            p3 = exp(M, p1, Y)
            pts = convert.(T, [p1, p2, p3])
            test_manifold(
                M,
                pts,
                test_exp_log=true,
                test_injectivity_radius=false,
                test_project_tangent=true,
                test_project_point=true,
                test_default_vector_transport=false,
                test_vee_hat=false,
                retraction_methods=[PolarRetraction(), QRRetraction()],
                inverse_retraction_methods=[
                    PolarInverseRetraction(),
                    QRInverseRetraction(),
                ],
                exp_log_atol_multiplier=10.0^3,
                is_tangent_atol_multiplier=20.0,
                projection_atol_multiplier=10.0,
                retraction_atol_multiplier=10.0,
                test_inplace=true,
                test_rand_point=true,
            )

            @testset "inner/norm" begin
                X1 = inverse_retract(M, pts[1], pts[2], PolarInverseRetraction())
                X2 = inverse_retract(M, pts[1], pts[3], PolarInverseRetraction())

                @test real(inner(M, pts[1], X1, X2)) ≈ real(inner(M, pts[1], X2, X1))
                @test imag(inner(M, pts[1], X1, X2)) ≈ -imag(inner(M, pts[1], X2, X1))
                @test isapprox(imag(inner(M, pts[1], X1, X1)), 0; atol=1e-30)

                @test norm(M, pts[1], X1) isa Real
                @test norm(M, pts[1], X1) ≈ sqrt(inner(M, pts[1], X1, X1))
            end
        end
    end

    @testset "Complex and conjugate" begin
        G = Grassmann(3, 1, ℂ)
        p = reshape([im, 0.0, 0.0], 3, 1)
        @test is_point(G, p)
        X = reshape([-0.5; 0.5; 0], 3, 1)
        @test_throws ManifoldDomainError is_vector(G, p, X, true)
        Y = project(G, p, X)
        @test is_vector(G, p, Y)
    end

    @testset "Projector representation" begin
        M = Grassmann(3, 2)
        p = ProjectorPoint([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0])
        X = ProjectorTVector([0.0 0.0 1.0; 0.0 0.0 1.0; 1.0 1.0 0.0])
        pS = StiefelPoint([1.0 0.0; 0.0 1.0; 0.0 0.0])
        Xs = StiefelTVector([0.0 1.0; -1.0 0.0; 0.0 0.0])
        @test representation_size(M, p) == (3, 3)

        q = embed(M, p)
        @test q == p.value
        q2 = similar(q)
        embed!(M, q2, p)
        @test q2 == p.value
        q3 = similar(q)
        embed!(M, q3, p.value)
        @test q3 == p.value
        Y = embed(M, p, X)
        @test Y == X.value
        Y2 = similar(Y)
        embed!(M, Y2, p, X)
        @test Y2 == X.value
        Y3 = similar(Y)
        embed!(M, Y3, p.value, X.value)
        @test Y3 == X.value

        pSe = similar(pS.value)
        embed!(M, pSe, pS)
        @test pSe == pS.value
        Xse = similar(Xs.value)
        embed!(M, Xse, pS, Xs)
        Xse == Xs.value

        p2 = ProjectorPoint(similar(p.value))
        pC = pS.value * pS.value'
        canonical_project!(M, p2, pS)
        @test p2.value == pC
        p3 = ProjectorPoint(similar(p.value))
        canonical_project!(M, p2, pS.value)
        @test p2.value == pC
        p3 = canonical_project(M, pS)
        @test p3.value == pC

        Y = ProjectorTVector(similar(X.value))
        Yc = Xs.value * pS.value' + pS.value * Xs.value'
        differential_canonical_project!(M, Y, pS, Xs)
        @test Y.value == Yc
        Y2 = ProjectorTVector(similar(X.value))
        differential_canonical_project!(M, Y2, pS.value, Xs.value)
        @test Y2.value == Yc
        Y3 = differential_canonical_project(M, pS, Xs)
        @test Y3.value == Yc
        Y4 = differential_canonical_project(M, pS.value, Xs.value)
        @test Y4.value == Yc

        @test horizontal_lift(Stiefel(3, 2), pS.value, X) == X.value[:, 1:2]

        exppx = exp(X.value * p.value - p.value * X.value)
        qc = exppx * p.value / exppx
        q = exp(M, p, X)
        @test qc == q.value

        d = -X
        edppd = exp(d.value * p.value - p.value * d.value)
        Yc2 = edppd * X.value / edppd
        Xp = parallel_transport_direction(M, p, X, d)
        @test Xp.value == Yc2
    end

    @testset "is_point & convert & show" begin
        M = Grassmann(3, 2)
        p = StiefelPoint([1.0 0.0; 0.0 1.0; 0.0 0.0])
        X = StiefelTVector([0.0 1.0; -1.0 0.0; 0.0 0.0])
        @test is_point(M, p, true)
        @test is_vector(M, p, X, true)
        @test repr(p) == "StiefelPoint($(p.value))"
        @test repr(X) == "StiefelTVector($(X.value))"
        M2 = Stiefel(3, 2)
        @test is_point(M2, p, true)
        @test is_vector(M2, p, X, true)

        p2 = convert(ProjectorPoint, p)
        @test is_point(M, p2, true)
        p3 = convert(ProjectorPoint, p.value)
        @test p2.value == p3.value
        X2 = ProjectorTVector([0.0 0.0 1.0; 0.0 0.0 1.0; 1.0 1.0 0.0])
        @test is_vector(M, p2, X2)
        @test repr(p2) == "ProjectorPoint($(p2.value))"
        @test repr(X2) == "ProjectorTVector($(X2.value))"

        # rank just 1
        pF1 = ProjectorPoint([1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
        @test_throws DomainError is_point(M, pF1, true)
        # not equal to its square
        pF2 = ProjectorPoint([1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 1.0 0.0])
        @test_throws DomainError is_point(M, pF2, true)
        # not symmetric
        pF3 = ProjectorPoint([0.0 1.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0])
        @test_throws DomainError is_point(M, pF3, true)

        # not symmetric
        XF1 = ProjectorTVector([0.0 0.0 1.0; 0.0 0.0 1.0; 1.0 0.0 0.0])
        @test_throws DomainError is_vector(M, p2, XF1, true)
        # XF2 is not p2*XF2 + XF2*p2
        XF2 = ProjectorTVector(ones(3, 3))
        @test_throws DomainError is_vector(M, p2, XF2, true)

        # embed for Stiefel with its point
        M2 = Stiefel(3, 2)
        q = embed(M2, p)
        @test q == p.value
        q2 = similar(q)
        embed!(M2, q2, p)
        q3 = similar(q)
        @test q2 == p.value
        embed!(M2, q3, p.value)
        @test q3 == p.value
        Y = embed(M2, p, X)
        @test Y == X.value
        Y2 = similar(Y)
        embed!(M2, Y2, p, X)
        @test Y2 == X.value
        Y3 = similar(Y)
        embed!(M2, Y3, p.value, X.value)
        @test Y3 == X.value
    end

    @testset "small distance tests" begin
        n, k = 5, 3
        @testset for fT in (Float32, Float64), T in (fT, Complex{fT})
            𝔽 = T isa Complex ? ℂ : ℝ
            M = Grassmann(n, k, 𝔽)
            U = Unitary(k, 𝔽)
            rT = real(T)
            atol = rtol = sqrt(eps(rT))
            @testset for t in (zero(rT), eps(rT)^(1 // 4) / 8, eps(rT)^(1 // 4)),
                z in (I, rand(U))

                p = project(M, randn(T, representation_size(M)))
                X = project(M, p, randn(T, representation_size(M)))
                X ./= norm(M, p, X)
                project!(M, X, p, X)
                @test distance(M, p, exp(M, p, t * X) * z) ≈ t atol = atol rtol = rtol
            end
        end
    end
end
