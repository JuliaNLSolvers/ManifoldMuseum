include("../utils.jl")

using RecursiveArrayTools: ArrayPartition

@testset "Product manifold" begin
    @test_throws MethodError ProductManifold()
    M1 = Sphere(2)
    M2 = Euclidean(2)
    @test (@inferred ProductManifold(M1, M2)) isa ProductManifold
    Mse = ProductManifold(M1, M2)
    @test Mse == M1 × M2
    @test !is_flat(Mse)
    @test Mse == ProductManifold(M1) × M2
    @test Mse == ProductManifold(M1) × ProductManifold(M2)
    @test Mse == M1 × ProductManifold(M2)
    @test Mse[1] == M1
    @test Mse[2] == M2
    @test injectivity_radius(Mse) ≈ π
    @test injectivity_radius(
        Mse,
        ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
    ) ≈ π
    @test injectivity_radius(Mse, ExponentialRetraction()) ≈ π
    @test injectivity_radius(
        Mse,
        ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0]),
        ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
    ) ≈ π
    @test injectivity_radius(
        Mse,
        ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0]),
        ExponentialRetraction(),
    ) ≈ π
    @test is_default_metric(Mse, ProductMetric())

    @test Manifolds.number_of_components(Mse) == 2
    # test that arrays are not points
    @test_throws DomainError is_point(Mse, [1, 2]; error=:error)
    @test check_point(Mse, [1, 2]) isa DomainError
    @test_throws DomainError is_vector(Mse, 1, [1, 2]; error=:error, check_base_point=false)
    @test check_vector(Mse, 1, [1, 2]; check_base_point=false) isa DomainError
    #default fallbacks for check_size, Product not working with Arrays
    @test Manifolds.check_size(Mse, zeros(2)) isa DomainError
    @test Manifolds.check_size(Mse, zeros(2), zeros(3)) isa DomainError
    types = [Vector{Float64}]
    TEST_FLOAT32 && push!(types, Vector{Float32})
    TEST_STATIC_SIZED && push!(types, MVector{5,Float64})

    retraction_methods = [
        ProductRetraction(ExponentialRetraction(), ExponentialRetraction()),
        ExponentialRetraction(),
    ]
    inverse_retraction_methods = [
        InverseProductRetraction(
            LogarithmicInverseRetraction(),
            LogarithmicInverseRetraction(),
        ),
        LogarithmicInverseRetraction(),
    ]

    @testset "get_component, set_component!, getindex and setindex!" begin
        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        @test get_component(Mse, p1, 1) == p1.x[1]
        @test get_component(Mse, p1, Val(1)) == p1.x[1]
        @test p1[Mse, 1] == p1.x[1]
        @test p1[Mse, Val(1)] == p1.x[1]
        @test p1[Mse, 1] isa Vector
        @test p1[Mse, Val(1)] isa Vector
        p2 = [10.0, 12.0]
        set_component!(Mse, p1, p2, 2)
        @test get_component(Mse, p1, 2) == p2
        p1[Mse, 2] = 2 * p2
        @test p1[Mse, 2] == 2 * p2
        p3 = [11.0, 15.0]
        set_component!(Mse, p1, p3, Val(2))
        @test get_component(Mse, p1, Val(2)) == p3
        p1[Mse, Val(2)] = 2 * p3
        @test p1[Mse, Val(2)] == 2 * p3

        p1ap = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        @test get_component(Mse, p1ap, 1) == p1ap.x[1]
        @test get_component(Mse, p1ap, Val(1)) == p1ap.x[1]
        @test p1ap[Mse, 1] == p1ap.x[1]
        @test p1ap[Mse, Val(1)] == p1ap.x[1]
        @test p1ap[Mse, 1] isa Vector
        @test p1ap[Mse, Val(1)] isa Vector
        set_component!(Mse, p1ap, p2, 2)
        @test get_component(Mse, p1ap, 2) == p2
        p1ap[Mse, 2] = 2 * p2
        @test p1ap[Mse, 2] == 2 * p2
        p3 = [11.0, 15.0]
        set_component!(Mse, p1ap, p3, Val(2))
        @test get_component(Mse, p1ap, Val(2)) == p3
        p1ap[Mse, Val(2)] = 2 * p3
        @test p1ap[Mse, Val(2)] == 2 * p3

        p1c = copy(p1)
        p1c.x[1][1] = -123.0
        @test p1c.x[1][1] == -123.0
        @test p1.x[1][1] == 0.0
        copyto!(p1c, p1)
        @test p1c.x[1][1] == 0.0

        p1c.x[1][1] = -123.0
        copyto!(p1ap, p1c)
        @test p1ap.x[1][1] == -123.0
    end

    @testset "some ArrayPartition functions" begin
        p = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        q = allocate(p)
        @test q.x[1] isa Vector
        p = ArrayPartition([[0.0, 1.0, 0.0]], [0.0, 0.0])
        q = allocate(p, Int)
        @test q.x[1] isa Vector{Vector{Int}}
    end

    @testset "allocate on PowerManifold of ProductManifold" begin
        p = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        q = allocate([p])
        @test q[1] isa ArrayPartition
        @test q[1].x[1] isa Vector

        p = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        q = allocate([p])
        @test q[1] isa ArrayPartition
        @test q[1].x[1] isa Vector
    end

    @testset "Broadcasting" begin
        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0, 1.0])
        p2 = ArrayPartition([3.0, 4.0, 5.0], [2.0, 5.0])
        br_result = p1 .+ 2.0 .* p2
        @test br_result isa ArrayPartition
        @test br_result.x[1] ≈ [6.0, 9.0, 10.0]
        @test br_result.x[2] ≈ [4.0, 11.0]

        br_result .= 2.0 .* p1 .+ p2
        @test br_result.x[1] ≈ [3.0, 6.0, 5.0]
        @test br_result.x[2] ≈ [2.0, 7.0]

        br_result .= p1
        @test br_result.x[1] ≈ [0.0, 1.0, 0.0]
        @test br_result.x[2] ≈ [0.0, 1.0]

        @test axes(p1) == (Base.OneTo(5),)

        # errors
        p3 = ArrayPartition([3.0, 4.0, 5.0], [2.0, 5.0], [3.0, 2.0])
        @test_throws DimensionMismatch p1 .+ p3
        @test_throws DimensionMismatch p1 .= p3
    end

    @testset "CompositeManifoldError" begin
        Mpr = Sphere(2) × Sphere(2)
        p1 = [1.0, 0.0, 0.0]
        p2 = [0.0, 1.0, 0.0]
        X1 = [0.0, 1.0, 0.2]
        X2 = [1.0, 0.0, 0.2]
        p = ArrayPartition(p1, p2)
        X = ArrayPartition(X1, X2)
        pf = ArrayPartition(p1, X1)
        Xf = ArrayPartition(X1, p2)
        @test is_point(Mpr, p; error=:error)
        @test_throws CompositeManifoldError is_point(Mpr, X; error=:error)
        @test_throws ComponentManifoldError is_vector(Mpr, pf, X; error=:error)
        @test_throws ComponentManifoldError is_vector(Mpr, p, Xf; error=:error)
    end

    @testset "arithmetic" begin
        Mee = ProductManifold(Euclidean(3), Euclidean(2))
        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0, 1.0])
        p2 = ArrayPartition([1.0, 2.0, 0.0], [2.0, 3.0])

        @test isapprox(Mee, p1 + p2, ArrayPartition([1.0, 3.0, 0.0], [2.0, 4.0]))
        @test isapprox(Mee, p1 - p2, ArrayPartition([-1.0, -1.0, 0.0], [-2.0, -2.0]))
        @test isapprox(Mee, -p1, ArrayPartition([0.0, -1.0, 0.0], [0.0, -1.0]))
        @test isapprox(Mee, p1 * 2, ArrayPartition([0.0, 2.0, 0.0], [0.0, 2.0]))
        @test isapprox(Mee, 2 * p1, ArrayPartition([0.0, 2.0, 0.0], [0.0, 2.0]))
        @test isapprox(Mee, p1 / 2, ArrayPartition([0.0, 0.5, 0.0], [0.0, 0.5]))
    end

    @testset "Show methods" begin
        Mse2 = ProductManifold(M1, M1, M2, M2)
        @test sprint(show, Mse2) == "ProductManifold($(M1), $(M1), $(M2), $(M2))"
        withenv("LINES" => 10, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", ProductManifold(M1)) ==
                  "ProductManifold with 1 submanifold:\n $(M1)"
            @test sprint(show, "text/plain", Mse2) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n $(M1)\n $(M2)\n $(M2)"
            return nothing
        end
        withenv("LINES" => 7, "COLUMNS" => 100) do
            @test sprint(show, "text/plain", Mse2) ==
                  "ProductManifold with 4 submanifolds:\n $(M1)\n ⋮\n $(M2)"
            return nothing
        end

        @test sprint(show, "text/plain", ProductManifold(Mse, Mse)) == """
        ProductManifold with 2 submanifolds:
         ProductManifold(Sphere(2, ℝ), Euclidean(2; field=ℝ))
         ProductManifold(Sphere(2, ℝ), Euclidean(2; field=ℝ))"""
    end

    M3 = Rotations(2)
    Mser = ProductManifold(M1, M2, M3)

    @test submanifold(Mser, 2) == M2
    @test (@inferred submanifold(Mser, Val((1, 3)))) == M1 × M3
    @test submanifold(Mser, 2:3) == M2 × M3
    @test submanifold(Mser, [1, 3]) == M1 × M3

    pts_sphere = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    pts_r2 = [[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]]
    angles = (0.0, π / 2, 2π / 3)
    pts_rot = [[cos(ϕ) sin(ϕ); -sin(ϕ) cos(ϕ)] for ϕ in angles]
    pts = [ArrayPartition(p[1], p[2], p[3]) for p in zip(pts_sphere, pts_r2, pts_rot)]
    test_manifold(
        Mser,
        pts,
        test_injectivity_radius=false,
        is_tangent_atol_multiplier=1,
        exp_log_atol_multiplier=1,
        test_inplace=true,
        test_rand_point=true,
        test_rand_tvector=true,
    )

    @testset "product vector transport" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0])
        q = ArrayPartition([0.0, 1.0, 0.0], [2.0, 0.0])
        X = log(Mse, p, q)
        m = ProductVectorTransport(ParallelTransport(), ParallelTransport())
        Y = vector_transport_to(Mse, p, X, q, m)
        Z = -log(Mse, q, p)
        @test isapprox(Mse, q, Y, Z)
    end

    @testset "Implicit product vector transport" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0])
        q = ArrayPartition([0.0, 1.0, 0.0], [2.0, 0.0])
        X = log(Mse, p, q)
        for m in [ParallelTransport(), SchildsLadderTransport(), PoleLadderTransport()]
            Y = vector_transport_to(Mse, p, X, q, m)
            Z1 = vector_transport_to(
                Mse.manifolds[1],
                submanifold_component.([p, X, q], Ref(1))...,
                m,
            )
            Z2 = vector_transport_to(
                Mse.manifolds[2],
                submanifold_component.([p, X, q], Ref(2))...,
                m,
            )
            Z = ArrayPartition(Z1, Z2)
            @test isapprox(Mse, q, Y, Z)
            Y2 = allocate(Mse, Y)
            vector_transport_to!(Mse, Y2, p, X, q, m)
            @test isapprox(Mse, q, Y2, Z)
        end
        for m in [ParallelTransport(), SchildsLadderTransport(), PoleLadderTransport()]
            Y = vector_transport_direction(Mse, p, X, X, m)
            Z1 = vector_transport_direction(
                Mse.manifolds[1],
                submanifold_component.([p, X, X], Ref(1))...,
                m,
            )
            Z2 = vector_transport_direction(
                Mse.manifolds[2],
                submanifold_component.([p, X, X], Ref(2))...,
                m,
            )
            Z = ArrayPartition(Z1, Z2)
            @test isapprox(Mse, q, Y, Z)
        end
    end
    @testset "Parallel transport" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0])
        q = ArrayPartition([0.0, 1.0, 0.0], [2.0, 0.0])
        X = log(Mse, p, q)
        # to
        Y = parallel_transport_to(Mse, p, X, q)
        Z1 = parallel_transport_to(
            Mse.manifolds[1],
            submanifold_component.([p, X, q], Ref(1))...,
        )
        Z2 = parallel_transport_to(
            Mse.manifolds[2],
            submanifold_component.([p, X, q], Ref(2))...,
        )
        Z = ArrayPartition(Z1, Z2)
        @test isapprox(Mse, q, Y, Z)
        Ym = allocate(Y)
        parallel_transport_to!(Mse, Ym, p, X, q)
        @test isapprox(Mse, q, Y, Z)

        # direction
        Y = parallel_transport_direction(Mse, p, X, X)
        Z1 = parallel_transport_direction(
            Mse.manifolds[1],
            submanifold_component.([p, X, X], Ref(1))...,
        )
        Z2 = parallel_transport_direction(
            Mse.manifolds[2],
            submanifold_component.([p, X, X], Ref(2))...,
        )
        Z = ArrayPartition(Z1, Z2)
        @test isapprox(Mse, q, Y, Z)
        Ym = allocate(Y)
        parallel_transport_direction!(Mse, Ym, p, X, X)
        @test isapprox(Mse, q, Ym, Z)
    end

    @testset "ArrayPartition" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0])
        @test submanifold_component(Mse, p, 1) === p.x[1]
        @test submanifold_component(Mse, p, Val(1)) === p.x[1]
        @test submanifold_component(p, 1) === p.x[1]
        @test submanifold_component(p, Val(1)) === p.x[1]
        @test submanifold_components(Mse, p) === p.x
        @test submanifold_components(p) === p.x
    end

    @testset "manifold tests (static size)" begin
        Ts = SizedVector{3,Float64}
        Tr2 = SizedVector{2,Float64}
        pts_sphere = [
            convert(Ts, [1.0, 0.0, 0.0]),
            convert(Ts, [0.0, 1.0, 0.0]),
            convert(Ts, [0.0, 0.0, 1.0]),
        ]
        pts_r2 =
            [convert(Tr2, [0.0, 0.0]), convert(Tr2, [1.0, 0.0]), convert(Tr2, [0.0, 0.1])]

        pts = [ArrayPartition(p[1], p[2]) for p in zip(pts_sphere, pts_r2)]
        basis_types = (
            DefaultOrthonormalBasis(),
            ProjectedOrthonormalBasis(:svd),
            get_basis(Mse, pts[1], DefaultOrthonormalBasis()),
            DiagonalizingOrthonormalBasis(
                ArrayPartition(SizedVector{3}([0.0, 1.0, 0.0]), SizedVector{2}([1.0, 0.0])),
            ),
        )
        distr_M1 = Manifolds.uniform_distribution(M1, pts_sphere[1])
        distr_M2 = Manifolds.projected_distribution(
            M2,
            Distributions.MvNormal(zero(pts_r2[1]), 1.0 * I),
        )
        distr_tv_M1 = Manifolds.normal_tvector_distribution(M1, pts_sphere[1], 1.0)
        distr_tv_M2 = Manifolds.normal_tvector_distribution(M2, pts_r2[1], 1.0)
        @test injectivity_radius(Mse, pts[1]) ≈ π
        @test injectivity_radius(Mse) ≈ π
        @test injectivity_radius(Mse, pts[1], ExponentialRetraction()) ≈ π
        @test injectivity_radius(Mse, ExponentialRetraction()) ≈ π

        @test ManifoldsBase.allocate_coordinates(
            Mse,
            pts[1],
            Float64,
            number_of_coordinates(Mse, DefaultOrthogonalBasis()),
        ) isa Vector{Float64}

        Y = allocate(pts[1])
        inverse_retract!(Mse, Y, pts[1], pts[2], default_inverse_retraction_method(Mse))
        @test isapprox(
            Mse,
            pts[1],
            Y,
            inverse_retract(Mse, pts[1], pts[2], default_inverse_retraction_method(Mse)),
        )

        test_manifold(
            Mse,
            pts;
            point_distributions=[Manifolds.ProductPointDistribution(distr_M1, distr_M2)],
            tvector_distributions=[
                Manifolds.ProductFVectorDistribution(distr_tv_M1, distr_tv_M2),
            ],
            test_injectivity_radius=true,
            test_musical_isomorphisms=true,
            musical_isomorphism_bases=[DefaultOrthonormalBasis()],
            test_tangent_vector_broadcasting=true,
            test_project_tangent=true,
            test_project_point=true,
            test_mutating_rand=true,
            retraction_methods=retraction_methods,
            inverse_retraction_methods=inverse_retraction_methods,
            test_riesz_representer=true,
            test_default_vector_transport=true,
            test_rand_point=true,
            test_rand_tvector=true,
            vector_transport_methods=[
                ProductVectorTransport(ParallelTransport(), ParallelTransport()),
                ProductVectorTransport(SchildsLadderTransport(), SchildsLadderTransport()),
                ProductVectorTransport(PoleLadderTransport(), PoleLadderTransport()),
            ],
            basis_types_vecs=(basis_types[1], basis_types[3], basis_types[4]),
            basis_types_to_from=basis_types,
            is_tangent_atol_multiplier=1,
            exp_log_atol_multiplier=1,
        )
        @test number_eltype(pts[1]) === Float64

        @test (@inferred ManifoldsBase._get_vector_cache_broadcast(pts[1])) === Val(false)
    end

    @testset "vee/hat" begin
        M1 = Rotations(3)
        M2 = Euclidean(3)
        M = M1 × M2

        e = Matrix{Float64}(I, 3, 3)
        p = ArrayPartition(exp(M1, e, hat(M1, e, [1.0, 2.0, 3.0])), [1.0, 2.0, 3.0])
        X = [0.1, 0.2, 0.3, -1.0, 2.0, -3.0]

        Xc = hat(M, p, X)
        X2 = vee(M, p, Xc)
        @test isapprox(X, X2)
    end

    @testset "get_coordinates" begin
        # make sure `get_coordinates` does not return an `ArrayPartition`
        p1 = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        X1 = ArrayPartition([1.0, 0.0, -1.0], [1.0, 0.0])
        Tp1Mse = TangentSpace(Mse, p1)
        c = get_coordinates(Tp1Mse, p1, X1, DefaultOrthonormalBasis())
        @test c isa Vector

        p1ap = ArrayPartition([0.0, 1.0, 0.0], [0.0, 0.0])
        X1ap = ArrayPartition([1.0, 0.0, -1.0], [1.0, 0.0])
        Tp1apMse = TangentSpace(Mse, p1ap)
        cap = get_coordinates(Tp1apMse, p1ap, X1ap, DefaultOrthonormalBasis())
        @test cap isa Vector
    end

    @testset "Basis printing" begin
        p = ArrayPartition([1.0, 0.0, 0.0], [1.0, 0.0])
        B = DefaultOrthonormalBasis()
        Bc = get_basis(Mse, p, B)
        Bc_components_s = sprint.(show, "text/plain", Bc.data.parts)
        @test sprint(show, "text/plain", Bc) == """
        $(typeof(B)) for a product manifold
        Basis for component 1:
        $(Bc_components_s[1])
        Basis for component 2:
        $(Bc_components_s[2])
        """
    end

    @testset "Basis-related errors" begin
        a = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0])
        B = CachedBasis(DefaultOrthonormalBasis(), ProductBasisData(([],)))
        @test_throws AssertionError get_vector!(
            Mse,
            a,
            ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0]),
            [1.0, 2.0, 3.0, 4.0, 5.0], # this is one element too long, hence assertionerror
            B,
        )
        @test_throws MethodError get_vector!(
            Mse,
            a,
            ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0]),
            [1.0, 2.0, 3.0, 4.0],
            B, # empty elements yield a submanifold MethodError
        )
    end

    @testset "allocation promotion" begin
        M2c = Euclidean(2; field=ℂ)
        Msec = ProductManifold(M1, M2c)
        @test Manifolds.allocation_promotion_function(Msec, get_vector, ()) === complex
        @test Manifolds.allocation_promotion_function(Mse, get_vector, ()) === identity
    end

    @testset "empty allocation" begin
        p = allocate_result(Mse, uniform_distribution)
        @test isa(p, ArrayPartition)
        @test size(p[Mse, 1]) == (3,)
        @test size(p[Mse, 2]) == (2,)
    end

    @testset "Uniform distribution" begin
        Mss = ProductManifold(Sphere(2), Sphere(2))
        p = rand(uniform_distribution(Mss))
        @test is_point(Mss, p)
        @test is_point(Mss, rand(uniform_distribution(Mss, p)))
    end

    @testset "Atlas & Induced Basis" begin
        M = ProductManifold(Euclidean(2), Euclidean(2))
        p = ArrayPartition(zeros(2), ones(2))
        X = ArrayPartition(ones(2), 2 .* ones(2))
        A = RetractionAtlas()
        a = get_parameters(M, A, p, p)
        p2 = get_point(M, A, p, a)
        @test all(submanifold_components(p2) .== submanifold_components(p))
    end

    @testset "metric conversion" begin
        M = SymmetricPositiveDefinite(3)
        N = ProductManifold(M, M)
        e = EuclideanMetric()
        p = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1]
        q = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1]
        P = ArrayPartition(p, q)
        X = ArrayPartition(log(M, p, q), log(M, q, p))
        Y = change_metric(N, e, P, X)
        Yc = ArrayPartition(
            change_metric(M, e, p, log(M, p, q)),
            change_metric(M, e, q, log(M, q, p)),
        )
        @test norm(N, P, Y - Yc) ≈ 0
        Z = change_representer(N, e, P, X)
        Zc = ArrayPartition(
            change_representer(M, e, p, log(M, p, q)),
            change_representer(M, e, q, log(M, q, p)),
        )
        @test norm(N, P, Z - Zc) ≈ 0
    end

    @testset "default retraction, inverse retraction and VT" begin
        Mstb = ProductManifold(M1, TangentBundle(M1))
        T_p_ap = ArrayPartition{
            Float64,
            Tuple{
                Matrix{Float64},
                ArrayPartition{Float64,Tuple{Matrix{Float64},Matrix{Float64}}},
            },
        }
        @test Manifolds.default_retraction_method(Mstb) === ProductRetraction(
            ExponentialRetraction(),
            Manifolds.FiberBundleProductRetraction(),
        )
        @test Manifolds.default_retraction_method(Mstb, T_p_ap) === ProductRetraction(
            ExponentialRetraction(),
            Manifolds.FiberBundleProductRetraction(),
        )

        @test Manifolds.default_inverse_retraction_method(Mstb) ===
              Manifolds.InverseProductRetraction(
            LogarithmicInverseRetraction(),
            Manifolds.FiberBundleInverseProductRetraction(),
        )
        @test Manifolds.default_inverse_retraction_method(Mstb, T_p_ap) ===
              Manifolds.InverseProductRetraction(
            LogarithmicInverseRetraction(),
            Manifolds.FiberBundleInverseProductRetraction(),
        )

        @test Manifolds.default_vector_transport_method(Mstb) === ProductVectorTransport(
            ParallelTransport(),
            Manifolds.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
        @test Manifolds.default_vector_transport_method(Mstb, T_p_ap) ===
              ProductVectorTransport(
            ParallelTransport(),
            Manifolds.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
        @test Manifolds.default_vector_transport_method(Mstb, T_p_ap) ===
              ProductVectorTransport(
            ParallelTransport(),
            Manifolds.FiberBundleProductVectorTransport(
                ParallelTransport(),
                ParallelTransport(),
            ),
        )
    end

    @testset "Riemann tensor" begin
        p = ArrayPartition([0.0, 1.0, 0.0], [2.0, 3.0])
        X = ArrayPartition([1.0, 0.0, 0.0], [2.0, 3.0])
        Y = ArrayPartition([0.0, 0.0, 3.0], [-2.0, 3.0])
        Z = ArrayPartition([-1.0, 0.0, 2.0], [2.0, -3.0])
        Xresult = ArrayPartition([6.0, 0.0, 3.0], [0.0, 0.0])
        @test isapprox(riemann_tensor(Mse, p, X, Y, Z), Xresult)
        Xresult2 = allocate(Xresult)
        riemann_tensor!(Mse, Xresult2, p, X, Y, Z)
        @test isapprox(Xresult2, Xresult)
    end

    @testset "ManifoldDiff" begin
        p = ArrayPartition([0.0, 1.0, 0.0], [2.0, 3.0])
        q = ArrayPartition([1.0, 0.0, 0.0], [-2.0, 3.0])
        X = ArrayPartition([1.0, 0.0, 0.0], [2.0, 3.0])
        # ManifoldDiff
        @test ManifoldDiff.adjoint_Jacobi_field(
            Mse,
            p,
            q,
            0.5,
            X,
            ManifoldDiff.βdifferential_shortest_geodesic_startpoint,
        ) == ArrayPartition([0.5, 0.0, 0.0], [1.0, 1.5])
        X2 = allocate(X)
        ManifoldDiff.adjoint_Jacobi_field!(
            Mse,
            X2,
            p,
            q,
            0.5,
            X,
            ManifoldDiff.βdifferential_shortest_geodesic_startpoint,
        )
        @test X2 == ArrayPartition([0.5, 0.0, 0.0], [1.0, 1.5])
        @test ManifoldDiff.jacobi_field(
            Mse,
            p,
            q,
            0.5,
            X,
            ManifoldDiff.βdifferential_shortest_geodesic_startpoint,
        ) == ArrayPartition([0.3535533905932738, -0.35355339059327373, 0.0], [1.0, 1.5])
        X2 = allocate(X)
        ManifoldDiff.jacobi_field!(
            Mse,
            X2,
            p,
            q,
            0.5,
            X,
            ManifoldDiff.βdifferential_shortest_geodesic_startpoint,
        )
        @test X2 ==
              ArrayPartition([0.3535533905932738, -0.35355339059327373, 0.0], [1.0, 1.5])
    end

    @testset "Hessian conversion" begin
        M = Sphere(2)
        N = M × M
        p = ArrayPartition([1.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        q = 1 / sqrt(2) * ArrayPartition([1.0, 1.0, 0.0], [0.0, 1.0, 1.0])
        q = 1 / sqrt(2) * ArrayPartition([0.0, 1.0, 1.0], [1.0, 1.0, 0.0])
        r = 1 / sqrt(3) * ArrayPartition([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
        X = log(M, p, q)
        Y = log(M, p, r)
        Z = -X
        H1 = riemannian_Hessian(N, p, Y, Z, X)
        H2 = ArrayPartition(
            [riemannian_Hessian(M, p[i, :], Y[i, :], Z[i, :], X[i, :]) for i in 1:2]...,
        )
        @test H1 == H2
        V = ArrayPartition([0.2, 0.0, 0.0], [0.0, 0.0, 0.3])
        W1 = Weingarten(N, p, X, V)
        W2 = ArrayPartition([Weingarten(M, p[i, :], X[i, :], V[i, :]) for i in 1:2]...)
        @test W1 == W2
    end
    @testset "Manifold volume" begin
        MS2 = Sphere(2)
        MS3 = Sphere(3)
        PM = ProductManifold(MS2, MS3)
        @test manifold_volume(PM) ≈ manifold_volume(MS2) * manifold_volume(MS3)
        p1 = [-0.9171596991960276, 0.39792260844341604, -0.02181017790481868]
        p2 = [
            -0.5427653626654726,
            5.420303965772687e-5,
            -0.8302022885580579,
            -0.12716099333369416,
        ]
        X1 = [-0.35333565579879633, -0.7896159441709865, 0.45204526334685574]
        X2 = [
            -0.33940201562492356,
            0.8092470417550779,
            0.18290591742514573,
            0.2548785571950708,
        ]
        @test volume_density(PM, ArrayPartition(p1, p2), ArrayPartition(X1, X2)) ≈
              volume_density(MS2, p1, X1) * volume_density(MS3, p2, X2)
    end
end
