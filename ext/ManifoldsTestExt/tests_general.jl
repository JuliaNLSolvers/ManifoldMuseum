"""
    find_eps(x...)

Find an appropriate tolerance for given points or tangent vectors, or their types.
"""
find_eps(x...) = find_eps(Base.promote_type(map(number_eltype, x)...))
find_eps(x::Type{TN}) where {TN<:Number} = eps(real(TN))
find_eps(x) = find_eps(number_eltype(x))

function test_manifold(
    M::AbstractManifold;
    points=[rand(M) for _ in 1:3],
    tangent_vectors=[rand(M; vector_at=p) for p in points],
    features=ManifoldFeatures(M),
    expectations=ManifoldExpectations(),
    non_points=[],
    non_tangent_vectors=[],
    # Errors should be in expectations as well
)
    length(points) ≥ 3 || error("$(length(points)) are not enough points, 3 required")
    length(tangent_vectors) == length(points) || error(
        "$(length(tangent_vectors)) has to be the same as number of points $(length(points)).",
    )
    Test.@testset "$(rpad("Testing the Manifold $M",60," "))" begin
        #
        # Manifold functions in alphabetical order – check ManifoldsBase for the mainlist
        #
        # exp – the exponential map
        if has_feature_expectations(features, expectations, :exp)
            in_place = get(features.properties, :inplace, true)
            atol = get(expectations.tolerances, :exp_atol, 0)
            test_exp(
                M,
                Tuple(points),
                Tuple(tangent_vectors);
                in_place=in_place,
                in_place_self=get(features.properties, :inplaceself, in_place),
                atol=atol,
                rtol=get(expectations.tolerances, :exp_rtol, atol > 0 ? 0 : eps),
            )
        end
        # is_point - verify points on a manifold
        if has_feature_expectations(features, expectations, :is_point)
            atol = get(expectations.tolerances, :is_point_atol, 0)
            e = get(expectations.errors, :is_point, DomainError)
            test_is_point(
                M,
                Tuple(points),
                Tuple(non_points);
                errors=e isa AbstractVector ? e : fill(e, length(non_points)),
                atol=atol,
            )
        end
        # is_vector – verify tangent vectors
        if has_feature_expectations(features, expectations, :is_vector)
            atol = get(expectations.tolerances, :is_vector_atol, 0)
            e = get(expectations.errors, :is_vector, DomainError)
            test_is_vector(
                M,
                Tuple(points),
                Tuple(tangent_vectors),
                Tuple(non_tangent_vectors),
                Tuple(non_points);
                atol=atol,
                errors=e isa AbstractVector ? e : fill(e, length(non_points)),
            )
        end
        # manifold_dimension – the dimension on the manifold
        if has_feature_expectations(features, expectations, :manifold_dimension)
            d = expectations.values[:manifold_dimension]
            Test.@testset "Manifold dimension" begin
                Test.@test manifold_dimension(M) ≥ 0
                Test.@test manifold_dimension(M) == d
            end
        end
        # representation_size – the dimension on the manifold
        if has_feature_expectations(features, expectations, :representation_size)
            r = expectations.values[:representation_size]
            Test.@testset "representation size" begin
                Test.@test representation_size(M) == r
            end
        end
    end
    #
end

#
# The old method for reference
#
"""
    test_manifold(
        M::AbstractManifold,
        pts::AbstractVector;
        tests = DEFAULT_TESTS;
        tolerances=DEFAULT_TOLERANCES,
        args,
    )

Test general properties of manifold `M`, given at least three different points
that lie on it (contained in `pts`), where `tests` is a Dictionary specifying which tests to run.

The following Keyword arguments can be used to set properties of the tests, that are actually run.


# Arguments
* `basis_has_specialized_diagonalizing_get = false`: if true, assumes that
    `DiagonalizingOrthonormalBasis` given in `basis_types` has
    [`get_coordinates`](@ref) and [`get_vector`](@ref) that work without caching.
* `basis_types_to_from = ()`: basis types that will be tested based on
    [`get_coordinates`](@ref) and [`get_vector`](@ref).
* `basis_types_vecs = ()` : basis types that will be tested based on `get_vectors`
* `inverse_retraction_method – (`inverse_retraction_method(M)`) default inverse retrcation,
  can be deactivated using `nothing`
* `default_retraction_method – (`default_retraction_method(M)`) default retraction to use,
  can be deactivated using `nothing`
* `exp_log_atol_multiplier = 0`: change absolute tolerance of exp/log tests
    (0 use default, i.e. deactivate atol and use rtol).
* `exp_log_rtol_multiplier = 1`: change the relative tolerance of exp/log tests
    (1 use default). This is deactivated if the `exp_log_atol_multiplier` is nonzero.
- `expected_dimension_type = Integer`: expected type of value returned by
    `manifold_dimension`.
- `inverse_retraction_methods = []`: inverse retraction methods that will be tested.
- `is_mutating = true`: whether mutating variants of functions should be tested.
- `is_tangent_atol_multiplier = 0`: determines atol of `is_vector` checks.
- `mid_point12 = shortest_geodesic(M, pts[1], pts[2], 0.5) if `exp` and `log` are available or a midpoint is explicitly set
- `point_distributions = []` : point distributions to test.
- `rand_tvector_atol_multiplier = 0` : chage absolute tolerance in testing random vectors
    (0 use default, i.e. deactivate atol and use rtol) random tangent vectors are tangent
    vectors.
- `retraction_atol_multiplier = 0`: change absolute tolerance of (inverse) retraction tests
    (0 use default, i.e. deactivate atol and use rtol).
- `retraction_rtol_multiplier = 1`: change the relative tolerance of (inverse) retraction
    tests (1 use default). This is deactivated if the `exp_log_atol_multiplier` is nonzero.
- `retraction_methods = []`: retraction methods that will be tested.
- `test_atlases = []`: Vector or tuple of atlases that should be tested.
- `test_injectivity_radius = true`: whether implementation of [`injectivity_radius`](@ref)
    should be tested.
- `test_inplace = false` : if true check if inplace variants work if they are activated,
   e.g. check that `exp!(M, p, p, X)` work if `exp`tests are activated
   This in general requires `is_mutating` to be true.
- `test_is_tangent`: if true check that the `inverse_retraction_method`
    actually returns valid tangent vectors.
- `test_musical_isomorphisms = false` : test musical isomorphisms.
- `test_mutating_rand = false` : test the mutating random function for points on manifolds.
- `test_project_point = false`: test projections onto the manifold.
- `test_project_tangent = false` : test projections on tangent spaces.
- `test_representation_size = true` : test repersentation size of points/tvectprs.
- `test_tangent_vector_broadcasting = true` : test boradcasting operators on TangentSpace.
- `test_vector_spaces = true` : test Vector bundle of this manifold.
- `test_default_vector_transport = false` : test the default vector transport (usually
   parallel transport).
- `test_vee_hat = false`: test [`vee`](@ref) and [`hat`](@ref) functions.
- `tvector_distributions = []` : tangent vector distributions to test.
- `vector_transport_methods = []`: vector transport methods that should be tested.
- `vector_transport_inverse_retractions = [inverse_retraction_method for _ in 1:length(vector_transport_methods)]``
  inverse retractions to use with the vector transport method (especially the differentiated ones)
- `vector_transport_to = [ true for _ in 1:length(vector_transport_methods)]`: whether
   to check the `to` variant of vector transport
- `vector_transport_direction = [ true for _ in 1:length(vector_transport_methods)]`: whether
   to check the `direction` variant of vector transport
"""
function test_manifold(
    M::AbstractManifold,
    pts::AbstractVector;
    test_features::Dict{Symbol,Bool}=Dict{Symbol,Bool}(),
    tol_features::Dict{Symbol,<:Real}=Dict{Symbol,Float64}(),
    test_functions::Dict{Function,Bool}=Dict{Function,Bool}(),
    tol_functions::Dict{Function,<:Real}=Dict{Function,Float64}(),
    #specific setups
    expected_dimension_type=Integer,
    inverse_retraction_method=default_inverse_retraction_method(M),
    retraction_method=default_retraction_method(M),
    #
    # Old kwargs...
    #
    basis_has_specialized_diagonalizing_get=false,
    basis_types_to_from=(),
    basis_types_vecs=(),
    # Tolerances
    exp_log_atol_multiplier=0,
    exp_log_rtol_multiplier=1,
    inverse_retraction_methods=[],
    is_mutating=true,
    is_tangent_atol_multiplier=0,
    musical_isomorphism_bases=[],
    point_distributions=[],
    projection_atol_multiplier=0,
    rand_tvector_atol_multiplier=0,
    retraction_atol_multiplier=0,
    retraction_methods=[],
    retraction_rtol_multiplier=1,
    test_atlases=(),
    test_is_tangent=true,
    test_injectivity_radius=true,
    test_inplace=false,
    test_musical_isomorphisms=false,
    test_mutating_rand=false,
    parallel_transport=false,
    parallel_transport_to=parallel_transport,
    parallel_transport_along=parallel_transport,
    parallel_transport_direction=parallel_transport,
    test_inner=true,
    test_norm=true,
    test_project_point=false,
    test_project_tangent=false,
    test_rand_point=false,
    test_rand_tvector=false,
    test_riesz_representer=false,
    test_tangent_vector_broadcasting=true,
    test_default_vector_transport=false,
    test_vector_spaces=true,
    test_vee_hat=false,
    tvector_distributions=[],
    vector_transport_methods=[],
    vector_transport_inverse_retractions=[
        inverse_retraction_method for _ in 1:length(vector_transport_methods)
    ],
    vector_transport_retractions=[
        retraction_method for _ in 1:length(vector_transport_methods)
    ],
    test_vector_transport_to=[true for _ in 1:length(vector_transport_methods)],
    test_vector_transport_direction=[true for _ in 1:length(vector_transport_methods)],
    mid_point12=(do_test(test_functions, exp) && do_test(test_functions, log)) ?
                shortest_geodesic(M, pts[1], pts[2], 0.5) : nothing,
)
    #
    # Test calls for the functions
    #
    Test.@testset "$(repr(M)) Tests" begin
        do_test(test_functions, manifold_dimension) &&
            Test.@testset "Manifold dimension                                " begin
                Test.@test isa(manifold_dimension(M), expected_dimension_type)
                Test.@test manifold_dimension(M) ≥ 0
                Test.@test manifold_dimension(M) == vector_space_dimension(
                    Manifolds.VectorBundleFibers(Manifolds.TangentSpace, M),
                )
                Test.@test manifold_dimension(M) == vector_space_dimension(
                    Manifolds.VectorBundleFibers(Manifolds.CotangentSpace, M),
                )
            end
        do_test(test_functions, representation_size) &&
            Test.@testset "representation_size                               " begin
                rs = Manifolds.representation_size(M)
                Test.@test isa(rs, Tuple)
                for s in rs
                    Test.@test s > 0
                end
                for fiber in (Manifolds.TangentSpace, Manifolds.CotangentSpace)
                    rs = Manifolds.representation_size(
                        Manifolds.VectorBundleFibers(fiber, M),
                    )
                    Test.@test isa(rs, Tuple)
                    for s in rs
                        Test.@test s > 0
                    end
                end
            end
        do_test(test_functions, is_point) &&
            Test.@testset "is_point                                         " begin
                for pt in pts
                    atol = get_tolerance(tol_functions, is_point; eps=find_eps(pt))
                    Test.@test is_point(M, pt; atol=atol)
                    Test.@test check_point(M, pt; atol=atol) === nothing
                end
            end
    end
    test_injectivity_radius && Test.@testset "injectivity radius" begin
        Test.@test injectivity_radius(M, pts[1]) > 0
        Test.@test injectivity_radius(M, pts[1]) ≥ injectivity_radius(M)
        for rm in retraction_methods
            Test.@test injectivity_radius(M, rm) > 0
            Test.@test injectivity_radius(M, pts[1], rm) ≥ injectivity_radius(M, rm)
            Test.@test injectivity_radius(M, pts[1], rm) ≤ injectivity_radius(M, pts[1])
        end
    end

    test_is_tangent && Test.@testset "is_vector" begin
        for (p, X) in zip(pts, tv)
            atol = is_tangent_atol_multiplier * find_eps(p)
            if !(check_vector(M, p, X; atol=atol) === nothing)
                print(check_vector(M, p, X; atol=atol))
            end
            Test.@test is_vector(M, p, X; atol=atol)
            Test.@test check_vector(M, p, X; atol=atol) === nothing
        end
    end

    do_test(test_functions, exp) &&
        !do_test(test_functions, :ExpLog) &&
        test_exp(
            M,
            Tuple(pts),
            Tuple(tv);
            in_place=is_mutating,
            in_place_self=test_inplace,
            atol_multiplier=exp_log_atol_multiplier,
            rtol_multiplier=exp_log_rtol_multiplier,
        )
    do_test(test_functions, log) &&
        !do_test(test_features, :ExpLog) &&
        test_log(
            M,
            Tuple(pts),
            Tuple(tv);
            in_place=is_mutating,
            atol_multiplier=exp_log_atol_multiplier,
            rtol_multiplier=exp_log_rtol_multiplier,
        )
    do_test(test_features, :ExpLog) && test_explog(
        M,
        Tuple(pts),
        Tuple(tv);
        in_place=is_mutating,
        in_place_self=test_inplace,
        atol_multiplier=exp_log_atol_multiplier,
        rtol_multiplier=exp_log_rtol_multiplier,
    )
    do_test(test_features, :ExpLog) &&
        Test.@testset "Interplay distance, inner and norm" begin
            X1 = log(M, pts[1], pts[2])
            if test_norm
                Test.@test distance(M, pts[1], pts[2]) ≈ norm(M, pts[1], X1)
            end

            if test_inner
                X3 = log(M, pts[1], pts[3])
                Test.@test inner(M, pts[1], X1, X3) ≈ conj(inner(M, pts[1], X3, X1))
                Test.@test inner(M, pts[1], X1, X1) ≈ real(inner(M, pts[1], X1, X1))

                Test.@test norm(M, pts[1], X1) ≈ sqrt(inner(M, pts[1], X1, X1))
            end
            if test_norm
                Test.@test norm(M, pts[1], X1) isa Real
            end
        end

    parallel_transport && test_parallel_transport(
        M,
        pts;
        along=parallel_transport_along,
        to=parallel_transport_to,
        direction=parallel_transport_direction,
        mutating=is_mutating,
    )
    Test.@testset "retraction tests" begin
        for retr_method in retraction_methods
            test_retr(
                M,
                Tuple(pts),
                Tuple(tv),
                retr_method;
                in_place=is_mutating,
                atol_multiplier=retraction_atol_multiplier,
                rtol_multiplier=retraction_rtol_multiplier,
            )
        end
    end
    Test.@testset "inverse retraction" begin
        for inv_retr_method in inverse_retraction_methods
            test_inv_retr(
                M,
                Tuple(pts),
                Tuple(tv),
                inv_retr_method;
                in_place=is_mutating,
                atol_multiplier=retraction_atol_multiplier,
                rtol_multiplier=retraction_rtol_multiplier,
            )
        end
    end

    Test.@testset "atlases" begin
        if !isempty(test_atlases)
            Test.@test get_default_atlas(M) isa AbstractAtlas{ℝ}
        end
        for A in test_atlases
            i = get_chart_index(M, A, pts[1])
            a = get_parameters(M, A, i, pts[1])
            Test.@test isa(a, AbstractVector)
            Test.@test length(a) == manifold_dimension(M)
            Test.@test isapprox(M, pts[1], get_point(M, A, i, a))
            if is_mutating
                get_parameters!(M, a, A, i, pts[2])
                Test.@test a ≈ get_parameters(M, A, i, pts[2])

                q = allocate(pts[1])
                get_point!(M, q, A, i, a)
                Test.@test isapprox(M, pts[2], q)
            end
        end
    end

    test_riesz_representer && Test.@testset "RieszRepresenterCotangentVector" begin
        rrcv = flat(M, pts[1], tv[1])
        Test.@test rrcv isa RieszRepresenterCotangentVector
        Test.@test rrcv.p === pts[1]
        Test.@test rrcv.X === tv[1]
        basis = dual_basis(M, pts[1], basis_types_to_from[1])
        coords = get_coordinates(M, pts[1], rrcv, basis)
        rrcv2 = get_vector(M, pts[1], coords, basis)
        Test.@test isapprox(M, pts[1], rrcv.X, rrcv2.X)
    end

    test_vector_spaces && Test.@testset "vector spaces" begin
        for p in pts
            X = zero_vector(M, p)
            mts = Manifolds.VectorBundleFibers(Manifolds.TangentSpace, M)
            Test.@test isapprox(M, p, X, zero_vector(mts, p))
            if is_mutating
                zero_vector!(mts, X, p)
                Test.@test isapprox(M, p, X, zero_vector(M, p))
            end
        end
    end

    Test.@testset "basic linear algebra in tangent space" begin
        for (p, X) in zip(pts, tv)
            Test.@test isapprox(M, p, 0 * X, zero_vector(M, p); atol=find_eps(pts[1]))
            Test.@test isapprox(M, p, 2 * X, X + X)
            Test.@test isapprox(M, p, 0 * X, X - X)
            Test.@test isapprox(M, p, (-1) * X, -X)
        end
    end

    test_tangent_vector_broadcasting &&
        Test.@testset "broadcasted linear algebra in tangent space" begin
            for (p, X) in zip(pts, tv)
                Test.@test isapprox(M, p, 3 * X, 2 .* X .+ X)
                Test.@test isapprox(M, p, -X, X .- 2 .* X)
                Test.@test isapprox(M, p, -X, .-X)
                if (isa(X, AbstractArray))
                    Y = allocate(X)
                    Y .= 2 .* X .+ X
                else
                    Y = 2 * X + X
                end
                Test.@test isapprox(M, p, Y, 3 * X)
            end
        end

    test_project_tangent && Test.@testset "project tangent" begin
        for (p, X) in zip(pts, tv)
            atol = find_eps(p) * projection_atol_multiplier
            X_emb = embed(M, p, X)
            Test.@test isapprox(M, p, X, project(M, p, X_emb); atol=atol)
            if is_mutating
                X2 = allocate(X)
                project!(M, X2, p, X_emb)
            else
                X2 = project(M, p, X_emb)
            end
            Test.@test isapprox(M, p, X2, X; atol=atol)
        end
    end

    test_project_point && Test.@testset "project point test" begin
        for p in pts
            atol = find_eps(p) * projection_atol_multiplier
            p_emb = embed(M, p)
            Test.@test isapprox(M, p, project(M, p_emb); atol=atol)
            if is_mutating
                p2 = allocate(p)
                project!(M, p2, p_emb)
            else
                p2 = project(M, p_emb)
            end
            Test.@test isapprox(M, p2, p; atol=atol)
        end
    end

    !(retraction_method === nothing || inverse_retraction_method === nothing) &&
        Test.@testset "vector transport" begin
            tvatol = is_tangent_atol_multiplier * find_eps(pts[1])
            X1 = inverse_retract(M, pts[1], pts[2], inverse_retraction_method)
            X2 = inverse_retract(M, pts[1], pts[3], inverse_retraction_method)
            pts32 = retract(M, pts[1], X2, retraction_method)
            test_default_vector_transport && Test.@testset "default vector transport" begin
                v1t1 = vector_transport_to(M, pts[1], X1, pts32)
                v1t2 = vector_transport_direction(M, pts[1], X1, X2)
                Test.@test is_vector(M, pts32, v1t1; atol=tvatol)
                Test.@test is_vector(M, pts32, v1t2; atol=tvatol)
                Test.@test isapprox(M, pts32, v1t1, v1t2)
                Test.@test isapprox(
                    M,
                    pts[1],
                    vector_transport_to(M, pts[1], X1, pts[1]),
                    X1,
                )

                is_mutating && Test.@testset "mutating variants" begin
                    v1t1_m = allocate(v1t1)
                    v1t2_m = allocate(v1t2)
                    vector_transport_to!(M, v1t1_m, pts[1], X1, pts32)
                    vector_transport_direction!(M, v1t2_m, pts[1], X1, X2)
                    Test.@test isapprox(M, pts32, v1t1, v1t1_m)
                    Test.@test isapprox(M, pts32, v1t2, v1t2_m)
                end
            end

            for (vtm, test_to, test_dir, rtr_m, irtr_m) in zip(
                vector_transport_methods,
                test_vector_transport_to,
                test_vector_transport_direction,
                vector_transport_retractions,
                vector_transport_inverse_retractions,
            )
                Test.@testset "vector transport method $(vtm)" begin
                    tvatol = is_tangent_atol_multiplier * find_eps(pts[1])
                    X1 = inverse_retract(M, pts[1], pts[2], irtr_m)
                    X2 = inverse_retract(M, pts[1], pts[3], irtr_m)
                    pts32 = retract(M, pts[1], X2, rtr_m)
                    test_to && (v1t1 = vector_transport_to(M, pts[1], X1, pts32, vtm))
                    test_dir && (v1t2 = vector_transport_direction(M, pts[1], X1, X2, vtm))
                    test_to && Test.@test is_vector(M, pts32, v1t1, true; atol=tvatol)
                    test_dir && Test.@test is_vector(M, pts32, v1t2, true; atol=tvatol)
                    (test_to && test_dir) &&
                        Test.@test isapprox(M, pts32, v1t1, v1t2, atol=tvatol)
                    test_to && Test.@test isapprox(
                        M,
                        pts[1],
                        vector_transport_to(M, pts[1], X1, pts[1], vtm),
                        X1;
                        atol=tvatol,
                    )
                    test_dir && Test.@test isapprox(
                        M,
                        pts[1],
                        vector_transport_direction(
                            M,
                            pts[1],
                            X1,
                            zero_vector(M, pts[1]),
                            vtm,
                        ),
                        X1;
                        atol=tvatol,
                    )

                    is_mutating && Test.@testset "mutating variants" begin
                        if test_to
                            v1t1_m = allocate(v1t1)
                            vector_transport_to!(M, v1t1_m, pts[1], X1, pts32, vtm)
                            Test.@test isapprox(M, pts32, v1t1, v1t1_m; atol=tvatol)
                            test_inplace &&
                                Test.@testset "inplace test for vector_transport_to!" begin
                                    X1a = copy(M, pts[1], X1)
                                    Xt = vector_transport_to(M, pts[1], X1, pts32, vtm)
                                    vector_transport_to!(M, X1a, pts[1], X1a, pts32, vtm)
                                    Test.@test isapprox(M, pts[1], X1a, Xt; atol=tvatol)
                                end
                        end
                        if test_dir
                            v1t2_m = allocate(v1t2)
                            vector_transport_direction!(M, v1t2_m, pts[1], X1, X2, vtm)
                            Test.@test isapprox(M, pts32, v1t2, v1t2_m; atol=tvatol)
                            test_inplace &&
                                Test.@testset "inplace test for vector_transport_direction!" begin
                                    X1a = copy(M, pts[1], X1)
                                    X2a = copy(M, pts[1], X2)
                                    Xt = vector_transport_direction(M, pts[1], X1, X2, vtm)
                                    vector_transport_direction!(
                                        M,
                                        X1a,
                                        pts[1],
                                        X1a,
                                        X2,
                                        vtm,
                                    )
                                    vector_transport_direction!(
                                        M,
                                        X2a,
                                        pts[1],
                                        X1,
                                        X2a,
                                        vtm,
                                    )
                                    Test.@test isapprox(M, pts[1], X1a, Xt; atol=tvatol)
                                    Test.@test isapprox(M, pts[1], X2a, Xt; atol=tvatol)
                                end
                        end
                    end
                end
            end
        end

    for btype in basis_types_vecs
        Test.@testset "Basis support for $(btype)" begin
            p = pts[1]
            b = get_basis(M, p, btype)
            Test.@test isa(b, CachedBasis)
            bvectors = get_vectors(M, p, b)
            N = length(bvectors)

            # test orthonormality
            for i in 1:N
                Test.@test norm(M, p, bvectors[i]) ≈ 1
                for j in (i + 1):N
                    Test.@test real(inner(M, p, bvectors[i], bvectors[j])) ≈ 0 atol =
                        sqrt(find_eps(p))
                end
            end
            if isa(btype, ProjectedOrthonormalBasis)
                # check projection idempotency
                for i in 1:N
                    Test.@test norm(M, p, bvectors[i]) ≈ 1
                    for j in (i + 1):N
                        Test.@test real(inner(M, p, bvectors[i], bvectors[j])) ≈ 0 atol =
                            sqrt(find_eps(p))
                    end
                end
                # check projection idempotency
                for i in 1:N
                    Test.@test isapprox(M, p, project(M, p, bvectors[i]), bvectors[i])
                end
            end
            if !isa(btype, ProjectedOrthonormalBasis) && (
                basis_has_specialized_diagonalizing_get ||
                !isa(btype, DiagonalizingOrthonormalBasis)
            )
                X1 = inverse_retract(M, p, pts[2], inverse_retraction_method)
                Xb = get_coordinates(M, p, X1, btype)

                Test.@test get_coordinates(M, p, X1, b) ≈ Xb
                Test.@test isapprox(
                    M,
                    p,
                    get_vector(M, p, Xb, b),
                    get_vector(M, p, Xb, btype),
                )
            end
        end
    end

    for btype in (basis_types_to_from..., basis_types_vecs...)
        p = pts[1]
        N = number_of_coordinates(M, btype)
        if !isa(btype, ProjectedOrthonormalBasis) && (
            basis_has_specialized_diagonalizing_get ||
            !isa(btype, DiagonalizingOrthonormalBasis)
        )
            X1 = inverse_retract(M, p, pts[2], inverse_retraction_method)

            Xb = get_coordinates(M, p, X1, btype)
            #Test.@test isa(Xb, AbstractVector{<:Real})
            Test.@test length(Xb) == N
            Xbi = get_vector(M, p, Xb, btype)
            Test.@test isapprox(M, p, X1, Xbi)

            Xs = [[ifelse(i == j, 1, 0) for j in 1:N] for i in 1:N]
            Xs_invs = [get_vector(M, p, Xu, btype) for Xu in Xs]
            # check orthonormality of inverse representation
            for i in 1:N
                Test.@test norm(M, p, Xs_invs[i]) ≈ 1 atol = sqrt(find_eps(p))
                for j in (i + 1):N
                    Test.@test real(inner(M, p, Xs_invs[i], Xs_invs[j])) ≈ 0 atol =
                        sqrt(find_eps(p))
                end
            end

            if is_mutating
                Xb_s = allocate(Xb)
                Test.@test get_coordinates!(M, Xb_s, p, X1, btype) === Xb_s
                Test.@test isapprox(Xb_s, Xb; atol=find_eps(p))

                Xbi_s = allocate(Xbi)
                Test.@test get_vector!(M, Xbi_s, p, Xb, btype) === Xbi_s
                Test.@test isapprox(M, p, X1, Xbi_s)
            end
        end
    end

    test_vee_hat && Test.@testset "vee and hat" begin
        p = pts[1]
        q = pts[2]
        X = inverse_retract(M, p, q, inverse_retraction_method)
        Y = vee(M, p, X)
        Test.@test length(Y) == number_of_coordinates(M, ManifoldsBase.VeeOrthogonalBasis())
        Test.@test isapprox(M, p, X, hat(M, p, Y))
        Y2 = allocate(Y)
        vee_ret = vee!(M, Y2, p, X)
        Test.@test vee_ret === Y2
        Test.@test isapprox(Y, Y2)
        X2 = allocate(X)
        hat_ret = hat!(M, X2, p, Y)
        Test.@test hat_ret === X2
        Test.@test isapprox(M, p, X2, X)
    end

    mid_point12 !== nothing && Test.@testset "midpoint" begin
        epsp1p2 = find_eps(pts[1], pts[2])
        atolp1p2 = exp_log_atol_multiplier * epsp1p2
        rtolp1p2 =
            exp_log_atol_multiplier == 0.0 ? sqrt(epsp1p2) * exp_log_rtol_multiplier : 0
        mp = mid_point(M, pts[1], pts[2])
        Test.@test isapprox(M, mp, mid_point12; atol=atolp1p2, rtol=rtolp1p2)
        if is_mutating
            mpm = allocate(mp)
            mid_point!(M, mpm, pts[1], pts[2])
            Test.@test isapprox(M, mpm, mid_point12; atol=atolp1p2, rtol=rtolp1p2)
            test_inplace && Test.@testset "inplace test for midpoint!" begin
                p1 = copy(M, pts[1])
                p2 = copy(M, pts[2])
                p3 = mid_point(M, pts[1], pts[2])
                mid_point!(M, p1, p1, pts[2])
                mid_point!(M, p2, pts[1], p2)
                Test.@test isapprox(M, p3, p1)
                Test.@test isapprox(M, p3, p2)
            end
        end
    end

    test_musical_isomorphisms && Test.@testset "Musical isomorphisms" begin
        if inverse_retraction_method !== nothing
            tv_m = inverse_retract(M, pts[1], pts[2], inverse_retraction_method)
        else
            tv_m = zero_vector(M, pts[1])
        end
        ctv_m = flat(M, pts[1], tv_m)
        Test.@test ctv_m(tv_m) ≈ norm(M, pts[1], tv_m)^2
        tv_m_back = sharp(M, pts[1], ctv_m)
        Test.@test isapprox(M, pts[1], tv_m, tv_m_back)

        if is_mutating
            ctv_m_s = allocate(ctv_m)
            flat!(M, ctv_m_s, pts[1], tv_m)
            Test.@test ctv_m_s(tv_m) ≈ ctv_m(tv_m)
            tv_m_s_back = allocate(tv_m_back)
            sharp!(M, tv_m_s_back, pts[1], ctv_m_s)
            Test.@test isapprox(M, pts[1], tv_m, tv_m_s_back)
        end

        for basis in musical_isomorphism_bases
            tv_m_f = ManifoldsBase.TFVector(get_coordinates(M, pts[1], tv_m, basis), basis)
            ctv_m_f = flat(M, pts[1], tv_m_f)
            Test.@test isa(ctv_m_f, CoTFVector)
            tv_m_f_back = sharp(M, pts[1], ctv_m_f)
            Test.@test isapprox(tv_m_f.data, tv_m_f_back.data)
        end
    end

    Test.@testset "number_eltype" begin
        for (p, X) in zip(pts, tv)
            Test.@test number_eltype(X) == number_eltype(p)
            p = retract(M, p, X, retraction_method)
            Test.@test number_eltype(p) == number_eltype(p)
        end
    end

    is_mutating && Test.@testset "copyto!" begin
        for (p, X) in zip(pts, tv)
            p2 = allocate(p)
            copyto!(p2, p)
            Test.@test isapprox(M, p2, p)

            X2 = allocate(X)
            if inverse_retraction_method === nothing
                X3 = zero_vector(M, p)
                copyto!(X2, X3)
                Test.@test isapprox(M, p, X2, zero_vector(M, p))
            else
                q = retract(M, p, X, retraction_method)
                X3 = inverse_retract(M, p, q, inverse_retraction_method)
                copyto!(X2, X3)
                Test.@test isapprox(
                    M,
                    p,
                    X2,
                    inverse_retract(M, p, q, inverse_retraction_method),
                )
            end
        end
    end

    is_mutating && Test.@testset "point distributions" begin
        for p in pts
            prand = allocate(p)
            for pd in point_distributions
                Test.@test Manifolds.support(pd) isa Manifolds.MPointSupport{typeof(M)}
                for _ in 1:10
                    Test.@test is_point(M, rand(pd))
                    if test_mutating_rand
                        rand!(pd, prand)
                        Test.@test is_point(M, prand)
                    end
                end
            end
        end
    end

    test_rand_point && Test.@testset "Base.rand point generation" begin
        rng_a = MersenneTwister(123)
        rng_b = MersenneTwister(123)
        Test.@test is_point(M, rand(M), true)
        # ensure that the RNG source is actually used
        Test.@test rand(rng_a, M) == rand(rng_b, M)
        # generation of multiple points
        Test.@test all(p -> is_point(M, p, true), rand(M, 3))
        Test.@test all(p -> is_point(M, p, true), rand(rng_a, M, 3))

        if test_inplace && is_mutating
            rng_a = MersenneTwister(123)
            rng_b = MersenneTwister(123)

            p = allocate(pts[1])
            rand!(M, p)
            Test.@test is_point(M, p, true)
            p = allocate(pts[1])
            rand!(rng_a, M, p)
            Test.@test is_point(M, p, true)
            # ensure that the RNG source is actually used
            q = allocate(pts[1])
            rand!(rng_b, M, q)
            Test.@test p == q
        end
    end

    test_rand_tvector && Test.@testset "Base.rand tangent vector generation" begin
        p = pts[1]
        rng_a = MersenneTwister(123)
        rng_b = MersenneTwister(123)
        randX = rand(M; vector_at=p)
        atol = rand_tvector_atol_multiplier * find_eps(randX)
        Test.@test is_vector(M, p, randX, true; atol=atol)
        # ensure that the RNG source is actually used
        Test.@test rand(rng_a, M; vector_at=p) == rand(rng_b, M; vector_at=p)
        # generation of multiple tangent vectors
        Test.@test all(X -> is_vector(M, p, X, true; atol=atol), rand(M, 3; vector_at=p))
        Test.@test all(
            X -> is_vector(M, p, X, true; atol=atol),
            rand(rng_a, M, 3; vector_at=p),
        )

        if test_inplace && is_mutating
            rng_a = MersenneTwister(123)
            rng_b = MersenneTwister(123)

            X = allocate(tv[1])
            rand!(M, X; vector_at=p)
            Test.@test is_vector(M, p, X, true; atol=atol)
            X = allocate(tv[1])
            rand!(rng_a, M, X; vector_at=p)
            Test.@test is_point(M, p, true)
            # ensure that the RNG source is actually used
            Y = allocate(tv[1])
            rand!(rng_b, M, Y; vector_at=p)
            Test.@test X == Y
        end
    end

    Test.@testset "tangent vector distributions" begin
        for tvd in tvector_distributions
            supp = Manifolds.support(tvd)
            Test.@test supp isa Manifolds.FVectorSupport{TangentBundleFibers{typeof(M)}}
            for _ in 1:10
                randtv = rand(tvd)
                atol = rand_tvector_atol_multiplier * find_eps(randtv)
                Test.@test is_vector(M, supp.point, randtv, true; atol=atol)
            end
        end
    end
    return nothing
end

"""
    test_parallel_transport(M, pts; along=false, to=true, diretion=true)

Generic tests for parallel transport on `M`given at least two pointsin `P`.

The single functions to transport `along` (a curve), `to` (a point) or (towards a) `direction`
are sub-tests that can be activated by the keywords arguemnts

!!! Note
Since the interface to specify curves is not yet provided, the along keyword does not have an effect yet
"""
function test_parallel_transport(
    M::AbstractManifold,
    P,
    Ξ=inverse_retract.(
        Ref(M),
        P[1:(end - 1)],
        P[2:end],
        Ref(default_inverse_retraction_method(M)),
    );
    along=false,
    to=true,
    direction=true,
    mutating=true,
)
    length(P) < 2 &&
        error("The Parallel Transport test set requires at least 2 points in P")
    Test.@testset "Parallel Transport" begin
        along && @warn "A test for parallel transport along test not yet implemented"
        Test.@testset "To (a point)" begin # even with to =false this displays no
            if to
                for i in 1:(length(P) - 1)
                    p = P[i]
                    q = P[i + 1]
                    X = Ξ[i]
                    Y1 = parallel_transport_to(M, p, X, q)
                    if mutating
                        Y2 = similar(X)
                        parallel_transport_to!(M, Y2, p, X, q)
                        # test that mutating and allocating to the same
                        Test.@test isapprox(M, q, Y1, Y2)
                        parallel_transport_to!(M, Y2, q, Y1, p)
                        # Test that transporting there and back again yields the identity
                        Test.@test isapprox(M, q, X, Y2)
                        parallel_transport_to!(M, Y1, q, Y1, p)
                        # Test that inplace does not have side effects
                    else
                        Y1 = parallel_transport_to(M, q, Y1, p)
                    end
                    Test.@test isapprox(M, q, X, Y1)
                end
            end
        end
        Test.@testset "(Tangent Vector) Direction" begin
            if direction
                for i in 1:(length(P) - 1)
                    p = P[i]
                    X = Ξ[i]
                    Y1 = parallel_transport_direction(M, p, X, X)
                    q = exp(M, p, X)
                    if mutating
                        Y2 = similar(X)
                        parallel_transport_direction!(M, Y2, p, X, X)
                        # test that mutating and allocating to the same
                        Test.@test isapprox(M, q, Y1, Y2)
                    end
                    # Test that Y is a tangent vector at q
                    Test.@test is_vector(M, p, Y1, true)
                end
            end
        end
    end
end