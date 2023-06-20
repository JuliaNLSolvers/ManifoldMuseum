using Manifolds, Test, JLD2
using ManifoldsBase: ManifoldDomainError
using Manifolds: test_manifold, has_feature_expectations
# Generate or load Test Scenario – can be set here or globally by env
generate_test = false;
config_file = (@__DIR__) * "/config/sphere.jld2"
generate = get(ENV, "TEST_MANIFOLD_GENERATE_TESTS", generate_test)

M = Sphere(2)
M2 = ArraySphere(2, 2)
if generate || !isfile(config_file)
    # Generate (semi-automatically) and save setup
    using Manifolds: find_manifold_functions, ManifoldFeatures, ManifoldExpectations
    features = ManifoldFeatures(M)
    expectations = ManifoldExpectations(
        values=Dict(
            :manifold_dimension => 2,
            :repr_manifold => "Sphere(2, ℝ)",
            :representation_size => (3,),
        ),
        tolerances=Dict(:exp_atol => 1e-9),
        errors=Dict(
            :is_point => DomainError,
            :is_vector => [DomainError, ManifoldDomainError],
        ),
    )
    expectations2 = Manifolds.ManifoldExpectations(
        values=Dict(
            :manifold_dimension => 3,
            :repr_manifold => "ArraySphere(2, 2; field = ℝ)",
            :representation_size => (2, 2),
        ),
        tolerances=Dict(:exp_atol => 1e-9),
        errors=Dict(
            :is_point => :DomainError,
            :is_vector => [DomainError, ManifoldDomainError],
        ),
    )
    jldsave(config_file; features, expectations, expectations2)
    @warn "Configuration for the Sphere regenerated. This should not be actve by default."
else
    file = jldopen(config_file)
    features = file["features"]
    expectations = file["expectations"]
    expectations2 = file["expectations2"]
    close(file)
end

#
# (classical / vector) Sphere
ps = [[1.0, 0.0, 0.0], 1 / sqrt(2) .* [1.0, 1.0, 0.0], 1 / sqrt(2) .* [1.0, 0.0, 1.0]]
nps = [[2.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]]
Xs = [1 / sqrt(2) .* [0.0, 1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
nXs = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]

test_manifold(
    M;
    points=ps,
    tangent_vectors=Xs,
    features=features,
    expectations=expectations,
    non_points=nps,
    non_tangent_vectors=nXs,
)

#
# ArraySphere
p2s = [[1.0 0.0; 0.0 0.0], [1/sqrt(2) 0.0; 0.0 1/sqrt(2)], [0.0 1/sqrt(2); 1/sqrt(2) 0.0]]
X2s = [1 / sqrt(2) .* [0.0 0.0; 1.0 1.0], [0.0 0.0; 1.0 0.0], [1.0 0.0; 0.0 0.0]]
np2s = [[1.0, 0.0, 0.0, 0.0], [2.0 0.0; 0.0 0.0]]
nX2s = [[2.0 0.0; 0.0 0.0], [0.0, 1.0, 0.0, 0.0]]

test_manifold(
    M2;
    points=p2s,
    tangent_vectors=X2s,
    features=features,
    expectations=expectations2,
    non_points=np2s,
    non_tangent_vectors=nX2s,
)