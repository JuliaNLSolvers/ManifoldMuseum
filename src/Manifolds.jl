"""
`Manifolds.jl` provides a library of manifolds aiming for an easy-to-use and fast implementation.
"""
module Manifolds

import Base:
    ^,
    angle,
    copyto!,
    convert,
    foreach,
    identity,
    in,
    inv,
    isempty,
    length,
    ndims,
    show,
    showerror,
    size,
    transpose
import ManifoldsBase:
    @trait_function,
    _access_nested,
    _get_basis,
    _injectivity_radius,
    _inverse_retract,
    _inverse_retract!,
    _read,
    _retract,
    _retract!,
    _write,
    active_traits,
    allocate,
    allocate_coordinates,
    allocate_result,
    allocate_result_type,
    allocation_promotion_function,
    array_value,
    base_manifold,
    change_basis,
    change_basis!,
    change_metric,
    change_metric!,
    change_representer,
    change_representer!,
    check_point,
    check_size,
    check_vector,
    copy,
    copyto!,
    default_approximation_method,
    default_inverse_retraction_method,
    default_retraction_method,
    default_vector_transport_method,
    decorated_manifold,
    distance,
    dual_basis,
    embed,
    embed!,
    exp,
    exp!,
    fiber_dimension,
    get_basis,
    get_basis_default,
    get_basis_diagonalizing,
    get_basis_orthogonal,
    get_basis_orthonormal,
    get_basis_vee,
    get_component,
    get_coordinates,
    get_coordinates!,
    get_coordinates_diagonalizing,
    get_coordinates_diagonalizing!,
    get_coordinates_orthogonal,
    get_coordinates_orthonormal,
    get_coordinates_orthogonal!,
    get_coordinates_orthonormal!,
    get_coordinates_vee!,
    get_embedding,
    get_iterator,
    get_vector,
    get_vector!,
    get_vector_diagonalizing,
    get_vector_diagonalizing!,
    get_vector_orthogonal,
    get_vector_orthonormal,
    get_vector_orthogonal!,
    get_vector_orthonormal!,
    get_vectors,
    gram_schmidt,
    hat,
    hat!,
    injectivity_radius,
    _injectivity_radius,
    injectivity_radius_exp,
    inner,
    isapprox,
    _isapprox,
    is_flat,
    is_point,
    is_vector,
    inverse_retract,
    inverse_retract!,
    _inverse_retract!,
    inverse_retract_cayley!,
    inverse_retract_embedded!,
    inverse_retract_pade!,
    inverse_retract_polar!,
    inverse_retract_project!,
    inverse_retract_qr!,
    inverse_retract_shooting!,
    inverse_retract_softmax!,
    log,
    log!,
    manifold_dimension,
    mid_point,
    mid_point!,
    norm,
    number_eltype,
    number_of_coordinates,
    parallel_transport_along,
    parallel_transport_along!,
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
    parent_trait,
    power_dimensions,
    project,
    project!,
    rand,
    rand!,
    representation_size,
    retract,
    retract!,
    _retract,
    retract!,
    retract_cayley!,
    retract_exp_ode!,
    retract_pade!,
    retract_polar!,
    retract_project!,
    retract_qr!,
    retract_sasaki!,
    retract_softmax!,
    riemann_tensor,
    riemann_tensor!,
    sectional_curvature,
    sectional_curvature_max,
    sectional_curvature_min,
    set_component!,
    submanifold,
    submanifold_component,
    submanifold_components,
    vector_space_dimension,
    vector_transport_along,           # just specified in Euclidean - the next 5 as well
    vector_transport_along!,
    vector_transport_along_diff!,     # For consistency these are imported, but for now not
    vector_transport_along_project!,  # overwritten with new definitons.
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_direction_diff!,
    vector_transport_to,
    vector_transport_to!,
    vector_transport_to_diff!,
    vector_transport_to_project!, # some overwrite layer 2
    _vector_transport_direction,
    _vector_transport_direction!,
    _vector_transport_to,
    _vector_transport_to!,
    vee,
    vee!,
    Weingarten,
    Weingarten!,
    zero_vector,
    zero_vector!
import ManifoldDiff:
    adjoint_Jacobi_field,
    adjoint_Jacobi_field!,
    diagonalizing_projectors,
    jacobi_field,
    jacobi_field!,
    riemannian_gradient,
    riemannian_gradient!,
    riemannian_Hessian,
    riemannian_Hessian!

import Statistics: mean, mean!, median, median!, cov, var
import StatsBase: mean_and_var

using Base.Iterators: repeated
using Distributions
using Einsum: @einsum
using HybridArrays
using InteractiveUtils: subtypes
using Kronecker
using Graphs
using LinearAlgebra
using ManifoldsBase:
    @next_trait_function,
    ℝ,
    ℂ,
    ℍ,
    AbstractApproximationMethod,
    AbstractBasis,
    AbstractDecoratorManifold,
    AbstractInverseRetractionMethod,
    AbstractManifold,
    AbstractManifoldPoint,
    AbstractMetric,
    AbstractNumbers,
    AbstractOrthogonalBasis,
    AbstractOrthonormalBasis,
    AbstractPowerManifold,
    AbstractPowerRepresentation,
    AbstractRetractionMethod,
    AbstractTrait,
    AbstractVectorTransportMethod,
    AbstractLinearVectorTransportMethod,
    ApproximateInverseRetraction,
    ApproximateRetraction,
    CachedBasis,
    CayleyRetraction,
    CayleyInverseRetraction,
    ComplexNumbers,
    ComponentManifoldError,
    CompositeManifoldError,
    CotangentSpace,
    CotangentSpaceType,
    CoTFVector,
    CoTVector,
    CyclicProximalPointEstimation,
    DefaultBasis,
    DefaultOrthogonalBasis,
    DefaultOrthonormalBasis,
    DefaultOrDiagonalizingBasis,
    DiagonalizingBasisData,
    DiagonalizingOrthonormalBasis,
    DifferentiatedRetractionVectorTransport,
    EfficientEstimator,
    EmbeddedManifold,
    EmptyTrait,
    EuclideanMetric,
    ExponentialRetraction,
    ExtrinsicEstimation,
    Fiber,
    FiberType,
    FVector,
    GeodesicInterpolation,
    GeodesicInterpolationWithinRadius,
    GradientDescentEstimation,
    InverseProductRetraction,
    IsIsometricEmbeddedManifold,
    IsEmbeddedManifold,
    IsEmbeddedSubmanifold,
    IsExplicitDecorator,
    LogarithmicInverseRetraction,
    ManifoldDomainError,
    ManifoldsBase,
    NestedPowerRepresentation,
    NestedReplacingPowerRepresentation,
    TraitList,
    NLSolveInverseRetraction,
    ODEExponentialRetraction,
    OutOfInjectivityRadiusError,
    PadeRetraction,
    PadeInverseRetraction,
    ParallelTransport,
    PolarInverseRetraction,
    PolarRetraction,
    PoleLadderTransport,
    PowerBasisData,
    PowerManifold,
    PowerManifoldNested,
    PowerManifoldNestedReplacing,
    ProductBasisData,
    ProductManifold,
    ProductMetric,
    ProductRetraction,
    ProductVectorTransport,
    ProjectedOrthonormalBasis,
    ProjectionInverseRetraction,
    ProjectionRetraction,
    ProjectionTransport,
    QuaternionNumbers,
    QRInverseRetraction,
    QRRetraction,
    RealNumbers,
    RiemannianMetric,
    SasakiRetraction,
    ScaledVectorTransport,
    SchildsLadderTransport,
    ShootingInverseRetraction,
    SoftmaxRetraction,
    SoftmaxInverseRetraction,
    TangentSpace,
    TangentSpaceType,
    TCoTSpaceType,
    TFVector,
    TVector,
    TypeParameter,
    ValidationCoTVector,
    ValidationManifold,
    ValidationMPoint,
    ValidationTVector,
    VectorSpaceFiber,
    VectorSpaceType,
    VeeOrthogonalBasis,
    WeiszfeldEstimation,
    @invoke_maker,
    _euclidean_basis_vector,
    combine_allocation_promotion_functions,
    geodesic,
    geodesic!,
    get_parameter,
    merge_traits,
    next_trait,
    number_of_components,
    number_system,
    real_dimension,
    rep_size_to_colons,
    shortest_geodesic,
    shortest_geodesic!,
    size_to_tuple,
    trait,
    wrap_type_parameter,
    ziptuples
using ManifoldDiff: ManifoldDiff
using ManifoldDiff:
    default_differential_backend,
    _derivative,
    _derivative!,
    differential,
    differential!,
    gradient,
    gradient!,
    _gradient,
    _gradient!,
    hessian,
    _hessian,
    jacobian,
    _jacobian,
    _jacobian!,
    set_default_differential_backend!
using ManifoldDiff:
    AbstractDiffBackend,
    AbstractRiemannianDiffBackend,
    CoprojectorOntoVector,
    ExplicitEmbeddedBackend,
    IdentityProjector,
    NoneDiffBackend,
    ProjectorOntoVector,
    RiemannianProjectionBackend,
    TangentDiffBackend
using Markdown: @doc_str
using MatrixEquations: lyapc, sylvc
using Quaternions: Quaternions
using Random
using RecursiveArrayTools: ArrayPartition
using Requires
using SimpleWeightedGraphs: AbstractSimpleWeightedGraph, get_weight
using SpecialFunctions
using StaticArrays
using Statistics
using StatsBase
using StatsBase: AbstractWeights

include("utils/helpers.jl")

include("product_representations.jl")

include("manifold_fallbacks.jl")

# Main Meta Manifolds
include("manifolds/ConnectionManifold.jl")
include("manifolds/MetricManifold.jl")
include("manifolds/QuotientManifold.jl")
include("manifolds/Fiber.jl")
include("manifolds/FiberBundle.jl")
include("manifolds/VectorFiber.jl")
include("manifolds/VectorBundle.jl")
include("groups/group.jl")

# Features I: Which are extended on Meta Manifolds
include("distributions.jl")
include("projected_distribution.jl")
include("statistics.jl")

# Meta Manifolds II: Products
include("manifolds/ProductManifold.jl")

METAMANIFOLDS = [
    AbstractManifold,
    AbstractDecoratorManifold,
    AbstractPowerManifold,
    PowerManifoldNested,
    PowerManifoldNestedReplacing,
    ProductManifold,
    TangentSpace,
    ValidationManifold,
    VectorBundle,
]

# Features II: That require metas
include("atlases.jl")
include("differentiation/ode_callback.jl")
include("cotangent_space.jl")

# Meta Manifolds II: Power Manifolds
include("manifolds/PowerManifold.jl")
include("manifolds/GraphManifold.jl")

#
# Manifolds
#
include("manifolds/Euclidean.jl")
include("manifolds/Lorentz.jl")

include("manifolds/CenteredMatrices.jl")
include("manifolds/CholeskySpace.jl")
include("manifolds/Circle.jl")
include("manifolds/Elliptope.jl")
include("manifolds/EmbeddedTorus.jl")
include("manifolds/FixedRankMatrices.jl")
include("manifolds/Flag.jl")
include("manifolds/FlagOrthogonal.jl")
include("manifolds/FlagStiefel.jl")
include("manifolds/GeneralizedGrassmann.jl")
include("manifolds/GeneralizedStiefel.jl")
include("manifolds/Hyperbolic.jl")
include("manifolds/Hyperrectangle.jl")
include("manifolds/MultinomialDoublyStochastic.jl")
include("manifolds/MultinomialSymmetric.jl")
include("manifolds/MultinomialSymmetricPositiveDefinite.jl")
include("manifolds/PositiveNumbers.jl")
include("manifolds/ProjectiveSpace.jl")
include("manifolds/SkewHermitian.jl")
include("manifolds/Spectrahedron.jl")
include("manifolds/Stiefel.jl")
include("manifolds/StiefelEuclideanMetric.jl")
include("manifolds/StiefelCanonicalMetric.jl")
include("manifolds/StiefelSubmersionMetric.jl")
include("manifolds/Sphere.jl")
include("manifolds/SphereSymmetricMatrices.jl")
include("manifolds/Symmetric.jl")
include("manifolds/SymmetricPositiveDefinite.jl")
include("manifolds/SPDFixedDeterminant.jl")
include("manifolds/SymmetricPositiveDefiniteBuresWasserstein.jl")
include("manifolds/SymmetricPositiveDefiniteGeneralizedBuresWasserstein.jl")
include("manifolds/SymmetricPositiveDefiniteAffineInvariant.jl")
include("manifolds/SymmetricPositiveDefiniteLogCholesky.jl")
include("manifolds/SymmetricPositiveDefiniteLogEuclidean.jl")
include("manifolds/SymmetricPositiveSemidefiniteFixedRank.jl")
include("manifolds/Tucker.jl")
#
include("manifolds/ProbabilitySimplex.jl")
include("manifolds/ProbabilitySimplexEuclideanMetric.jl")
include("manifolds/GeneralUnitaryMatrices.jl")
include("manifolds/Unitary.jl")
include("manifolds/Rotations.jl")
include("manifolds/Orthogonal.jl")

# shape spaces require Sphere
include("manifolds/KendallsPreShapeSpace.jl")
include("manifolds/KendallsShapeSpace.jl")

# Introduce the quotient, Grassmann, only after Stiefel
include("manifolds/Grassmann.jl")

# Introduce Symplectic and so on manifolds only after Grassmann
# Since that defines the StiefelPoint, StiefelTVector
include("manifolds/Symplectic.jl")
include("manifolds/Hamiltonian.jl") # Hamiltonian requires symplectic
include("manifolds/SymplecticStiefel.jl")
include("manifolds/SymplecticGrassmann.jl") # Requires SymplecticStiefel

# Product or power based manifolds
include("manifolds/Torus.jl")
include("manifolds/Multinomial.jl")
include("manifolds/Oblique.jl")
include("manifolds/EssentialManifold.jl")

#
# Group Manifolds
include("groups/GroupManifold.jl")

# a) generics
include("groups/addition_operation.jl")
include("groups/multiplication_operation.jl")
include("groups/connections.jl")
include("groups/metric.jl")
include("groups/group_action.jl")
include("groups/group_operation_action.jl")
include("groups/validation_group.jl")
include("groups/product_group.jl")
include("groups/semidirect_product_group.jl")
include("groups/power_group.jl")

# generic group - commopn (special) unitary/orthogonal functions
include("groups/general_unitary_groups.jl")
# Special Group Manifolds
include("groups/general_linear.jl")
include("groups/special_linear.jl")
include("groups/translation_group.jl")
include("groups/unitary.jl")
include("groups/orthogonal.jl")
include("groups/special_unitary.jl")
include("groups/special_orthogonal.jl")
include("groups/circle_group.jl")
include("groups/heisenberg.jl")

include("groups/translation_action.jl")
include("groups/rotation_action.jl")

include("groups/special_euclidean.jl")

include("groups/rotation_translation_action.jl")

# final utilities
include("trait_recursion_breaking.jl")
include("utils/features.jl")

@doc raw"""
    Base.in(p, M::AbstractManifold; kwargs...)
    p ∈ M

Check, whether a point `p` is a valid point (i.e. in) a [`AbstractManifold`](@extref `ManifoldsBase.AbstractManifold`)  `M`.
This method employs [`is_point`](@extref `ManifoldsBase.is_point-Tuple{AbstractManifold, Any, Bool}`) deactivating the error throwing option.
"""
Base.in(p, M::AbstractManifold; kwargs...) = is_point(M, p, false; kwargs...)

@doc raw"""
    Base.in(p, TpM::TangentSpace; kwargs...)
    X ∈ TangentSpace(M, p)

Check whether `X` is a tangent vector from (in) the tangent space ``T_p\mathcal M``, i.e.
the [`TangentSpace`](@extref `ManifoldsBase.TangentSpace`)
at `p` on the [`AbstractManifold`](@extref `ManifoldsBase.AbstractManifold`)  `M`.
This method uses [`is_vector`](@extref `ManifoldsBase.is_vector-Tuple{AbstractManifold, Any, Any, Bool, Bool}`) deactivating the error throw option.
"""
function Base.in(X, TpM::TangentSpace; kwargs...)
    return is_vector(base_manifold(TpM), TpM.point, X, false; kwargs...)
end

@doc raw"""
    manifold_volume(M::AbstractManifold)

Volume of manifold `M` defined through integration of Riemannian volume element in a chart.
Note that for many manifolds there is no universal agreement over the exact ranges over
which the integration should happen. For details see [BoyaSudarshanTilma:2003](@cite).
"""
manifold_volume(::AbstractManifold)

@doc raw"""
    volume_density(M::AbstractManifold, p, X)

Volume density function of manifold `M`, i.e. determinant of the differential of exponential map
`exp(M, p, X)`. Determinant can be understood as computed in a basis, from the matrix
of the linear operator said differential corresponds to. Details are available in Section 4.1
of [ChevallierLiLuDunson:2022](@cite).

Note that volume density is well-defined only for `X` for which `exp(M, p, X)` is injective.
"""
volume_density(::AbstractManifold, p, X)

# functions populated with methods by extensions

function solve_chart_log_bvp end
function estimate_distance_from_bvp end

function solve_chart_exp_ode end
function solve_chart_parallel_transport_ode end

function find_eps end
function test_parallel_transport end
function test_manifold end
function test_group end
function test_action end

# end of functions populated with methods by extensions

function __init__()
    @static if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f === solve_exp_ode
                print(io, "\nDid you forget to load OrdinaryDiffEq? For example: ")
                printstyled(io, "`using OrdinaryDiffEq`", color=:cyan)
            end
        end
    end

    @static if !isdefined(Base, :get_extension)
        @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
            @require DiffEqCallbacks = "459566f4-90b8-5000-8ac3-15dfb0a30def" begin
                include("../ext/ManifoldsOrdinaryDiffEqDiffEqCallbacksExt.jl")
            end
        end

        @require OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed" begin
            include("../ext/ManifoldsOrdinaryDiffEqExt.jl")
        end

        @require BoundaryValueDiffEq = "764a87c0-6b3e-53db-9096-fe964310641d" begin
            include("../ext/ManifoldsBoundaryValueDiffEqExt.jl")
        end

        @require NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56" begin
            include("../ext/ManifoldsNLsolveExt.jl")
        end

        @require RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01" begin
            @require Colors = "5ae59095-9a9b-59fe-a467-6f913c188581" begin
                include("../ext/ManifoldsRecipesBaseExt.jl")
            end
        end

        @require Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40" begin
            include("../ext/ManifoldsTestExt/ManifoldsTestExt.jl")
        end
    end

    return nothing
end

include("deprecated.jl")

export test_manifold
export test_group, test_action

# Abstract main types
export CoTVector, AbstractManifold, AbstractManifoldPoint, TVector
# Manifolds
export AbstractSphere, AbstractProjectiveSpace
export Euclidean,
    ArrayProjectiveSpace,
    ArraySphere,
    CenteredMatrices,
    CholeskySpace,
    Circle,
    Elliptope,
    EssentialManifold,
    FixedRankMatrices,
    Flag,
    GeneralizedGrassmann,
    GeneralizedStiefel,
    Grassmann,
    HamiltonianMatrices,
    HeisenbergGroup,
    Hyperbolic,
    Hyperrectangle,
    KendallsPreShapeSpace,
    KendallsShapeSpace,
    Lorentz,
    MultinomialDoubleStochastic,
    MultinomialMatrices,
    MultinomialSymmetric,
    MultinomialSymmetricPositiveDefinite,
    Oblique,
    OrthogonalMatrices,
    PositiveArrays,
    PositiveMatrices,
    PositiveNumbers,
    PositiveVectors,
    ProbabilitySimplex,
    ProjectiveSpace,
    Rotations,
    SkewHermitianMatrices,
    SkewSymmetricMatrices,
    Spectrahedron,
    Sphere,
    SphereSymmetricMatrices,
    Stiefel,
    SymmetricMatrices,
    SymmetricPositiveDefinite,
    SPDFixedDeterminant,
    SymmetricPositiveSemidefiniteFixedRank,
    Symplectic,
    SymplecticGrassmann,
    SymplecticMatrices,
    SymplecticStiefel,
    SymplecticElement,
    Torus,
    Tucker,
    UnitaryMatrices
# Point representation types
export HyperboloidPoint,
    PoincareBallPoint,
    PoincareHalfSpacePoint,
    SVDMPoint,
    TuckerPoint,
    StiefelPoint,
    ProjectorPoint,
    SPDPoint
# Tangent vector representation types
export HyperboloidTVector,
    PoincareBallTVector,
    PoincareHalfSpaceTVector,
    TuckerTVector,
    UMVTVector,
    ProjectorTVector,
    StiefelTVector
export AbstractNumbers, ℝ, ℂ, ℍ
export Hamiltonian
# decorator manifolds
export AbstractDecoratorManifold
export IsIsometricEmbeddedManifold, IsEmbeddedManifold, IsEmbeddedSubmanifold
export IsDefaultMetric, IsDefaultConnection, IsMetricManifold, IsConnectionManifold
export ValidationManifold, ValidationMPoint, ValidationTVector, ValidationCoTVector
export Fiber, FiberBundle, CotangentBundle, CotangentSpace, FVector
export AbstractPowerManifold,
    AbstractPowerRepresentation,
    ArrayPowerRepresentation,
    NestedPowerRepresentation,
    NestedReplacingPowerRepresentation,
    PowerManifold,
    QuotientManifold
export ProductManifold, EmbeddedManifold
export GraphManifold, GraphManifoldType, VertexManifold, EdgeManifold
export ArrayPartition
export ProjectedPointDistribution, TangentBundle
export TangentSpace, VectorSpaceFiber, VectorSpaceType, VectorBundle
export AbstractVectorTransportMethod,
    DifferentiatedRetractionVectorTransport, ParallelTransport, ProjectedPointDistribution
export PoleLadderTransport, SchildsLadderTransport
export ProductVectorTransport
export AbstractAffineConnection, ConnectionManifold, LeviCivitaConnection
export AbstractCartanSchoutenConnection,
    CartanSchoutenMinus, CartanSchoutenPlus, CartanSchoutenZero
export MetricManifold
# Metric types
export AbstractMetric,
    AffineInvariantMetric,
    BuresWassersteinMetric,
    CanonicalMetric,
    EuclideanMetric,
    ExtendedSymplecticMetric,
    FisherRaoMetric,
    GeneralizedBuresWassersteinMetric,
    LogCholeskyMetric,
    LogEuclideanMetric,
    LorentzMetric,
    MinkowskiMetric,
    PowerMetric,
    ProductMetric,
    RealSymplecticMetric,
    RiemannianMetric,
    StiefelSubmersionMetric
export AbstractAtlas, RetractionAtlas
# Vector transport types
export AbstractVectorTransportMethod, ParallelTransport, ProjectionTransport
# Retraction types
export AbstractRetractionMethod,
    CayleyRetraction,
    ExponentialRetraction,
    QRRetraction,
    PolarRetraction,
    ProjectionRetraction,
    SoftmaxRetraction,
    ODEExponentialRetraction,
    OrthographicRetraction,
    PadeRetraction,
    ProductRetraction,
    SasakiRetraction
# Inverse Retraction types
export AbstractInverseRetractionMethod,
    ApproximateInverseRetraction,
    ApproximateLogarithmicMap,
    CayleyInverseRetraction,
    LogarithmicInverseRetraction,
    QRInverseRetraction,
    OrthographicInverseRetraction,
    PolarInverseRetraction,
    ProjectionInverseRetraction,
    ShootingInverseRetraction,
    SoftmaxInverseRetraction
# Estimation methods for median and mean
export AbstractApproximationMethod,
    GradientDescentEstimation,
    CyclicProximalPointEstimation,
    EfficientEstimator,
    GeodesicInterpolation,
    GeodesicInterpolationWithinRadius,
    ExtrinsicEstimation
# Tangent space bases
export CachedBasis,
    DefaultBasis,
    DefaultOrthogonalBasis,
    DefaultOrthonormalBasis,
    DiagonalizingOrthonormalBasis,
    InducedBasis,
    ProjectedOrthonormalBasis
# Errors on Manifolds
export ComponentManifoldError, CompositeManifoldError, ManifoldDomainError
# Functions on Manifolds
export ×,
    ^,
    action_side,
    allocate,
    allocate_result,
    base_manifold,
    bundle_projection,
    canonical_project,
    canonical_project!,
    change_basis,
    change_basis!,
    change_metric,
    change_metric!,
    change_representer,
    change_representer!,
    check_chart_switch,
    check_point,
    check_vector,
    christoffel_symbols_first,
    christoffel_symbols_second,
    christoffel_symbols_second_jacobian,
    convert,
    complex_dot,
    decorated_manifold,
    default_approximation_method,
    default_inverse_retraction_method,
    default_retraction_method,
    default_vector_transport_method,
    det_local_metric,
    differential_canonical_project,
    differential_canonical_project!,
    distance,
    dual_basis,
    einstein_tensor,
    embed,
    embed!,
    exp,
    exp!,
    flat,
    flat!,
    gaussian_curvature,
    geodesic,
    geodesic!,
    get_default_atlas,
    get_component,
    get_embedding,
    get_orbit_action,
    get_total_space,
    hat,
    hat!,
    horizontal_component,
    horizontal_component!,
    horizontal_lift,
    horizontal_lift!,
    identity_element,
    identity_element!,
    induced_basis,
    incident_log,
    injectivity_radius,
    inner,
    inverse_local_metric,
    inverse_retract,
    inverse_retract!,
    isapprox,
    is_default_connection,
    is_default_metric,
    is_flat,
    is_group_manifold,
    is_hamiltonian,
    is_identity,
    is_point,
    is_vector,
    kurtosis,
    local_metric,
    local_metric_jacobian,
    log,
    log!,
    log_local_metric_density,
    manifold_dimension,
    manifold_volume,
    metric,
    mean,
    mean!,
    mean_and_var,
    mean_and_std,
    median,
    median!,
    mid_point,
    mid_point!,
    minkowski_metric,
    moment,
    norm,
    normal_tvector_distribution,
    number_eltype,
    number_of_coordinates,
    one,
    power_dimensions,
    parallel_transport_along,
    parallel_transport_along!,
    parallel_transport_direction,
    parallel_transport_direction!,
    parallel_transport_to,
    parallel_transport_to!,
    project,
    project!,
    projected_distribution,
    rand,
    rand!,
    real_dimension,
    ricci_curvature,
    ricci_tensor,
    representation_size,
    retract,
    retract!,
    riemannian_gradient,
    riemannian_gradient!,
    riemannian_Hessian,
    riemannian_Hessian!,
    riemann_tensor,
    riemann_tensor!,
    sectional_curvature,
    sectional_curvature_matrix,
    sectional_curvature_max,
    sectional_curvature_min,
    set_component!,
    sharp,
    sharp!,
    shortest_geodesic,
    shortest_geodesic!,
    skewness,
    std,
    sym_rem,
    symplectic_inverse,
    symplectic_inverse!,
    symplectic_inverse_times,
    symplectic_inverse_times!,
    submanifold,
    submanifold_component,
    submanifold_components,
    uniform_distribution,
    var,
    vector_space_dimension,
    vector_transport_along,
    vector_transport_along!,
    vector_transport_direction,
    vector_transport_direction!,
    vector_transport_to,
    vector_transport_to!,
    vee,
    vee!,
    vertical_component,
    vertical_component!,
    volume_density,
    Weingarten,
    Weingarten!,
    zero_vector,
    zero_vector!
# Lie group types & functions
export AbstractGroupAction,
    AbstractGroupOperation,
    ActionDirection,
    AdditionOperation,
    CircleGroup,
    GeneralLinear,
    GroupManifold,
    GroupOperationAction,
    Identity,
    LeftAction,
    LeftInvariantMetric,
    LeftSide,
    MultiplicationOperation,
    Orthogonal,
    PowerGroup,
    ProductGroup,
    ProductOperation,
    RealCircleGroup,
    RightAction,
    RightInvariantMetric,
    RightSide,
    RotationAction,
    RotationTranslationAction,
    RotationTranslationActionOnVector,
    SemidirectProductGroup,
    SpecialEuclidean,
    SpecialLinear,
    SpecialOrthogonal,
    SpecialUnitary,
    TranslationGroup,
    TranslationAction,
    Unitary
export AbstractInvarianceTrait
export IsMetricManifold, IsConnectionManifold
export IsGroupManifold,
    HasLeftInvariantMetric, HasRightInvariantMetric, HasBiinvariantMetric
export adjoint_action,
    adjoint_action!,
    adjoint_apply_diff_group,
    adjoint_apply_diff_group!,
    adjoint_inv_diff,
    adjoint_inv_diff!,
    affine_matrix,
    apply,
    apply!,
    apply_diff,
    apply_diff!,
    apply_diff_group,
    apply_diff_group!,
    base_group,
    center_of_orbit,
    has_approx_invariant_metric,
    compose,
    compose!,
    direction,
    direction_and_side,
    exp_lie,
    exp_lie!,
    group_manifold,
    geodesic,
    get_coordinates_lie,
    get_coordinates_lie!,
    get_coordinates_orthogonal,
    get_coordinates_orthonormal,
    get_coordinates_orthogonal!,
    get_coordinates_orthonormal!,
    get_coordinates_vee!,
    get_vector_diagonalizing!,
    get_vector_lie,
    get_vector_lie!,
    get_vector_orthogonal,
    get_vector_orthonormal,
    has_biinvariant_metric,
    has_invariant_metric,
    identity_element,
    identity_element!,
    inv,
    inv!,
    inv_diff,
    inv_diff!,
    inverse_apply,
    inverse_apply!,
    inverse_apply_diff,
    inverse_apply_diff!,
    inverse_translate,
    inverse_translate!,
    inverse_translate_diff,
    inverse_translate_diff!,
    lie_bracket,
    lie_bracket!,
    log_lie,
    log_lie!,
    optimal_alignment,
    optimal_alignment!,
    screw_matrix,
    switch_direction,
    switch_side,
    translate,
    translate!,
    translate_diff,
    translate_diff!
# Orthonormal bases
export AbstractBasis,
    AbstractOrthonormalBasis,
    DefaultOrthonormalBasis,
    DiagonalizingOrthonormalBasis,
    ProjectedOrthonormalBasis,
    CachedBasis,
    DiagonalizingBasisData,
    ProductBasisData,
    PowerBasisData
export OutOfInjectivityRadiusError
export get_basis,
    get_coordinates, get_coordinates!, get_vector, get_vector!, get_vectors, number_system

# atlases and charts
export get_point, get_point!, get_parameters, get_parameters!

end # module
