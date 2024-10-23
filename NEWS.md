# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.10.5] – unreleased

### Changed

* rewrote the `CONTRIBUTING.md` and adapt it to todays links and references.
* Deprecate direction and side parameter in group exponential/logarithm retractions

## [0.10.4] - 2024-10-20

### Added

* `uniform_distribution` now has an error hint explaining what has to be done to make it work.
* `Euclidean` now follows the new `has_components` function from `ManifoldsBase.jl` (0.15.18)
  and can handle also the `r`-norms now.
- Union type `MatrixGroup`
- Columnwise group action with arbitrary matrix groups
- `uniform_distribution` now has an error hint explaining what has to be done to make it work.
- `lie_bracket` is exactly zero on orthogonal Lie algebra in 2D

## [0.10.3] - 2024-10-04

### Changed

* **Mildly breaking**: the number system parameter now corresponds to the coefficients standing in front of basis vectors in a linear combination instead of components of a vector. For example, `DefaultOrthonormalBasis() == DefaultOrthonormalBasis(ℝ)` of `Euclidean(3, field=ℂ)` now has 6 vectors, and `DefaultOrthonormalBasis(ℂ)` of the same manifold has 3 basis vectors.

### Fixed

* Fixed `solve_exp_ode` only returning the starting position ([#744](https://github.com/JuliaManifolds/Manifolds.jl/issues/744))
* Fixed documentation of `solve_exp_ode` function signature ([#740](https://github.com/JuliaManifolds/Manifolds.jl/issues/740))

## [0.10.2] - 2024-09-24

### Added

* `GroupManifold` can now be called with two arguments, the third one defaulting to `LeftInvariantRepresentation`.

### Changed

* fixes a few typographical errors.

## [0.10.1] – 2024-08-29

### Changed

* `identity_element` on `ProductManifold` without `RecursiveArrayTools.jl` now prints a useful error message.

## [0.10.0] – 2024-08-24

### Changed

* `Distributions.jl`, `RecursiveArrayTools.jl` and `HybridArrays.jl` were moved to weak dependencies to reduce load time and improve extensibility.
* `translate_diff`, `inv_diff` and thus `apply_diff_group`, are available for all the groups with invariant tangent vector storage.
* `SpecialEuclidean` group now has a different default tangent vector representation, the left-invariant one; to get the old representation pass `vectors=HybridTangentRepresentation()` to the constructor of `SpecialEuclidean`.
* `adjoint_action` takes a direction argument; by default it is `LeftAction`.
* `adjoint_action!` is the necessary method to implement in groups with left-invariant tangent vector representation.
* Fixed a few typos in the doc string of the SPD fixed determinant description.
* Random point on the `MultinomialSymmetricPositiveDefinite` manifold was improved to make it more robust.

### Added

* Introduced `exp_inv` and `log_inv` based on `exp_lie` and `log_lie`. They are invariant to the group operation.
* A tutorial about usage of group-related functionality.

### Removed

* Deprecated bindings:
  * `ExtrinsicEstimation()` (should be replaced with `ExtrinsicEstimation(EfficientEstimator())`),
  * `Symplectic` (renamed to `SymplecticMatrices`),
  * `SymplecticMatrix` (renamed to `SymplecticElement`).
  * `AbstractEstimationMethod` (renamed to `AbstractApproximationMethod`).
  * `VectorBundleVectorTransport` (renamed to `FiberBundleProductVectorTransport`).
  * `rand` on `SymplecticMatrices` and `SymplecticStiefel` no longer accepts `hamiltonian_norm` as an alias for `σ`.
  * `mean!` and `median!` no longer accept `extrinsic_method` (should be replaced with `e = ExtrinsicEstimation(extrinsic_method)`).
* As a result of making `Distributions.jl` and `RecursiveArrayTools.jl` weak dependencies the following symbols are no longer exported from `Manifolds.jl`. Essential functionality is still available but distribution-related features may change in the future without a breaking release.
  * `ArrayPartition` (`RecursiveArrayTools.jl` needs to be explicitly imported),
  * `ProjectedPointDistribution` (not exported),
  * `normal_tvector_distribution` (not exported),
  * `projected_distribution` (not exported),
  * `uniform_distribution` (not exported).
* Ability to create non-real `SymplecticStiefel` and `SymplecticGrassmann` manifolds; essential functionality was missing so it was removed until a more developed version is developed.

## [0.9.20] – 2024-06-17

### Added

* implemented parallel transport on the Grassmann manifold with respect to Stiefel representation

### Changed

* since now all exp/log/parallel transport are available for all representations of `Grassmann`,
  these are now also set as defaults, since they are more exact.

## [0.9.19] – 2024-06-12

### Changed

* Updated `Project.toml` compatibility entries.
* Updated CI for Julia 1.11-beta

### Fixed

* a few typos in the doc string of the SPD fixed determinant description
* several other typographical errors throughout the documentation

## [0.9.18] – 2024-05-07

### Added

* added the injectivity radius for the Stiefel manifold with Euclidean metric

## [0.9.17] – 2024-04-23

### Added

* `Hyperrectangle` manifold with boundary.

## [0.9.16] – 2024-04-01

### Changed

* `NonlinearSolve.jl` and `PythonCall.jl` are no longer an upper bounded dependency (bugs were fixed).

## [0.9.15] – 2024-03-24

### Added

* using `DocumenterInterLinks` for links to other Julia packages documentation.
* Implementation of `sectional_curvature`, `sectional_curvature_min` and `sectional_curvature_max` for several manifolds.
* `sectional_curvature_matrix` function and a tutorial on coordinate-free curvature.

### Changed

* `default_vector_transport_method` for `GeneralUnitaryMatrices` other than `Rotations` was changed to `ProjectionTransport`.

### Fixed

* typographical errors in tutorials/working-in-charts.jl.
* several typographical errors in the docs
* unifies to use two backticks ``` `` ``` for math instead of ` $ ` further in the docs

## [0.9.14] – 2024-01-31

### Added

* `rand` on `UnitaryMatrices`
* `rand` on arbitrary `GroupManifold`s and manifolds with `IsGroupManifold` trait
  generating points and elements from the Lie algebra, respectively

## [0.9.13] – 2024-01-24

### Added

* added the real symplectic Grassmann manifold `SymplecticGrassmann`
* Introduce the manifold of `HamiltonianMatrices` and a wrapper for `Hamiltonian` matrices
* introduce `rand(:HamiltonianMatrices)`
* extend `rand` to also `rand!` for `HamiltonianMatrices`, `SymplecticMatrices` and `SymplecticStiefel`
* implement `riemannian_gradient` conversion for `SymplecticMatrices` and `SymplecticGrassmann`
* the new manifold of `MultinomialSymmetricPositiveDefinite` matrices
* `rand!` for `MultinomialDoublyStochastic` and `MultinomialSymmetric`

### Deprecated

* Rename `Symplectic` to `SimplecticMatrices` in order to have a `Symplectic` wrapper for such matrices as well in the future for the next breaking change.
* Rename `SymplecticMatrix` to `SymplecticElement` to clarify that it is the special matrix ``J_{2n}`` and not an arbitrary symplectic matrix.

### Fixed

* a bug that cause `project` for tangent vectors to return wrong results on `MultinomialDoublyStochastic`

## [0.9.12] – 2024-01-21

### Fixed

* Fixed `var` on `TranslationGroup`.

## [0.9.11] – 2023-12-27

### Fixed

* Fixed mixed array index number in-place `parallel_transport_to!` on zero-index `Euclidean`, on Julia 1.6.

## [0.9.10] – 2023-12-27

### Added

* Compatibility with `RecursiveArrayTools` v3.

### Fixed

* Fixed mixed array index number in-place `parallel_transport_to!` on real `Circle`, on Julia 1.6.
* Violations of MD004 lint rule in this file.

## [0.9.9] – 2023-12-25

### Fixed

* introduced a nonzero `atol` for all point and vector checks that compre to zero.
  This makes those checks a bit more relaxed by default and resolves [#630](https://github.com/JuliaManifolds/Manifolds.jl/issues/630).
* `default_estimation_method(M, f)` is deprecated, use `default_approximation_method(M, f)` for your specific method `f` on the manifold `M`.
* `AbstractEstimationMethod` is deprecated, use `AbstractApproximationMethod` instead.

## [0.9.8] - 2023-11-17

### Fixed

* Improved distribution of random vector generation for rotation matrices and complex circle.

## [0.9.7] – 2023-11-14

### Fixed

* Fixed `is_flat` for `CholeskySpace` and `SymmetricPositiveDefinite` with `LogCholeskyMetric` [https://github.com/JuliaManifolds/Manifolds.jl/issues/684](https://github.com/JuliaManifolds/Manifolds.jl/issues/684).

## [0.9.6] - 2023-11-09

### Fixed

* Fixed real coefficient basis for complex circle (an issue exposed by [https://github.com/JuliaManifolds/ManifoldsBase.jl/pull/173](https://github.com/JuliaManifolds/ManifoldsBase.jl/pull/173)).
* Fixed `VeeOrthogonalBasis` test for non-real manifolds.

## [0.9.5] - 2023-11-08

### Changed

* `identity_element` now returns a complex matrix for unitary group.
* `number_of_coordinates` is now exported.

## [0.9.4] - 2023-11-06

### Added

* Functions `inv_diff`, `inv_diff!`, `adjoint_inv_diff` and `adjoint_inv_diff!` that correspond to differentials and pullbacks of group inversion.
* Julia 1.10-rc CI workflow.

### Changed

* Documentation project files are marked as compatible with `BoundaryValueDiffEq` v5.

### Fixed

* Fixed issue with incorrect implementation of `apply_diff_group` in `GroupOperationAction` with left backward and right forward action [#669](https://github.com/JuliaManifolds/Manifolds.jl/issues/669).

## [0.9.3] - 2023-10-28

### Added

* Support for `BoundaryValueDiffEq` v5.

## [0.9.2] - 2023-10-27

### Added

* `rand(G; vector_at=Identity(G))` now works for translation, special orthogonal and special Euclidean groups `G` (issue [#665](https://github.com/JuliaManifolds/Manifolds.jl/issues/665)).
* `get_embedding` now works for `GeneralUnitaryMultiplicationGroup`.
* Github action that checks for NEWS.md changes.

## [0.9.1] - 2023-10-25

### Added

* a new retraction and its inverse for the fixed Rank Manifolds, the orthographic retraction.

## [0.9.0] - 2023-10-24

### Added

* Vector bundles are generalized to fiber bundles. Old `BundleFibers` functionality was reworked to better match mathematical abstractions. Fiber bundle functionality is experimental and minor changes may happen without a breaking release, with the exception of `TangentBundle` which is considered to be stable.
* `RotationTranslationAction` is introduced.

### Changed

* Sizes of all manifolds can now be either encoded in type or stored in a field to avoid over-specialization.
  The default is set to store the size in type parameter (except for `PowerManifold` and its variants), replicating the previous behavior.
  For field storage, pass the `parameter=:field` keyword argument to manifold constructor.
  For example statically sized `CenteredMatrices{m,n}` is now `CenteredMatrices{TypeParameter{Tuple{m,n}}}`, whereas the type of special Euclidean group with field-stored size is `CenteredMatrices{Tuple{Int,Int}}`. Similar change applies to:
  * `CenteredMatrices{m,n}`,
  * `CholeskySpace{N}`,
  * `Elliptope{N,K}`,
  * `Euclidean`,
  * `FixedRankMatrices{m,n,k}`,
  * `KendallsPreShapeSpace{n,k}`,
  * `KendallsShapeSpace{n,k}`,
  * `GeneralLinear{n}`,
  * `GeneralUnitaryMultiplicationGroup{n}`,
  * `GeneralizedGrassmann{n,k}`,
  * `GeneralizedStiefel{n,k}`,
  * `Grassmann{n,k}`,
  * `Heisenberg{n}`,
  * `Hyperbolic{n}`,
  * `MultinomialMatrices{N,M}`,
  * `MultinomialDoublyStochastic{n}`,
  * `MultinomialSymmetric{n}`,
  * `Orthogonal{n}`,
  * `PowerManifold`,
  * `PositiveArrays`,
  * `PositiveMatrices`,
  * `PositiveNumbers`,
  * `ProbabilitySimplex{n}`,
  * `SPDFixedDeterminant{n}`,
  * `SpecialLinear{n}`,
  * `SpecialOrthogonal{n}`,
  * `SpecialUnitary{n}`,
  * `SpecialEuclidean{n}`,
  * `SpecialEuclideanManifold{n}`,
  * `Spectrahedron{n,k}`,
  * `SphereSymmetricMatrices{N}`,
  * `Stiefel{n,k}`,
  * `SymmetricMatrices{N}`,
  * `SymmetricPositiveDefinite{n}`,
  * `SymmetricPositiveSemidefiniteFixedRank{n,k}`,
  * `Symplectic{n}`,
  * `SymplecticStiefel{n,k}`,
  * `TranslationGroup`,
  * `Tucker`.

  For example

  ```{julia}
  function Base.show(io::IO, ::CenteredMatrices{m,n}) where {m,n}
      return print(io, "CenteredMatrices($m, $n)")
  end
  ```

  needs to be replaced with

  ```{julia}
  function Base.show(io::IO, ::CenteredMatrices{TypeParameter{Tuple{m,n}}}) where {m,n}
      return print(io, "CenteredMatrices($m, $n)")
  end
  ```

  for statically-sized groups and

  ```{julia}
  function Base.show(io::IO, M::CenteredMatrices{Tuple{Int,Int}})
      m, n = get_parameter(M.size)
      return print(io, "CenteredMatrices($m, $n; parameter=:field)")
  end
  ```

  for groups with size stored in field. Alternatively, you can use a single generic method like this:

  ```{julia}
  function Base.show(io::IO, M::CenteredMatrices{T}) where {T}
      m, n = get_parameter(M)
      if T <: TypeParameter
          return print(io, "CenteredMatrices($m, $n)")
      else
          return print(io, "CenteredMatrices($m, $n; parameter=:field)")
      end
  end
  ```

* Argument order for type aliases `RotationActionOnVector` and `RotationTranslationActionOnVector`: most often dispatched on argument is now first.
* A more consistent handling of action direction was introduced. 4-valued `ActionDirection` was split into 2-valued `ActionDirection` (either left or right action) and `GroupActionSide` (action acting from the left or right side). See [https://github.com/JuliaManifolds/Manifolds.jl/issues/637](https://github.com/JuliaManifolds/Manifolds.jl/issues/637) for a design discussion.

### Removed

* `ProductRepr` is removed; please use `ArrayPartition` instead.
* Default methods throwing "not implemented" `ErrorException` for some group-related operations. Standard `MethodError` is now thrown instead.
* `LinearAffineMetric` was deprecated in a previous release and the symbol is now removed.
  Please use `AffineInvariantMetric` instead.
