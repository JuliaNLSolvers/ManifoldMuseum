@doc raw"""
    BurresWassertseinMetric <: AbstractMetric

The Bures Wasserstein metric for symmetric positive definite matrices[^MalagoMontruccioPistone2018].

[^MalagoMontruccioPistone2018]:
    > Malagò, L., Montrucchio, L., Pistone, G.:
    > _Wasserstein Riemannian geometry of Gaussian densities_.
    > Information Geometry, 1, pp. 137–179, 2018.
    > doi: [10.1007/s41884-018-0014-4](https://doi.org/10.1007/s41884-018-0014-4)
"""
struct BuresWassersteinMetric <: RiemannianMetric end

@doc raw"""
    change_representer(M::MetricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric}, E::EuclideanMetric, p, X)

Given a tangent vector ``X ∈ T_p\mathcal M`` representing a linear function on the tangent
space at `p` with respect to the [`EuclideanMetric`](@ref) `g_E`,
this is turned into the representer with respect to the (default) metric,
the [`BuresWassersteinMetric`](@ref) on the [`SymmetricPositiveDefinite`](@ref) `M`.

To be precise we are looking for ``Z∈T_p\mathcal P(n)`` such that for all ``Y∈T_p\mathcal P(n)```
it holds

```math
⟨X,Y⟩ = \operatorname{tr}(XY) = ⟨Z,Y⟩_{\mathrm{BW}}
```
for all ``Y`` and hence we get
``Z``= 2(A+A^{\mathrm{T}})`` with ``A=Xp``.
"""
change_representer(
    ::MetricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric},
    ::EuclideanMetric,
    p,
    X,
)

function change_representer!(
    ::MetricManifold{ℝ,<:SymmetricPositiveDefinite,BuresWassersteinMetric},
    Y,
    ::EuclideanMetric,
    p,
    X,
)
    Y .= 2 .* (X * p + p' * X')
    return Y
end

@doc raw"""
    distance(::MatricManifold{SymmetricPositiveDefinite,BuresWassersteinMetric}, p, q)

Compute the distance with respect to the [`BuresWassersteinMetric`](@ref) on [`SymmetricPositiveDefinite`](@ref) matrices, i.e.

```math
    d(p.q) = \bigl( \operatorname{tr}(p) + \operatorname{tr}(q) + 2\operatorname{tr}(p^{\frac{1}{2}qp^{\frac{1}[2}) \bigr)^\frac{1}{2},
```

where the last trace can be simplified (by rotating the matrix products in the trace) to ``\operatorname{tr}(pq)``.
"""
function distance(
    ::MetricManifold{ℝ,<:SymmetricPositiveDefinite,BuresWassersteinMetric},
    p,
    q,
)
    return sqrt(tr(p) + tr(q) + 2 * sqrt(tr(q * p)))
end

@doc raw"""
    exp(::MatricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric}, p, X)

Compute the exponential map on [`SymmetricPositiveDefinite`](@ref) with respect to
the [`BuresWassersteinMetric`](@ref) given by

```math
    \exp_p(X) = p+X+\mathcal L_p(X)pL_p(X)
```

where ``q=L_p(X)`` denotes the lyaponov operator, i.e. it solves ``pq + qp = X``.
"""
exp(::MetricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric}, p, X)

function exp!(
    ::MetricManifold{ℝ,<:SymmetricPositiveDefinite,BuresWassersteinMetric},
    q,
    p,
    X,
)
    q .= lyap(p, -X) #lyap solves qp+pq-X=0
    q .= p .+ X .+ q * p * q
    return q
end

@doc raw"""
    inner(::MetricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric}, p, X, Y)

Compute the inner product [`SymmetricPositiveDefinite`](@ref) with respect to
the [`BuresWassersteinMetric`](@ref) given by

```math
    ⟨X,Y⟩_{\mathrm{BW}} = \frac{1}{2}\operatorname{tr}( \mathcal L_p(X)Y)
```

where ``q=L_p(X)`` denotes the lyaponov operator, i.e. it solves ``pq + qp = X``.
"""
function inner(
    ::MetricManifold{ℝ,<:SymmetricPositiveDefinite,BuresWassersteinMetric},
    p,
    X,
    Y,
)
    return 1 / 2 * tr(lyap(p, -X) * Y)
end

@doc raw"""
    log(::MatricManifold{SymmetricPositiveDefinite,BuresWassersteinMetric}, p, q)

Compute the logarithmic map on [`SymmetricPositiveDefinite`](@ref) with respect to
the [`BuresWassersteinMetric`](@ref) given by

```math
    \loq_p(q) = (pq)^{\frac{1}{2}} + (qp)^{\frac{1}{2}} - 2 p
```

where ``q=L_p(X)`` denotes the lyaponov operator, i.e. it solves ``pq + qp = X``.
"""
log(::MetricManifold{ℝ,SymmetricPositiveDefinite,BuresWassersteinMetric}, p, q)

function log!(
    ::MetricManifold{ℝ,<:SymmetricPositiveDefinite,BuresWassersteinMetric},
    X,
    p,
    q,
)
    X .= Symmetric(sqrt(p * q) + sqrt(q * p)) - 2 * p
    return X
end