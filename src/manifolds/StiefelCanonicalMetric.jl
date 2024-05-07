@doc raw"""
    CanonicalMetric <: AbstractMetric

The Canonical Metric refers to a metric for the [`Stiefel`](@ref)
manifold, see[EdelmanAriasSmith:1998](@cite).

"""
struct CanonicalMetric <: RiemannianMetric end

"""
    ApproximateLogarithmicMap <: ApproximateInverseRetraction

An approximate implementation of the logarithmic map, which is an [`inverse_retract`](@ref)ion.
See [`inverse_retract(::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, ::Any, ::Any, ::ApproximateLogarithmicMap)`](@ref) for a use case.

# Fields

* `max_iterations` – maximal number of iterations used in the approximation
* `tolerance` – a tolerance used as a stopping criterion

"""
struct ApproximateLogarithmicMap{T} <: ApproximateInverseRetraction
    max_iterations::Int
    tolerance::T
end

function distance(M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, q, p)
    return norm(M, p, log(M, p, q))
end

@doc raw"""
    q = exp(M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, p, X)
    exp!(M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, q, p, X)

Compute the exponential map on the [`Stiefel`](@ref)`(n, k)` manifold with respect to the [`CanonicalMetric`](@ref).

First, decompose The tangent vector ``X`` into its horizontal and vertical component with
respect to ``p``, i.e.

```math
X = pp^{\mathrm{T}}X + (I_n-pp^{\mathrm{T}})X,
```
where ``I_n`` is the ``n×n`` identity matrix.
We introduce ``A=p^{\mathrm{T}}X`` and ``QR = (I_n-pp^{\mathrm{T}})X`` the `qr` decomposition
of the vertical component. Then using the matrix exponential ``\operatorname{Exp}`` we introduce ``B`` and ``C`` as

```math
\begin{pmatrix}
B\\C
\end{pmatrix}
\coloneqq
\operatorname{Exp}\left(
\begin{pmatrix}
A & -R^{\mathrm{T}}\\ R & 0
\end{pmatrix}
\right)
\begin{pmatrix}I_k\\0\end{pmatrix}
```

the exponential map reads

```math
q = \exp_p X = pC + QB.
```
For more details, see [EdelmanAriasSmith:1998](@cite)[Zimmermann:2017](@cite).
"""
exp(::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, ::Any...)

function exp!(M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, q, p, X)
    n, k = get_parameter(M.manifold.size)
    A = p' * X
    n == k && return mul!(q, p, exp(A))
    QR = qr(X - p * A)
    BC_ext = exp([A -QR.R'; QR.R 0*I])
    @views begin # COV_EXCL_LINE
        mul!(q, p, BC_ext[1:k, 1:k])
        mul!(q, Matrix(QR.Q), BC_ext[(k + 1):(2 * k), 1:k], true, true)
    end
    return q
end

@doc raw"""
    inner(M::MetricManifold{ℝ, Stiefel{<:Any,ℝ}, X, CanonicalMetric}, p, X, Y)

Compute the inner product on the [`Stiefel`](@ref) manifold with respect to the
[`CanonicalMetric`](@ref). The formula reads

```math
g_p(X,Y) = \operatorname{tr}\bigl( X^{\mathrm{T}}(I_n - \frac{1}{2}pp^{\mathrm{T}})Y \bigr).
```
"""
function inner(M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric}, p, X, Y)
    n, k = get_parameter(M.manifold.size)
    T = Base.promote_eltype(p, X, Y)
    if n == k
        return T(dot(X, Y)) / 2
    else
        return T(dot(X, Y)) - T(dot(p'X, p'Y)) / 2
    end
end

@doc raw"""
    X = inverse_retract(M::MetricManifold{ℝ, Stiefel{<:Any,ℝ}, CanonicalMetric}, p, q, a::ApproximateLogarithmicMap)
    inverse_retract!(M::MetricManifold{ℝ, Stiefel{<:Any,ℝ}, X, CanonicalMetric}, p, q, a::ApproximateLogarithmicMap)

Compute an approximation to the logarithmic map on the [`Stiefel`](@ref)`(n, k)` manifold with respect to the [`CanonicalMetric`](@ref)
using a matrix-algebraic based approach to an iterative inversion of the formula of the
[`exp`](@ref exp(::MetricManifold{ℝ, Stiefel{<:Any,ℝ}, CanonicalMetric}, ::Any...)).

The algorithm is derived in [Zimmermann:2017](@cite) and it uses the `max_iterations` and the `tolerance` field
from the [`ApproximateLogarithmicMap`](@ref).
"""
inverse_retract(
    ::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric},
    ::Any,
    ::Any,
    ::ApproximateLogarithmicMap,
)

function log(
    M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric},
    p,
    q;
    maxiter::Int=10000,
    tolerance=1e-9,
)
    X = allocate_result(M, log, p, q)
    inverse_retract!(M, X, p, q, ApproximateLogarithmicMap(maxiter, tolerance))
    return X
end

function log!(
    M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric},
    X,
    p,
    q;
    maxiter::Int=10000,
    tolerance=1e-9,
)
    inverse_retract!(M, X, p, q, ApproximateLogarithmicMap(maxiter, tolerance))
    return X
end

function inverse_retract!(
    M::MetricManifold{ℝ,<:Stiefel{<:Any,ℝ},CanonicalMetric},
    X,
    p,
    q,
    a::ApproximateLogarithmicMap,
)
    n, k = get_parameter(M.manifold.size)
    qfact = stiefel_factorization(p, q)
    V = allocate(qfact.Z, Size(2k, 2k))
    LV = allocate(V)
    Zcompl = qr(qfact.Z).Q[1:(2k), (k + 1):(2k)]
    @views begin # COV_EXCL_LINE
        Vpcols = V[1:(2k), (k + 1):(2k)] #second half of the columns
        B = LV[(k + 1):(2k), 1:k]
        C = LV[(k + 1):(2k), (k + 1):(2k)]
        copyto!(V[1:(2k), 1:k], qfact.Z)
        F = svd(Zcompl[(k + 1):(2k), 1:k]) # preprocessing: Procrustes
    end
    new_Vpcols = allocate(Vpcols)
    mul!(new_Vpcols, Zcompl, F.U)
    mul!(Vpcols, new_Vpcols, F.V')
    S = allocate(B)
    for _ in 1:(a.max_iterations)
        log_safe!(LV, V)
        norm(C) ≤ a.tolerance && break
        copyto!(S, I)
        mul!(S, B, B', -1 // 12, 1 // 2)
        Γ = lyap(S, C)
        expΓ = exp(Γ)
        mul!(new_Vpcols, Vpcols, expΓ)
        copyto!(Vpcols, new_Vpcols)
    end
    @views mul!(X, qfact.U, LV[1:(2k), 1:k])
    return X
end

@doc raw"""
    Y = riemannian_Hessian(M::MetricManifold{ℝ, Stiefel, CanonicalMetric}, p, G, H, X)
    riemannian_Hessian!(M::MetricManifold{ℝ, Stiefel, CanonicalMetric}, Y, p, G, H, X)

Compute the Riemannian Hessian ``\operatorname{Hess} f(p)[X]`` given the
Euclidean gradient ``∇ f(\tilde p)`` in `G` and the Euclidean Hessian ``∇^2 f(\tilde p)[\tilde X]`` in `H`,
where ``\tilde p, \tilde X`` are the representations of ``p,X`` in the embedding,.

Here, we adopt Eq. (5.6) [Nguyen:2023](@cite), for the [`CanonicalMetric`](@ref)
``α_0=1, α_1=\frac{1}{2}`` in their formula. The formula reads

```math
    \operatorname{Hess}f(p)[X]
    =
    \operatorname{proj}_{T_p\mathcal M}\Bigl(
        ∇^2f(p)[X] - \frac{1}{2} X \bigl( (∇f(p))^{\mathrm{H}}p + p^{\mathrm{H}}∇f(p)\bigr)
        - \frac{1}{2} \bigl( P ∇f(p) p^{\mathrm{H}} + p ∇f(p))^{\mathrm{H}} P)X
    \Bigr),
```
where ``P = I-pp^{\mathrm{H}}``.
"""
riemannian_Hessian(M::MetricManifold{𝔽,Stiefel,CanonicalMetric}, p, G, H, X) where {𝔽}

function riemannian_Hessian!(
    M::MetricManifold{𝔽,<:Stiefel{<:Any,𝔽},CanonicalMetric},
    Y,
    p,
    G,
    H,
    X,
) where {𝔽}
    Gp = symmetrize(G' * p)
    Z = symmetrize((I - p * p') * G * p')
    project!(M, Y, p, H - X * Gp - Z * X)
    return Y
end
