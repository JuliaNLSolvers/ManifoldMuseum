@doc raw"""
    AbstractAffineConnection

Abstract type for affine connections on a manifold.
"""
abstract type AbstractAffineConnection end

"""
    LeviCivitaConnection

The [Levi-Civita connection](https://en.wikipedia.org/wiki/Levi-Civita_connection) of a Riemannian manifold.
"""
struct LeviCivitaConnection <: AbstractAffineConnection end

struct MetricDecoratorType <: AbstractDecoratorType end

"""
    AbstractConnectionManifold{𝔽,M<:AbstractManifold{𝔽},G<:AbstractAffineConnection} <: AbstractDecoratorManifold{𝔽}

Equip an [`AbstractManifold`](@ref) explicitly with an [`AbstractAffineConnection`](@ref) `G`.

`AbstractConnectionManifold` is defined by values of [`christoffel_symbols_second`](@ref),
which is used for a default implementation of [`exp`](@ref), [`log`](@ref) and
[`vector_transport_to`](@ref). Closed-form formulae for particular connection manifolds may
be explicitly implemented when available.

An overview of basic properties of affine connection manifolds can be found in [^Pennec2020].

[^Pennec2020]:
    > X. Pennec and M. Lorenzi, “5 - Beyond Riemannian geometry: The affine connection
    > setting for transformation groups,” in Riemannian Geometric Statistics in Medical Image
    > Analysis, X. Pennec, S. Sommer, and T. Fletcher, Eds. Academic Press, 2020, pp. 169–229.
    > doi: 10.1016/B978-0-12-814725-2.00012-1.
"""
abstract type AbstractConnectionManifold{𝔽} <:
              AbstractDecoratorManifold{𝔽,MetricDecoratorType} end

"""
    connection(M::AbstractManifold)

Get the connection (an object of a subtype of [`AbstractAffineConnection`](@ref))
of [`AbstractManifold`](@ref) `M`.
"""
connection(::AbstractManifold)

"""
    ConnectionManifold(M, C)

Decorate the [`AbstractManifold`](@ref) `M` with [`AbstractAffineConnection`](@ref) `C`.
"""
struct ConnectionManifold{𝔽,M<:AbstractManifold{𝔽},C<:AbstractAffineConnection} <:
       AbstractConnectionManifold{𝔽}
    manifold::M
    connection::C
end

@doc raw"""
    christoffel_symbols_second(
        M::AbstractManifold,
        p,
        B::AbstractBasis;
        backend::AbstractDiffBackend = diff_backend(),
    )

Compute the Christoffel symbols of the second kind in local coordinates of basis `B`.
For affine connection manifold the Christoffel symbols need to be explicitly implemented
while, for a [`MetricManifold`](@ref) they are computed as (in Einstein summation convention)

````math
Γ^{l}_{ij} = g^{kl} Γ_{ijk},
````

where ``Γ_{ijk}`` are the Christoffel symbols of the first kind
(see [`christoffel_symbols_first`](@ref)), and ``g^{kl}`` is the inverse of the local
representation of the metric tensor. The dimensions of the resulting multi-dimensional array
are ordered ``(l,i,j)``.
"""
christoffel_symbols_second(::AbstractManifold, ::Any, ::AbstractBasis)

@decorator_transparent_signature christoffel_symbols_second(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis;
    kwargs...,
)

@doc raw"""
    christoffel_symbols_second_jacobian(
        M::AbstractManifold,
        p,
        B::AbstractBasis;
        backend::AbstractDiffBackend = diff_backend(),
    )

Get partial derivatives of the Christoffel symbols of the second kind
for manifold `M` at `p` with respect to the coordinates of `B`, i.e.

```math
\frac{∂}{∂ p^l} Γ^{k}_{ij} = Γ^{k}_{ij,l}.
```

The dimensions of the resulting multi-dimensional array are ordered ``(i,j,k,l)``.
"""
christoffel_symbols_second_jacobian(::AbstractManifold, ::Any, B::AbstractBasis)
function christoffel_symbols_second_jacobian(
    M::AbstractManifold,
    p,
    B::AbstractBasis;
    backend::AbstractDiffBackend=diff_backend(),
)
    n = size(p, 1)
    ∂Γ = reshape(
        _jacobian(q -> christoffel_symbols_second(M, q, B; backend=backend), p, backend),
        n,
        n,
        n,
        n,
    )
    return ∂Γ
end
@decorator_transparent_signature christoffel_symbols_second_jacobian(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis;
    kwargs...,
)

"""
    connection(M::ConnectionManifold)

Return the connection associated with [`ConnectionManifold`](@ref) `M`.
"""
connection(M::ConnectionManifold) = M.connection

Base.copyto!(M::AbstractConnectionManifold, q, p) = copyto!(M.manifold, q, p)
Base.copyto!(M::AbstractConnectionManifold, Y, p, X) = copyto!(M.manifold, Y, p, X)

@doc raw"""
    exp(M::AbstractConnectionManifold, p, X)

Compute the exponential map on the [`AbstractConnectionManifold`](@ref) `M` equipped with
corresponding affine connection.

If `M` is a [`MetricManifold`](@ref) with a metric that was declared the default metric
using [`is_default_metric`](@ref), this method falls back to `exp(M, p, X)`.

Otherwise it numerically integrates the underlying ODE, see [`solve_exp_ode`](@ref).
Currently, the numerical integration is only accurate when using a single
coordinate chart that covers the entire manifold. This excludes coordinates
in an embedded space.
"""
exp(::AbstractConnectionManifold, ::Any...)

@decorator_transparent_fallback function exp!(M::AbstractConnectionManifold, q, p, X)
    tspan = (0.0, 1.0)
    A = get_default_atlas(M)
    i = get_chart_index(M, A, p)
    B = induced_basis(M, A, i, TangentSpace)
    sol = solve_exp_ode(M, p, X, tspan, B; dense=false, saveat=[1.0])
    return copyto!(q, sol)
end

"""
    gaussian_curvature(M::AbstractManifold, p, B::AbstractBasis; backend::AbstractDiffBackend = diff_backend())

Compute the Gaussian curvature of the manifold `M` at the point `p` using basis `B`.
This is equal to half of the scalar Ricci curvature, see [`ricci_curvature`](@ref).
"""
gaussian_curvature(::AbstractManifold, ::Any, ::AbstractBasis)
function gaussian_curvature(M::AbstractManifold, p, B::AbstractBasis; kwargs...)
    return ricci_curvature(M, p, B; kwargs...) / 2
end
@decorator_transparent_signature gaussian_curvature(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis;
    kwargs...,
)

function injectivity_radius(M::AbstractConnectionManifold, p)
    return injectivity_radius(base_manifold(M), p)
end
function injectivity_radius(M::AbstractConnectionManifold, m::AbstractRetractionMethod)
    return injectivity_radius(base_manifold(M), m)
end
function injectivity_radius(M::AbstractConnectionManifold, m::ExponentialRetraction)
    return injectivity_radius(base_manifold(M), m)
end
function injectivity_radius(M::AbstractConnectionManifold, p, m::AbstractRetractionMethod)
    return injectivity_radius(base_manifold(M), p, m)
end
function injectivity_radius(M::AbstractConnectionManifold, p, m::ExponentialRetraction)
    return injectivity_radius(base_manifold(M), p, m)
end

"""
    ricci_tensor(M::AbstractManifold, p, B::AbstractBasis; backend::AbstractDiffBackend = diff_backend())

Compute the Ricci tensor, also known as the Ricci curvature tensor,
of the manifold `M` at the point `p` using basis `B`,
see [`https://en.wikipedia.org/wiki/Ricci_curvature#Introduction_and_local_definition`](https://en.wikipedia.org/wiki/Ricci_curvature#Introduction_and_local_definition).
"""
ricci_tensor(::AbstractManifold, ::Any, ::AbstractBasis)
function ricci_tensor(M::AbstractManifold, p, B::AbstractBasis; kwargs...)
    R = riemann_tensor(M, p, B; kwargs...)
    n = size(R, 1)
    Ric = allocate(R, Size(n, n))
    @einsum Ric[i, j] = R[l, i, l, j]
    return Ric
end
@decorator_transparent_signature ricci_tensor(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis;
    kwargs...,
)

@doc raw"""
    riemann_tensor(M::AbstractManifold, p, B::AbstractBasis; backend::AbstractDiffBackend=diff_backend())

Compute the Riemann tensor ``R^l_{ijk}``, also known as the Riemann curvature
tensor, at the point `p` in local coordinates defined by `B`. The dimensions of the
resulting multi-dimensional array are ordered ``(l,i,j,k)``.

The function uses the coordinate expression involving the second Christoffel symbol,
see [`https://en.wikipedia.org/wiki/Riemann_curvature_tensor#Coordinate_expression`](https://en.wikipedia.org/wiki/Riemann_curvature_tensor#Coordinate_expression)
for details.

# See also

[`christoffel_symbols_second`](@ref), [`christoffel_symbols_second_jacobian`](@ref)
"""
riemann_tensor(::AbstractManifold, ::Any, ::AbstractBasis)
function riemann_tensor(
    M::AbstractManifold,
    p,
    B::AbstractBasis;
    backend::AbstractDiffBackend=diff_backend(),
)
    n = size(p, 1)
    Γ = christoffel_symbols_second(M, p, B; backend=backend)
    ∂Γ = christoffel_symbols_second_jacobian(M, p, B; backend=backend) ./ n
    R = allocate(∂Γ, Size(n, n, n, n))
    @einsum R[l, i, j, k] =
        ∂Γ[l, i, k, j] - ∂Γ[l, i, j, k] + Γ[s, i, k] * Γ[l, s, j] - Γ[s, i, j] * Γ[l, s, k]
    return R
end
@decorator_transparent_signature riemann_tensor(
    M::AbstractDecoratorManifold,
    p,
    B::AbstractBasis;
    kwargs...,
)

@doc raw"""
    solve_exp_ode(
        M::AbstractConnectionManifold,
        p,
        X,
        tspan,
        B::AbstractBasis;
        backend::AbstractDiffBackend = diff_backend(),
        solver = AutoVern9(Rodas5()),
        kwargs...,
    )

Approximate the exponential map on the manifold over the provided timespan
assuming the default connection of the given manifold by solving the ordinary differential
equation

```math
\frac{d^2}{dt^2} p^k + Γ^k_{ij} \frac{d}{dt} p_i \frac{d}{dt} p_j = 0,
```

where ``Γ^k_{ij}`` are the Christoffel symbols of the second kind, and
the Einstein summation convention is assumed. The arguments `tspan` and
`solver` follow the `OrdinaryDiffEq` conventions. `kwargs...` specify keyword
arguments that will be passed to `OrdinaryDiffEq.solve`.

Currently, the numerical integration is only accurate when using a single
coordinate chart that covers the entire manifold. This excludes coordinates
in an embedded space.

!!! note
    This function only works when
    [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) is loaded with
    ```julia
    using OrdinaryDiffEq
    ```
"""
function solve_exp_ode(M, p, X, tspan, B::AbstractBasis; kwargs...)
    return error(
        "solve_exp_ode not implemented on $(typeof(M)) for point $(typeof(p)), vector $(typeof(X)), and timespan $(typeof(tspan)). For a suitable default, enter `using OrdinaryDiffEq` on Julia 1.1 or greater.",
    )
end

#
# Introduce transparency
# (a) new functions & other parents
for f in [
    christoffel_symbols_second_jacobian,
    exp,
    gaussian_curvature,
    get_coordinates,
    get_vector,
    log,
    mean,
    median,
    project,
    ricci_tensor,
    riemann_tensor,
    vector_transport_along,
    vector_transport_direction,
    vector_transport_direction!, #since it has a default using _to!
    vector_transport_to,
]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractConnectionManifold,
                args...,
            )
                return Val(:parent)
            end
        end,
    )
end

# (b) changes / intransparencies.
for f in [
    christoffel_symbols_second, # this is basic for connection manifolds but not for metric manifolds
    exp!,
    get_coordinates!,
    get_vector!,
    get_basis,
    inner,
    inverse_retract!,
    log!,
    mean!,
    median!,
    norm,
    project!,
    retract!,
    vector_transport_along!,
    vector_transport_to!,
]
    eval(
        quote
            function decorator_transparent_dispatch(
                ::typeof($f),
                ::AbstractConnectionManifold,
                args...,
            )
                return Val(:intransparent)
            end
        end,
    )
end
# (c) special cases
function decorator_transparent_dispatch(
    ::typeof(exp!),
    M::AbstractConnectionManifold,
    q,
    p,
    X,
    t,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(inverse_retract!),
    M::AbstractConnectionManifold,
    X,
    p,
    q,
    m::LogarithmicInverseRetraction,
)
    return Val(:parent)
end
function decorator_transparent_dispatch(
    ::typeof(retract!),
    M::AbstractConnectionManifold,
    q,
    p,
    X,
    m::ExponentialRetraction,
)
    return Val(:parent)
end
