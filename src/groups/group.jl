_dep_warn_group = """
GroupManifold functionality will move to its own package `LieGroups.jl`.
"""

@doc raw"""
    AbstractGroupOperation

Abstract type for smooth binary operations ``∘`` on elements of a Lie group ``\mathcal{G}``:
```math
∘ : \mathcal{G} × \mathcal{G} → \mathcal{G}
```
An operation can be either defined for a specific group manifold over
number system `𝔽` or in general, by defining for an operation `Op` the following methods:

    identity_element!(::AbstractDecoratorManifold, q, q)
    inv!(::AbstractDecoratorManifold, q, p)
    _compose!(::AbstractDecoratorManifold, x, p, q)

Note that a manifold is connected with an operation by wrapping it with a decorator,
[`AbstractDecoratorManifold`](@extref `ManifoldsBase.AbstractDecoratorManifold`)
using the [`IsGroupManifold`](@ref) to specify the operation.
For a concrete case the concrete wrapper [`GroupManifold`](@ref) can be used.
"""
abstract type AbstractGroupOperation end

"""
    abstract type AbstractGroupVectorRepresentation end

An abstract supertype for indicating representation of tangent vectors on a group manifold.
The most common representations are [`LeftInvariantRepresentation`](@ref),
[`TangentVectorRepresentation`](@ref) and [`HybridTangentRepresentation`](@ref).
"""
abstract type AbstractGroupVectorRepresentation end

"""
    TangentVectorRepresentation

Specify that tangent vectors in a group are stored in a non-invariant way, corresponding to
the storage implied by the underlying manifold.
"""
struct TangentVectorRepresentation <: AbstractGroupVectorRepresentation end

"""
    LeftInvariantRepresentation

Specify that tangent vectors in a group are stored in Lie algebra using left-invariant
representation.
"""
struct LeftInvariantRepresentation <: AbstractGroupVectorRepresentation end

"""
    IsGroupManifold{O<:AbstractGroupOperation} <: AbstractTrait

A trait to declare an [`AbstractManifold`](@extref `ManifoldsBase.AbstractManifold`)  as a manifold with group structure
with operation of type `O`.

Using this trait you can turn a manifold that you implement _implicitly_ into a Lie group.
If you wish to decorate an existing manifold with one (or different) [`AbstractGroupAction`](@ref)s,
see [`GroupManifold`](@ref).

# Constructor

    IsGroupManifold(op::AbstractGroupOperation, vectors::AbstractGroupVectorRepresentation)
"""
struct IsGroupManifold{O<:AbstractGroupOperation,VR<:AbstractGroupVectorRepresentation} <:
       AbstractTrait
    op::O
    vectors::VR
end

"""
    AbstractInvarianceTrait <: AbstractTrait

A common supertype for anz [`AbstractTrait`](@extref `ManifoldsBase.AbstractTrait`) related to metric invariance
"""
abstract type AbstractInvarianceTrait <: AbstractTrait end

"""
    HasLeftInvariantMetric <: AbstractInvarianceTrait

Specify that the default metric functions for the left-invariant metric on a [`GroupManifold`](@ref)
are to be used.
"""
struct HasLeftInvariantMetric <: AbstractInvarianceTrait end

direction_and_side(::HasLeftInvariantMetric) = LeftForwardAction()
direction_and_side(::Type{HasLeftInvariantMetric}) = LeftForwardAction()

"""
    HasRightInvariantMetric <: AbstractInvarianceTrait

Specify that the default metric functions for the right-invariant metric on a [`GroupManifold`](@ref)
are to be used.
"""
struct HasRightInvariantMetric <: AbstractInvarianceTrait end

direction_and_side(::HasRightInvariantMetric) = RightBackwardAction()
direction_and_side(::Type{HasRightInvariantMetric}) = RightBackwardAction()

"""
    HasBiinvariantMetric <: AbstractInvarianceTrait

Specify that the default metric functions for the bi-invariant metric on a [`GroupManifold`](@ref)
are to be used.
"""
struct HasBiinvariantMetric <: AbstractInvarianceTrait end
function parent_trait(::HasBiinvariantMetric)
    return ManifoldsBase.TraitList(HasLeftInvariantMetric(), HasRightInvariantMetric())
end

"""
    is_group_manifold(G::GroupManifold)
    is_group_manifold(G::AbstractManifold, o::AbstractGroupOperation)

returns whether an [`AbstractDecoratorManifold`](@extref `ManifoldsBase.AbstractDecoratorManifold`)
is a group manifold with [`AbstractGroupOperation`](@ref) `o`.
For a [`GroupManifold`](@ref) `G` this checks whether the right operations is stored within `G`.
"""
is_group_manifold(::AbstractManifold, ::AbstractGroupOperation) = false

@trait_function is_group_manifold(M::AbstractDecoratorManifold, op::AbstractGroupOperation)
function is_group_manifold(
    ::TraitList{<:IsGroupManifold{<:O}},
    ::AbstractDecoratorManifold,
    ::O,
) where {O<:AbstractGroupOperation}
    return true
end
@trait_function is_group_manifold(M::AbstractDecoratorManifold)
is_group_manifold(::AbstractManifold) = false
function is_group_manifold(
    t::TraitList{<:IsGroupManifold{<:AbstractGroupOperation}},
    M::AbstractDecoratorManifold,
)
    return is_group_manifold(M, t.head.op)
end

base_group(M::MetricManifold) = decorated_manifold(M)
base_group(M::ConnectionManifold) = decorated_manifold(M)
base_group(M::AbstractDecoratorManifold) = M

"""
    ActionDirection

Direction of action on a manifold, either [`LeftAction`](@ref) or [`RightAction`](@ref).
"""
abstract type ActionDirection end

@doc raw"""
    LeftAction()

Left action of a group on a manifold. For a forward action ``α: G×X → X`` it is characterized by
```math
α(g, α(h, x)) = α(gh, x)
```
for all ``g, h ∈ G`` and ``x ∈ X``.
"""
struct LeftAction <: ActionDirection end

"""
    RightAction()

Right action of a group on a manifold. For a forward action ``α: G×X → X`` it is characterized by
```math
α(g, α(h, x)) = α(hg, x)
```
for all ``g, h ∈ G`` and ``x ∈ X``.

Note that a right action may act from either left or right side in an expression.
"""
struct RightAction <: ActionDirection end

"""
    GroupActionSide

Side of action on a manifold, either [`LeftSide`](@ref) or [`RightSide`](@ref).
"""
abstract type GroupActionSide end

"""
    LeftSide()

An action of a group on a manifold that acts from the left side, i.e. ``α: G×X → X``.
"""
struct LeftSide <: GroupActionSide end

"""
    RightSide()

An action of a group on a manifold that acts from the right side, i.e. ``α: X×G → X``.
"""
struct RightSide <: GroupActionSide end

"""
    switch_direction(::ActionDirection)

Returns type of action between left and right.
This function does not affect side of action, see [`switch_side`](@ref).
"""
switch_direction(::ActionDirection)
switch_direction(::LeftAction) = RightAction()
switch_direction(::RightAction) = LeftAction()

"""
    switch_side(::GroupActionSide)

Returns side of action between left and right.
This function does not affect the action being left or right, see [`switch_direction`](@ref).
"""
switch_side(::GroupActionSide)
switch_side(::LeftSide) = RightSide()
switch_side(::RightSide) = LeftSide()

const ActionDirectionAndSide = Tuple{ActionDirection,GroupActionSide}

const LeftForwardAction = Tuple{LeftAction,LeftSide}
const LeftBackwardAction = Tuple{LeftAction,RightSide}
const RightForwardAction = Tuple{RightAction,LeftSide}
const RightBackwardAction = Tuple{RightAction,RightSide}

LeftForwardAction() = (LeftAction(), LeftSide())
LeftBackwardAction() = (LeftAction(), RightSide())
RightForwardAction() = (RightAction(), LeftSide())
RightBackwardAction() = (RightAction(), RightSide())

@doc raw"""
    Identity{O<:AbstractGroupOperation}

Represent the group identity element ``e ∈ \mathcal{G}`` on a Lie group ``\mathcal G``
with [`AbstractGroupOperation`](@ref) of type `O`.

Similar to the philosophy that points are agnostic of their group at hand, the identity
does not store the group `g` it belongs to. However it depends on the type of the [`AbstractGroupOperation`](@ref) used.

See also [`identity_element`](@ref) on how to obtain the corresponding [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`) or array representation.

# Constructors

    Identity(G::AbstractDecoratorManifold{𝔽})
    Identity(o::O)
    Identity(::Type{O})

create the identity of the corresponding subtype `O<:`[`AbstractGroupOperation`](@ref)
"""
struct Identity{O<:AbstractGroupOperation} end

@trait_function Identity(M::AbstractDecoratorManifold)
function Identity(
    ::TraitList{<:IsGroupManifold{O}},
    ::AbstractDecoratorManifold,
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`Identity` will move and keep its name.", :Idenity)
    return Identity{O}()
end
Identity(::O) where {O<:AbstractGroupOperation} = Identity(O)
Identity(::Type{O}) where {O<:AbstractGroupOperation} = Identity{O}()

# To ensure allocate_result_type works in general if identity appears in the tuple
number_eltype(::Identity) = Bool

@doc raw"""
    identity_element(G::AbstractDecoratorManifold)

Return a point representation of the [`Identity`](@ref) on the [`IsGroupManifold`](@ref) `G`.
By default this representation is the default array or number representation.
It should return the corresponding default representation of ``e`` as a point on `G` if
points are not represented by arrays.
"""
identity_element(G::AbstractDecoratorManifold)
@trait_function identity_element(G::AbstractDecoratorManifold)
function identity_element(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold)
    Base.depwarn(
        _dep_warn_group * "\n`identity_element` will move and keep its name.",
        :identity_element,
    )
    BG = base_group(G)
    q = allocate_result(BG, identity_element)
    return identity_element!(BG, q)
end

@trait_function identity_element!(G::AbstractDecoratorManifold, p)

function allocate_result(G::AbstractDecoratorManifold, f::typeof(identity_element))
    apf = allocation_promotion_function(G, f, ())
    rs = representation_size(G)
    return ManifoldsBase.allocate_result_array(G, f, apf(Float64), rs)
end

@doc raw"""
    identity_element(G::AbstractDecoratorManifold, p)

Return a point representation of the [`Identity`](@ref) on the [`IsGroupManifold`](@ref) `G`,
where `p` indicates the type to represent the identity.
"""
identity_element(G::AbstractDecoratorManifold, p)
@trait_function identity_element(G::AbstractDecoratorManifold, p)
function identity_element(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, p)
    BG = base_group(G)
    q = allocate_result(BG, identity_element, p)
    return identity_element!(BG, q)
end

Base.adjoint(e::Identity) = e

function check_size(
    ::TraitList{<:IsGroupManifold{O}},
    M::AbstractDecoratorManifold,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    return nothing
end
function check_size(::EmptyTrait, M::AbstractDecoratorManifold, e::Identity)
    return DomainError(0, "$M seems to not be a group manifold with $e.")
end

@doc raw"""
    is_identity(G::AbstractDecoratorManifold, q; kwargs)

Check whether `q` is the identity on the [`IsGroupManifold`](@ref) `G`, i.e. it is either
the [`Identity`](@ref)`{O}` with the corresponding [`AbstractGroupOperation`](@ref) `O`, or
(approximately) the correct point representation.
"""
is_identity(G::AbstractDecoratorManifold, q)
@trait_function is_identity(G::AbstractDecoratorManifold, q; kwargs...)
function is_identity(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    q;
    kwargs...,
)
    Base.depwarn(
        _dep_warn_group * "\n`is_identity` will move and keep its name.",
        :is_identity,
    )
    BG = base_group(G)
    return isapprox(BG, identity_element(BG), q; kwargs...)
end
function is_identity(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    ::Identity{O};
    kwargs...,
) where {O<:AbstractGroupOperation}
    return true
end
function is_identity(
    ::TraitList{<:IsGroupManifold},
    ::AbstractDecoratorManifold,
    ::Identity;
    kwargs...,
)
    return false
end

@inline function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p::Identity{O},
    q;
    kwargs...,
) where {O<:AbstractGroupOperation}
    return is_identity(G, q; kwargs...)
end
@inline function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p,
    q::Identity{O};
    kwargs...,
) where {O<:AbstractGroupOperation}
    BG = base_group(G)
    return is_identity(BG, p; kwargs...)
end
function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p::Identity{O},
    q::Identity{O};
    kwargs...,
) where {O<:AbstractGroupOperation}
    return true
end
function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p::Identity{O},
    q::Identity;
    kwargs...,
) where {O<:AbstractGroupOperation}
    return false
end
function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p::Identity,
    q::Identity{O};
    kwargs...,
) where {O<:AbstractGroupOperation}
    return false
end

@inline function isapprox(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p::Identity{O},
    X,
    Y;
    kwargs...,
) where {O<:AbstractGroupOperation}
    BG = base_group(G)
    return isapprox(BG, identity_element(BG), X, Y; kwargs...)
end
function isapprox(
    ::TraitList{<:IsGroupManifold},
    ::AbstractDecoratorManifold,
    ::Identity,
    ::Identity;
    kwargs...,
)
    return false
end

function Base.show(io::IO, ::Identity{O}) where {O<:AbstractGroupOperation}
    return print(io, "Identity($O)")
end

function is_point(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    e::Identity;
    error::Symbol=:none,
    kwargs...,
)
    ie = is_identity(G, e; kwargs...)
    if !ie
        s = "The provided identity is not a point on $G."
        (error === :error) && throw(DomainError(e, s))
        (error === :info) && @info s
        (error === :warn) && @warn s
    end
    return ie
end

function is_vector(
    t::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    e::Identity,
    X,
    cbp::Bool=true;
    error::Symbol=:none,
    kwargs...,
)
    if cbp
        # pass te down so this throws an error if error=:error
        # if error is not `:error` and is_point was false -> return false, otherwise continue
        (!is_point(G, e; error=error, kwargs...)) && return false
    end
    return is_vector(
        next_trait(t),
        G,
        identity_element(G),
        X,
        false;
        error=error,
        kwargs...,
    )
end

@doc raw"""
    adjoint_action(G::AbstractDecoratorManifold, p, X, dir=LeftAction())

Adjoint action of the element `p` of the Lie group `G` on the element `X`
of the corresponding Lie algebra.

If `dir` is `LeftAction()`, it is defined as the differential of the group automorphism ``Ψ_p(q) = pqp⁻¹`` at
the identity of `G`.

The formula reads
````math
\operatorname{Ad}_p(X) = dΨ_p(e)[X]
````
where ``e`` is the identity element of `G`.

If `dir` is `RightAction()`, then the formula is
````math
\operatorname{Ad}_p(X) = dΨ_{p^{-1}}(e)[X]
````

Note that the adjoint representation of a Lie group isn't generally faithful.
Notably the adjoint representation of SO(2) is trivial.
"""
adjoint_action(G::AbstractDecoratorManifold, p, X, dir)
@trait_function adjoint_action(G::AbstractDecoratorManifold, p, Xₑ, dir)
@trait_function adjoint_action!(G::AbstractDecoratorManifold, Y, p, Xₑ, dir)
function adjoint_action(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    Xₑ,
    dir,
)
    BG = base_group(G)
    Y = allocate_result(BG, adjoint_action, Xₑ, p)
    return adjoint_action!(BG, Y, p, Xₑ, dir)
end
function adjoint_action(::AbstractDecoratorManifold, ::Identity, Xₑ, ::LeftAction)
    return Xₑ
end
function adjoint_action(::AbstractDecoratorManifold, ::Identity, Xₑ, ::RightAction)
    return Xₑ
end
# backward compatibility
function adjoint_action(G::AbstractDecoratorManifold, p, X)
    return adjoint_action(G, p, X, LeftAction())
end
function adjoint_action!(G::AbstractDecoratorManifold, Y, p, X)
    return adjoint_action!(G, Y, p, X, LeftAction())
end
# fall back method: the right action is defined from the left action
function adjoint_action!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    Y,
    p,
    X,
    ::RightAction,
)
    BG = base_group(G)
    return adjoint_action!(BG, Y, inv(BG, p), X, LeftAction())
end

@doc raw"""
    adjoint_inv_diff(G::AbstractDecoratorManifold, p, X)

Compute the value of pullback of inverse ``p^{-1} ∈ \mathcal{G}`` of an element
``p ∈ \mathcal{G}`` at tangent vector `X` at ``p^{-1}``. The result is a tangent vector
at ``p``.
"""
adjoint_inv_diff(G::AbstractDecoratorManifold, p)

@trait_function adjoint_inv_diff(G::AbstractDecoratorManifold, p, X)
function adjoint_inv_diff(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    X,
)
    Y = allocate_result(G, inv_diff, X, p)
    return adjoint_inv_diff!(G, Y, p, X)
end

@trait_function adjoint_inv_diff!(G::AbstractDecoratorManifold, Y, p, X)

"""
    adjoint_matrix(G::AbstractManifold, p, B::AbstractBasis=DefaultOrthonormalBasis())

Compute the adjoint matrix related to conjugation of vectors by element `p` of Lie group `G`
for basis `B`. It is the matrix ``A`` such that for each element `X` of the Lie algebra
with coefficients ``c`` in basis `B`, ``Ac`` is the vector of coefficients of `X` conjugated
by `p` in basis `B`.
"""
function adjoint_matrix(G::AbstractManifold, p, B::AbstractBasis=DefaultOrthonormalBasis())
    J = allocate_jacobian(G, G, adjoint_matrix, p)
    return adjoint_matrix!(G, J, p, B)
end

function adjoint_matrix!(
    G::AbstractManifold,
    J,
    p,
    B::AbstractBasis=DefaultOrthonormalBasis(),
)
    Bb = get_basis(G, p, B)
    Vs = get_vectors(G, p, Bb)
    for i in eachindex(Vs)
        get_coordinates!(G, view(J, :, i), p, adjoint_action(G, p, Vs[i]), B)
    end
    return J
end

function ManifoldDiff.differential_exp_argument_lie_approx!(
    M::AbstractManifold,
    Z,
    p,
    X,
    Y;
    n=20,
)
    Base.depwarn(
        _dep_warn_group *
        "\n`differential_exp_argument_lie_approx` will be removed, cf `differential_exp_argument` instead.",
        :differential_exp_argument_lie_approx,
    )
    tmp = copy(M, p, Y)
    a = -1.0
    zero_vector!(M, Z, p)
    for k in 0:n
        a *= -1 // (k + 1)
        Z .+= a .* tmp
        if k < n
            copyto!(tmp, lie_bracket(M, X, tmp))
        end
    end
    q = exp(M, p, X)
    translate_diff!(M, Z, q, Identity(M), Z)
    return Z
end

@doc raw"""
    inv(G::AbstractDecoratorManifold, p)

Inverse ``p^{-1} ∈ \mathcal{G}`` of an element ``p ∈ \mathcal{G}``, such that
``p \circ p^{-1} = p^{-1} \circ p = e ∈ \mathcal{G}``, where ``e`` is the [`Identity`](@ref)
element of ``\mathcal{G}``.
"""
inv(::AbstractDecoratorManifold, ::Any...)
@trait_function Base.inv(G::AbstractDecoratorManifold, p)
function Base.inv(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, p)
    Base.depwarn(_dep_warn_group * "\n`inv` is moved and keeps its name.", :inv)
    q = allocate_result(G, inv, p)
    BG = base_group(G)
    return inv!(BG, q, p)
end

function Base.inv(
    ::TraitList{<:IsGroupManifold{O}},
    ::AbstractDecoratorManifold,
    e::Identity{O},
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`inv` is moved and keeps its name.", :inv)
    return e
end

@trait_function inv!(G::AbstractDecoratorManifold, q, p)

function inv!(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    q,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    BG = base_group(G)
    return identity_element!(BG, q)
end
function inv!(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    ::Identity{O},
    e::Identity{O},
) where {O<:AbstractGroupOperation}
    return e
end

@doc raw"""
    inv_diff(G::AbstractDecoratorManifold, p, X)

Compute the value of differential of inverse ``p^{-1} ∈ \mathcal{G}`` of an element
``p ∈ \mathcal{G}`` at tangent vector `X` at `p`. The result is a tangent vector at ``p^{-1}``.

*Note*: the default implementation of `inv_diff` and `inv_diff!`
assumes that the tangent vector ``X`` is stored at
the point ``p ∈ \mathcal{G}`` as the vector ``Y ∈ \mathfrak{g}``
 where ``X = pY``.
"""
inv_diff(G::AbstractDecoratorManifold, p)

@trait_function inv_diff(G::AbstractDecoratorManifold, p, X)
function inv_diff(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, p, X)
    Base.depwarn(_dep_warn_group * "\n`inv_diff` is renamed to `diff_inv`.", :inv_diff)
    return -adjoint_action(base_group(G), p, X)
end

@trait_function inv_diff!(G::AbstractDecoratorManifold, Y, p, X)

function inv_diff!(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, Y, p, X)
    adjoint_action!(G, Y, p, X)
    Y .*= -1
    return Y
end

function Base.copyto!(
    ::TraitList{<:IsGroupManifold{O}},
    ::AbstractDecoratorManifold,
    e::Identity{O},
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    return e
end
function Base.copyto!(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    p,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    BG = base_group(G)
    return identity_element!(BG, p)
end

@doc raw"""
    compose(G::AbstractDecoratorManifold, p, q)

Compose elements ``p,q ∈ \mathcal{G}`` using the group operation ``p \circ q``.

For implementing composition on a new group manifold, please overload `_compose`
instead so that methods with [`Identity`](@ref) arguments are not ambiguous.
"""
compose(::AbstractDecoratorManifold, ::Any...)

@trait_function compose(G::AbstractDecoratorManifold, p, q)
function compose(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, p, q)
    Base.depwarn(_dep_warn_group * "\n`compose` is moved and keeps its name.", :compose)
    return _compose(base_group(G), p, q)
end
function compose(
    ::AbstractDecoratorManifold,
    ::Identity{O},
    p,
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`compose` is moved and keeps its name.", :compose)
    return p
end
function compose(
    ::AbstractDecoratorManifold,
    p,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`compose` is moved and keeps its name.", :compose)
    return p
end
function compose(
    ::AbstractDecoratorManifold,
    e::Identity{O},
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`compose` is moved and keeps its name.", :compose)
    return e
end

function _compose(G::AbstractDecoratorManifold, p, q)
    x = allocate_result(G, compose, p, q)
    return _compose!(G, x, p, q)
end

@trait_function compose!(M::AbstractDecoratorManifold, x, p, q)

function compose!(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, x, q, p)
    return _compose!(base_group(G), x, q, p)
end
function compose!(
    G::AbstractDecoratorManifold,
    q,
    p,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    return copyto!(G, q, p)
end
function compose!(
    G::AbstractDecoratorManifold,
    q,
    ::Identity{O},
    p,
) where {O<:AbstractGroupOperation}
    return copyto!(G, q, p)
end
function compose!(
    G::AbstractDecoratorManifold,
    q,
    ::Identity{O},
    e::Identity{O},
) where {O<:AbstractGroupOperation}
    return identity_element!(G, q)
end
function compose!(
    ::AbstractDecoratorManifold,
    e::Identity{O},
    ::Identity{O},
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    return e
end

Base.transpose(e::Identity) = e

@trait_function hat(M::AbstractDecoratorManifold, e::Identity, X)
@trait_function hat!(M::AbstractDecoratorManifold, Y, e::Identity, X)

@doc raw"""
    hat(M::AbstractDecoratorManifold{𝔽,O}, ::Identity{O}, Xⁱ) where {𝔽,O<:AbstractGroupOperation}

Given a basis ``e_i`` on the tangent space at a the [`Identity`](@ref) and tangent
component vector ``X^i``, compute the equivalent vector representation
``X=X^i e_i**, where Einstein summation notation is used:

````math
∧ : X^i ↦ X^i e_i
````

For array manifolds, this converts a vector representation of the tangent
vector to an array representation. The [`vee`](@ref) map is the `hat` map's
inverse.
"""
function hat(
    ::TraitList{<:IsGroupManifold{O}},
    M::AbstractDecoratorManifold,
    ::Identity{O},
    X,
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`hat` is moved and keeps its name.", :hat)
    return get_vector_lie(M, X, VeeOrthogonalBasis())
end
function hat!(
    ::TraitList{<:IsGroupManifold{O}},
    M::AbstractDecoratorManifold,
    Y,
    ::Identity{O},
    X,
) where {O<:AbstractGroupOperation}
    return get_vector_lie!(M, Y, X, VeeOrthogonalBasis())
end
function hat(M::AbstractManifold, e::Identity, ::Any)
    return throw(ErrorException("On $M there exists no identity $e"))
end
function hat!(M::AbstractManifold, c, e::Identity, X)
    return throw(ErrorException("On $M there exists no identity $e"))
end

@trait_function vee(M::AbstractDecoratorManifold, e::Identity, X)
@trait_function vee!(M::AbstractDecoratorManifold, Y, e::Identity, X)

@doc raw"""
    vee(M::AbstractManifold, p, X)

Given a basis ``e_i`` on the tangent space at a point `p` and tangent
vector `X`, compute the vector components ``X^i``, such that ``X = X^i e_i``, where
Einstein summation notation is used:

````math
\vee : X^i e_i ↦ X^i
````

For array manifolds, this converts an array representation of the tangent
vector to a vector representation. The [`hat`](@ref) map is the `vee` map's
inverse.
"""
function vee(
    ::TraitList{<:IsGroupManifold{O}},
    M::AbstractDecoratorManifold,
    ::Identity{O},
    X,
) where {O<:AbstractGroupOperation}
    Base.depwarn(_dep_warn_group * "\n`vee` is moved and keeps its name.", :vee)
    return get_coordinates_lie(M, X, VeeOrthogonalBasis())
end
function vee!(
    ::TraitList{<:IsGroupManifold{O}},
    M::AbstractDecoratorManifold,
    Y,
    ::Identity{O},
    X,
) where {O<:AbstractGroupOperation}
    return get_coordinates_lie!(M, Y, X, VeeOrthogonalBasis())
end
function vee(M::AbstractManifold, e::Identity, X)
    return throw(ErrorException("On $M there exists no identity $e"))
end
function vee!(M::AbstractManifold, c, e::Identity, X)
    return throw(ErrorException("On $M there exists no identity $e"))
end

"""
    lie_bracket(G::AbstractDecoratorManifold, X, Y)

Lie bracket between elements `X` and `Y` of the Lie algebra corresponding to
the Lie group `G`, cf. [`IsGroupManifold`](@ref).

This can be used to compute the adjoint representation of a Lie algebra.
Note that this representation isn't generally faithful. Notably the adjoint
representation of 𝔰𝔬(2) is trivial.
"""
lie_bracket(G::AbstractDecoratorManifold, X, Y)
@trait_function lie_bracket(M::AbstractDecoratorManifold, X, Y)

@trait_function lie_bracket!(M::AbstractDecoratorManifold, Z, X, Y)

_action_order(BG::AbstractDecoratorManifold, p, q, ::LeftForwardAction) = (p, q)
_action_order(BG::AbstractDecoratorManifold, p, q, ::LeftBackwardAction) = (q, inv(BG, p))
_action_order(BG::AbstractDecoratorManifold, p, q, ::RightForwardAction) = (inv(BG, p), q)
_action_order(BG::AbstractDecoratorManifold, p, q, ::RightBackwardAction) = (q, p)

@doc raw"""
    translate(G::AbstractDecoratorManifold, p, q, conv::ActionDirectionAndSide=LeftForwardAction()])

Translate group element ``q`` by ``p`` with the translation ``τ_p`` with the specified
`conv`ention, either
- left forward ``τ_p(q) = p \circ q``
- left backward ``τ_p(q) = q \circ p^{-1}``
- right backward ``τ_p(q) = q \circ p``
- right forward ``τ_p(q) = p^{-1} \circ q``
"""
translate(::AbstractDecoratorManifold, ::Any...)
@trait_function translate(
    G::AbstractDecoratorManifold,
    p,
    q,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function translate(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    q,
    conv::ActionDirectionAndSide,
)
    Base.depwarn(
        _dep_warn_group * "\n`translate` is discontinued – use `compose` instead.",
        :translate,
    )
    BG = base_group(G)
    return compose(BG, _action_order(BG, p, q, conv)...)
end

@trait_function translate!(
    G::AbstractDecoratorManifold,
    X,
    p,
    q,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function translate!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    X,
    p,
    q,
    conv::ActionDirectionAndSide,
)
    BG = base_group(G)
    return compose!(BG, X, _action_order(BG, p, q, conv)...)
end

@doc raw"""
    inverse_translate(G::AbstractDecoratorManifold, p, q, conv::ActionDirectionAndSide=LeftForwardAction())

Inverse translate group element ``q`` by ``p`` with the translation ``τ_p^{-1}``
with the specified `conv`ention, either left forward (``L_p^{-1}``), left backward
(``R'_p^{-1}``), right backward (``R_p^{-1}``) or right forward (``L'_p^{-1}``), defined as
```math
\begin{aligned}
L_p^{-1} &: q ↦ p^{-1} \circ q\\
L'_p^{-1} &: q ↦ p \circ q\\
R_p^{-1} &: q ↦ q \circ p^{-1}\\
R'_p^{-1} &: q ↦ q \circ p.
\end{aligned}
"""
inverse_translate(::AbstractDecoratorManifold, ::Any...)
@trait_function inverse_translate(
    G::AbstractDecoratorManifold,
    p,
    q,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function inverse_translate(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    q,
    conv::ActionDirectionAndSide,
)
    Base.depwarn(
        _dep_warn_group *
        "\n`inverse_translate` is discontinued – use `compose` with the `inv` instead.",
        :inverse_translate,
    )
    BG = base_group(G)
    return translate(BG, inv(BG, p), q, conv)
end

@trait_function inverse_translate!(
    G::AbstractDecoratorManifold,
    X,
    p,
    q,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function inverse_translate!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    X,
    p,
    q,
    conv::ActionDirectionAndSide,
)
    BG = base_group(G)
    return translate!(BG, X, inv(BG, p), q, conv)
end

@doc raw"""
    translate_diff(G::AbstractDecoratorManifold, p, q, X, conv::ActionDirectionAndSide=LeftForwardAction())

For group elements ``p, q ∈ \mathcal{G}`` and tangent vector ``X ∈ T_q \mathcal{G}``, compute
the action of the differential of the translation ``τ_p`` by ``p`` on ``X``, with the specified
left or right `conv`ention. The differential transports vectors:
```math
(\mathrm{d}τ_p)_q : T_q \mathcal{G} → T_{τ_p q} \mathcal{G}\\
```

*Note*: the default implementation of `translate_diff` and `translate_diff!`
assumes that a tangent vector ``X`` at
a point ``q ∈ \mathcal{G}`` is stored  as the vector ``Y ∈ \mathfrak{g}``
 where ``X = qY``.
The implementation at `q = Identity` is independent of the storage choice.
"""
translate_diff(::AbstractDecoratorManifold, ::Any...)
@trait_function translate_diff(
    G::AbstractDecoratorManifold,
    p,
    q,
    X,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function translate_diff(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    Base.depwarn(
        _dep_warn_group *
        "\n`translate_diff` is discontinued – use `diff_right_compose` and `diff_left_compose`, respectively instead.",
        :translate_diff,
    )
    Y = allocate_result(G, translate_diff, X, p, q)
    BG = base_group(G)
    translate_diff!(BG, Y, p, q, X, conv)
    return Y
end
@trait_function translate_diff!(
    G::AbstractDecoratorManifold,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)

function translate_diff!(
    ::TraitList{<:IsGroupManifold{<:AbstractGroupOperation,LeftInvariantRepresentation}},
    G::AbstractDecoratorManifold,
    Y,
    ::Any,
    ::Any,
    X,
    ::LeftForwardAction,
)
    return copyto!(G, Y, X)
end
function translate_diff!(
    ::TraitList{<:IsGroupManifold{<:AbstractGroupOperation,LeftInvariantRepresentation}},
    G::AbstractDecoratorManifold,
    Y,
    ::Any,
    ::Any,
    X,
    ::RightForwardAction,
)
    return copyto!(G, Y, X)
end
function translate_diff!(
    ::TraitList{<:IsGroupManifold{<:AbstractGroupOperation,LeftInvariantRepresentation}},
    G::AbstractDecoratorManifold,
    Y,
    p,
    ::Any,
    X,
    ::LeftBackwardAction,
)
    return adjoint_action!(G, Y, p, X, LeftAction())
end
function translate_diff!(
    ::TraitList{<:IsGroupManifold{<:AbstractGroupOperation,LeftInvariantRepresentation}},
    G::AbstractDecoratorManifold,
    Y,
    p,
    ::Any,
    X,
    ::RightBackwardAction,
)
    return adjoint_action!(G, Y, p, X, RightAction())
end

translate_diff(::AbstractDecoratorManifold, ::Identity, q, X, ::LeftForwardAction) = X
translate_diff(::AbstractDecoratorManifold, ::Identity, q, X, ::RightForwardAction) = X
translate_diff(::AbstractDecoratorManifold, ::Identity, q, X, ::LeftBackwardAction) = X
translate_diff(::AbstractDecoratorManifold, ::Identity, q, X, ::RightBackwardAction) = X

@doc raw"""
    inverse_translate_diff(G::AbstractDecoratorManifold, p, q, X, conv::ActionDirectionAndSide=LeftForwardAction())

For group elements ``p, q ∈ \mathcal{G}`` and tangent vector ``X ∈ T_q \mathcal{G}``, compute
the action on ``X`` of the differential of the inverse translation ``τ_p`` by ``p``, with the
specified left or right `conv`ention. The differential transports vectors:
```math
(\mathrm{d}τ_p^{-1})_q : T_q \mathcal{G} → T_{τ_p^{-1} q} \mathcal{G}\\
```
"""
inverse_translate_diff(::AbstractDecoratorManifold, ::Any...)
@trait_function inverse_translate_diff(
    G::AbstractDecoratorManifold,
    p,
    q,
    X,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function inverse_translate_diff(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    Base.depwarn(
        _dep_warn_group *
        "\n`inverse_translate_diff` is discontinued – use `diff_right_compose` and `diff_left_compose`, respectively, together with `inv` instead.",
        :inverse_translate_diff,
    )
    BG = base_group(G)
    return translate_diff(BG, inv(BG, p), q, X, conv)
end

@trait_function inverse_translate_diff!(
    G::AbstractDecoratorManifold,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide=LeftForwardAction(),
)
function inverse_translate_diff!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    BG = base_group(G)
    return translate_diff!(BG, Y, inv(BG, p), q, X, conv)
end

"""
    log_inv(G::AbstractManifold, p, q)

Compute logarithmic map on a Lie group `G` invariant to group operation. For groups with a
bi-invariant metric or a Cartan-Schouten connection, this is the same as `log` but for
other groups it may differ.
"""
function log_inv(G::AbstractManifold, p, q)
    Base.depwarn(
        _dep_warn_group *
        "\n`log_inv` is discontinued – use `log` on the Lie group instead.",
        :log_inv,
    )
    BG = base_group(G)
    return log_lie(BG, compose(BG, inv(BG, p), q))
end
function log_inv!(G::AbstractManifold, X, p, q)
    x = allocate_result(G, inv, p)
    BG = base_group(G)
    inv!(BG, x, p)
    compose!(BG, x, x, q)
    log_lie!(BG, X, x)
    return X
end

"""
    exp_inv(G::AbstractManifold, p, X, t::Number=1)

Compute exponential map on a Lie group `G` invariant to group operation. For groups with a
bi-invariant metric or a Cartan-Schouten connection, this is the same as `exp` but for
other groups it may differ.
"""
function exp_inv(G::AbstractManifold, p, X, t::Number=1)
    Base.depwarn(
        _dep_warn_group *
        "\n`exp_inv` is discontinued – use `exp` on the Lie group instead.",
        :log_inv,
    )
    BG = base_group(G)
    return compose(BG, p, exp_lie(BG, t * X))
end
function exp_inv!(G::AbstractManifold, q, p, X)
    BG = base_group(G)
    exp_lie!(BG, q, X)
    compose!(BG, q, p, q)
    return q
end

@doc raw"""
    exp_lie(G, X)
    exp_lie!(G, q, X)

Compute the group exponential of the Lie algebra element `X`. It is equivalent to the
exponential map defined by the [`CartanSchoutenMinus`](@ref) connection.

Given an element ``X ∈ 𝔤 = T_e \mathcal{G}``, where ``e`` is the [`Identity`](@ref) element of
the group ``\mathcal{G}``, and ``𝔤`` is its Lie algebra, the group exponential is the map

````math
\exp : 𝔤 → \mathcal{G},
````
such that for ``t,s ∈ ℝ``, ``γ(t) = \exp (t X)`` defines a one-parameter subgroup with the
following properties. Note that one-parameter subgroups are commutative (see [Suhubi:2013](@cite),
section 3.5), even if the Lie group itself is not commutative.

````math
\begin{aligned}
γ(t) &= γ(-t)^{-1}\\
γ(t + s) &= γ(t) \circ γ(s) = γ(s) \circ γ(t)\\
γ(0) &= e\\
\lim_{t → 0} \frac{d}{dt} γ(t) &= X.
\end{aligned}
````

!!! note
    In general, the group exponential map is distinct from the Riemannian exponential map
    [`exp`](@ref).

For example for the [`MultiplicationOperation`](@ref) and either `Number` or `AbstractMatrix`
the Lie exponential is the numeric/matrix exponential.

````math
\exp X = \operatorname{Exp} X = \sum_{n=0}^∞ \frac{1}{n!} X^n.
````

Since this function also depends on the group operation, make sure to implement
the corresponding trait version `exp_lie(::TraitList{<:IsGroupManifold}, G, X)`.
"""
exp_lie(G::AbstractManifold, X)
@trait_function exp_lie(M::AbstractDecoratorManifold, X)
function exp_lie(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, X)
    Base.depwarn(
        _dep_warn_group *
        "\n`exp_lie`` is discontinued – use `exp(G, p, X)` on the Lie group at `p` the `Identity` instead.",
        :exp_lie,
    )
    BG = base_group(G)
    q = allocate_result(BG, exp_lie, X)
    return exp_lie!(BG, q, X)
end

@trait_function exp_lie!(M::AbstractDecoratorManifold, q, X)

@doc raw"""
    log_lie(G, q)
    log_lie!(G, X, q)

Compute the Lie group logarithm of the Lie group element `q`. It is equivalent to the
logarithmic map defined by the [`CartanSchoutenMinus`](@ref) connection.

Given an element ``q ∈ \mathcal{G}``, compute the right inverse of the group exponential map
[`exp_lie`](@ref), that is, the element ``\log q = X ∈ 𝔤 = T_e \mathcal{G}``, such that
``q = \exp X``

!!! note
    In general, the group logarithm map is distinct from the Riemannian logarithm map
    [`log`](@ref).

    For matrix Lie groups this is equal to the (matrix) logarithm:

````math
\log q = \operatorname{Log} q = \sum_{n=1}^∞ \frac{(-1)^{n+1}}{n} (q - e)^n,
````

where ``e`` here is the [`Identity`](@ref) element, that is, ``1`` for numeric ``q`` or the
identity matrix ``I_m`` for matrix ``q ∈ ℝ^{m×m}``.

Since this function also depends on the group operation, make sure to implement
either
* `_log_lie(G, q)` and `_log_lie!(G, X, q)` for the points not being the [`Identity`](@ref)
* the trait version `log_lie(::TraitList{<:IsGroupManifold}, G, e)`, `log_lie(::TraitList{<:IsGroupManifold}, G, X, e)` for own implementations of the identity case.
"""
log_lie(::AbstractDecoratorManifold, q)
@trait_function log_lie(G::AbstractDecoratorManifold, q)
function log_lie(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, q)
    Base.depwarn(
        _dep_warn_group *
        "\n`log_lie`` is discontinued – use `log(G, p, q)` on the Lie group at `p` the `Identity` instead.",
        :exp_lie,
    )
    BG = base_group(G)
    return _log_lie(BG, q)
end
function log_lie(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    BG = base_group(G)
    return zero_vector(BG, identity_element(BG))
end
# though identity was taken care of – as usual restart decorator dispatch
function _log_lie(G::AbstractDecoratorManifold, q)
    BG = base_group(G)
    X = allocate_result(BG, log_lie, q)
    return log_lie!(BG, X, q)
end

@trait_function log_lie!(G::AbstractDecoratorManifold, X, q)
function log_lie!(::TraitList{<:IsGroupManifold}, G::AbstractDecoratorManifold, X, q)
    BG = base_group(G)
    return _log_lie!(BG, X, q)
end
function log_lie!(
    ::TraitList{<:IsGroupManifold{O}},
    G::AbstractDecoratorManifold,
    X,
    ::Identity{O},
) where {O<:AbstractGroupOperation}
    return zero_vector!(G, X, identity_element(G))
end

"""
    GroupExponentialRetraction{D<:ActionDirectionAndSide} <: AbstractRetractionMethod

Retraction using the group exponential [`exp_lie`](@ref) "translated" to any point on the
manifold.

For more details, see
[`retract`](@ref retract(::GroupManifold, p, X, ::GroupExponentialRetraction)).

# Constructor

    GroupExponentialRetraction(conv::ActionDirectionAndSide = LeftAction())
"""
struct GroupExponentialRetraction{D<:ActionDirectionAndSide} <: AbstractRetractionMethod end

function GroupExponentialRetraction(conv::ActionDirectionAndSide=LeftForwardAction())
    return GroupExponentialRetraction{typeof(conv)}()
end

"""
    GroupLogarithmicInverseRetraction{D<:ActionDirectionAndSide} <: AbstractInverseRetractionMethod

Retraction using the group logarithm [`log_lie`](@ref) "translated" to any point on the
manifold.

For more details, see
[`inverse_retract`](@ref inverse_retract(::GroupManifold, p, q ::GroupLogarithmicInverseRetraction)).

# Constructor

    GroupLogarithmicInverseRetraction(conv::ActionDirectionAndSide = LeftForwardAction())
"""
struct GroupLogarithmicInverseRetraction{D<:ActionDirectionAndSide} <:
       AbstractInverseRetractionMethod end

function GroupLogarithmicInverseRetraction(conv::ActionDirectionAndSide=LeftForwardAction())
    return GroupLogarithmicInverseRetraction{typeof(conv)}()
end

direction_and_side(::GroupExponentialRetraction{D}) where {D} = D()
direction_and_side(::GroupLogarithmicInverseRetraction{D}) where {D} = D()

@doc raw"""
    retract(
        G::AbstractDecoratorManifold,
        p,
        X,
        method::GroupExponentialRetraction,
    )

Compute the retraction using the group exponential [`exp_lie`](@ref) "translated" to any
point on the manifold.
With a group translation ([`translate`](@ref)) ``τ_p`` in a specified direction, the
retraction is

````math
\operatorname{retr}_p = τ_p \circ \exp \circ (\mathrm{d}τ_p^{-1})_p,
````

where ``\exp`` is the group exponential ([`exp_lie`](@ref)), and ``(\mathrm{d}τ_p^{-1})_p`` is
the action of the differential of inverse translation ``τ_p^{-1}`` evaluated at ``p`` (see
[`inverse_translate_diff`](@ref)).
"""
function retract(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    X,
    method::GroupExponentialRetraction,
)
    conv = direction_and_side(method)
    Xₑ = inverse_translate_diff(G, p, p, X, conv)
    pinvq = exp_lie(G, Xₑ)
    q = translate(G, p, pinvq, conv)
    return q
end
function retract(
    tl::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    X,
    t::Number,
    method::GroupExponentialRetraction,
)
    return retract(tl, G, p, t * X, method)
end

function retract!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    q,
    p,
    X,
    method::GroupExponentialRetraction,
)
    conv = direction_and_side(method)
    Xₑ = inverse_translate_diff(G, p, p, X, conv)
    pinvq = exp_lie(G, Xₑ)
    return translate!(G, q, p, pinvq, conv)
end
function retract!(
    tl::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    q,
    p,
    X,
    t::Number,
    method::GroupExponentialRetraction,
)
    return retract!(tl, G, q, p, t * X, method)
end

@doc raw"""
    inverse_retract(
        G::AbstractDecoratorManifold,
        p,
        X,
        method::GroupLogarithmicInverseRetraction,
    )

Compute the inverse retraction using the group logarithm [`log_lie`](@ref) "translated"
to any point on the manifold.
With a group translation ([`translate`](@ref)) ``τ_p`` in a specified direction, the
retraction is

````math
\operatorname{retr}_p^{-1} = (\mathrm{d}τ_p)_e \circ \log \circ τ_p^{-1},
````

where ``\log`` is the group logarithm ([`log_lie`](@ref)), and ``(\mathrm{d}τ_p)_e`` is the
action of the differential of translation ``τ_p`` evaluated at the identity element ``e``
(see [`translate_diff`](@ref)).
"""
function inverse_retract(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    p,
    q,
    method::GroupLogarithmicInverseRetraction,
)
    conv = direction_and_side(method)
    pinvq = inverse_translate(G, p, q, conv)
    Xₑ = log_lie(G, pinvq)
    return translate_diff(G, p, Identity(G), Xₑ, conv)
end

function inverse_retract!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractDecoratorManifold,
    X,
    p,
    q,
    method::GroupLogarithmicInverseRetraction,
)
    conv = direction_and_side(method)
    pinvq = inverse_translate(G, p, q, conv)
    Xₑ = log_lie(G, pinvq)
    return translate_diff!(G, X, p, Identity(G), Xₑ, conv)
end

@trait_function get_vector_lie(G::AbstractManifold, X, B::AbstractBasis)
@trait_function get_vector_lie!(G::AbstractManifold, Y, X, B::AbstractBasis)

@doc raw"""
    get_vector_lie(G::AbstractDecoratorManifold, a, B::AbstractBasis)

Reconstruct a tangent vector from the Lie algebra of `G` from coordinates `a` of a basis `B`.
This is similar to calling [`get_vector`](@ref) at the `p=`[`Identity`](@ref)`(G)`.
"""
function get_vector_lie(
    ::TraitList{<:IsGroupManifold},
    G::AbstractManifold,
    X,
    B::AbstractBasis,
)
    return get_vector(base_manifold(G), identity_element(G), X, B)
end
function get_vector_lie!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractManifold,
    Y,
    X,
    B::AbstractBasis,
)
    return get_vector!(base_manifold(G), Y, identity_element(G), X, B)
end

@trait_function get_coordinates_lie(G::AbstractManifold, X, B::AbstractBasis)
@trait_function get_coordinates_lie!(G::AbstractManifold, a, X, B::AbstractBasis)

@doc raw"""
    get_coordinates_lie(G::AbstractManifold, X, B::AbstractBasis)

Get the coordinates of an element `X` from the Lie algebra og `G` with respect to a basis `B`.
This is similar to calling [`get_coordinates`](@ref) at the `p=`[`Identity`](@ref)`(G)`.
"""
function get_coordinates_lie(
    ::TraitList{<:IsGroupManifold},
    G::AbstractManifold,
    X,
    B::AbstractBasis,
)
    return get_coordinates(base_manifold(G), identity_element(G), X, B)
end
function get_coordinates_lie!(
    ::TraitList{<:IsGroupManifold},
    G::AbstractManifold,
    a,
    X,
    B::AbstractBasis,
)
    return get_coordinates!(base_manifold(G), a, identity_element(G), X, B)
end

function ManifoldsBase._pick_basic_allocation_argument(
    ::AbstractManifold,
    f,
    ::Identity,
    x...,
)
    return x[1]
end
