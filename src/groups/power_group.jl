
const PowerGroup{𝔽,M<:AbstractManifold{𝔽},TPR<:AbstractPowerRepresentation} =
    GroupManifold{𝔽,<:AbstractPowerManifold{𝔽,M,TPR},ProductOperation}

"""
    PowerGroupNested

Alias to [`PowerGroup`](@ref) with [`NestedPowerRepresentation`](@extref `ManifoldsBase.NestedPowerRepresentation`)
representation.
"""
const PowerGroupNested{𝔽,M<:AbstractManifold{𝔽}} = PowerGroup{𝔽,M,NestedPowerRepresentation}

"""
    PowerGroupNestedReplacing

Alias to [`PowerGroup`](@ref) with [`NestedReplacingPowerRepresentation`](@extref `ManifoldsBase.NestedReplacingPowerRepresentation`)
representation.
"""
const PowerGroupNestedReplacing{𝔽,M<:AbstractManifold{𝔽}} =
    PowerGroup{𝔽,M,NestedReplacingPowerRepresentation}

"""
    PowerGroup{𝔽,T} <: GroupManifold{𝔽,<:AbstractPowerManifold{𝔽,M,RPT},ProductOperation}

Decorate a power manifold with a [`ProductOperation`](@ref).

Constituent manifold of the power manifold must also have a [`IsGroupManifold`](@ref) or
a decorated instance of one. This type is mostly useful for equipping the direct product of
group manifolds with an [`Identity`](@ref) element.

# Constructor

    PowerGroup(manifold::AbstractPowerManifold)
"""
function PowerGroup(manifold::AbstractPowerManifold)
    Base.depwarn(
        _dep_warn_group *
        "\n `PowerGroup` will be named to `PowerLiegroup` and its argument will be a Lie group instead of a manifold.",
        :PowerGroup,
    )
    if !is_group_manifold(manifold.manifold)
        error("All powered manifold must be or decorate a group.")
    end
    op = ProductOperation()
    return GroupManifold(manifold, op, vector_representation(manifold.manifold))
end

function ManifoldsBase._access_nested(
    M::PowerManifold,
    ::Identity{ProductOperation},
    i::Tuple,
)
    return Identity(M.manifold)
end

# lower level methods are added instead of top level ones to not have to deal
# with `Identity` disambiguation

_compose!(G::PowerGroup, x, p, q) = _compose!(G.manifold, x, p, q)
function _compose!(M::AbstractPowerManifold, x, p, q)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        compose!(
            N,
            _write(M, rep_size, x, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
        )
    end
    return x
end
function _compose!(M::PowerManifoldNestedReplacing, x, p, q)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        x[i...] = compose(N, _read(M, rep_size, p, i), _read(M, rep_size, q, i))
    end
    return x
end

@inline function active_traits(f, M::PowerGroup, args...)
    if is_metric_function(f)
        #pass to manifold by default - but keep Group Decorator for the retraction
        return merge_traits(IsGroupManifold(M.op, M.vectors), IsExplicitDecorator())
    else
        return merge_traits(
            IsGroupManifold(M.op, M.vectors),
            active_traits(f, M.manifold, args...),
            IsExplicitDecorator(),
        )
    end
end

function adjoint_inv_diff!(G::PowerGroup, Y, p, X)
    GM = G.manifold
    rep_size = representation_size(GM.manifold)
    for i in get_iterator(GM)
        adjoint_inv_diff!(
            GM.manifold,
            _write(GM, rep_size, Y, i),
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, X, i),
        )
    end
    return Y
end
function adjoint_inv_diff!(G::PowerGroupNestedReplacing, Y, p, X)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        Y[i...] = adjoint_inv_diff(N, _read(GM, rep_size, p, i), _read(GM, rep_size, X, i))
    end
    return Y
end

function identity_element!(G::PowerGroup, p)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        identity_element!(N, _write(GM, rep_size, p, i))
    end
    return p
end
function identity_element!(G::PowerGroupNestedReplacing, p)
    GM = G.manifold
    N = GM.manifold
    for i in get_iterator(GM)
        p[i...] = identity_element(N)
    end
    return p
end

function is_identity(G::PowerGroup, p::Identity{ProductOperation}; kwargs...)
    return true
end
function is_identity(G::PowerGroup, p::Identity; kwargs...)
    return false
end
function is_identity(G::PowerGroup, p; kwargs...)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        if !is_identity(N, _read(GM, rep_size, p, i))
            return false
        end
    end
    return true
end

inv!(G::PowerGroup, q, ::Identity{ProductOperation}) = identity_element!(G, q)
function inv!(G::PowerGroupNestedReplacing, q, ::Identity{ProductOperation})
    return identity_element!(G, q)
end
function inv!(G::PowerGroup, q, p)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        inv!(N, _write(GM, rep_size, q, i), _read(GM, rep_size, p, i))
    end
    return q
end
function inv!(G::PowerGroupNestedReplacing, q, p)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        q[i...] = inv(N, _read(GM, rep_size, p, i))
    end
    return q
end
inv!(::PowerGroup, q::Identity{ProductOperation}, ::Identity{ProductOperation}) = q
function inv!(
    ::PowerGroupNestedReplacing,
    q::Identity{ProductOperation},
    ::Identity{ProductOperation},
)
    return q
end

function inv_diff!(G::PowerGroup, Y, p, X)
    GM = G.manifold
    rep_size = representation_size(GM.manifold)
    for i in get_iterator(GM)
        inv_diff!(
            GM.manifold,
            _write(GM, rep_size, Y, i),
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, X, i),
        )
    end
    return Y
end
function inv_diff!(G::PowerGroupNestedReplacing, Y, p, X)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        Y[i...] = inv_diff(N, _read(GM, rep_size, p, i), _read(GM, rep_size, X, i))
    end
    return Y
end

function translate!(G::PowerGroup, x, p, q, conv::ActionDirectionAndSide)
    return translate!(G.manifold, x, p, q, conv)
end
function translate!(M::AbstractPowerManifold, x, p, q, conv::ActionDirectionAndSide)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        translate!(
            N,
            _write(M, rep_size, x, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
            conv,
        )
    end
    return x
end
function translate!(M::PowerManifoldNestedReplacing, x, p, q, conv::ActionDirectionAndSide)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        x[i...] = translate(N, _read(M, rep_size, p, i), _read(M, rep_size, q, i), conv)
    end
    return x
end

function inverse_translate!(G::PowerGroup, x, p, q, conv::ActionDirectionAndSide)
    return inverse_translate!(G.manifold, x, p, q, conv)
end
function inverse_translate!(M::AbstractPowerManifold, x, p, q, conv::ActionDirectionAndSide)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        inverse_translate!(
            N,
            _write(M, rep_size, x, i),
            _read(M, rep_size, p, i),
            _read(M, rep_size, q, i),
            conv,
        )
    end
    return x
end
function inverse_translate!(
    M::PowerManifoldNestedReplacing,
    x,
    p,
    q,
    conv::ActionDirectionAndSide,
)
    N = M.manifold
    rep_size = representation_size(N)
    for i in get_iterator(M)
        x[i...] =
            inverse_translate(N, _read(M, rep_size, p, i), _read(M, rep_size, q, i), conv)
    end
    return x
end

function _common_power_adjoint_action!(G, Y, p, X, conv)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        adjoint_action!(
            N,
            _write(GM, rep_size, Y, i),
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, X, i),
            conv,
        )
    end
    return Y
end
function adjoint_action!(G::PowerGroup, Y, p, X, conv::LeftAction)
    return _common_power_adjoint_action!(G, Y, p, X, conv)
end
function adjoint_action!(G::PowerGroup, Y, p, X, conv::RightAction)
    return _common_power_adjoint_action!(G, Y, p, X, conv)
end

function _common_power_replacing_adjoint_action!(G, Y, p, X, conv)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        Y[i...] =
            adjoint_action(N, _read(GM, rep_size, p, i), _read(GM, rep_size, X, i), conv)
    end
    return Y
end
function adjoint_action!(G::PowerGroupNestedReplacing, Y, p, X, conv::LeftAction)
    return _common_power_replacing_adjoint_action!(G, Y, p, X, conv)
end
function adjoint_action!(G::PowerGroupNestedReplacing, Y, p, X, conv::RightAction)
    return _common_power_replacing_adjoint_action!(G, Y, p, X, conv)
end

function translate_diff!(G::PowerGroup, Y, p, q, X, conv::LeftForwardAction)
    return _common_power_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(G::PowerGroup, Y, p, q, X, conv::RightForwardAction)
    return _common_power_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(G::PowerGroup, Y, p, q, X, conv::LeftBackwardAction)
    return _common_power_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(G::PowerGroup, Y, p, q, X, conv::RightBackwardAction)
    return _common_power_translate_diff!(G, Y, p, q, X, conv)
end

function _common_power_translate_diff!(
    G::PowerGroup,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        translate_diff!(
            N,
            _write(GM, rep_size, Y, i),
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, q, i),
            _read(GM, rep_size, X, i),
            conv,
        )
    end
    return Y
end

function translate_diff!(G::PowerGroupNestedReplacing, Y, p, q, X, conv::LeftForwardAction)
    return _common_power_replacing_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(G::PowerGroupNestedReplacing, Y, p, q, X, conv::RightForwardAction)
    return _common_power_replacing_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(G::PowerGroupNestedReplacing, Y, p, q, X, conv::LeftBackwardAction)
    return _common_power_replacing_translate_diff!(G, Y, p, q, X, conv)
end
function translate_diff!(
    G::PowerGroupNestedReplacing,
    Y,
    p,
    q,
    X,
    conv::RightBackwardAction,
)
    return _common_power_replacing_translate_diff!(G, Y, p, q, X, conv)
end
function _common_power_replacing_translate_diff!(
    G::PowerGroupNestedReplacing,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        Y[i...] = translate_diff(
            N,
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, q, i),
            _read(GM, rep_size, X, i),
            conv,
        )
    end
    return Y
end

function inverse_translate_diff!(G::PowerGroup, Y, p, q, X, conv::ActionDirectionAndSide)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        inverse_translate_diff!(
            N,
            _write(GM, rep_size, Y, i),
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, q, i),
            _read(GM, rep_size, X, i),
            conv,
        )
    end
    return Y
end
function inverse_translate_diff!(
    G::PowerGroupNestedReplacing,
    Y,
    p,
    q,
    X,
    conv::ActionDirectionAndSide,
)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        Y[i...] = inverse_translate_diff(
            N,
            _read(GM, rep_size, p, i),
            _read(GM, rep_size, q, i),
            _read(GM, rep_size, X, i),
            conv,
        )
    end
    return Y
end

function exp_lie!(G::PowerGroup, q, X)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        exp_lie!(N, _write(GM, rep_size, q, i), _read(GM, rep_size, X, i))
    end
    return q
end
function exp_lie!(G::PowerGroupNestedReplacing, q, X)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        q[i...] = exp_lie(N, _read(GM, rep_size, X, i))
    end
    return q
end

# on this meta level we first pass down before we resolve identity.
function log_lie!(G::PowerGroup, X, q)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        log_lie!(N, _write(GM, rep_size, X, i), _read(GM, rep_size, q, i))
    end
    return X
end
function log_lie!(G::PowerGroupNestedReplacing, X, q)
    GM = G.manifold
    N = GM.manifold
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        X[i...] = log_lie(N, _read(GM, rep_size, q, i))
    end
    return X
end

#overwrite identity case to avoid allocating identity too early.
function log_lie!(G::PowerGroup, X, ::Identity{ProductOperation})
    GM = G.manifold
    N = GM.manifold
    qi = Identity(N)
    rep_size = representation_size(N)
    for i in get_iterator(GM)
        log_lie!(N, _write(GM, rep_size, X, i), qi)
    end
    return X
end
function log_lie!(G::PowerGroupNestedReplacing, X, ::Identity{ProductOperation})
    GM = G.manifold
    N = GM.manifold
    qi = Identity(N)
    for i in get_iterator(GM)
        X[i...] = log_lie(N, qi)
    end
    return X
end

_filter_out_identities() = ()
_filter_out_identities(x, y...) = (x, _filter_out_identities(y...)...)
_filter_out_identities(::Identity, y...) = _filter_out_identities(y...)

function allocate_result(M::PowerGroup, f, x...)
    return allocate_result(M.manifold, f, _filter_out_identities(x...)...)
end
function allocate_result(M::PowerGroupNestedReplacing, f, x...)
    return allocate_result(M.manifold, f, _filter_out_identities(x...)...)
end
function allocate_result(M::PowerGroup, f::typeof(identity_element))
    return allocate_result(M.manifold, f)
end
function allocate_result(M::PowerGroupNestedReplacing, f::typeof(identity_element))
    return [allocate_result(M.manifold.manifold, f) for _ in get_iterator(M.manifold)]
end

@inline function _read(
    M::AbstractPowerManifold,
    ::Tuple,
    ::Identity{ProductOperation},
    ::Int,
)
    return Identity(M.manifold)
end
