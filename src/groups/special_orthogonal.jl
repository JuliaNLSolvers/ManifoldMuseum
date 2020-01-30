@doc doc"""
    SpecialOrthogonal{n} <: GroupManifold{Rotations{n},MultiplicationOperation}

Special orthogonal group $\mathrm{SO}(n)$ represented by rotation matrices.

# Constructor
    SpecialOrthogonal(n)
"""
const SpecialOrthogonal{n} = GroupManifold{Rotations{n},MultiplicationOperation}

has_invariant_metric(::SpecialOrthogonal, ::ActionDirection) = Val(true)

is_default_metric(::MetricManifold{<:SpecialOrthogonal,EuclideanMetric}) = Val(true)

SpecialOrthogonal(n) = SpecialOrthogonal{n}(Rotations(n), MultiplicationOperation())

show(io::IO, ::SpecialOrthogonal{n}) where {n} = print(io, "SpecialOrthogonal($(n))")

inv(::SpecialOrthogonal, x) = transpose(x)

inverse_translate(G::SpecialOrthogonal, x, y, conv::LeftAction) = inv(G, x) * y
inverse_translate(G::SpecialOrthogonal, x, y, conv::RightAction) = y * inv(G, x)

translate_diff(::SpecialOrthogonal, x, y, v, ::LeftAction) = v
translate_diff(G::SpecialOrthogonal, x, y, v, ::RightAction) = inv(G, x) * v * x

function translate_diff!(G::SpecialOrthogonal, vout, x, y, v, conv::ActionDirection)
    return copyto!(vout, translate_diff(G, x, y, v, conv))
end

function inverse_translate_diff(G::SpecialOrthogonal, x, y, v, conv::ActionDirection)
    return translate_diff(G, inv(G, x), y, v, conv)
end

function inverse_translate_diff!(G::SpecialOrthogonal, vout, x, y, v, conv::ActionDirection)
    return copyto!(vout, inverse_translate_diff(G, x, y, v, conv))
end

group_exp!(G::SpecialOrthogonal, y, v) = exp!(G, y, Identity(G), v)

group_log!(G::SpecialOrthogonal, v, y) = log!(G, v, Identity(G), y)

function allocate_result(
    ::GT,
    ::typeof(exp),
    ::Identity{GT},
    v,
) where {n,GT<:SpecialOrthogonal{n}}
    return allocate(v)
end
function allocate_result(
    ::GT,
    ::typeof(log),
    ::Identity{GT},
    y,
) where {n,GT<:SpecialOrthogonal{n}}
    return allocate(y)
end
