"""
    process_manifolds(args, warn=true)

Given a vector of strings of manifold types, check whether they exist and add them
to an array of types.
"""
function process_manifolds(args::Vector{String}, warn::Bool=true)
    manifolds = Type[]
    for s in args
        m = Symbol(s)
        if hasproperty(Manifolds, m)
            M = getproperty(Manifolds, m)
            if M <: AbstractManifold
                push!(manifolds, M)
            else
                warn && (@warn "`$s` seems not to be a manifold.")
            end
        else
            warn && (@warn "`$s` seems not to be defined in `Manifolds` namespace.")
        end
    end
    return manifolds
end
"""
    process_functionds(args, warn=true)

Given a vector of strings of functions types, check whether they exist and add them
to an array of Functions.
"""
function process_functions(args::Vector{String}, warn::Bool=true)
    functions = Function[]
    for s in args
        m = Symbol(s)
        if hasproperty(Manifolds, m)
            f = getproperty(Manifolds, m)
            if f isa Function
                push!(functions, f)
            else
                warn && (@warn "`$s` seems not to be a function.")
            end
        else
            warn && (@warn "`$s` seems not to be defined in `Main` namespace.")
        end
    end
    return functions
end

"""
    collect_args(start, args)

Given a vector of arguments, and a start, this function collects all
arguments at start, start+1,... until an argument starting with `--` appears.
It returns a tupe (j, candidates)
of the collected arguments as `candidates` and `j` indicating the first argment
not collected, i.e. the one starting again with `--`.
"""
function collect_args(start, args)
    candidates = String[]
    # collect manifolds
    i = start
    while i <= length(args)
        startswith(args[i], "--") && break
        push!(candidates, args[i])
        i = i + 1
    end
    return (i, candidates)
end
