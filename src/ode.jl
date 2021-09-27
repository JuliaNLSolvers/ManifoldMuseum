function solve_exp_ode(
    M::MetricManifold,
    p,
    X,
    tspan,
    B::InducedBasis{𝔽,TangentSpaceType};
    solver=AutoVern9(Rodas5()),
    backend=diff_backend(),
    kwargs...,
) where {𝔽}
    pc = get_parameters(M, B.A, B.i, p)
    Xc = get_coordinates(M, p, X, B)
    n = length(pc)
    iX = SVector{n}(1:n)
    ip = SVector{n}((n + 1):(2n))
    u0 = allocate(pc, 2n)
    u0[iX] .= Xc
    u0[ip] .= pc

    function exp_problem(u, params, t)
        M = params[1]
        dp = u[iX]
        p = u[ip]
        ddp = allocate(u, Size(n))
        du = allocate(u)
        Γ = christoffel_symbols_second(M, pc, B; backend=backend)
        @einsum ddp[k] = -Γ[k, i, j] * dp[i] * dp[j]
        du[iX] .= ddp
        du[ip] .= dp
        return Base.convert(typeof(u), du)
    end

    params = (M,)
    prob = ODEProblem(exp_problem, u0, tspan, params)
    sol = solve(prob, solver; kwargs...)
    return get_point(M, B.A, B.i, sol.u[1][(n + 1):end])
end
