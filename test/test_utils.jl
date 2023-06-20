@testset "utils test" begin
    Random.seed!(42)
    @testset "usinc_from_cos" begin
        @test Manifolds.usinc_from_cos(-1) == 0
        @test Manifolds.usinc_from_cos(-1.0) == 0.0
    end
    @testset "log_safe!" begin
        n = 8
        Q = qr(randn(n, n)).Q
        A1 = Matrix(Hermitian(Q * Diagonal(rand(n)) * Q'))
        @test exp(Manifolds.log_safe!(similar(A1), A1)) ≈ A1 atol = 1e-6
        A1_fail = Matrix(Hermitian(Q * Diagonal([-1; rand(n - 1)]) * Q'))
        @test_throws DomainError Manifolds.log_safe!(similar(A1_fail), A1_fail)

        T = triu!(randn(n, n))
        T[diagind(T)] .= rand.()
        @test exp(Manifolds.log_safe!(similar(T), T)) ≈ T atol = 1e-6
        T_fail = copy(T)
        T_fail[1] = -1
        @test_throws DomainError Manifolds.log_safe!(similar(T_fail), T_fail)

        A2 = Q * T * Q'
        @test exp(Manifolds.log_safe!(similar(A2), A2)) ≈ A2 atol = 1e-6
        A2_fail = Q * T_fail * Q'
        @test_throws DomainError Manifolds.log_safe!(similar(A2_fail), A2_fail)

        A3 = exp(SizedMatrix{n,n}(randn(n, n)))
        @test A3 isa SizedMatrix
        @test exp(Manifolds.log_safe!(similar(A3), A3)) ≈ A3 atol = 1e-6
        @test exp(Manifolds.log_safe(A3)) ≈ A3 atol = 1e-6

        A3_fail = Float64[1 2; 3 1]
        @test_throws DomainError Manifolds.log_safe!(similar(A3_fail), A3_fail)

        A4 = randn(ComplexF64, n, n)
        @test exp(Manifolds.log_safe!(similar(A4), A4)) ≈ A4 atol = 1e-6
    end
    @testset "isnormal" begin
        @test !Manifolds.isnormal([1.0 2.0; 3.0 4.0])
        @test !Manifolds.isnormal(complex.(reshape(1:4, 2, 2), reshape(5:8, 2, 2)))

        # diagonal
        @test Manifolds.isnormal(diagm(randn(5)))
        @test Manifolds.isnormal(diagm(randn(ComplexF64, 5)))
        @test Manifolds.isnormal(Diagonal(randn(5)))
        @test Manifolds.isnormal(Diagonal(randn(ComplexF64, 5)))

        # symmetric/hermitian
        @test Manifolds.isnormal(Symmetric(randn(3, 3)))
        @test Manifolds.isnormal(Hermitian(randn(3, 3)))
        @test Manifolds.isnormal(Hermitian(randn(ComplexF64, 3, 3)))
        x = Matrix(Symmetric(randn(3, 3)))
        x[3, 1] += eps()
        @test !Manifolds.isnormal(x)
        @test Manifolds.isnormal(x; atol=sqrt(eps()))

        # skew-symmetric/skew-hermitian
        skew(x) = x - x'
        @test Manifolds.isnormal(skew(randn(3, 3)))
        @test Manifolds.isnormal(skew(randn(ComplexF64, 3, 3)))

        # orthogonal/unitary
        @test Manifolds.isnormal(Matrix(qr(randn(3, 3)).Q); atol=sqrt(eps()))
        @test Manifolds.isnormal(Matrix(qr(randn(ComplexF64, 3, 3)).Q); atol=sqrt(eps()))
    end
    @testset "realify/unrealify!" begin
        # round trip real
        x = randn(3, 3)
        @test Manifolds.realify(x, ℝ) === x
        @test Manifolds.unrealify!(similar(x), x, ℝ) == x

        # round trip complex
        x2 = randn(ComplexF64, 3, 3)
        x2r = Manifolds.realify(x2, ℂ)
        @test eltype(x2r) <: Real
        @test size(x2r) == (6, 6)
        x2c = Manifolds.unrealify!(similar(x2), x2r, ℂ)
        @test x2c ≈ x2

        # matrix multiplication is preserved
        x3 = randn(ComplexF64, 3, 3)
        x3r = Manifolds.realify(x3, ℂ)
        @test x2 * x3 ≈ Manifolds.unrealify!(similar(x2), x2r * x3r, ℂ)
    end
    @testset "allocation" begin
        @test allocate([1 2; 3 4], Float64, Size(3, 3)) isa Matrix{Float64}
        @test allocate(SA[1 2; 3 4], Float64, Size(3, 3)) isa MMatrix{3,3,Float64}
        @test allocate(SA[1 2; 3 4], Size(3, 3)) isa MMatrix{3,3,Int}
        @test Manifolds.quat_promote(Float64) === Quaternions.QuaternionF64
        @test Manifolds.quat_promote(Float32) === Quaternions.QuaternionF32
        @test Manifolds.quat_promote(QuaternionF64) === Quaternions.QuaternionF64
        @test Manifolds.quat_promote(QuaternionF32) === Quaternions.QuaternionF32
    end
    @testset "eigen_safe" begin
        @test Manifolds.eigen_safe(SA[1.0 0.0; 0.0 1.0]) isa
              Eigen{Float64,Float64,<:SizedMatrix{2,2},<:SizedVector{2}}
    end
    @testset "max_eps" begin
        x64 = randn(Float64, 2)
        x32 = randn(Float32, 2)
        z32 = randn(ComplexF32, 2)
        xi = rand(0:1, 2)
        @test Manifolds.max_eps(x64, x64) == eps()
        @test Manifolds.max_eps(x64, x32) == eps(Float32)
        @test Manifolds.max_eps(x32, x64) == eps(Float32)
        @test Manifolds.max_eps(xi, xi) == 0
        @test Manifolds.max_eps(xi, x64) == eps()
        @test Manifolds.max_eps(xi, x32) == eps(Float32)
        @test Manifolds.max_eps(xi, z32) == eps(Float32)
        @test Manifolds.max_eps(xi, x64, x32, z32) == eps(Float32)
    end
    @test Manifolds.is_metric_function(flat)
    @test Manifolds.is_metric_function(sharp)
end