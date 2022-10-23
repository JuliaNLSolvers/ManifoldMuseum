
@doc raw"""
based on
.. [1] Awad H. Al-Mohy and Nicholas J. Higham (2009)
           Computing the Frechet Derivative of the Matrix Exponential,
           with an application to Condition Number Estimation.
           SIAM Journal On Matrix Analysis and Applications.,
           30 (4). pp. 1639-1657. ISSN 1095-7162
"""

ell_table_61 = (
        nothing,
        # 1
        2.11e-8,
        3.56e-4,
        1.08e-2,
        6.49e-2,
        2.00e-1,
        4.37e-1,
        7.83e-1,
        1.23e0,
        1.78e0,
        2.42e0,
        # 11
        3.13e0,
        3.90e0,
        4.74e0,
        5.63e0,
        6.56e0,
        7.52e0,
        8.53e0,
        9.56e0,
        1.06e1,
        1.17e1,
)

@doc raw"""
    _diff_pade3(A, E)

Compute expm using 3-term pade approximant
"""
@inline function _diff_pade3(A, E)
    b = (120., 60., 12., 1.)
    A2 = A * A
    M2 = A * E + E*A
    U = A * (b[4]*A2 + UniformScaling(b[2]))
    V = b[3]*A2 + UniformScaling(b[1])
    Lu = A * (b[3]*M2) + E * (b[3]*A2 + UniformScaling(b[1]))
    Lv = b[3] .* M2
    return U, V, Lu, Lv
end        

@doc raw"""
    _diff_pade5(A, E)

Compute expm using 5-term pade approximant
"""
@inline function _diff_pade5(A, E)
    b = (30240., 15120., 3360., 420., 30., 1.)
    A2 = A * A
    M2 = A * E + E * A
    A4 = A2 * A2
    M4 = A2 * M2 + M2 * A2
    U = A * (b[6]*A4 + b[4]*A2 + UniformScaling(b[2]))
    V = b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    Lu = (A * (b[6]*M4 + b[4]*M2) +
            E * (b[6]*A4 + b[4]*A2 + UniformScaling(b[2])))
    Lv = b[5]*M4 + b[3]*M2
    return U, V, Lu, Lv
end

@doc raw"""
    _diff_pade7(A, E)

Compute expm using 7-term pade approximant
"""
@inline function _diff_pade7(A, E)
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    A2 = A * A
    M2 = A * E + E * A
    A4 = A2 * A2
    M4 = A2 * M2 + M2 * A2
    A6 = A2 * A4
    M6 = A4 * M2 + M4 * A2
    U = A * (b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2]))
    V = b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    Lu = (A*(b[8]*M6 + b[6]*M4 + b[4]*M2) +
            E*(b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])))
    Lv = b[7]*M6 + b[5]*M4 + b[3]*M2
    return U, V, Lu, Lv
end

@doc raw"""
    _diff_pade9(A, E)

Compute expm using 9-term pade approximant
"""
@inline function _diff_pade9(A, E)
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
            2162160., 110880., 3960., 90., 1.)
    A2 = A * A
    M2 = A * E + E * A
    A4 = A2 * A2
    M4 = A2 * M2 + M2 * A2
    A6 = A2 * A4
    M6 = A4 * M2 + M4 * A2
    A8 = A4 * A4
    M8 = A4 * M4 + M4 * A4
    U = A * (b[10]*A8 + b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2]))
    V = b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    Lu = (A *(b[10]*M8 + b[8]*M6 + b[6]*M4 + b[4]*M2) +
            E * (b[10]*A8 + b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])))
    Lv = b[9]*M8 + b[7]*M6 + b[5]*M4 + b[3]*M2
    return U, V, Lu, Lv
end

@doc raw"""
    _diff_pade13(A, E)
Compute expm using 13-term pade approximant
"""
@inline function _diff_pade13(A, E)
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
         1187353796428800., 129060195264000., 10559470521600.,
         670442572800., 33522128640., 1323241920., 40840800., 960960.,
         16380., 182., 1.)
    
    A2 = A * A
        
    M2 = A * E + E * A
    A4 = A2 * A2
    M4 = A2 * M2 + M2 * A2
    A6 = A2 * A4
    M6 = A4 * M2 + M4 * A2
    W1 = b[14]*A6 + b[12]*A4 + b[10]*A2
    W2 = b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])
    Z1 = b[13]*A6 + b[11]*A4 + b[9]*A2
    Z2 = b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    W = A6 * W1 + W2
    U = A * W
    V = A6 * Z1 + Z2
    Lw1 = b[14]*M6 + b[12]*M4 + b[10]*M2
    Lw2 = b[8]*M6 + b[6]*M4 + b[4]*M2
    Lz1 = b[13]*M6 + b[11]*M4 + b[9]*M2
    Lz2 = b[7]*M6 + b[5]*M4 + b[3]*M2
    Lw = A6 * Lw1 + M6 * W1 + Lw2
    Lu = A * Lw + E * W
    Lv = A6 * Lz1 + M6 * Z1 + Lz2
    return U, V, Lu, Lv
end    

@doc raw"""
    expmm_frechet(A, E)
Compute frechet derivative of expm(A) in direction E using algorithm 6.4 of @ref
"""
function expm_frechet(A, E)
    n = size(A, 1)
    s = nothing    
    A_norm_1 = maximum(sum(abs.(A), dims=1))
    m_pade_pairs = (
        (3, _diff_pade3),
        (5, _diff_pade5),
        (7, _diff_pade7),
        (9, _diff_pade9)
    )
    
    for m_pade in m_pade_pairs
        m, pade = m_pade
        if A_norm_1 <= ell_table_61[m]
            U, V, Lu, Lv = pade(A, E)
            s = 0
            break
        end            
    end
    if isnothing(s)
        # scaling
        s = max(0, Int(ceil(log2(A_norm_1 / ell_table_61[14]))))
        # pade order 13
        U, V, Lu, Lv = _diff_pade13((2.0^-s) * A, (2.0^-s) * E)
    end
    # factor once and solve twice    
    lu_piv = lu(-U + V)
    eA = lu_piv \ (U + V)
    eAf = lu_piv \ (Lu + Lv + (Lu - Lv)* eA)

    # squaring
    for k in 1:s
        eAf = eA * eAf + eAf * eA
        eA = eA * eA
    end

    return eA, eAf
end

@doc raw"""
    _diff_pade3!(buff, A, E)
Compute expm using 3-term pade approximant, with buff used as temporary storage
The returns U, V, Lu, Lv are stored in the first four blocks of buff
"""
@inline function _diff_pade3!(buff, A, E)
    b = (120., 60., 12., 1.)
    k = size(A)[1]
    
    @views begin
        U = buff[1:k, 1:end]
        V = buff[k+1:2*k, 1:end]
        Lu = buff[2*k+1:3*k, 1:end]
        Lv = buff[3*k+1:4*k, 1:end]
        
        A2 = buff[4*k+1:5*k, 1:end]
        M2 = buff[5*k+1:6*k, 1:end]
    end
    
    mul!(A2, A, A)
    mul!(M2, A, E)
    mul!(M2, E, A, 1, 1)
    mul!(U, A, b[2])
    mul!(U, A, A2, b[4], 1)
    V .= b[3]*A2 + UniformScaling(b[1])
    mul!(Lu, E, V)
    mul!(Lu, A, M2, b[3], 1)
    mul!(Lv,  b[3], M2)
end        

@doc raw"""
    _diff_pade5!(buff, A, E)
Compute expm using 5-term pade approximant, with buff used as temporary storage
The returns U, V, Lu, Lv are stored in the first four blocks of buff
"""
@inline function _diff_pade5!(buff, A, E)
    b = (30240., 15120., 3360., 420., 30., 1.)
    k = size(A)[1]    
    @views begin
        U = buff[1:k, 1:end]
        V = buff[k+1:2*k, 1:end]
        Lu = buff[2*k+1:3*k, 1:end]
        Lv = buff[3*k+1:4*k, 1:end]
        
        A2 = buff[4*k+1:5*k, 1:end]
        M2 = buff[5*k+1:6*k, 1:end]

        A4 = buff[6*k+1:7*k, 1:end]
        M4 = buff[7*k+1:8*k, 1:end]
    end
    
    mul!(A2, A, A)
    mul!(M2, A, E)
    mul!(M2, E, A, 1, 1)
    
    mul!(A4, A2, A2)
    mul!(M4, A2, M2)
    mul!(M4, M2, A2, 1, 1)
    Z = b[6]*A4 + b[4]*A2 + UniformScaling(b[2])
    mul!(U, A, Z)
    V .= b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    mul!(Lu, E, Z)
    mul!(Lu, A, M4, b[6], 1)
    mul!(Lu, A, M2, b[4], 1)
    
    Lv .= b[5]*M4 + b[3]*M2
end

@doc raw"""
    _diff_pade7!(buff, A, E)
Compute expm using 7-term pade approximant, with buff used as temporary storage
The returns U, V, Lu, Lv are stored in the first four blocks of buff
"""
@inline function _diff_pade7!(buff, A, E)
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    k = size(A)[1]
    @views begin
        U = buff[1:k, 1:end]
        V = buff[k+1:2*k, 1:end]
        Lu = buff[2*k+1:3*k, 1:end]
        Lv = buff[3*k+1:4*k, 1:end]
        
        A2 = buff[4*k+1:5*k, 1:end]
        M2 = buff[5*k+1:6*k, 1:end]

        A4 = buff[6*k+1:7*k, 1:end]
        M4 = buff[7*k+1:8*k, 1:end]

        A6 = buff[8*k+1:9*k, 1:end]
        M6 = buff[9*k+1:10*k, 1:end]        
    end    
    
    mul!(A2, A, A)
    mulsym!(M2, A,  E)
    mul!(A4, A2, A2)
    mulsym!(M4, A2, M2)
    
    mul!(A6,  A2, A4)
    mul!(M6, A4, M2)
    mul!(M6, M4, A2, 1, 1)
    
    Z= b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])
    mul!(U, A, Z)
    V .= b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    mul!(Lu, E, Z)
    mul!(Lu, A, b[8]*M6 + b[6]*M4 + b[4]*M2, 1, 1)
    Lv .= b[7]*M6 + b[5]*M4 + b[3]*M2
end

@doc raw"""
    _diff_pade9!(buff, A, E)
Compute expm using 9-term pade approximant, with buff used as temporary storage
The returns U, V, Lu, Lv are stored in the first four blocks of buff
"""
@inline function _diff_pade9!(buff, A, E)
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
         2162160., 110880., 3960., 90., 1.)

    k = size(A)[1]
    @views begin
        U = buff[1:k, 1:end]
        V = buff[k+1:2*k, 1:end]
        Lu = buff[2*k+1:3*k, 1:end]
        Lv = buff[3*k+1:4*k, 1:end]
        
        A2 = buff[4*k+1:5*k, 1:end]
        M2 = buff[5*k+1:6*k, 1:end]

        A4 = buff[6*k+1:7*k, 1:end]
        M4 = buff[7*k+1:8*k, 1:end]

        A6 = buff[8*k+1:9*k, 1:end]
        M6 = buff[9*k+1:10*k, 1:end]

        A8 = buff[10*k+1:11*k, 1:end]
        M8 = buff[11*k+1:12*k, 1:end]        
    end    
    
    mul!(A2, A, A)
    mulsym!(M2, A, E)
    mul!(A4, A2, A2)
    mulsym!(M4, A2, M2)
    mul!(A6, A2, A4)
    mul!(M6, A4, M2)
    mul!(M6, M4, A2, 1, 1)
    mul!(A8, A4, A4)
    mulsym!(M8, A4, M4)
    Z = b[10]*A8 + b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])
    mul!(U, A, Z)
    V .= b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])
    
    mul!(Lu, E, Z)
    mul!(Lu, A, b[10]*M8 + b[8]*M6 + b[6]*M4 + b[4]*M2, 1, 1)
    Lv .= b[9]*M8 + b[7]*M6 + b[5]*M4 + b[3]*M2
end

@doc raw"""
    _diff_pade13!(buff, A, E)
Compute expm using 13-term pade approximant, with buff used as temporary storage
The returns U, V, Lu, Lv are stored in the first four blocks of buff
"""
@inline function _diff_pade13!(buff, A, E)
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
         1187353796428800., 129060195264000., 10559470521600.,
         670442572800., 33522128640., 1323241920., 40840800., 960960.,
         16380., 182., 1.)
    k = size(A)[1]

    @views begin
        U = buff[1:k, 1:end]
        V = buff[k+1:2*k, 1:end]
        Lu = buff[2*k+1:3*k, 1:end]
        Lv = buff[3*k+1:4*k, 1:end]
        
        A2 = buff[4*k+1:5*k, 1:end]
        M2 = buff[5*k+1:6*k, 1:end]
        
        A4 = buff[6*k+1:7*k, 1:end]
        M4 = buff[7*k+1:8*k, 1:end]
        
        A6 = buff[8*k+1:9*k, 1:end]
        M6 = buff[9*k+1:10*k, 1:end]
        
        W1 = buff[10*k+1:11*k, 1:end]
        Z1 = buff[11*k+1:12*k, 1:end]
        W = buff[12*k+1:13*k, 1:end]
        Lw1 = buff[13*k+1:14*k, 1:end]
        Lz1 = buff[14*k+1:15*k, 1:end]        
        Lw = buff[15*k+1:16*k, 1:end]
    end
    
    mul!(A2, A, A)
    mul!(M2, A, E)
    mul!(M2, E, A, 1, 1)
    mul!(A4, A2, A2)
    mul!(M4, A2, M2)
    mul!(M4, M2, A2, 1, 1)
    mul!(A6, A2, A4)
    mul!(M6, A4, M2)
    mul!(M6, M4, A2, 1, 1)
    mul!(W1, b[14], A6)
    W1 .+= b[12].*A4
    W1 .+= b[10].*A2
    
    mul!(Z1, b[13], A6)
    Z1 .+= b[11].*A4
    Z1 .+= b[9].*A2
    
    mul!(W, A6, W1)
    W .+= b[8]*A6 + b[6]*A4 + b[4]*A2 + UniformScaling(b[2])
    
    mul!(U, A,  W)
    mul!(V, A6, Z1)
    V .+= b[7]*A6 + b[5]*A4 + b[3]*A2 + UniformScaling(b[1])

    Lw1 .= b[14]*M6 + b[12]*M4 + b[10]*M2
    mul!(Lz1, b[13], M6)
    Lz1 .+= b[11]*M4 + b[9]*M2
    
    mul!(Lw, A6, Lw1)
    mul!(Lw, M6, W1, 1, 1)
    Lw .+= b[8]*M6 + b[6]*M4 + b[4]*M2
    
    mul!(Lu, A, Lw)
    mul!(Lu, E, W, 1, 1)
    mul!(Lv, A6, Lz1)
    mul!(Lv, M6,  Z1, 1, 1)
    Lv .+=  b[7]*M6 + b[5]*M4 + b[3]*M2
end    


@doc raw"""
    expmm_frechet!(buff, A, E)
Compute frechet derivative of expm(A) in direction E using algorithm 6.4 of @ref
buff is a matrix of size 16*k times k
the returns, eA = exp(A), eAf = dexp(A, E) are stored in the first two blocks
the remaining blocks are used as temporary storage
"""
function expm_frechet!(buff, A, E)
    n = size(A, 1)
    s = nothing
    A_norm_1 = maximum(sum(abs.(A), dims=1))
    k = size(A)[1]    
    m_pade_pairs = (
        (3, _diff_pade3!),
        (5, _diff_pade5!),
        (7, _diff_pade7!),
        (9, _diff_pade9!)
    )
    
    for m_pade in m_pade_pairs
        m, pade = m_pade
        if A_norm_1 <= ell_table_61[m]
            U, V, Lu, Lv = pade(buff, A, E)
            s = 0
            break
        end            
    end
    if isnothing(s)
        # scaling
        s = max(0, Int(ceil(log2(A_norm_1 / ell_table_61[14]))))
        # pade order 13
        _diff_pade13!(buff, (2.0^-s) * A, (2.0^-s) * E)
    end
    buff[4*k+1:8*k, :] .= buff[1:4*k, :]
    @views begin
        eA = buff[1:k, :]
        eAf = buff[k+1:2*k, :]
        tmp = buff[2*k+1:3*k, :]
        
        U = buff[4*k+1:5*k, :]
        V = buff[5*k+1:6*k, :]
        Lu = buff[6*k+1:7*k, :]
        Lv = buff[7*k+1:8*k, :]
    end
    
    # factor once and solve twice    
    lu_piv = lu(-U + V)
    broadcast!(+, eA, U, V)
    ldiv!(lu_piv, eA)
    broadcast!(-, tmp, Lu, Lv)
    mul!(eAf, tmp, eA)
    eAf .+= Lu .+ Lv
    ldiv!(lu_piv, eAf)

    # squaring
    for k in 1:s
        mulsym!(tmp, eA, eAf)
        eAf .= tmp
        mul!(tmp, eA, eA)
        eA .= tmp
    end
end

@doc raw"""
    mulsym!(C, A, E)
Compute C = A*E + E*A by mul    
"""
@inline function mulsym!(C, A, B)
    mul!(C, A, B)
    mul!(C, B, A, 1, 1)
end    

