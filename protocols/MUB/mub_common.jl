module MUB_Common

using LinearAlgebra
using Ket # For mutually unbiased bases

export cond_prob, HAgB

function cond_prob(; d, x, y, a, b, q)
    psi = vec(Matrix(I, d, d))/sqrt(d)
    rho = (1 - q)*psi*psi' + q*Matrix(I, d^2, d^2)/d^2
    Bs = Ket.mub(d)
    v1 = Bs[x][:, a]
    v2 = Bs[y][:, b]
    return tr(kron(v1*v1', v2*v2') * rho) |> real
end

function HAgB(d::Int; q::Number)
    HAB = 0.0
    for a in 1:d, b in 1:d
        p = cond_prob(; d, x=1, y=1, a, b, q)
        if p != 0
            HAB += -p*log2(p)
        end
    end

    HB = 0.0
    for b in 1:d
        p = 0.0
        for a in 1:d
            p += cond_prob(; d, x=1, y=1, a, b, q)
        end
        if p != 0
            HB += -p*log2(p)
        end
    end

    return HAB - HB
end

end # module MUB_Common
