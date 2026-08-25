module BB84_SDI

export finite_key_rate, finite_key_rate_opt

using LinearAlgebra
using JuMP, Dualization, MosekTools
using Optim

# const Ent_BFF = include("../common/ent_BFF.jl").Ent_BFF
# using .Ent_BFF

const Ent_KS = include("../common/ent_KS.jl").Ent_KS
using .Ent_KS

const FiniteSize = include("finite_size.jl").FiniteSize

"""
Conditional probability before binning Bob's measurement outcomes.
"""
function cond_prob_full(; x, y, a, b, q, η_A, η_B)
    phi = [1, 0, 0, 1]/sqrt(2)
    rho = (1 - q)*phi*phi' + q*Matrix(I, 4, 4)/4
    id = [1 0; 0 1]
    Ms = [
        [η_A*[1 0; 0 0], η_A*[0 0; 0 1], (1-η_A)*id],
        [η_A*0.5*[1 1; 1 1], η_A*0.5*[1 -1; -1 1], (1-η_A)*id],
    ]
    Ns = [
        [η_B*[1 0; 0 0], η_B*[0 0; 0 1], (1-η_B)*id],
        [η_B*0.5*[1 1; 1 1], η_B*0.5*[1 -1; -1 1], (1-η_B)*id],
    ]
    return tr(kron(Ms[x][a], Ns[y][b])*rho)
end

"""
Conditional probability after binning Bob's non-detection outcome to b=0.
"""
function cond_prob(; x, y, a, b, q, η_A, η_B)
    p(; x, y, a, b) = cond_prob_full(; x, y, a, b, q, η_A, η_B)
    if b == 1
        return p(; x, y, a, b=1) + p(; x, y, a, b=3)
    elseif b == 2
        return p(; x, y, a, b=2)
    end
end

function p_hon(q, η_A, η_B)
    p_err = 0.5*(cond_prob(; x=2, y=2, a=1, b=2, q, η_A, η_B) + cond_prob(; x=2, y=2, a=2, b=1, q, η_A, η_B))
    p_loss = 0.5*sum(cond_prob(; x=2, y=2, a=3, b, q, η_A, η_B) for b in 1:2)
    p_rem = 1 - p_err - p_loss
    @assert p_rem >= 0
    return [p_err, p_loss, p_rem]
end

function vN_entropy_lossy(p, γ)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)

    f = 0.5*(1-γ) # Sifting factor
    Ms = [
        [
            [1 0 0; 0 0 0; 0 0 0],
            [0 0 0; 0 1 0; 0 0 0],
        ],
        [
            0.5*[1 1 0; 1 1 0; 0 0 0],
            0.5*[1 -1 0; -1 1 0; 0 0 0],
        ],
    ]

    function constraints(Ns)
        cs = []
        id = Matrix(I, 3, 3)

        # Upper-bound on QBER in complementary basis
        c = p[1]*id - 0.5*(kron(Ms[2][1], Ns[2][2]) + kron(Ms[2][2], Ns[2][1]))
        push!(cs, c)

        # Upper-bound on losses
        N_tot = Ns[1][1] + Ns[1][2]
        c = p[2]*id - 0.5*kron(id - (Ms[1][1] + Ms[1][2]), N_tot) 
        push!(cs, c)

        return cs
    end

    return vN_entropy(; Ms=f*Ms[1], dY=2, dB=2, objective=objective_A, constraints, optimizer)
end

function vN_entropy_lossy2(grad, γ)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)

    constraints(_) = []

    f = 0.5*(1-γ) # Sifting factor
    Ms = [
        [
            [1 0 0; 0 0 0; 0 0 0],
            [0 0 0; 0 1 0; 0 0 0],
        ],
        [
            0.5*[1 1 0; 1 1 0; 0 0 0],
            0.5*[1 -1 0; -1 1 0; 0 0 0],
        ],
    ]

    function offset(Ns)
        id = Matrix(I, 3, 3)
        N_tot = Ns[2][1] + Ns[2][2]

        delta = zeros(3, 3)
        delta += grad[1]*0.5*(kron(Ms[2][1], Ns[2][2]) + kron(Ms[2][2], Ns[2][1]))
        delta += grad[2]*0.5*kron(id - (Ms[2][1] + Ms[2][2]), N_tot)
        return -delta
    end

    h_CA, _  = vN_entropy(; Ms=f*Ms[1], dY=2, dB=2, objective=objective_A, offset, constraints, optimizer)
    return h_CA
end

# function vN_entropy_lossy_full_stats(q, η_A, η_B)
#     optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
#     optimizer = dual_optimizer(optimizer)
# 
#     function constraints(; Ms, Ns)
#         cs = []
#         id = Matrix(I, 3, 3)
#         for x in 1:2, y in 1:2
#             for a in 1:2, b in 1:2
#                 c = kron(Ms[x][a], Ns[y][b]) - id*cond_prob(; x, y, a, b, q, η_A, η_B)
#                 push!(cs, c)
#                 push!(cs, -c)
#             end
#         end
# 
#         return cs
#     end
# 
#     Ms = [
#         [
#             [1 0 0; 0 0 0; 0 0 0],
#             [0 0 0; 0 1 0; 0 0 0],
#         ],
#         [
#             0.5*[1 1 0; 1 1 0; 0 0 0],
#             0.5*[1 -1 0; -1 1 0; 0 0 0],
#         ],
#     ]
# 
#     return vN_entropy(; Ms, dY=2, dB=2, objective=objective_A, constraints, optimizer)
# end

function HAgB(q, η_A, η_B)
    HAB = 0
    for a in 1:2, b in 1:3
        p = cond_prob_full(; x=1, y=1, a, b, q, η_A, η_B)
        if iszero(p)
            continue
        end
        HAB += -p*log2(p)
    end

    HB = 0
    for b in 1:3
        p = 0
        for a in 1:2
            p += cond_prob_full(; x=1, y=1, a, b, q, η_A, η_B)
        end
        if iszero(p)
            continue
        end
        HB += -p*log2(p)
    end

    return HAB - HB
end

function find_delta(ε_EV_cmp, n)
    # b = -2*Delta/(3*n)*log2(1/eps_EV_cmp)
    # c = -6*var/(3*n)*log2(1/eps_EV_cmp)
    # delta = (-b + sqrt(b^2 - 4*c))/2
    # return delta
    dC = 3
    return sqrt(1/(2n)*log(2dC/ε_EV_cmp))
end

function EC_size(HAgB, n, γ, ε_KV_cmp)
    dS = 3
    return n*0.5*(1-γ)*HAgB + 2*sqrt(n)*sqrt(1 + 2*log2(2/ε_KV_cmp))*log2(1+2*dS) + 2*log2(2/ε_KV_cmp)
end

function finite_key_rate(p_hon::Vector, HAgB; n, α=1+1/√n, γ=0.1)
    h_CA_hon, neg_grad = vN_entropy_lossy(p_hon, γ)
    grad = -neg_grad
    r1 = finite_key_rate(p_hon, h_CA_hon, grad, HAgB; n, α, γ)
    # r2 = finite_key_rate(p_hon, grad, HAgB; n, α, γ)
    # @show r1, r2
    return r1
end

function finite_key_rate(p_hon::Vector, grad::Vector, HAgB; n, α=1+1/√n, γ=0.1)
    h = vN_entropy_lossy2(grad, γ)
    h_CA_hon = h + dot([grad; 0], p_hon)
    grad = grad

    return finite_key_rate(p_hon, h_CA_hon, grad, HAgB; n, α, γ)
end

function finite_key_rate(p_hon::Vector, h_CA_hon::Number, grad::Vector, HAgB; n, α=1+1/√n, γ)
    ε_KV_cmp = ε_EV_cmp = 0.5*1e-2
    λ_EC = EC_size(HAgB, n, γ, ε_KV_cmp)

    δ = find_delta(ε_EV_cmp, n)
    max_grad = max(grad...)
    h_CA = h_CA_hon - 1/γ*dot(max_grad*[1, 1] - grad, [δ, δ])

    # Crossover min-tradeoff function
    f = FiniteSize.MinTradeoffData(
        h_CA,
        max(grad...) - min(grad...),
        1/γ*(max(grad...) - min(grad...))^2,
    )

    ε_snd = 1e-12
    ε_KV = ε_PA = 0.1*ε_snd
    H_max = n*(γ*(0.5 - p_hon[2]) + δ) # We have that only two values of C are possible. Hence log|C|=1.

    l = FiniteSize.l_key(; n, α, H_max, λ_EC, ε_snd, ε_KV, ε_PA, f)
    return l/n
end

function finite_key_rate_opt(p_hon, HAgB; n, α=1+1/√n, γ=0.1)
    # Optimize the min-tradeoff function
    function f(grad)
        return -finite_key_rate(p_hon, grad, HAgB; n, α)
    end
    _, grad0 = vN_entropy_lossy(p_hon, γ)
    res = optimize(f, grad0)
    return -Optim.minimum(res)
end

function test()
    q = 0.01
    η_A = 0.9
    η_B = 0.8
    γ = 0.1

    p = p_hon(q, η_A, η_B)
    H = HAgB(q, η_A, η_B)
    @show vN_entropy_lossy(p, γ)[1] - (1-γ)*0.5*HAgB(q, η_A, η_B) - γ*(0.5 - p[2])
    @show finite_key_rate(p, H; n=Int(1e12), γ)
    # @show finite_key_rate_opt(p, H; n=Int(1e12))
end

function asymptotic_key_rate(p_hon, HAgB)
    γ = 0.1
    ent, _ = vN_entropy_lossy(p_hon, γ)
    return ent - (1 - γ)*0.5*HAgB - γ*(0.5 - p_hon[2])
end

end # module BB84_SDI

