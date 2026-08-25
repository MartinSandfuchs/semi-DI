module FiniteSize

using JuMP

struct MinTradeoffData
    h_CA::Float64
    delta::Float64
    var::Float64
end

function l_key(; n, α, H_max, λ_EC, ε_snd, ε_KV, ε_PA, f::MinTradeoffData)
    ε_s = 0.5*(ε_snd - ε_KV - ε_PA)
    @assert ε_s > 0

    h_min = h_GEAT(; n, α, eps=0.25*ε_s, p_Ω=(ε_snd - ε_KV), f)
    return n*h_min - H_max - 2*log2(1/ε_PA) - λ_EC - ceil(log2(1/ε_KV)) - 2*log2(2/(ε_s/4)^2)
end

function h_GEAT(; n, α, eps, p_Ω, f::MinTradeoffData)
    # dA = 3*(3*3) # Why did I choose that?
    dA = 3

    V = log2(2*dA^2 + 1) + sqrt(2 + f.var)
    g = 2*log2(sqrt(2)/eps)
    K_exp = 2*log2(dA) + f.delta

    K1 = (2-α)^3/(6*(3 - 2α)^3*log(2))*2^((α-1)/(2-α)*K_exp)
    K2 = (K_exp*log(2) + log(1 + exp(2)/2^K_exp))^3
    K = K1*K2

    return f.h_CA - log(2)/2*(α-1)/(2-α)*V^2 - (g + α*log2(1/p_Ω))/(n*(α - 1)) - ((α-1)/(2-α))^2*K
end

end # module FiniteSize
