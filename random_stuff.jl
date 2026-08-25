module RandomStuff

using LinearAlgebra

function rand_state(d)
    M = rand(d, d)
    P = M'*M
    return P/tr(P)
end

function compute_q(rho_1, rho_2)
    rho_E = rho_1 + rho_2

    function helper(rho_x, rho_E)
        tmp = rho_E^(-1/2)
        v = eigvals(tmp*rho_x*tmp)
        v_max = max(v...)
        # @show eigvals(v_max*rho_E - rho_x)
        return v_max
    end

    q = max(helper(rho_1, rho_E), helper(rho_2, rho_E))
    return q
end

function vN_entropy(rho_1, rho_2)
    rho_E = rho_1 + rho_2

    HXE = -tr(rho_1*log(rho_1) + rho_2*log(rho_2))
    HE = -tr(rho_E*log(rho_E))
    return (HXE - HE)/log(2)
end

h2(q) = -q*log2(q) - (1-q)*log2(1-q)

function test(d)
    rho_1 = 0.5*rand_state(d)
    rho_2 = 0.5*rand_state(d)

    q = compute_q(rho_1, rho_2)
    HXgE = vN_entropy(rho_1, rho_2)

    @assert HXgE >= -log2(q) - 1e-4
    if !(HXgE >= h2(q) - 1e-5)
        @show q
        @show HXgE
        @show h2(q)
        @assert false "Found a violation of the inequality"
    end
end

function test(d, n)
    for _ in 1:n
        test(d)
    end
end

test(3, 1_000_000)


end # module RandomStuff
