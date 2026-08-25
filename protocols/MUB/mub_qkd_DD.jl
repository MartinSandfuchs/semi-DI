module MUB_QKD_DD

export vN_entropy, cond_prob

using FastGaussQuadrature, DataStructures, LinearAlgebra
using JuMP, MosekTools, Dualization

const MatrixNPA = include("../../NPA-MP/src/MatrixNPA.jl").MatrixNPA
using .MatrixNPA

const T = ComplexF64
const dA = 3

function bases()::Vector{Vector{Vector{ComplexF64}}}
    ξ = exp(2*π*im/3)
    b = []
    z = ComplexF64[1, 0, 0]
    o = ComplexF64[0, 1, 0]
    t = ComplexF64[0, 0, 1]

    push!(b, [z, o, t])
    push!(b, [(z + o + t)/√(3), (z + ξ^2*o + ξ*t)/√(3), (z + ξ*o + ξ^2*t)/√(3)])
    push!(b, [(z + o + ξ*t)/√(3), (z + ξ^2*o + ξ^2*t)/√(3), (z + ξ*o + t)/√(3)])
    push!(b, [(z + o + ξ^2*t)/√(3), (z + ξ^2*o + t)/√(3), (z + ξ*o + ξ*t)/√(3)])

    return b
end

function get_quad(m::Int)
	ts, ws = gaussradau(m)
    return 0.5*(1 .- reverse(ts)), 0.5*reverse(ws)
end

function cond_prob(; x, y, a, b, q)
    z = [1, 0, 0]
    o = [0, 1, 0]
    t = [0, 0, 1]
    psi = (kron(z, z) + kron(o, o) + kron(t, t))/sqrt(3)
    rho = (1 - q)*psi*psi' + q*Matrix(I, 3^2, 3^2)/3^2
    Bs = bases()
    v1 = Bs[x][a]
    v2 = Bs[y][b]
    return tr(kron(v1*v1', v2*v2') * rho)
end

function vN_entropy(q::Number)
    Bs = bases()
    Ms = map(vec -> map(v -> v*v', vec), Bs)
    Ns = map(vec -> map(v -> v*v', vec), Bs)

    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    optimizer = dual_optimizer(optimizer)
    model = Model(optimizer)

    m = 12
    ts, ws = get_quad(m)
    id_B = Matrix(I, dA, dA)

    @variable(model, sigma[1:dA^2, 1:dA^2] in HermitianPSDCone())
    @constraint(model, real(tr(sigma)) == 1)

    obj = 0
    for i in 1:m
        t_i = ts[i]
        w_i = ws[i]
        c_i = w_i/(t_i*log(2))
        for a in 1:dA
            zeta_ia = @variable(model, [1:dA^2, 1:dA^2], set=ComplexPlane())
            eta_ia = @variable(model, [1:dA^2, 1:dA^2] in HermitianPSDCone())
            theta_ia = @variable(model, [1:dA^2, 1:dA^2] in HermitianPSDCone())

            Gamma_1 = [sigma zeta_ia; zeta_ia' eta_ia]
            @constraint(model, Hermitian(Gamma_1) in HermitianPSDCone())

            Gamma_2 = [sigma zeta_ia'; zeta_ia theta_ia]
            @constraint(model, Hermitian(Gamma_2) in HermitianPSDCone())

            obj_i = tr(kron(Ms[1][a], id_B)*sigma) + tr(kron(Ms[1][a], id_B)*(zeta_ia + zeta_ia' + (1 - t_i)*eta_ia) + t_i*theta_ia)
            obj += c_i*obj_i
        end
    end

    # Add measurement constraints
    for x in 1:4, y in 1:4, a in 1:3, b in 1:3
        p = cond_prob(; x, y, a, b, q)
        mmt = kron(Ms[x][a], Ns[y][b])
        @constraint(model, real(tr(sigma*mmt) - p) >= 0)
        @constraint(model, real(p - tr(sigma*mmt)) >= 0)
    end

    @objective(model, Min, real(obj))

    optimize!(model)
    ent = objective_value(model)
    return ent
end

end # module MUB_QKD_DD

