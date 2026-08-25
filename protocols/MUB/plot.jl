module Plot

using JLD2
using CairoMakie, LaTeXStrings

const MUB_QKD_SDI = include("mub_qkd_SDI.jl").MUB_QKD_SDI
using .MUB_QKD_SDI
const MUB_QKD_DD = include("mub_qkd_DD.jl").MUB_QKD_DD
const MUB_SDP = include("mub_sdp_mothsara.jl").MUB_SDP

function key_rate(objective, q)
    return vN_entropy_MUB(objective, q) - HAgB(q)
end

function key_rate_DD(q)
    h1 = MUB_QKD_DD.vN_entropy(q)
    h2 = HAgB(q)
    return h1 - h2
end

fname = "MUB_rates.jld2"
if isfile(fname)
    data = load(fname)
    qs = data["qs"]
    rates_DD = data["rates_DD"]
    rates_A = data["rates_A"]
    rates_B = data["rates_B"]
    rates_SDP = data["rates_SDP"]
    rates_SDP_big = data["rates_SDP_big"]
else
    qs = LinRange(1e-4, 0.25, 30) |> collect
    rates_DD = key_rate_DD.(qs)
    rates_A = key_rate.(objective_A, qs)
    rates_B = key_rate.(objective_B, qs)
    rates_SDP = [MUB_SDP.key_rate(3; q) for q in qs]
    rates_SDP_big = [MUB_SDP.key_rate(7; q) for q in qs]
    jldsave(fname; qs, rates_A, rates_B, rates_DD, rates_SDP, rates_SDP_big)
end

# function fixup()
#     data = load(fname)
#     qs = data["qs"]
#     rates_DD = data["rates_DD"]
#     rates_A = data["rates_A"]
#     rates_B = data["rates_B"]
# 
#     rates_SDP = [MUB_SDP.key_rate(3; q) for q in qs]
#     rates_SDP_big = [MUB_SDP.key_rate(7; q) for q in qs]
#     jldsave(fname; qs, rates_A, rates_B, rates_DD, rates_SDP, rates_SDP_big)
# end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        xlabel=L"\text{Depolarizing probability } q", ylabel=L"\text{key rate}",
        # title=L"\text{Mutually unbiased bases protocol}",
    )
    xlims!(ax, 0, last(qs))
    ylims!(ax, 0, log2(3))

    lines!(ax, qs, rates_A; linestyle=:dash)
    scatter!(ax, qs, rates_A; label=L"Alice trusted, $d=3$, NPA-MP")

    lines!(ax, qs, rates_B; linestyle=:dash)
    scatter!(ax, qs, rates_B; label=L"Bob trusted, $d=3$, NPA-MP")

    lines!(ax, qs, rates_SDP; linestyle=:dot)
    scatter!(ax, qs, rates_SDP; label=L"Alice trusted, $d=3$, SDP", marker=:diamond)

    lines!(ax, qs, rates_SDP_big; linestyle=:dot)
    scatter!(ax, qs, rates_SDP_big; label=L"Alice trusted, $d=7$, SDP", marker=:diamond)

    h_bin(q) = -q*log2(q) - (1-q)*log2(1-q)
    q2 = collect(qs)
    rates_DD2 = log2(3) .- (1-1/9)*q2*log2(3^2-1) - h_bin.(1 .- q2 + q2/3^2)
    # lines!(ax, qs, rates_DD; linestyle=:dash, label="Device dependent")
    lines!(ax, qs, rates_DD2; linestyle=:dash, label=L"Device dependent, $d=3$")

    axislegend(ax; position=:rt)

    save("plot_mub_qkd_semiDI.pdf", fig)
    return fig
end

end # module Plot
