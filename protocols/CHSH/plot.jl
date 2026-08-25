module Plot

using JLD2
using CairoMakie, LaTeXStrings

const CHSH_QKD_SDI = include("CHSH_SDI.jl").CHSH_QKD_SDI
using .CHSH_QKD_SDI

function vN_entropy_CHSH_DI(q)
    w_max = (2 + sqrt(2))/4
    w = w_max*(1 - q) + q*0.5
    if w <= 0.75
        return 0.0
    end
    h2(p) = -p*log2(p) - (1 - p)*log2(1 - p)
    return 1 - h2(1/2 + √((8*w - 4)^2/4 - 1)/2)
end

fname = "CHSH_rates2.jld2"
if isfile(fname)
    data = load(fname)
    qs = data["qs"]
    ent_A = data["ent_A"]
    ent_B = data["ent_B"]
    ent_DI = data["ent_DI"]
else
    qs = LinRange(1e-5, 0.17, 30)
    ent_A = vN_entropy_CHSH.(objective_A, qs)
    ent_B = vN_entropy_CHSH.(objective_B, qs)
    ent_DI = vN_entropy_CHSH_DI.(qs)
    jldsave(fname; qs, ent_A, ent_B, ent_DI)
end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        xlabel=L"\text{Depolarizing probability } q",
        ylabel=L"\text{key rate}",
        # title="CHSH (SDI)",
    )
    xlims!(ax, 0, last(qs))
    ylims!(ax, 0, 1)

    h_bin(q) = -q*log2(q) - (1-q)*log2(1-q)
    ec = h_bin.(qs/2)

    lines!(ax, qs, ent_DI - ec; label="Fully DI", color="#2a2")

    lines!(ax, qs, ent_A - ec; linestyle=:dash)
    scatter!(ax, qs, ent_A - ec; label="Alice trusted")

    lines!(ax, qs, ent_B - ec; linestyle=:dash)
    scatter!(ax, qs, ent_B - ec; label="Bob trusted")

    axislegend(ax; position=:rt)

    save("plot_CHSH_qkd_semiDI.pdf", fig)
    return fig
end

end # module Plot
