module Plot

using JLD2
using CairoMakie, LaTeXStrings

const I3322_QKD_SDI = include("I3322_qkd_SDI.jl").I3322_QKD_SDI
using .I3322_QKD_SDI

fname = "I3322_rates.jld2"
if isfile(fname)
    data = load(fname)
    qs = data["qs"]
    ent_A = data["ent_A"]
    ent_B = data["ent_B"]
else
    qs = LinRange(1e-4, 0.25, 20)
    ent_A = vN_entropy_I3322.(objective_A, qs)
    ent_B = vN_entropy_I3322.(objective_B, qs)
    jldsave(fname; qs, ent_A, ent_B)
end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        xlabel="Depolarizing probability", ylabel="key rate",
        title="I3322 (SDI)",
    )
    xlims!(ax, 0, last(qs))
    ylims!(ax, 0, 1)

    h_bin(q) = -q*log2(q) - (1-q)*log2(1-q)
    ec = h_bin.(qs/2)

    lines!(ax, qs, ent_A - ec; linestyle=:dash)
    scatter!(ax, qs, ent_A - ec; label="Alice trusted")

    lines!(ax, qs, ent_B - ec; linestyle=:dash)
    scatter!(ax, qs, ent_B - ec; label="Bob trusted")

    axislegend(ax; position=:rt)

    save("plot_I3322_qkd_semiDI.pdf", fig)
    return fig
end

end # module Plot
