module Plot

using JLD2
using CairoMakie, LaTeXStrings

const BB84_SDI = include("BB84_SDI.jl").BB84_SDI
using .BB84_SDI

function key_rate(objective, q)
    return vN_entropy_BB84(objective, q) - HAgB(q)
end

fname = "BB84_rates.jld2"
if isfile(fname)
    data = load(fname)
    qs = data["qs"]
    rates_A = data["rates_A"]
    rates_B = data["rates_B"]
else
    qs = LinRange(1e-4, 0.25, 20) |> collect
    rates_A = key_rate.(objective_A, qs)
    rates_B = key_rate.(objective_B, qs)
    jldsave(fname; qs, rates_A, rates_B)
end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        xlabel="Depolarizing probability", ylabel="key rate",
        title="BB84 (SDI)",
    )
    ylims!(ax, 0, 1.0)
    xlims!(ax, 0, last(qs))

    lines!(ax, qs, rates_A; linestyle=:dash)
    scatter!(ax, qs, rates_A; label="Alice trusted")

    lines!(ax, qs, rates_B; linestyle=:dash)
    scatter!(ax, qs, rates_B; label="Bob trusted")

    h_bin(q) = -q*log2(q) - (1-q)*log2(1-q)
    q2 = collect(qs)
    rates_DD = 1 .- 2*h_bin.(q2/2)
    lines!(ax, qs, rates_DD; linestyle=:dash, label="Device dependent")

    axislegend(ax; position=:rt)

    save("plot_BB84_semiDI.pdf", fig)
    return fig
end

end # module Plot
