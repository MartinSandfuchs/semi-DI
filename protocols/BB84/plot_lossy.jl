module Plot

using CairoMakie, LaTeXStrings
using JLD2

const BB84_SDI = include("BB84_lossy_SDI.jl").BB84_SDI
using .BB84_SDI

fname = "BB84_lossy_rates.jld2"
if isfile(fname)
    data = load(fname)
    etas = data["etas"]
    rates_all = data["rates_all"]
else
    etas = 2 .^ LinRange(0, -10, 20)
    qs = [0.0, 0.05, 0.1, 0.15, 0.2]
    rates_all = Dict((q, rate_lossy.(q, etas)) for q in qs)
    jldsave(fname; etas, rates_all)
end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        title="BB84 with losses, Alice trusted",
        xlabel="Loss (dB)", ylabel="key rate",
        yscale=log10,
    )
    loss_dB = -10*log10.(etas)

    for (q, rates) in sort(rates_all)
        pos_rates = rates[rates .> 0]
        pos_loss = loss_dB[rates .> 0]
        if !isempty(pos_rates)
            lines!(ax, pos_loss, pos_rates, label="q = $(100*q)%")
        end
    end
    axislegend(ax; position=:rt)

    save("plot_BB84_lossy.pdf", fig)
    return fig
end

end # module Plot
