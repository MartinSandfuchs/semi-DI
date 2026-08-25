module Plot

using CairoMakie, LaTeXStrings
using JLD2

function h_bin(q)
    if iszero(q)
        return 0.0
    else
        return -q*log2(q) - (1-q)*log2(1-q)
    end
end

function make_figure()
    fig = Figure(; size=(700, 550))
    # fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        xlabel=L"\text{Depolarizing probability } q",
        ylabel=L"\text{key rate}",
    )
    ylims!(ax, 0, log2(3))
    xlims!(ax, 0, nothing)

    colors = Makie.wong_colors()

    # # protocols = ["BB84", "MUB", "I3322", "CHSH"]
    # protocols = ["BB84", "MUB", "I3322"]
    # for (i, protocol) in enumerate(protocols)
    #     data = load("$protocol/$protocol_rates.jld2")
    #     qs = data["qs"]
    #     rates_A = data["rates_A"]
    #     rates_B = data["rates_B"]
    #     color = colors[i]
    #     scatter!(ax, qs, rates_A; label="BB84, A trusted", color, marker=:diamond)
    #     lines!(ax, qs, rates_A; color, linestyle=:dash)
    #     scatter!(ax, qs, rates_B; label="BB84, B trusted", color, marker=:x)
    #     lines!(ax, qs, rates_B; color, linestyle=:dash)
    # end

    # Plot BB84
    data = load("BB84/BB84_rates.jld2")
    qs = data["qs"]
    rates_A = data["rates_A"]
    rates_B = data["rates_B"]
    color = colors[1]
    scatter!(ax, qs, rates_A; label=L"\text{BB84, A trusted}", color, marker=:diamond)
    lines!(ax, qs, rates_A; color, linestyle=:dash)
    scatter!(ax, qs, rates_B; label=L"\text{BB84, B trusted}", color, marker=:x)
    lines!(ax, qs, rates_B; color, linestyle=:dash)

    # Plot MUB
    data = load("MUB/MUB_rates.jld2")
    qs = data["qs"]
    rates_A = data["rates_A"]
    rates_B = data["rates_B"]
    color = colors[2]
    scatter!(ax, qs, rates_A; label=L"\text{MUB, A trusted}", color, marker=:diamond)
    lines!(ax, qs, rates_A; color, linestyle=:dash)
    scatter!(ax, qs, rates_B; label=L"\text{MUB, B trusted}", color, marker=:x)
    lines!(ax, qs, rates_B; color, linestyle=:dash)

    # Plot I3322
    data = load("I3322/I3322_rates.jld2")
    qs = data["qs"]
    ent_A = data["ent_A"]
    ent_B = data["ent_B"]
    ec = h_bin.(qs/2)
    color = colors[3]
    scatter!(ax, qs, ent_A - ec; label=L"I_{3322},\text{ A trusted}", color, marker=:diamond)
    lines!(ax, qs, ent_A - ec; color, linestyle=:dash)
    scatter!(ax, qs, ent_B - ec; label=L"I_{3322},\text{ B trusted}", color, marker=:x)
    lines!(ax, qs, ent_B - ec; color, linestyle=:dash)

    # Plot CHSH
    data = load("CHSH/CHSH_rates2.jld2")
    qs = data["qs"]
    ent_A = data["ent_A"]
    ent_B = data["ent_B"]
    ec = h_bin.(qs/2)
    color = colors[4]
    scatter!(ax, qs, ent_A - ec; label=L"\text{CHSH, A trusted}", color, marker=:diamond)
    lines!(ax, qs, ent_A - ec; color, linestyle=:dash)
    scatter!(ax, qs, ent_B - ec; label=L"\text{CHSH, B trusted}", color, marker=:x)
    lines!(ax, qs, ent_B - ec; color, linestyle=:dash)

    axislegend(ax; position=:rt)

    save("plot_all.pdf", fig)
    return fig
end

end # module Plot
