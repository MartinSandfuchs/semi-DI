module Plot

using JLD2
using CairoMakie
using LinearAlgebra
using LaTeXStrings
using Printf

const BB84_SDI = include("BB84_lossy_SDI.jl").BB84_SDI
using .BB84_SDI

fname = "BB84_lossy_finite_rates.jld2"
if isfile(fname)
    data = load(fname)
    etas = data["etas"]
    rates_all = data["rates_all"]
else
    etas = 10 .^ LinRange(0, -2.6, 50)
    η_B = 0.8
    qs = [0.01, 0.05]
    ns = [1e9, 1e10, 1e11]

    rates_all = Dict{Tuple{Float64, Float64}, Vector{Float64}}()
    for q in qs
        for n in ns
            @show q, n
            GC.gc()

            rates = zeros(length(etas))
            Threads.@threads for i in 1:length(etas)
                p_hon = BB84_SDI.p_hon(q, etas[i], η_B)
                HAgB = BB84_SDI.HAgB(q, etas[i], η_B)
                rate = finite_key_rate_opt(p_hon, HAgB; n)
                # rate = finite_key_rate(p_hon, HAgB; n)
                rates[i] = rate

                # p_hon = BB84_SDI.p_hon(q, etas[i], η_B)
                # if i == length(etas)
                #     p_foot = p_hon
                # else
                #     p_foot = BB84_SDI.p_hon(q, etas[i + 1], η_B)
                # end
                # rate = finite_key_rate(p_hon, p_foot; n)
                # push!(rates, rate)
            end
            rates_all[(q, n)] = rates
        end
        rates_asymptotic = zeros(length(etas))
        for i in 1:length(etas)
            p_hon = BB84_SDI.p_hon(q, etas[i], η_B)
            HAgB = BB84_SDI.HAgB(q, etas[i], η_B)
            rate = BB84_SDI.asymptotic_key_rate(p_hon, HAgB)
            rates_asymptotic[i] = rate
        end
        rates_all[(q, Inf)] = rates_asymptotic
    end
    jldsave(fname; etas, rates_all)
end

function make_figure()
    fig = Figure()
    ax = fig[1, 1] = Axis(
        fig;
        # title = "Qubit BB84 with losses",
        xlabel = L"\text{loss } \eta_A \text{ (dB)}",
        ylabel = L"\text{key rate}",
        yscale = log10,
    )
    ylims!(1e-7, 10)
    xlims!(0, 25)

    loss_dB = -10*log10.(etas)

    rates_all_vec = collect(rates_all)
    rates_all_sorted = sort(rates_all_vec; by=x -> reverse(first(x)))
    for ((q, n), rates) in rates_all_sorted
        if n < Inf
            rates = max.(rates, floatmin(Float64))
            q_str = @sprintf "%.6g"  100*q
            lines!(ax, loss_dB, rates; label=L"q=%$(q_str)\%,\,n=10^{%$(log10(n) |> Int)}")
        else
            rates = max.(rates, floatmin(Float64))
            q_str = @sprintf "%.6g"  100*q
            lines!(ax, loss_dB, rates; label=L"q=%$(q_str)\%,\,\text{asymptotic}", linestyle=:dash)
        end
    end

    axislegend(ax; position=:rt, nbanks=2)

    save("plot_BB84_lossy_finite.pdf", fig)
    return fig
end

end # module Plot
