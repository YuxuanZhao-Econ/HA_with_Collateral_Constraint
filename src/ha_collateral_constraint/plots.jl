"""Render a compact results table in the notebook."""
function markdown_table(headers, rows; digits=6)
    showval(x) = x isa AbstractFloat ? @sprintf("%.*g", digits, x) : string(x)
    lines = ["| "*join(headers, " | ")*" |",
             "| "*join(fill("---", length(headers)), " | ")*" |"]
    append!(lines, ["| "*join(showval.(collect(row)), " | ")*" |" for row in rows])
    join(lines, "\n")
end

showtable(headers, rows; digits=6) =
    display(MIME"text/html"(), Markdown.parse(markdown_table(headers, rows; digits)))

"""Plot the simulated mean assets, prices, constrained share, and average asset CDF."""
function plot_simulation(m, sim; periods=1:min(500, sim.T))
    fig = Plots.plot(layout=(2,2), size=(1100,720))
    Plots.plot!(fig[1], periods, sim.B[periods];
                label="Actual B", ylabel="Mean assets", xlabel="Date")
    Plots.plot!(fig[2], periods, sim.price[periods];
                label="Clearing price", ylabel="Nontradable price", xlabel="Date")
    Plots.plot!(fig[2], periods, sim.ppred[periods]; label="Forecast", linestyle=:dash)
    Plots.plot!(fig[3], periods, 100 .* sim.binding[periods];
                label=false, ylabel="Constrained households (%)", xlabel="Date")
    Plots.plot!(fig[4], m.bf, cumsum(sim.Faverage);
                label=false, ylabel="Average asset CDF", xlabel="Inherited assets b",
                xlims=(-1.15,-0.35), ylims=(0,1))
    fig
end

"""
Plot saving and tradable consumption at the first saved simulation date with state z.
Hold its full distribution, mean assets, and clearing price fixed; vary b and epsilon.
"""
function plot_household_policies(result, sim, z; assets=range(-1.10,-0.35,length=220),
                                  verbose=true)
    m = result.m
    snapshot_index = findfirst(x -> x.z == z, sim.snapshots)
    isnothing(snapshot_index) && error("No saved simulation snapshot with z = $z")
    ss = sim.snapshots[snapshot_index]
    B = dot(m.bf, ss.F)
    isapprox(B, ss.B; atol=1e-10) || error("Snapshot assets and mean are inconsistent")
    market = clear_market(m, result.hh, result.forecast, ss.F, z)
    policies = current_policies(m, market.M, z, market.price, assets)
    if verbose
        @printf("Fixed policy state: z = %d, B = %.5f, pN = %.5f, constrained share = %.2f%%, simulated date = %d\n",
                z, B, market.price, 100*market.binding, ss.k)
    end

    fig = Plots.plot(layout=(1,2), size=(1100,400))
    colors = [:dodgerblue, :darkorange]
    for e in 1:m.ne
        label = "epsilon = $(m.e[e])"
        color = colors[mod1(e, length(colors))]
        Plots.plot!(fig[1], assets, policies.g[:,e]; label, color,
                    xlabel="Inherited assets b", ylabel="Next assets b'")
        Plots.hline!(fig[1], [collateral(m,z,m.e[e],market.price)];
                     label=false, color, linestyle=:dot)
        Plots.plot!(fig[2], assets, policies.cT[:,e]; label, color,
                    xlabel="Inherited assets b", ylabel="Tradable consumption")
    end
    Plots.plot!(fig[1], assets, assets;
                label="45-degree line", color=:gray, linestyle=:dash)
    Plots.plot!(fig[1]; title="z = $z, B = $(round(B,digits=3))", titlefontsize=10)
    Plots.plot!(fig[2]; title="Constrained share: $(round(100*market.binding,digits=1))%",
                titlefontsize=10)
    fig
end

"""Plot the fitted H or implied log price by z, with observation counts and visited bands."""
function plot_aggregate_rule(result, which)
    which in (:price, :bond) || throw(ArgumentError("which must be :price or :bond"))
    m = result.m; forecast = result.forecast; sim = result.simulation
    ncol = ceil(Int, sqrt(m.nz)); nrow = cld(m.nz, ncol)
    fig = Plots.plot(layout=(nrow,ncol), size=(400*ncol,310*nrow), legend=false,
                     left_margin=10Plots.mm, bottom_margin=4Plots.mm,
                     top_margin=4Plots.mm, right_margin=2Plots.mm)
    for z in 1:m.nz
        visited_B = sim.B[sim.z .== z]
        if !isempty(visited_B)
            lo,hi = quantile(visited_B, [0.01,0.99])
            Plots.vspan!(fig[z], [lo,hi]; color=:lightgray, alpha=0.25, label=false)
        end
        # Sample every polynomial interval as well as the two linear tails.
        breaks = forecast.basis[z].breaks
        curve_B = sort!(unique!(vcat(collect(range(first(m.B),last(m.B);length=600)),
            [B for j in 1:length(breaks)-1 for B in range(breaks[j],breaks[j+1];length=12)])))
        values = which == :price ? [log(price_forecast(forecast,z,B)) for B in curve_B] :
                                   [bond_forecast(forecast,z,B) for B in curve_B]
        Plots.plot!(fig[z], curve_B, values; color=:steelblue, linewidth=2, label=false)
        Plots.vline!(fig[z], [first(breaks),last(breaks)];
                     color=:gray, alpha=0.5, linewidth=0.8, linestyle=:dot, label=false)
        Plots.plot!(fig[z];
                    title="z=$z | n=$(length(visited_B)) | yT=$(round(m.yT[z],digits=3)), yN=$(round(m.yN[z],digits=3))",
                    titlefontsize=8, xtickfontsize=7, ytickfontsize=7,
                    xlims=(first(m.B),last(m.B)),
                    xlabel=z>(nrow-1)*ncol ? "Mean assets B" : "",
                    ylabel=(z-1)%ncol==0 ?
                           (which == :price ? "log forecast price" : "Forecast B next") : "")
    end
    fig
end

"""Plot five impulse responses and the two constrained-household share levels."""
function plot_girf(girf)
    periods = 0:(girf.horizon-1)
    titles = ("Tradable endowment", "Nontradable endowment", "Nontradable price",
              "Tradable consumption", "Next mean assets", "Constrained share")
    units = ("Percent", "Percent", "Percent", "Percent", "Tradable units", "Households (%)")
    fig = Plots.plot(layout=(2,3), size=(1300,760), legend=false,
                     left_margin=8Plots.mm, bottom_margin=8Plots.mm,
                     top_margin=4Plots.mm, right_margin=4Plots.mm)
    for j in 1:5
        Plots.plot!(fig[j], periods, girf.response[j,:];
                    ribbon=1.645 .* girf.mc_se[j,:], fillalpha=0.18,
                    color=:steelblue, linewidth=2, title=titles[j], ylabel=units[j],
                    xlabel=j>3 ? "Periods after impact" : "")
        Plots.hline!(fig[j], [0.0]; color=:gray, linestyle=:dot, linewidth=1)
    end
    Plots.plot!(fig[6], periods, girf.share_base_pct;
                ribbon=1.645 .* girf.share_se_base, fillalpha=0.16,
                color=:gray, linewidth=2, label="Baseline", legend=:topright,
                title=titles[6], ylabel=units[6], xlabel="Periods after impact", ylims=(0,100))
    Plots.plot!(fig[6], periods, girf.share_shock_pct;
                ribbon=1.645 .* girf.share_se_shock, fillalpha=0.16,
                color=:steelblue, linewidth=2, label="Adverse shock")
    fig
end
