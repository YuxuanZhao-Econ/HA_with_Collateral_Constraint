using Plots

const DE_COLOR="#D98E04"
const SP_COLOR="#2878B5"

function plot_defaults()
    default(fontfamily="sans-serif",guidefontsize=11,tickfontsize=10,
        titlefontsize=12,legendfontsize=9,linewidth=2.2,
        background_color=:white,foreground_color=:black,
        gridalpha=0.15,framestyle=:box,margin=5*Plots.mm,dpi=140)
end

function two_policies(x,a,b;kwargs...)
    plt=plot(x,a;color=DE_COLOR,label="Decentralized",linestyle=:solid,kwargs...)
    plot!(plt,x,b;color=SP_COLOR,label="Planner",linestyle=:dash)
    plt
end

function policy_figure(m,de,sp;s=6)
    plot_defaults()
    a=two_policies(m.b,de.g[:,s],sp.g[:,s];
        xlabel="Current bond holdings b",ylabel="Next-period holdings b'",
        title="Borrowing decisions",legend=:topleft,xlims=(-1.02,-0.60),ylims=(-1.04,-0.58))
    plot!(a,m.b,m.b;color=:gray,linewidth=1,linestyle=:dot,label="45-degree line")
    b=two_policies(m.b,de.price[:,s],sp.price[:,s];
        xlabel="Current bond holdings b",ylabel="Nontradable price (T units)",
        title="Collateral price",legend=false,xlims=(-1.02,-0.60))
    plot(a,b;layout=(1,2),size=(1080,430),
        plot_title=@sprintf("Original shock state %d: yT = %.4f, yN = %.4f",s,m.yT[s],m.yN[s]))
end

function distribution_figure(pair)
    plot_defaults()
    lo=min(minimum(pair.de.b),minimum(pair.sp.b))
    hi=max(maximum(pair.de.b),maximum(pair.sp.b))
    edges=range(lo-0.005,hi+0.005,length=75)
    ymax=1.12max(histogram_peak(pair.de.b,edges),histogram_peak(pair.sp.b,edges))
    plt=histogram(pair.de.b;bins=edges,normalize=:pdf,color=DE_COLOR,
        fillalpha=0.30,linealpha=0.6,label="Decentralized",xlabel="Bond holdings b",
        ylabel="Probability density",title="Ergodic bond distribution",size=(880,480),ylims=(0,ymax))
    histogram!(plt,pair.sp.b;bins=edges,normalize=:pdf,color=SP_COLOR,
        fillalpha=0.30,linealpha=0.6,label="Planner")
    plot!(plt;plot_title="$(pair.T) retained periods; common aggregate shocks")
end

function instruments_figure(m,sp,wt;s=6)
    plot_defaults()
    limits=[-collateral_limit(sp.c[i,s],m.yT[s],m.yN[s],m.p) for i in 1:m.nb]
    margin=100 .* (1 .+ sp.g[:,s]./limits)
    a=plot(m.b,100wt.tax[:,s];color=SP_COLOR,label=false,
        xlabel="Current bond holdings b",ylabel="Percent",title="State-contingent debt tax",
        xlims=(-1.02,-0.60),ylims=(0,25))
    b=plot(m.b,margin;color=SP_COLOR,label=false,
        xlabel="Current bond holdings b",ylabel="Percent",title="Tightening of collateral margins",
        xlims=(-1.02,-0.60),ylims=(0,40))
    plot(a,b;layout=(1,2),size=(1080,430),plot_title="Policy instruments at original shock state 6")
end

function welfare_figure(m,wt;s=6)
    plot_defaults()
    plot(m.b,100wt.welfare[:,s];color=SP_COLOR,label=false,xlims=(-1.02,-0.6),
        xlabel="Current bond holdings b",ylabel="Percent of permanent consumption",
        title="Welfare gain from correcting the externality",size=(880,460),
        plot_title="Switch from DE to planner, including transition costs; shock state 6")
end

function crisis_distribution_figure(pair)
    plot_defaults()
    d=100 .* (pair.de.C[pair.crisis_de]./mean(pair.de.C).-1)
    s=100 .* (pair.sp.C[pair.crisis_sp]./mean(pair.sp.C).-1)
    edges=range(min(-40,floor(min(minimum(d),minimum(s))/5)*5),
        max(5,ceil(max(maximum(d),maximum(s))/5)*5),length=31)
    ymax=1.12max(histogram_peak(d,edges),histogram_peak(s,edges))
    plt=histogram(d;bins=edges,normalize=:pdf,color=DE_COLOR,fillalpha=0.3,
        label="Decentralized ($(length(d)) events)",xlabel="Consumption deviation from own long-run mean (%)",
        ylabel="Conditional probability density",title="Consumption during each economy's own crises",size=(900,490),ylims=(0,ymax))
    histogram!(plt,s;bins=edges,normalize=:pdf,color=SP_COLOR,fillalpha=0.3,
        label="Planner ($(length(s)) events)")
    plt
end

function histogram_peak(values,edges)
    counts=zeros(Int,length(edges)-1)
    for x in values
        first(edges)<=x<=last(edges) || error("Histogram would omit observations")
        j=clamp(searchsortedlast(edges,x),1,length(counts))
        counts[j]+=1
    end
    maximum(counts./(length(values).*diff(collect(edges))))
end

function event_figure(event,pair)
    plot_defaults()
    t=event.relative_time
    d,s=event.de,event.sp
    a=two_policies(t,100 .* (d.C./mean(pair.de.C).-1),100 .* (s.C./mean(pair.sp.C).-1);
        title="Consumption",ylabel="Deviation from own mean (%)",xlabel="Years from crisis",legend=:bottomright)
    b=two_policies(t,100d.CAY,100s.CAY;title="Current account / GDP",ylabel="Percent",xlabel="Years from crisis",legend=false)
    c=two_policies(t,-100 .* (d.P./mean(pair.de.P).-1),-100 .* (s.P./mean(pair.sp.P).-1);
        title="Real exchange rate depreciation",ylabel="Percent below own mean",xlabel="Years from crisis",legend=false)
    for pl in (a,b,c)
        vline!(pl,[0];color=:gray,linestyle=:dot,linewidth=1,label=false)
    end
    plot(a,b,c;layout=(1,3),size=(1250,430),plot_title="Median DE crisis: same debt at t = -2 and same shocks")
end

function save_figure(plt,root,name)
    folder=joinpath(root,"figures","bianchi2011")
    mkpath(folder)
    savefig(plt,joinpath(folder,name*".png"))
    savefig(plt,joinpath(folder,name*".svg"))
    plt
end
