"""Simulate by linear interpolation of the converged bond policy.

Binding uses a 1e-5 slack tolerance. Interpolation errors are reported in
validation; no new household optimization or collateral root is solved here.
"""
function policy_at(m,sol,b,s; binding_tol=1e-5)
    first(m.b)-1e-10<=b<=last(m.b)+1e-10 || error("Simulation left the bond grid: $b")
    bp=interp(m.b,view(sol.g,:,s),b)
    c=m.p.R*b+m.yT[s]-bp
    q=nt_price(c,m.yN[s],m.p)
    gap=bp+m.p.kappaT*m.yT[s]+m.p.kappaN*q*m.yN[s]
    (;bp,c,price=q,slack=gap,binding=gap<binding_tol,
      C=basket(c,m.yN[s],m.p),P=price_index(q,m.p),
      Y=m.yT[s]+q*m.yN[s],CA=bp-b,TB=m.yT[s]-c)
end

function shock_path(m,T;seed=2011)
    rng=Xoshiro(seed)
    pi=stationary_shocks(m)
    draw(prob,u)=min(searchsortedfirst(cumsum(prob),u),length(prob))
    s=draw(pi,rand(rng))
    path=Vector{Int}(undef,T)
    cdf=cumsum(m.Pi,dims=2)
    for t in 1:T
        path[t]=s
        s=min(searchsortedfirst(view(cdf,s,:),rand(rng)),m.ns)
    end
    path
end

function simulate(m,sol,states;burn=0,b0=-0.9443)
    T=length(states)
    keys=(:b,:bp,:c,:C,:price,:P,:Y,:CA,:TB,:slack)
    arrays=NamedTuple{keys}(Tuple(zeros(T) for _ in keys))
    binding=falses(T)
    b=b0
    for t in 1:T
        v=policy_at(m,sol,b,states[t])
        arrays.b[t]=b
        for key in keys[2:end]
            getproperty(arrays,key)[t]=getproperty(v,key)
        end
        binding[t]=v.binding
        b=v.bp
    end
    keep=burn+1:T
    trimmed=NamedTuple{keys}(Tuple(getproperty(arrays,k)[keep] for k in keys))
    merge(trimmed,(;binding=binding[keep],states=states[keep],
        CAY=trimmed.CA./trimmed.Y,TBY=trimmed.TB./trimmed.Y))
end

function simulate_pair(m,de,sp;T=80000,burn=5000,seed=2011)
    states=shock_path(m,T+burn;seed)
    d=simulate(m,de,states;burn)
    s=simulate(m,sp,states;burn)
    threshold=mean(d.CA)+std(d.CA)
    (;de=d,sp=s,crisis_de=d.binding .& (d.CA .> threshold),
      crisis_sp=s.binding .& (s.CA .> threshold),threshold,T,burn,seed)
end

function welfare_and_tax(m,de,sp)
    Vde=evaluate_policy(m,de.g)
    Vsp=evaluate_policy(m,sp.g)
    welfare=(Vsp./Vde).^(1/(1-m.p.sigma)).-1
    mu=[marginal_utility(sp.c[i,s],m.yN[s],m.p) for i in 1:m.nb,s in 1:m.ns]
    expected=m.p.beta*m.p.R*(mu*m.Pi')
    tax=zeros(size(mu))
    for s in 1:m.ns,i in 1:m.nb
        if !sp.binding[i,s]
            tax[i,s]=max(0.0,mu[i,s]/interp(m.b,view(expected,:,s),sp.g[i,s])-1)
        end
    end
    (;Vde,Vsp,welfare,tax)
end

function simulation_moments(sim)
    (;std_C=100std(sim.C)/mean(sim.C),std_P=100std(sim.P)/mean(sim.P),
      std_CAY=100std(sim.CAY),std_TBY=100std(sim.TBY),
      corr_C=cor(sim.C,sim.Y),corr_P=cor(sim.P,sim.Y),
      corr_CAY=cor(sim.CAY,sim.Y),corr_TBY=cor(sim.TBY,sim.Y),
      debt=mean(-100sim.b./sim.Y),binding=100mean(sim.binding))
end

function summarize(m,pair,wt)
    d,s=pair.de,pair.sp
    w=[interp(m.b,view(wt.welfare,:,d.states[t]),d.b[t]) for t in eachindex(d.b)]
    tau=[s.binding[t] ? 0.0 : interp(m.b,view(wt.tax,:,s.states[t]),s.b[t]) for t in eachindex(s.b)]
    (;crisis_de=100mean(pair.crisis_de),crisis_sp=100mean(pair.crisis_sp),
      welfare=100mean(w),tax=100mean(tau),de=simulation_moments(d),sp=simulation_moments(s),
      crisis_count_de=count(pair.crisis_de),crisis_count_sp=count(pair.crisis_sp),
      min_slack=min(minimum(d.slack),minimum(s.slack)),
      grid_edge_share_de=mean((d.bp.<first(m.b)+1e-7).|(d.bp.>last(m.b)-1e-7)),
      grid_edge_share_sp=mean((s.bp.<first(m.b)+1e-7).|(s.bp.>last(m.b)-1e-7)))
end

"""Paper Table 2 / footnote 12 counterfactual, at the median DE crisis.

Order DE crisis events by net outflows CA=b'-b (the crisis reversal measure).
Start SP at the SAME DE debt two periods before the selected event and use
the SAME shocks. This differs from conditioning an independent SP history.
"""
function event_study(m,de,sp,pair;before=2,after=3,event_quantile=0.5)
    candidates=[t for t in findall(pair.crisis_de) if t>before && t+after<=pair.T]
    isempty(candidates) && error("No crisis events available")
    ordered=sort(candidates;by=t->pair.de.CA[t])
    0<event_quantile<=1 || error("Event quantile must lie in (0,1]")
    tstar=ordered[ceil(Int,event_quantile*length(ordered))]
    idx=tstar-before:tstar+after
    states=pair.de.states[idx]
    b0=pair.de.b[first(idx)]
    d=simulate(m,de,states;b0)
    s=simulate(m,sp,states;b0)
    pos=before+1
    impact(x,baseline)=(;consumption=100(x.C[pos]/mean(baseline.C)-1),
        current_account=100x.CAY[pos],depreciation=-100(x.P[pos]/mean(baseline.P)-1))
    (;tstar,relative_time=collect(-before:after),de=d,sp=s,
      impact_de=impact(d,pair.de),impact_sp=impact(s,pair.sp),
      median_CA=pair.de.CA[tstar],crisis_events=length(candidates))
end

function comparison_rows(summary,event)
    [
      ("Crisis probability (%)",5.5,0.4,summary.crisis_de,summary.crisis_sp),
      ("Mean debt / GDP (%)",29.2,28.6,summary.de.debt,summary.sp.debt),
      ("Consumption std. (%)",5.9,5.3,summary.de.std_C,summary.sp.std_C),
      ("Real exchange rate std. (%)",7.5,3.4,summary.de.std_P,summary.sp.std_P),
      ("CA / GDP std. (pp)",2.8,0.6,summary.de.std_CAY,summary.sp.std_CAY),
      ("TB / GDP std. (pp)",2.9,0.6,summary.de.std_TBY,summary.sp.std_TBY),
      ("corr(consumption, GDP)",0.83,0.86,summary.de.corr_C,summary.sp.corr_C),
      ("corr(real exchange rate, GDP)",0.79,0.44,summary.de.corr_P,summary.sp.corr_P),
      ("corr(CA / GDP, GDP)",-0.76,-0.05,summary.de.corr_CAY,summary.sp.corr_CAY),
      ("corr(TB / GDP, GDP)",-0.77,-0.16,summary.de.corr_TBY,summary.sp.corr_TBY),
      ("Median-crisis consumption (%)",-16.7,-10.1,event.impact_de.consumption,event.impact_sp.consumption),
      ("Median-crisis CA / GDP (%)",7.8,0.0,event.impact_de.current_account,event.impact_sp.current_account),
      ("Median-crisis depreciation (%)",19.2,1.1,event.impact_de.depreciation,event.impact_sp.depreciation),
    ]
end

function markdown_table(headers,rows;digits=3)
    fmt(x)=x isa AbstractFloat ? @sprintf("%.*f",digits,x) : string(x)
    header="| "*join(headers," | ")*" |\n| "*join(fill("---",length(headers))," | ")*" |\n"
    header*join(["| "*join(fmt.(collect(row))," | ")*" |" for row in rows],"\n")
end

function write_result_csv(path,headers,rows)
    escape(x)="\""*replace(string(x),"\""=>"\"\"")*"\""
    open(path,"w") do io
        println(io,join(escape.(headers),","))
        for row in rows
            println(io,join(escape.(collect(row)),","))
        end
    end
    path
end

function write_run_metadata(path,m,pair,source)
    open(path,"w") do io
        println(io,"Bianchi (2011) Julia replication")
        println(io,"Julia version: ",VERSION)
        println(io,"Run timestamp (Unix seconds): ",time())
        println(io,"Bond grid: ",m.nb," points on [",first(m.b),", ",last(m.b),"]")
        println(io,"Joint income states: ",m.ns)
        println(io,"Parameters: ",m.p)
        println(io,"RNG: Xoshiro; seed: ",pair.seed)
        println(io,"Retained periods: ",pair.T,"; burn-in: ",pair.burn)
        println(io,"DE and SP solver tolerance: 2e-10")
        println(io,"DE outflow crisis threshold: ",pair.threshold)
        println(io,"Original proc_shock.mat SHA-256: ",source.sha256)
        for f in sort(filter(f->endswith(f,".jl"),readdir(@__DIR__)))
            println(io,"Source SHA-256 ",f,": ",bytes2hex(sha256(read(joinpath(@__DIR__,f)))))
        end
    end
    path
end
