function clear_market(m,hh,f,F,z;price_tol=2e-9)
    B=dot(m.bf,F)
    first(m.B)<=B<=last(m.B) || error("Realized mean assets outside aggregate grid: $B")
    Bpred=bond_forecast(f,z,B)
    first(m.B)<=Bpred<=last(m.B) || error("On-path B forecast outside aggregate grid: $Bpred")
    M=continuation(m,hh.lambda,z,Bpred)
    active=findall(>(0.0),F); assets=m.bf[active]; weights=F[active]
    previous_policy=Ref{Union{Nothing,Matrix{Float64}}}(nothing)
    function demand(price)
        pol=current_policies(m,M,z,price,assets;initial=previous_policy[])
        previous_policy[]=pol.g
        sum(weights[i]*m.pe[e]*pol.cN[i,e] for i in eachindex(active),e in 1:m.ne)-m.yN[z]
    end
    p0=price_forecast(f,z,B)
    # Feasibility of positive expenditure at the economic borrowing limit.
    pfloor=maximum((-m.p.R*b-(1+m.p.kappaT)*(m.yT[z]+e))/((1+m.p.kappaN)*m.yN[z])
                   for b in assets,e in m.e)
    # The individual grid must include the economic borrowing limit.
    pceiling=minimum((-first(m.b)-m.p.kappaT*(m.yT[z]+e))/(m.p.kappaN*m.yN[z]) for e in m.e)
    left=max(0.8p0,pfloor+0.01,0.02); right=min(1.2p0,pceiling-0.001)
    fl=demand(left); fr=demand(right)
    for _ in 1:12
        fl*fr<=0 && break
        left=max(left*0.8,pfloor+0.001,0.005)
        right=min(right*1.2,pceiling-0.0001)
        fl=demand(left); fr=demand(right)
    end
    fl*fr<=0 || error("Cannot bracket nontradable market: z=$z B=$B residuals=($fl,$fr)")
    price=Roots.find_zero(demand,(left,right),Roots.Brent();xatol=price_tol,xrtol=0.0)
    pol=current_policies(m,M,z,price,assets;initial=previous_policy[])
    Fnext=zeros(m.nf)
    for e in 1:m.ne,i in eachindex(active)
        j,w=bracket(m.bf,pol.g[i,e])
        mass=weights[i]*m.pe[e]
        Fnext[j]+=(1-w)*mass; Fnext[j+1]+=w*mass
    end
    cT=sum(weights[i]*m.pe[e]*pol.cT[i,e] for i in eachindex(active),e in 1:m.ne)
    cN=sum(weights[i]*m.pe[e]*pol.cN[i,e] for i in eachindex(active),e in 1:m.ne)
    Bnext=dot(m.bf,Fnext)
    binding=sum(weights[i]*m.pe[e]*(pol.g[i,e]-collateral(m,z,m.e[e],price)<1e-7)
                for i in eachindex(active),e in 1:m.ne)
    (;price,Fnext,B,Bnext,Bpred,cT,cN,binding,active,pol,M,
      market_error=cN-m.yN[z],resource_error=cT+Bnext-m.p.R*B-m.yT[z],
      variance=dot((m.bf.-B).^2,F))
end

function shock_path(m,T,seed)
    rng=Xoshiro(seed); z=argmin(abs.(m.yT.-1).+abs.(m.yN.-1)); path=zeros(Int,T)
    for t in 1:T
        path[t]=z
        z=min(searchsortedfirst(cumsum(m.Pi[z,:]),rand(rng)),m.nz)
    end
    path
end

function simulate(m,hh,f;T=2400,burn=400,seed=20260925,F0=nothing,keep=0,
                  stress=false,verbose=false)
    F=F0===nothing ? initial_distribution(m;B=mean(m.B)) : copy(F0)
    path=shock_path(m,T+burn,seed)
    B=zeros(T); Bnext=zeros(T); Bpred=zeros(T); price=zeros(T); ppred=zeros(T)
    cT=zeros(T); binding=zeros(T); variance=zeros(T)
    market=zeros(T); resource=zeros(T); mass=zeros(T); endpoint=zeros(T)
    snapshots=NamedTuple[]; Faverage=zeros(m.nf)
    stress_values=fill(-Inf,4); stress_snapshots=Dict{Int,NamedTuple}()
    for t in eachindex(path)
        v=clear_market(m,hh,f,F,path[t])
        if t>burn
            k=t-burn
            B[k]=v.B; Bnext[k]=v.Bnext; Bpred[k]=v.Bpred
            price[k]=v.price; ppred[k]=price_forecast(f,path[t],v.B)
            cT[k]=v.cT; binding[k]=v.binding; variance[k]=v.variance
            market[k]=v.market_error; resource[k]=v.resource_error
            mass[k]=sum(v.Fnext)-1; endpoint[k]=v.Fnext[1]+v.Fnext[end]
            Faverage .+= F./T
            if keep>0 && (k%max(1,div(T,keep))==0)
                push!(snapshots,(;F=copy(F),z=path[t],B=v.B,price=v.price,k))
            end
            if stress
                scores=(-v.price,v.binding,abs(v.price/ppred[k]-1),abs(v.Bnext-v.Bpred))
                for j in eachindex(scores)
                    if scores[j]>stress_values[j]
                        stress_values[j]=scores[j]
                        stress_snapshots[j]=(;F=copy(F),z=path[t],B=v.B,price=v.price,k)
                    end
                end
            end
        end
        F=v.Fnext
        verbose && t%1000==0 && @printf("simulation %d/%d B=%.5f\n",t,length(path),v.B)
    end
    for j in sort(collect(keys(stress_snapshots)))
        s=stress_snapshots[j]
        any(x->x.k==s.k,snapshots) || push!(snapshots,s)
    end
    (;B,Bnext,Bpred,price,ppred,cT,binding,variance,market,resource,mass,endpoint,
      z=path[burn+1:end],Ffinal=F,Faverage,snapshots,seed,T,burn)
end

function fit_forecasts(m,f,sim;ridge=0.002,smooth=0.001)
    logP=copy(f.logP); H=copy(f.H)
    counts=zeros(Int,m.nz)
    for z in 1:m.nz
        ids=findall(==(z),sim.z); counts[z]=length(ids)
        length(ids)>=12 || error("Too few observations in aggregate state $z: $(length(ids))")
        X=zeros(length(ids),m.nB)
        for (row,t) in enumerate(ids)
            j,w=bracket(m.B,sim.B[t]); X[row,j]=1-w; X[row,j+1]=w
        end
        D=zeros(m.nB-2,m.nB)
        for j in 1:m.nB-2
            D[j,j:j+2]=[1.0,-2.0,1.0]
        end
        # Local linear basis avoids polynomial extrapolation in rare states.
        # A disclosed roughness penalty and old-rule prior identify unused nodes.
        A=X'X+ridge*I+smooth*(D'D)
        logP[:,z]=A\(X'*log.(sim.price[ids])+ridge*f.logP[:,z])
        H[:,z]=A\(X'*sim.Bnext[ids]+ridge*f.H[:,z])
        lo,hi=forecast_price_bounds(m,z)
        logP[:,z]=clamp.(logP[:,z],lo,hi)
        H[:,z]=clamp.(H[:,z],first(m.B)+1e-5,last(m.B)-1e-5)
    end
    (;grid=f.grid,logP,H),counts
end

function forecast_distance(m,a,b)
    ep=maximum(abs(log(price_forecast(a,z,B)/price_forecast(b,z,B))) for z in 1:m.nz for B in m.B)
    eb=maximum(abs(bond_forecast(a,z,B)-bond_forecast(b,z,B)) for z in 1:m.nz for B in m.B)
    max(ep,eb)
end

function observed_forecast_distance(sim,a,b)
    maximum(max(abs(log(price_forecast(a,sim.z[t],sim.B[t])/price_forecast(b,sim.z[t],sim.B[t]))),
                abs(bond_forecast(a,sim.z[t],sim.B[t])-bond_forecast(b,sim.z[t],sim.B[t]))) for t in 1:sim.T)
end

function solve_equilibrium(m;T=2400,burn=400,seed=20260925,maxiter=35,
                           damping=0.25,tol=2e-4,initial=nothing,verbose=true)
    f=initial===nothing ? initial_forecast(m) : regrid_forecast(m,initial)
    history=NamedTuple[]; oldg=nothing; previous=Inf
    for it in 1:maxiter
        hhtol=it==1 ? 2e-8 : clamp(0.001previous,2e-8,2e-5)
        hh=solve_households(m,f;initial=oldg,tol=hhtol)
        sim=simulate(m,hh,f;T,burn,seed)
        proposed,counts=fit_forecasts(m,f,sim)
        update=observed_forecast_distance(sim,f,proposed)
        # A candidate stopping point is always rechecked with tight households.
        if (update<tol || it==maxiter) && hhtol>2e-8
            hh=solve_households(m,f;initial=hh.g,tol=2e-8)
            hhtol=2e-8
            sim=simulate(m,hh,f;T,burn,seed)
            proposed,counts=fit_forecasts(m,f,sim)
            update=observed_forecast_distance(sim,f,proposed)
        end
        global_update=forecast_distance(m,f,proposed)
        bpmax=maximum(abs.(sim.Bnext.-sim.Bpred))
        ppmax=maximum(abs.(sim.price./sim.ppred.-1))
        push!(history,(;iteration=it,hh_iterations=hh.iterations,hhtol,update,global_update,bpmax,ppmax,
                        meanB=mean(sim.B),binding=mean(sim.binding),min_count=minimum(counts)))
        verbose && @printf("KS %2d | HH %4d | update %.3e | B err %.3e | price err %.3e | mean B %.5f\n",
            it,hh.iterations,update,bpmax,ppmax,mean(sim.B))
        if update<tol || it==maxiter
            return (;m,forecast=f,hh,simulation=sim,history,converged=update<tol)
        end
        f=(;grid=f.grid,logP=(1-damping).*f.logP.+damping.*proposed.logP,
             H=(1-damping).*f.H.+damping.*proposed.H)
        oldg=hh.g
        previous=update
    end
end
