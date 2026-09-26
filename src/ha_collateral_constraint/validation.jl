function household_diagnostics(m,hh,f)
    maxinterior=0.0; inequality=0.0; slackmin=Inf; upper=0
    for z in 1:m.nz,j in 1:m.nB
        price=price_forecast(f,z,m.B[j]); st=consumption_at_price(price,m.p)
        M=continuation(m,hh.lambda,z,bond_forecast(f,z,m.B[j]))
        for e in 1:m.ne,i in 1:m.nb
            g=hh.g[i,e,j,z]
            c=(resources(m,m.b[i],m.e[e],z,price)-g)/st.d
            lhs=st.K*c^(-m.p.sigma); rhs=interpolate(m.b,M,g)
            slack=g-collateral(m,z,m.e[e],price); slackmin=min(slackmin,slack)
            residual=(lhs-rhs)/lhs
            if g>last(m.b)-1e-7
                upper+=1
            elseif slack<1e-7
                inequality=max(inequality,-residual)
            else
                maxinterior=max(maxinterior,abs(residual))
            end
        end
    end
    (;maxinterior,inequality,slackmin,upper,policy_update=hh.error,
      aggregate_grid_clips=hh.forecast_grid_clips)
end

function forecast_diagnostics(sim,f)
    berr=sim.Bnext.-sim.Bpred; perr=sim.price./sim.ppred.-1
    # Recursively use predicted B instead of resetting it to observed mean assets.
    pred=sim.B[1]; recursive=zeros(sim.T)
    for t in 1:sim.T
        pred=bond_forecast(f,sim.z[t],pred)
        recursive[t]=pred-sim.Bnext[t]
    end
    (;B_rmse=sqrt(mean(berr.^2)),B_max=maximum(abs.(berr)),
      price_rmse=sqrt(mean(perr.^2)),price_max=maximum(abs.(perr)),
      B_recursive_max=maximum(abs.(recursive)),
      R2_B=1-sum(berr.^2)/sum((sim.Bnext.-mean(sim.Bnext)).^2))
end

function path_diagnostics(m,hh,f,sim)
    maxinterior=0.0; inequality=0.0; slackmin=Inf; budget=0.0; rows=0
    for s in sim.snapshots
        v=clear_market(m,hh,f,s.F,s.z)
        st=consumption_at_price(v.price,m.p)
        for e in 1:m.ne,i in eachindex(v.active)
            s.F[v.active[i]]>1e-12 || continue
            b=m.bf[v.active[i]]; g=v.pol.g[i,e]; c=v.pol.cT[i,e]
            lhs=st.K*c^(-m.p.sigma); rhs=interpolate(m.b,v.M,g)
            slack=g-collateral(m,s.z,m.e[e],v.price)
            slackmin=min(slackmin,slack)
            budget=max(budget,abs(g+c+v.price*v.pol.cN[i,e]-resources(m,b,m.e[e],s.z,v.price)))
            if slack<1e-7
                inequality=max(inequality,(rhs-lhs)/lhs)
            else
                maxinterior=max(maxinterior,abs((lhs-rhs)/lhs))
            end
            rows+=1
        end
    end
    (;maxinterior,inequality,slackmin,budget,rows,
      market=maximum(abs.(sim.market)),resource=maximum(abs.(sim.resource)),
      mass=maximum(abs.(sim.mass)),endpoint=maximum(sim.endpoint))
end

function summarize(result)
    s=result.simulation; f=result.forecast
    (;converged=result.converged,iterations=length(result.history),
      meanB=mean(s.B),sdB=std(s.B),minB=minimum(s.B),maxB=maximum(s.B),
      mean_price=mean(s.price),binding_share=mean(s.binding),
      mean_asset_variance=mean(s.variance),forecast_diagnostics(s,f)...)
end

"""Independent Euler check using next-period clearing on the full next F.

This deliberately replaces both forecasted B' and forecasted future prices by
the actual distribution update and branch-by-branch next-period market clearing.
It measures economically relevant forecast/interpolation error, not merely the
residual of the household equation used by the solver.
"""
function equilibrium_euler_diagnostics(m,hh,f,sim;maxsnapshots=8,masscut=1e-8)
    interior=Float64[]; weights=Float64[]; ineq=0.0; bindingmax=0.0; periods=0
    for s in sim.snapshots[1:min(maxsnapshots,length(sim.snapshots))]
        v=clear_market(m,hh,f,s.F,s.z)
        keep=findall(i->s.F[v.active[i]]>masscut,eachindex(v.active))
        isempty(keep) && continue
        choices=vec(v.pol.g[keep,:]); future=zeros(length(choices))
        for zp in 1:m.nz
            nxt=clear_market(m,hh,f,v.Fnext,zp)
            individual=current_policies(m,nxt.M,zp,nxt.price,choices)
            st=consumption_at_price(nxt.price,m.p)
            for e in 1:m.ne,k in eachindex(choices)
                future[k]+=m.p.beta*m.p.R*m.Pi[s.z,zp]*m.pe[e]*
                    st.K*individual.cT[k,e]^(-m.p.sigma)
            end
        end
        st=consumption_at_price(v.price,m.p)
        for e in 1:m.ne,(ii,i) in enumerate(keep)
            k=ii+(e-1)*length(keep)
            lhs=st.K*v.pol.cT[i,e]^(-m.p.sigma)
            residual=(lhs-future[k])/lhs
            slack=v.pol.g[i,e]-collateral(m,s.z,m.e[e],v.price)
            if slack<1e-7
                ineq=max(ineq,-residual); bindingmax=max(bindingmax,abs(residual))
            else
                push!(interior,abs(residual)); push!(weights,s.F[v.active[i]]*m.pe[e])
            end
        end
        periods+=1
    end
    (;periods,points=length(interior),max_interior=maximum(interior;init=0.0),
      weighted_mean=isempty(weights) ? 0.0 : dot(interior,weights)/sum(weights),
      binding_inequality=ineq)
end

function forecast_guard_diagnostics(m,f,sim)
    active=0
    for t in 1:sim.T
        z=sim.z[t]; j,w=bracket(m.B,sim.B[t]); lo,hi=forecast_price_bounds(m,z)
        hit=false
        for (k,weight) in ((j,1-w),(j+1,w))
            weight>1e-10 || continue
            hit |= f.logP[k,z]<=lo+1e-6 || f.logP[k,z]>=hi-1e-6 ||
                   f.H[k,z]<=first(m.B)+2e-5 || f.H[k,z]>=last(m.B)-2e-5
        end
        active+=hit
    end
    (;share=active/sim.T,minimum_B=minimum(sim.B),maximum_B=maximum(sim.B),
      minimum_forecast_B=minimum(sim.Bpred),maximum_forecast_B=maximum(sim.Bpred))
end
