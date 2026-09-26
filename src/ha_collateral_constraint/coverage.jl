"""Expand the regression sample using actual histograms under every aggregate state.

Each date nearest a B grid point and 100 evenly spaced dates contributes one
market-clearing observation per aggregate state. All transitions in the reused
Bianchi chain are strictly positive. These rows receive unit regression weight;
they improve state coverage and are not used as observations for economic moments.
"""
function coverage_training(m,hh,f,sim;uniform_dates=100,selected_dates=nothing)
    all(m.Pi.>0) || error("Coverage experiment assumes positive aggregate transitions")
    length(sim.snapshots)==sim.T || error("Coverage fitting needs every retained histogram")
    stride=max(1,div(sim.T,uniform_dates))
    ids=selected_dates===nothing ?
        unique(vcat([argmin(abs.(sim.B.-B)) for B in m.B],collect(stride:stride:sim.T))) : selected_dates
    B=copy(sim.B); Bnext=copy(sim.Bnext); price=copy(sim.price); z=copy(sim.z)
    for t in ids,zp in 1:m.nz
        v=clear_market(m,hh,f,sim.snapshots[t].F,zp)
        push!(B,v.B); push!(Bnext,v.Bnext); push!(price,v.price); push!(z,zp)
    end
    (;B,Bnext,price,z,T=length(B),extra_rows=length(B)-sim.T,selected_dates=ids)
end

"""Coverage experiment with the same household and price-clearing blocks.

Only the fitting sample changes relative to solve_equilibrium. A separate entry
point keeps the original two-grid benchmark reproducible and directly comparable.
"""
function solve_covered_equilibrium(m,initial;T=8000,burn=500,seed=20260925,
                                   maxiter=65,damping=0.25,tol=2e-4,verbose=true)
    f=regrid_forecast(m,initial); oldg=nothing; previous=Inf
    history=NamedTuple[]; selected_dates=nothing
    for it in 1:maxiter
        hhtol=it==1 ? 2e-8 : clamp(0.001previous,2e-8,2e-5)
        hh=solve_households(m,f;initial=oldg,tol=hhtol)
        sim=simulate(m,hh,f;T,burn,seed,keep=T)
        training=coverage_training(m,hh,f,sim;selected_dates)
        selected_dates=training.selected_dates
        proposed,counts=fit_forecasts(m,f,training)
        update=observed_forecast_distance(training,f,proposed)
        if (update<tol || it==maxiter) && hhtol>2e-8
            hh=solve_households(m,f;initial=hh.g,tol=2e-8); hhtol=2e-8
            sim=simulate(m,hh,f;T,burn,seed,keep=T)
            training=coverage_training(m,hh,f,sim;selected_dates)
            proposed,counts=fit_forecasts(m,f,training)
            update=observed_forecast_distance(training,f,proposed)
        end
        global_update=forecast_distance(m,f,proposed)
        bpmax=maximum(abs.(sim.Bnext.-sim.Bpred))
        ppmax=maximum(abs.(sim.price./sim.ppred.-1))
        push!(history,(;iteration=it,hh_iterations=hh.iterations,hhtol,update,global_update,
                        bpmax,ppmax,meanB=mean(sim.B),binding=mean(sim.binding),
                        min_count=minimum(counts),extra_rows=training.extra_rows))
        verbose && @printf("Coverage %2d | HH %4d | update %.3e | B err %.3e | price err %.3e\n",
                           it,hh.iterations,update,bpmax,ppmax)
        if update<tol || it==maxiter
            # Histograms used to fit the rule are disposable; retain the ordinary
            # path and its average distribution for compact result caches.
            return (;m,forecast=f,hh,simulation=merge(sim,(snapshots=NamedTuple[],)),
                    history,selected_dates,converged=update<tol)
        end
        f=(;grid=f.grid,logP=(1-damping).*f.logP.+damping.*proposed.logP,
             H=(1-damping).*f.H.+damping.*proposed.H)
        oldg=hh.g; previous=update
    end
end
