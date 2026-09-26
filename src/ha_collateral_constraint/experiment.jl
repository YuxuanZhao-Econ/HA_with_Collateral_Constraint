using Serialization, SHA

"""Reproducible two-grid experiment, with source- and configuration-checked caches.

The coarse run starts from initial_forecast. The fine run starts from the coarse
forecast and resolves household policies on the finer grid. Validation is always
recomputed by the notebook, including when a solver cache is reused.
"""
function run_experiment(root;recompute=false,verbose=true)
    specs=(
        (name="baseline",nb=121,nB=41,nf=601,maxiter=140),
        (name="refined",nb=241,nB=81,nf=1001,maxiter=100),
    )
    solver_files=[joinpath(root,"src","ha_collateral_constraint",x)
                  for x in ("model.jl","households.jl","equilibrium.jl")]
    push!(solver_files,joinpath(root,"src","bianchi2011","shocks.jl"))
    source_sha=bytes2hex(sha256(join(read.(solver_files,String),"\n")))
    mkpath(joinpath(root,"tmp","ha_collateral_constraint"))
    results=NamedTuple[]
    for spec in specs
        settings=(;spec,source_sha,julia_version=string(VERSION),T=8000,burn=500,
                    seed=20260925,damping=0.25,tol=2e-4,Bmin=-1.10,Bmax=-0.30,risk=0.05)
        # Include the exact initial forecast, so changing a prior run invalidates
        # downstream caches even if its numerical configuration is unchanged.
        initial=isempty(results) ? nothing : results[end].result.forecast
        io=IOBuffer(); serialize(io,initial)
        key=(;settings,initial_sha=bytes2hex(sha256(take!(io))))
        cache=joinpath(root,"tmp","ha_collateral_constraint",spec.name*".jls")
        payload=!recompute && isfile(cache) ? deserialize(cache) : nothing
        cached=payload!==nothing && payload.key==key
        if cached
            result=payload.result
            verbose && println("Loaded verified solver cache: ",spec.name)
        else
            m=model(;nb=spec.nb,nB=spec.nB,nf=spec.nf,
                      Bmin=settings.Bmin,Bmax=settings.Bmax,risk=settings.risk)
            result=solve_equilibrium(m;T=settings.T,burn=settings.burn,seed=settings.seed,
                damping=settings.damping,tol=settings.tol,maxiter=spec.maxiter,initial,verbose)
            serialize(cache,(;key,result))
        end
        push!(results,(;name=spec.name,result,cached,key))
    end
    results
end

"""A third, mean-only run to address rare-state coverage and collateral kinks."""
function run_coverage_experiment(root,refined;recompute=false,verbose=true)
    base=refined.result
    B=sort(unique(round.(vcat(base.m.B,collect(range(-1.02,-0.82,length=101)));digits=12)))
    m=merge(base.m,(;B,nB=length(B)))
    files=[joinpath(root,"src","ha_collateral_constraint",x) for x in
           ("model.jl","households.jl","equilibrium.jl","coverage.jl")]
    push!(files,joinpath(root,"src","bianchi2011","shocks.jl"))
    source_sha=bytes2hex(sha256(join(read.(files,String),"\n")))
    io=IOBuffer(); serialize(io,base.forecast)
    key=(;source_sha,initial_sha=bytes2hex(sha256(take!(io))),julia_version=string(VERSION),
          T=8000,burn=500,seed=20260925,maxiter=65,damping=0.25,tol=2e-4,
          Bgrid=B,nb=m.nb,nf=m.nf,risk=m.e,parameters=m.p,uniform_dates=100)
    cache=joinpath(root,"tmp","ha_collateral_constraint","coverage.jls")
    payload=!recompute && isfile(cache) ? deserialize(cache) : nothing
    cached=payload!==nothing && payload.key==key
    if cached
        result=payload.result
        verbose && println("Loaded verified solver cache: coverage")
    else
        result=solve_covered_equilibrium(m,base.forecast;T=key.T,burn=key.burn,seed=key.seed,
            maxiter=key.maxiter,damping=key.damping,tol=key.tol,verbose)
        serialize(cache,(;key,result))
    end
    (;name="coverage",result,cached,key)
end
