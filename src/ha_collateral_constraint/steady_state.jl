"""Re-solve households with constant aggregate endowments and IID individual risk."""
function steady_state_households(m,price;initial=nothing,tol=2e-10,maxiter=6000)
    st=consumption_at_price(price,m.p)
    cash=[resources(m,b,e,1,price) for b in m.b,e in m.e]
    lower=[collateral(m,1,e,price) for e in m.e]
    g=initial===nothing ? repeat(m.b,1,m.ne) : copy(initial)
    for e in 1:m.ne,i in 1:m.nb
        g[i,e]=clamp(g[i,e],lower[e],min(last(m.b),cash[i,e]-1e-8))
    end
    proposal=similar(g)
    for it in 1:maxiter
        lambda=st.K.*((cash.-g)./st.d).^(-m.p.sigma)
        M=m.p.beta*m.p.R.*(lambda*m.pe)
        for e in 1:m.ne,i in 1:m.nb
            proposal[i,e]=current_choice(m,cash[i,e],lower[e],st,
                                         interpolate(m.b,M,g[i,e]))
        end
        err=maximum(abs.(proposal.-g))
        @. g=0.75g+0.25proposal
        if err<tol
            cT=(cash.-g)./st.d
            return (;g,cT,cN=st.q.*cT,lambda=st.K.*cT.^(-m.p.sigma),
                    iterations=it,error=err)
        end
    end
    error("Steady-state household iteration did not converge")
end

"""Find the invariant inherited-asset distribution using the two-node lottery."""
function steady_state_distribution(m,g;initial=nothing,tol=1e-12,maxiter=20000)
    saving=[interpolate(m.b,view(g,:,e),b) for b in m.bf,e in 1:m.ne]
    indices=zeros(Int,m.nf,m.ne); weights=zeros(m.nf,m.ne)
    for e in 1:m.ne,i in 1:m.nf
        indices[i,e],weights[i,e]=bracket(m.bf,saving[i,e])
    end
    F=initial===nothing ? initial_distribution(m;B=(first(m.bf)+last(m.bf))/2) : copy(initial)
    Fnext=similar(F)
    for it in 1:maxiter
        fill!(Fnext,0.0)
        for e in 1:m.ne,i in 1:m.nf
            j=indices[i,e]; w=weights[i,e]; mass=F[i]*m.pe[e]
            Fnext[j]+=(1-w)*mass; Fnext[j+1]+=w*mass
        end
        err=maximum(abs.(Fnext.-F))
        F,Fnext=Fnext,F
        if err<tol
            return (;F,saving,indices,weights,iterations=it,error=err)
        end
    end
    error("Steady-state distribution iteration did not converge")
end

"""
    solve_steady_state(m; yT=1.0, yN=1.0, ...)

Remove aggregate risk, retain the supplied preferences, asset grids and IID
income risk, and jointly solve the stationary policies, distribution and price.
No stochastic forecasting rule or aggregate-state transition is used.
"""
function solve_steady_state(m;yT=1.0,yN=1.0,household_tol=2e-10,
                            distribution_tol=1e-12,market_tol=1e-9,
                            household_maxiter=6000,distribution_maxiter=20000,
                            verbose=true)
    yT+minimum(m.e)>0 && yN>0 || throw(ArgumentError("Endowments must be positive"))
    abs(dot(m.e,m.pe))<1e-14 || throw(ArgumentError("Individual income must have zero mean"))
    # The single aggregate state is absorbing; households are solved on (b, epsilon).
    dm=(;p=m.p,b=m.b,bf=m.bf,e=m.e,pe=m.pe,nb=m.nb,nf=m.nf,ne=m.ne,
          yT=[Float64(yT)],yN=[Float64(yN)],Pi=ones(1,1),nz=1)
    pfloor=max(0.0,maximum((-m.p.R*b-(1+m.p.kappaT)*(yT+e))/
                           ((1+m.p.kappaN)*yN) for b in m.b,e in m.e))
    pceiling=m.p.kappaN>0 ? minimum((-first(m.b)-m.p.kappaT*(yT+e))/
                                   (m.p.kappaN*yN) for e in m.e) : Inf
    pfloor+2e-4<pceiling || error("Asset grid leaves no feasible steady-state price interval")
    p0=(1-m.p.omega)/m.p.omega*(yT/yN)^(1+m.p.eta)
    left=max(0.8p0,pfloor+1e-4); right=min(1.2p0,pceiling-1e-4)
    left<right || error("Cannot initialize the steady-state price bracket")
    previous_g=nothing; previous_F=nothing; evaluations=0
    function evaluate(price)
        hh=steady_state_households(dm,price;initial=previous_g,tol=household_tol,
                                    maxiter=household_maxiter)
        dist=steady_state_distribution(dm,hh.g;initial=previous_F,tol=distribution_tol,
                                       maxiter=distribution_maxiter)
        previous_g=hh.g; previous_F=dist.F; evaluations+=1
        st=consumption_at_price(price,m.p)
        cT=[(resources(dm,b,m.e[e],1,price)-dist.saving[i,e])/st.d
             for (i,b) in enumerate(m.bf),e in 1:m.ne]
        aggregate_cT=dot(dist.F,cT*m.pe)
        (;hh,dist,cT,cN=st.q.*cT,aggregate_cT,
          aggregate_cN=st.q*aggregate_cT,market_error=st.q*aggregate_cT-yN)
    end
    excess_demand(price)=evaluate(price).market_error
    fl=excess_demand(left); fr=excess_demand(right)
    for _ in 1:20
        fl*fr<=0 && break
        left=max(0.8left,pfloor+1e-4); right=min(1.2right,pceiling-1e-4)
        fl=excess_demand(left); fr=excess_demand(right)
    end
    fl*fr<=0 || error("Cannot bracket steady-state nontradable market: residuals=($fl,$fr)")
    price=Roots.find_zero(excess_demand,(left,right),Roots.Brent();xatol=1e-12,xrtol=0.0)
    v=evaluate(price); F=v.dist.F; B=dot(m.bf,F)
    Fnext=zeros(m.nf)
    for e in 1:m.ne,i in 1:m.nf
        j=v.dist.indices[i,e]; w=v.dist.weights[i,e]; mass=F[i]*m.pe[e]
        Fnext[j]+=(1-w)*mass; Fnext[j+1]+=w*mass
    end
    Bnext=dot(m.bf,Fnext)
    binding=0.0; upper_mass=0.0; interior=0.0; inequality=0.0
    M=m.p.beta*m.p.R.*(v.hh.lambda*m.pe)
    for e in 1:m.ne,i in 1:m.nb
        residual=(v.hh.lambda[i,e]-interpolate(m.b,M,v.hh.g[i,e]))/v.hh.lambda[i,e]
        slack=v.hh.g[i,e]-collateral(dm,1,m.e[e],price)
        if slack<1e-7
            inequality=max(inequality,-residual)
        elseif v.hh.g[i,e]<last(m.b)-1e-7
            interior=max(interior,abs(residual))
        end
    end
    for e in 1:m.ne,i in 1:m.nf
        mass=F[i]*m.pe[e]
        binding+=mass*(v.dist.saving[i,e]-collateral(dm,1,m.e[e],price)<1e-7)
        upper_mass+=mass*(v.dist.saving[i,e]>last(m.b)-1e-7)
    end
    budget_error=maximum(abs(v.dist.saving[i,e]+v.cT[i,e]+price*v.cN[i,e]-
                             resources(dm,m.bf[i],m.e[e],1,price))
                         for i in 1:m.nf,e in 1:m.ne)
    distribution_error=maximum(abs.(Fnext.-F))
    resource_error=v.aggregate_cT+Bnext-m.p.R*B-yT
    price_implied=(1-m.p.omega)/m.p.omega*((yT+(m.p.R-1)*B)/yN)^(1+m.p.eta)
    diagnostics=(;market_error=v.market_error,distribution_error,resource_error,
                   stationary_mean_error=Bnext-B,budget_error,mass_error=sum(F)-1,
                   euler_interior=interior,binding_inequality=inequality,upper_mass,
                   household_update=v.hh.error)
    converged=abs(v.market_error)<market_tol && distribution_error<distribution_tol &&
              v.hh.error<household_tol && abs(sum(F)-1)<1e-10
    result=(;m=dm,hh=v.hh,F,price,price_implied,B,Bnext,cT=v.aggregate_cT,cN=v.aggregate_cN,
              binding,variance=dot((m.bf.-B).^2,F),diagnostics,converged,
              price_evaluations=evaluations,distribution_iterations=v.dist.iterations)
    converged || error("Steady-state equilibrium did not meet its stopping conditions")
    verbose && @printf("Deterministic steady state: B = %.6f | pN = %.6f | constrained = %.2f%%\n",
                       B,price,100binding)
    result
end

"""Load a compatible deterministic equilibrium, or solve and save it locally."""
function load_or_solve_steady_state(m,cache_path;force_compute=false,verbose=true,solver_options...)
    root=normpath(joinpath(@__DIR__,"..",".."))
    files=[joinpath(@__DIR__,name) for name in ("model.jl","households.jl","steady_state.jl")]
    append!(files,[joinpath(root,"Project.toml"),joinpath(root,"Manifest.toml")])
    signature=(;format_version=1,model=(;m.p,m.b,m.bf,m.e,m.pe),
                 solver_options=(;solver_options...),julia_version=string(VERSION),
                 source_sha=bytes2hex(SHA.sha256(join(read.(files,String),"\n"))))
    if !force_compute && isfile(cache_path)
        try
            saved=Serialization.deserialize(cache_path)
            if saved.signature==signature && saved.result.converged
                verbose && println("Loaded deterministic steady state from ",cache_path)
                return saved.result
            end
            verbose && println("Saved deterministic steady state is incompatible; recomputing.")
        catch err
            @warn "Could not load deterministic steady state; recomputing" exception=(err,catch_backtrace())
        end
    end
    result=solve_steady_state(m;solver_options...,verbose)
    mkpath(dirname(cache_path))
    Serialization.serialize(cache_path,(;signature,result))
    verbose && println("Saved deterministic steady state to ",cache_path)
    result
end
