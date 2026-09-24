"""Recompute budget, collateral, and Euler/complementarity conditions."""
function diagnostics(m,sol)
    p=m.p
    mu=[marginal_utility(sol.c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
    expected=p.beta*p.R*(mu*m.Pi')
    residual=similar(mu)
    for s in 1:m.ns,i in 1:m.nb
        residual[i,s]=(mu[i,s]-interp(m.b,view(expected,:,s),sol.g[i,s]))/mu[i,s]
    end
    interior=(sol.slack.>1e-7).&(sol.g.>first(m.b)+1e-7).&(sol.g.<last(m.b)-1e-7)
    binding=sol.binding
    common=(;min_consumption=minimum(sol.c),min_collateral_slack=minimum(sol.slack),
      budget_error=maximum(abs.(sol.c+sol.g.-p.R.*m.b.-m.yT')),
      max_policy_update=sol.error,iterations=sol.iterations,
      grid_lower_choices=count(sol.g.<first(m.b)+1e-8),
      grid_upper_choices=count(sol.g.>last(m.b)-1e-8),
      psi_max=maximum([externality(sol.c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]))
    if sol.kind==:DE
        return merge(common,(;max_interior_euler=any(interior) ? maximum(abs.(residual[interior])) : NaN,
          binding_inequality_violation=any(binding) ? max(0.0,-minimum(residual[binding])) : 0.0))
    end
    psi=[externality(sol.c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
    lambda=mu+sol.multiplier.*psi
    expected_sp=p.beta*p.R*(lambda*m.Pi')
    planner_residual=similar(mu)
    for s in 1:m.ns,i in 1:m.nb
        planner_residual[i,s]=(lambda[i,s]-sol.multiplier[i,s]-
            interp(m.b,view(expected_sp,:,s),sol.g[i,s]))/mu[i,s]
    end
    away_from_grid=(sol.g.>first(m.b)+1e-7).&(sol.g.<last(m.b)-1e-7)
    merge(common,(;max_euler_residual=maximum(abs.(planner_residual[away_from_grid])),
        min_multiplier=minimum(sol.multiplier),
        max_complementarity=maximum(abs.(sol.multiplier.*sol.slack)./mu)))
end

"""Euler KKT check between grid points, recomputing future consumption and uT.
Unlike the iteration, this does not interpolate precomputed marginal utilities.
"""
function offgrid_diagnostics(m,de)
    errors=Float64[]; worst_inequality=0.0;worst_slack=0.0
    for s in 1:m.ns,i in 1:m.nb-1
        b=(m.b[i]+m.b[i+1])/2
        v=policy_at(m,de,b,s)
        lhs=marginal_utility(v.c,m.yN[s],m.p)
        rhs=0.0
        for sp in 1:m.ns
            nxt=policy_at(m,de,v.bp,sp)
            rhs+=m.Pi[s,sp]*marginal_utility(nxt.c,m.yN[sp],m.p)
        end
        gap=(lhs-m.p.beta*m.p.R*rhs)/lhs
        worst_slack=min(worst_slack,v.slack)
        if v.binding
            worst_inequality=max(worst_inequality,-gap)
        else
            push!(errors,abs(gap))
        end
    end
    (;maximum=maximum(errors),p99=quantile(errors,0.99),median=median(errors),
      binding_violation=worst_inequality,min_slack=worst_slack,points=m.ns*(m.nb-1))
end

function policy_distance(m,sol,mref,ref)
    maximum(abs(sol.g[i,s]-interp(mref.b,view(ref.g,:,s),m.b[i])) for i in 1:m.nb,s in 1:m.ns)
end

"""Numerical derivative checks of the two key analytical derivatives."""
function derivative_checks(m)
    p=m.p; h=1e-6; uerr=0.0;perr=0.0
    for n in m.yN,c in (0.4,0.8,1.2)
        fd=(utility(c+h,n,p)-utility(c-h,n,p))/(2h)
        uerr=max(uerr,abs(fd/marginal_utility(c,n,p)-1))
        pd=p.kappaN*n*(nt_price(c+h,n,p)-nt_price(c-h,n,p))/(2h)
        perr=max(perr,abs(pd/externality(c,n,p)-1))
    end
    (;marginal_utility_relative_error=uerr,collateral_derivative_relative_error=perr)
end
