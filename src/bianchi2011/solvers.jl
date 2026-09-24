"""DE time iteration: lagged prices, lagged bond choices, and Euler inversion.

The numbered blocks correspond to the notebook's algorithm. The collateral
limit is updated algebraically from the previous price; no boundary is solved.
"""
function solve_de(m; tol=1e-11, maxiter=50000, damping=0.2, initial=nothing, verbose=false)
    p=m.p
    g=initial===nothing ? repeat(m.b,1,m.ns) : copy(initial)
    cash=p.R.*m.b .+ m.yT'
    c=cash-g
    price=[nt_price(c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
    cnew=similar(c); gnew=similar(g)
    err=Inf
    for it in 1:maxiter
        # 1. Freeze expected marginal utility at LAST iteration's bond choice.
        uT=[marginal_utility(c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
        expected_grid=p.beta*p.R*(uT*m.Pi')
        for s in 1:m.ns,i in 1:m.nb
            emu=interp(m.b,view(expected_grid,:,s),g[i,s])

            # 2. Compute the temporary borrowing limit from LAST iteration's price.
            b_bound=clamp(-p.kappaT*m.yT[s]-p.kappaN*price[i,s]*m.yN[s],first(m.b),last(m.b))
            c_bound=cash[i,s]-b_bound
            c_min=cash[i,s]-last(m.b)

            # 3. Check the binding choice; otherwise invert current marginal utility.
            f(cc)=marginal_utility(cc,m.yN[s],p)-emu
            if f(c_bound)>=0
                cnew[i,s]=c_bound
            elseif f(c_min)<=0
                cnew[i,s]=c_min  # numerical saving-grid endpoint
            else
                cnew[i,s]=Roots.find_zero(f,(c_min,c_bound),Roots.Brent(); xatol=2e-14, xrtol=0.0)
            end
            gnew[i,s]=cash[i,s]-cnew[i,s]
        end

        # 4. Damp consumption/bonds, then update the equilibrium price.
        err=max(maximum(abs.(cnew-c)),maximum(abs.(gnew-g)))
        @. c=damping*cnew+(1-damping)*c
        @. g=cash-c
        for s in 1:m.ns,i in 1:m.nb
            price[i,s]=nt_price(c[i,s],m.yN[s],p)
        end
        verbose && it%100==0 && @printf("DE %5d | update %.3e\n",it,err)
        err<tol && return solution(m,g;kind=:DE,iterations=it,error=err)
    end
    error("DE time iteration failed to converge: $err")
end

"""Planner time iteration with the same price/bound update as DE.

Expected marginal value includes mu*Psi. Update the multiplier from the
previous marginal value minus its discounted expectation; interior mu is zero.
"""
function solve_planner(m; tol=1e-11,maxiter=100000,damping=0.005,initial=nothing,verbose=false)
    p=m.p
    g=initial===nothing ? copy(solve_de(m).g) : copy(initial)
    cash=p.R.*m.b .+ m.yT'
    c=cash-g
    price=[nt_price(c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
    multiplier=zeros(size(g))
    cnew=similar(c); gnew=similar(g); munew=similar(g)
    err=Inf
    for it in 1:maxiter
        # 1. Form the planner marginal value, then evaluate its expectation at old g.
        lambda=[marginal_utility(c[i,s],m.yN[s],p)+
                multiplier[i,s]*externality(c[i,s],m.yN[s],p) for i in 1:m.nb,s in 1:m.ns]
        expected_grid=p.beta*p.R*(lambda*m.Pi')
        for s in 1:m.ns,i in 1:m.nb
            emu=interp(m.b,view(expected_grid,:,s),g[i,s])

            # 2. Update the temporary bound using the old price, exactly as in DE.
            collateral_bound=-p.kappaT*m.yT[s]-p.kappaN*price[i,s]*m.yN[s]
            b_bound=clamp(collateral_bound,first(m.b),last(m.b))
            c_bound=cash[i,s]-b_bound
            c_min=cash[i,s]-last(m.b)

            # 3. Choose consumption and update the collateral multiplier.
            f(cc)=marginal_utility(cc,m.yN[s],p)-emu
            if f(c_bound)>=0
                cnew[i,s]=c_bound
                if first(m.b)<collateral_bound<last(m.b)
                    munew[i,s]=max(0.0,lambda[i,s]-emu)
                else
                    munew[i,s]=0.0  # a numerical grid bound is not collateral
                end
            elseif f(c_min)<=0
                cnew[i,s]=c_min
                munew[i,s]=0.0
            else
                cnew[i,s]=Roots.find_zero(f,(c_min,c_bound),Roots.Brent(); xatol=2e-14, xrtol=0.0)
                munew[i,s]=0.0
            end
            gnew[i,s]=cash[i,s]-cnew[i,s]
        end

        # 4. Damp consumption/bonds, retain the multiplier update, and recompute prices.
        err=max(maximum(abs.(cnew-c)),maximum(abs.(gnew-g)),maximum(abs.(munew-multiplier)))
        @. c=damping*cnew+(1-damping)*c
        @. g=cash-c
        multiplier .= munew
        for s in 1:m.ns,i in 1:m.nb
            price[i,s]=nt_price(c[i,s],m.yN[s],p)
        end
        verbose && it%500==0 && @printf("SP %5d | update %.3e\n",it,err)
        err<tol && return solution(m,g;kind=:SP,iterations=it,error=err,multiplier)
    end
    error("Planner time iteration failed to converge: $err")
end

"""Evaluate a fixed policy for welfare calculations; this does not choose policies."""
function evaluate_policy(m,g; tol=2e-11,maxiter=10000,initial=nothing)
    U=[utility(m.p.R*m.b[i]+m.yT[s]-g[i,s],m.yN[s],m.p) for i in 1:m.nb,s in 1:m.ns]
    V=initial===nothing ? U/(1-m.p.beta) : copy(initial)
    W=similar(V)
    for _ in 1:maxiter
        EV=V*m.Pi'
        for s in 1:m.ns,i in 1:m.nb
            W[i,s]=U[i,s]+m.p.beta*interp(m.b,view(EV,:,s),g[i,s])
        end
        err=maximum(abs.(W-V))
        V,W=W,V
        err<tol && return V
    end
    error("Policy evaluation failed")
end

