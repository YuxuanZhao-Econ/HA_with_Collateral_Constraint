"""Freeze expected marginal utility at the previous bond choice, as in solve_de.

CRRA and the optimal CES consumption ratio give an exact scalar Euler inverse.
This changes only the evaluation of the scalar equation, not the iteration map.
"""
@inline function current_choice(m,cash,lower,st,M)
    lower>=first(m.b)-1e-10 || error("Economic borrowing limit below individual grid: $lower")
    upper=min(last(m.b),cash-1e-12)
    lower<upper || error("No positive-consumption choice")
    cbound=(cash-lower)/st.d
    cmin=(cash-upper)/st.d
    c=clamp((st.K/M)^(1/m.p.sigma),cmin,cbound)
    cash-st.d*c
end

function marginal_values(m,f,g)
    lam=similar(g)
    for z in 1:m.nz,j in 1:m.nB
        price=price_forecast(f,z,m.B[j]); st=consumption_at_price(price,m.p)
        for e in 1:m.ne,i in 1:m.nb
            c=(resources(m,m.b[i],m.e[e],z,price)-g[i,e,j,z])/st.d
            c>0 || error("Nonpositive household consumption")
            lam[i,e,j,z]=st.K*c^(-m.p.sigma)
        end
    end
    lam
end

function continuation(m,lam,z,Bnext)
    # Forecast extensions outside the aggregate grid are constant. Counts are
    # reported; realized and validation paths must remain strictly inside it.
    j,w=bracket(m.B,clamp(Bnext,first(m.B),last(m.B)))
    M=zeros(m.nb)
    for zp in 1:m.nz,e in 1:m.ne,i in 1:m.nb
        M[i]+=m.p.beta*m.p.R*m.Pi[z,zp]*m.pe[e]*
            ((1-w)*lam[i,e,j,zp]+w*lam[i,e,j+1,zp])
    end
    M
end

function solve_households(m,f;initial=nothing,tol=2e-8,maxiter=6000,damping=0.25)
    g=zeros(m.nb,m.ne,m.nB,m.nz)
    for z in 1:m.nz,j in 1:m.nB,e in 1:m.ne,i in 1:m.nb
        price=price_forecast(f,z,m.B[j])
        lo=collateral(m,z,m.e[e],price)
        cash=resources(m,m.b[i],m.e[e],z,price)
        old=initial===nothing ? m.b[i] : initial[i,e,j,z]
        g[i,e,j,z]=clamp(old,lo,min(last(m.b),cash-1e-8))
    end
    gn=similar(g); err=Inf
    clipped=count(!(first(m.B)<=bond_forecast(f,z,B)<=last(m.B)) for z in 1:m.nz for B in m.B)
    for it in 1:maxiter
        lam=marginal_values(m,f,g)
        for z in 1:m.nz,j in 1:m.nB
            price=price_forecast(f,z,m.B[j]); st=consumption_at_price(price,m.p)
            M=continuation(m,lam,z,bond_forecast(f,z,m.B[j]))
            for e in 1:m.ne,i in 1:m.nb
                emu=interpolate(m.b,M,g[i,e,j,z])
                gn[i,e,j,z]=current_choice(m,resources(m,m.b[i],m.e[e],z,price),
                    collateral(m,z,m.e[e],price),st,emu)
            end
        end
        err=maximum(abs.(gn.-g))
        @. g=(1-damping)*g+damping*gn
        if err<tol
            return (;g,lambda=marginal_values(m,f,g),iterations=it,error=err,forecast_grid_clips=clipped)
        end
    end
    error("Household time iteration did not converge; residual=$err")
end

"""At a trial current price, solve today's choices with future lambda fixed.

The local iteration also freezes M at the previous bond choice. All household
choices are recomputed when a new price is tested; next B remains the forecast.
"""
function current_policies(m,M,z,price,assets;tol=5e-9,maxiter=1200,damping=0.45,initial=nothing)
    st=consumption_at_price(price,m.p)
    g=zeros(length(assets),m.ne)
    cash=similar(g); lo=[collateral(m,z,e,price) for e in m.e]
    for e in 1:m.ne,i in eachindex(assets)
        cash[i,e]=resources(m,assets[i],m.e[e],z,price)
        cash[i,e]>lo[e]+1e-10 || error("Trial price leaves occupied assets infeasible")
        g[i,e]=clamp(initial===nothing ? assets[i] : initial[i,e],lo[e],min(last(m.b),cash[i,e]-1e-8))
    end
    for it in 1:maxiter
        err=0.0
        for e in 1:m.ne,i in eachindex(assets)
            proposal=current_choice(m,cash[i,e],lo[e],st,interpolate(m.b,M,g[i,e]))
            err=max(err,abs(proposal-g[i,e]))
            g[i,e]=(1-damping)*g[i,e]+damping*proposal
        end
        if err<tol
            cT=(cash.-g)./st.d
            return (;g,cT,cN=st.q.*cT,iterations=it,error=err)
        end
    end
    error("Current-period policy iteration failed at price $price")
end
