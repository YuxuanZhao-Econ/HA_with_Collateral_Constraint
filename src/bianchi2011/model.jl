"""Unrounded baseline parameters in the supplied main.m, lines 23-28."""
parameters(; sigma=2.0, eta=1/0.83-1, beta=0.906, omega=0.307,
           kappaT=0.3235, kappaN=0.3235, r=0.04) =
    (; sigma, eta, beta, omega, kappaT, kappaN, r, R=1+r)

@inline basket(c, n, p) = (p.omega*c^(-p.eta)+(1-p.omega)*n^(-p.eta))^(-1/p.eta)
@inline function utility(c, n, p)
    C = basket(c, n, p)
    p.sigma == 1 ? log(C) : C^(1-p.sigma)/(1-p.sigma)
end
@inline marginal_utility(c, n, p) =
    p.omega * (p.omega*c^(-p.eta)+(1-p.omega)*n^(-p.eta))^((p.sigma-1)/p.eta-1) * c^(-p.eta-1)
@inline nt_price(c, n, p) = (1-p.omega)/p.omega*(c/n)^(1+p.eta)
@inline collateral_limit(c, t, n, p) = -p.kappaT*t-p.kappaN*nt_price(c,n,p)*n
@inline externality(c, n, p) = p.kappaN*(1+p.eta)*nt_price(c,n,p)*n/c
"""Domestic CES price index in tradable units; a fall is a real depreciation."""
@inline price_index(q, p) = (p.omega^(1/(1+p.eta)) +
    (1-p.omega)^(1/(1+p.eta))*q^(p.eta/(1+p.eta)))^((1+p.eta)/p.eta)

"""Linear interpolation on the uniform bond grid, with no extrapolation."""
@inline function interp(grid, values, x)
    n = length(grid)
    z = clamp((x-first(grid))/(last(grid)-first(grid))*(n-1), 0.0, Float64(n-1))
    i = min(floor(Int,z)+1,n-1)
    w = z-(i-1)
    @inbounds return (1-w)*values[i]+w*values[i+1]
end

function model(; nb=800, bmin=-1.02, bmax=-0.2, p=parameters())
    nb>=3 || error("Need at least three asset points")
    yT,yN,Pi=original_shocks()
    maximum(abs.(sum(Pi,dims=2).-1))<1e-12 || error("Invalid transition probabilities")
    minimum(Pi)>=0 || error("Negative transition probability")
    Pi=Pi./sum(Pi,dims=2)
    b=collect(range(bmin,bmax,length=nb))
    minimum(p.R*bi+t-bmax for bi in b, t in yT)>0 || error("Grid permits nonpositive consumption")
    (;p,b,yT,yN,Pi,nb,ns=length(yT))
end

function stationary_shocks(m; tol=1e-14)
    pi=fill(1/m.ns,m.ns)
    for _ in 1:10000
        nxt=m.Pi'*pi
        maximum(abs.(nxt-pi))<tol && return nxt/sum(nxt)
        pi=nxt
    end
    error("Shock invariant distribution did not converge")
end

function source_check(root)
    file=joinpath(root,"reference","MatlabCode_overborrowingAER_web","proc_shock.mat")
    actual=bytes2hex(sha256(read(file)))
    actual==SHOCK_FILE_SHA256 || error("Original shock file differs from the exported Julia data")
    (;file,sha256=actual,verified=true)
end

function solution(m,g; kind,iterations,error,multiplier=nothing)
    c=[m.p.R*m.b[i]+m.yT[s]-g[i,s] for i in 1:m.nb,s in 1:m.ns]
    price=[nt_price(c[i,s],m.yN[s],m.p) for i in 1:m.nb,s in 1:m.ns]
    slack=[g[i,s]+m.p.kappaT*m.yT[s]+m.p.kappaN*price[i,s]*m.yN[s] for i in 1:m.nb,s in 1:m.ns]
    (;g,c,price,slack,binding=slack.<1e-8,kind,iterations,error,multiplier)
end
