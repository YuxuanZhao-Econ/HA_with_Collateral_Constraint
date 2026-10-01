"""Bianchi structural parameters; idiosyncratic income risk is specified separately."""
parameters(;sigma=2.0,eta=1/0.83-1,beta=0.906,omega=0.307,
           kappaT=0.3235,kappaN=0.3235,r=0.04) =
    (;sigma,eta,beta,omega,kappaT,kappaN,r,R=1+r)

income_risk(scale=0.05) = scale==0 ? ([0.0],[1.0]) : ([-scale,scale],[0.5,0.5])

function model(;nb=121,nB=17,nf=401,bmin=-1.35,bmax=1.0,
               Bmin=-1.20,Bmax=-0.30,risk=0.05,p=parameters(),shocks=:original)
    yT,yN,Pi = if shocks == :original
        original_shocks()
    elseif shocks == :gh5
        gh_joint_shocks(5)
    elseif shocks == :toy
        ([0.98,1.02],[1.0,1.0],[0.9 0.1;0.1 0.9])
    else
        error("Unknown aggregate shock discretization: $shocks")
    end
    Pi=Pi./sum(Pi,dims=2)
    e,pe=income_risk(risk)
    minimum(yT)+minimum(e)>0 || error("Nonpositive tradable endowment")
    abs(dot(e,pe))<1e-14 || error("IID income must have zero mean")
    # Concentrate policy nodes around the debt range without truncating savers.
    nlow=round(Int,0.8nb)
    b=vcat(collect(range(bmin,-0.45,length=nlow)),
           collect(range(-0.45,bmax,length=nb-nlow+1))[2:end])
    bf=collect(range(bmin,bmax,length=nf))
    B=collect(range(Bmin,Bmax,length=nB))
    (;p,b,bf,B,e,pe,yT,yN,Pi,nb=length(b),nf,nB,ne=length(e),nz=length(yT))
end

@inline function bracket(grid,x)
    first(grid)-1e-9 <= x <= last(grid)+1e-9 || error("Interpolation outside grid: $x not in $(extrema(grid))")
    j=clamp(searchsortedlast(grid,x),1,length(grid)-1)
    w=(clamp(x,first(grid),last(grid))-grid[j])/(grid[j+1]-grid[j])
    j,w
end

@inline function interpolate(grid,v,x)
    j,w=bracket(grid,x)
    (1-w)*v[j]+w*v[j+1]
end

"""Static intratemporal allocation; lambda=K*cT^(-sigma)."""
@inline function consumption_at_price(price,p)
    q=((1-p.omega)/(p.omega*price))^(1/(1+p.eta))
    d=1+price*q
    a=abs(p.eta)<1e-12 ? q^(1-p.omega) :
      (p.omega+(1-p.omega)*q^(-p.eta))^(-1/p.eta)
    K=p.omega*a^(1+p.eta-p.sigma)
    (;q,d,K)
end

@inline collateral(m,z,e,price) = -m.p.kappaT*(m.yT[z]+e)-m.p.kappaN*price*m.yN[z]
@inline resources(m,b,e,z,price) = m.p.R*b+m.yT[z]+e+price*m.yN[z]

function initial_distribution(m;B=-0.8)
    F=zeros(m.nf)
    j,w=bracket(m.bf,B)
    F[j]=1-w; F[j+1]=w
    F
end
