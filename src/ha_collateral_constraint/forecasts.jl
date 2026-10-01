"""Fixed cubic B-spline basis with endpoint tangent lines outside its support."""
struct SplineBasis
    breaks::Vector{Float64}
    knots::Vector{Float64}
    polynomials::Array{Float64,3} # power, coefficient, interval
end

function spline_derivatives(knots, x)
    n = length(knots)-1
    v = zeros(4,n)
    for i in 1:n
        v[1,i] = knots[i] <= x < knots[i+1] ||
                 (x == last(knots) && knots[i] < x == knots[i+1])
    end
    for degree in 1:3
        w = zeros(4,n-degree)
        for i in axes(w,2)
            left = knots[i+degree]-knots[i]
            right = knots[i+degree+1]-knots[i+1]
            left > 0 && (w[1,i] += (x-knots[i])/left*v[1,i])
            right > 0 && (w[1,i] += (knots[i+degree+1]-x)/right*v[1,i+1])
            for derivative in 1:degree
                left > 0 && (w[derivative+1,i] += degree/left*v[derivative,i])
                right > 0 && (w[derivative+1,i] -= degree/right*v[derivative,i+1])
            end
        end
        v = w
    end
    v
end

function SplineBasis(breaks::AbstractVector)
    grid = Float64.(breaks)
    length(grid) >= 2 && all(diff(grid).>0) || error("Spline breaks must be strictly increasing")
    knots = vcat(fill(first(grid),4),grid[2:end-1],fill(last(grid),4))
    K = length(knots)-4
    polynomials = zeros(4,K,length(grid)-1)
    for j in 1:length(grid)-1
        polynomials[:,:,j] = spline_derivatives(knots,grid[j])./[1,1,2,6]
    end
    SplineBasis(grid,knots,polynomials)
end

"""Read fixed state-specific breaks, and verify the shock labels of their reference path."""
function load_forecast_breaks(m,path)
    rows = [split(line,',') for line in readlines(path)[2:end]]
    length(rows) == m.nz || error("Forecast knot file has the wrong number of states")
    breaks = Vector{Float64}[]
    for (z,row) in enumerate(rows)
        parse(Int,row[1]) == z || error("Unordered forecast states")
        isapprox(parse(Float64,row[2]),m.yT[z];atol=1e-12) &&
        isapprox(parse(Float64,row[3]),m.yN[z];atol=1e-12) || error("Forecast knots use different aggregate shocks")
        grid = parse.(Float64,row[5:end])
        length(grid)==18 && first(m.B)<=first(grid)<last(grid)<=last(m.B) ||
            error("Expected 18 breaks for a 20-coefficient spline within the B domain")
        all(diff(grid).>0) || error("Repeated forecast knots in state $z")
        push!(breaks,grid)
    end
    breaks
end

function spline_design(basis, points)
    K = size(basis.polynomials,2)
    X = zeros(length(points),K)
    grid = basis.breaks
    for (row,B) in enumerate(points)
        j = clamp(searchsortedlast(grid,B),1,length(grid)-1)
        c = view(basis.polynomials,:,:,j)
        x = clamp(B,first(grid),last(grid))-grid[j]
        tail = B-clamp(B,first(grid),last(grid))
        for k in 1:K
            X[row,k] = ((c[4,k]*x+c[3,k])*x+c[2,k])*x+c[1,k] +
                       tail*((3c[4,k]*x+2c[3,k])*x+c[2,k])
        end
    end
    X
end

"""Twenty coefficients per state; continuous cubic pieces and linear tails."""
struct Forecast{M}
    m::M
    basis::Vector{SplineBasis}
    theta::Matrix{Float64}
    coefficients::Vector{Matrix{Float64}}
end

function Forecast(m,basis,theta)
    coefficients = [zeros(4,length(b.breaks)-1) for b in basis]
    for z in 1:m.nz, j in axes(coefficients[z],2), r in 1:4
        coefficients[z][r,j] = dot(view(basis[z].polynomials,r,:,j),view(theta,:,z))
    end
    Forecast(m,basis,Matrix{Float64}(theta),coefficients)
end

@inline function bond_forecast(f::Forecast,z,B)
    grid = f.basis[z].breaks; c = f.coefficients[z]
    j = clamp(searchsortedlast(grid,B),1,size(c,2))
    x = clamp(B,first(grid),last(grid))-grid[j]
    value = ((c[4,j]*x+c[3,j])*x+c[2,j])*x+c[1,j]
    slope = (3c[4,j]*x+2c[3,j])*x+c[2,j]
    value+slope*(B-clamp(B,first(grid),last(grid)))
end

function price_forecast(f::Forecast,z,B)
    m=f.m; cT=m.p.R*B+m.yT[z]-bond_forecast(f,z,B)
    cT>0 || error("Nonpositive implied tradable consumption: z=$z B=$B cT=$cT")
    (1-m.p.omega)/m.p.omega*(cT/m.yN[z])^(1+m.p.eta)
end

function initial_forecast(m;forecast_breaks=nothing)
    grids = isnothing(forecast_breaks) ?
            [collect(range(first(m.B),last(m.B);length=18)) for _ in 1:m.nz] : forecast_breaks
    length(grids)==m.nz || error("One spline support is required for each state")
    basis = SplineBasis.(grids)
    all(size(b.polynomials,2)==20 for b in basis) || error("Expected 20 spline coefficients")
    center = mean(m.B)
    # Greville points reproduce the affine initial rule exactly in a cubic basis.
    theta = hcat([[center+0.85*(mean(b.knots[k+1:k+3])-center) for k in 1:20] for b in basis]...)
    f = Forecast(m,basis,theta)
    for z in 1:m.nz
        violation,_ = forecast_violation(m,f.basis[z],view(theta,:,z),z)
        violation<1e-9 || error("Affine initial forecast is infeasible in state $z")
    end
    f
end

function damp_forecast(a::Forecast,b::Forecast,damping)
    0<damping<=1 || error("Forecast damping must lie in (0,1]")
    all(a.basis[z].breaks==b.basis[z].breaks for z in 1:a.m.nz) || error("Spline knots changed during KS")
    Forecast(a.m,a.basis,(1-damping).*a.theta.+damping.*b.theta)
end

"""Four affine bounds: aggregate grid and current household feasibility at every B."""
function forecast_bounds(m,z,B)
    p=m.p
    pmin=max(1e-8,(-p.R*first(m.b)-(1+p.kappaT)*(m.yT[z]+minimum(m.e)))/
                       ((1+p.kappaN)*m.yN[z]))+1e-7
    pmax=(-first(m.b)-p.kappaT*(m.yT[z]+maximum(m.e)))/(p.kappaN*m.yN[z])-1e-7
    0<pmin<pmax || error("Infeasible forecast price bounds")
    cmin=m.yN[z]*(p.omega*pmin/(1-p.omega))^(1/(1+p.eta))
    cmax=m.yN[z]*(p.omega*pmax/(1-p.omega))^(1/(1+p.eta))
    max(first(m.B)+1e-5,p.R*B+m.yT[z]-cmax),
    min(last(m.B)-1e-5,p.R*B+m.yT[z]-cmin)
end

function quadratic_roots(a,b,c)
    if abs(a)<eps()*max(abs(b),abs(c),1.0)
        return abs(b)>eps()*max(abs(c),1.0) ? [-c/b] : Float64[]
    end
    disc=b*b-4a*c
    disc<0 && return Float64[]
    q=-0.5*(b+copysign(sqrt(disc),b))
    q==0 ? [-b/(2a)] : [q/a,c/q]
end

function forecast_extrema(m,basis,theta)
    points=vcat(first(m.B),last(m.B),basis.breaks)
    for j in 1:length(basis.breaks)-1
        c=view(basis.polynomials,:,:,j)*theta
        width=basis.breaks[j+1]-basis.breaks[j]
        for slope in (0.0,m.p.R), x in quadratic_roots(3c[4],2c[3],c[2]-slope)
            0<x<width && push!(points,basis.breaks[j]+x)
        end
    end
    sort!(unique!(points))
end

function forecast_violation(m,basis,theta,z)
    points=forecast_extrema(m,basis,theta)
    values=spline_design(basis,points)*theta
    violations=[max(lo-values[i],values[i]-hi) for (i,B) in enumerate(points)
                for (lo,hi) in (forecast_bounds(m,z,B),)]
    max(0.0,maximum(violations)),points[violations.>1e-10]
end

"""Solve the small strictly convex constrained least-squares problem in whitened coordinates.

The feasible old rule starts a primal active-set iteration. No clipping of fitted
values is applied: active constraints are part of the coefficient optimization.
"""
function constrained_coefficients(gram,raw,C,lower,upper,initial)
    factor=cholesky(Symmetric(gram))
    A=vcat(C,-C)/factor.U
    b=vcat(lower-C*raw,C*raw-upper)
    scales=sqrt.(sum(abs2,A;dims=2))[:,1]
    A./=scales; b./=scales
    y=factor.U*(initial-raw)
    minimum(A*y-b)>=-1e-9 || error("The old forecast must be feasible before refitting")
    active=Int[]
    for iteration in 1:3000
        if isempty(active)
            direction=-y; multipliers=Float64[]
        else
            decomposition=qr(Matrix(A[active,:]'))
            Q=Matrix(decomposition.Q)[:,1:length(active)]
            R=Matrix(decomposition.R)
            direction=-y+Q*(Q'*y)
            multipliers=R\(Q'*y)
        end
        if norm(direction)<=1e-10*max(1.0,norm(y))
            if isempty(active) || minimum(multipliers)>=-1e-10
                theta=raw+factor.U\y
                maximum(max.(lower-C*theta,C*theta-upper))<=1e-8 || error("Infeasible QP solution")
                return theta
            end
            deleteat!(active,argmin(multipliers))
            continue
        end
        slack=A*y-b; movement=A*direction
        step=1.0; blocking=0
        for i in eachindex(b)
            if movement[i]<-1e-12 && !(i in active)
                limit=max(0.0,slack[i])/(-movement[i])
                if limit<step
                    step=limit; blocking=i
                end
            end
        end
        y+=step*direction
        if blocking!=0
            # A blocking row must add an independent restriction in the metric
            # of the objective; redundant bounds can be skipped at roundoff.
            row=Vector(A[blocking,:])
            if !isempty(active)
                Q=Matrix(qr(Matrix(A[active,:]')).Q)[:,1:length(active)]
                row-=Q*(Q'*row)
            end
            norm(row)>1e-10 || error("Degenerate active forecast constraints")
            push!(active,blocking)
        end
    end
    error("Forecast constrained least squares did not converge")
end

function fit_forecast(m,f::Forecast,sim;ridge=0.002)
    ridge>0 || error("A positive old-rule stabilizer is required")
    theta=copy(f.theta); counts=zeros(Int,m.nz)
    for z in 1:m.nz
        ids=findall(==(z),sim.z); counts[z]=length(ids)
        counts[z]>=12 || error("Too few observations in aggregate state $z: $(counts[z])")
        basis=f.basis[z]
        X=spline_design(basis,sim.B[ids]); N=spline_design(basis,m.B)
        old=view(f.theta,:,z); center=mean(sim.Bnext[ids])
        gram=X'X+ridge*(N'N)
        rhs=X'*(sim.Bnext[ids].-center)+ridge*N'*(N*old.-center)
        raw=gram\rhs .+ center
        violation,_=forecast_violation(m,basis,raw,z)
        fitted=raw
        if violation>1e-10
            points=sort!(unique!(vcat(m.B,forecast_extrema(m,basis,raw))))
            for round in 1:12
                C=spline_design(basis,points)
                bounds=[forecast_bounds(m,z,B) for B in points]
                fitted=constrained_coefficients(gram,raw,C,first.(bounds),last.(bounds),old)
                violation,extra=forecast_violation(m,basis,fitted,z)
                violation<1e-9 && break
                points=sort!(unique!(vcat(points,extra)))
            end
        end
        violation,_=forecast_violation(m,basis,fitted,z)
        violation<1e-9 || error("Interior forecast constraints did not converge in state $z")
        theta[:,z]=fitted
    end
    Forecast(m,f.basis,theta),counts
end
