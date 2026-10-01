function _draw_state(weights, total, u)
    threshold = u*total
    cumulative = 0.0
    for k in eachindex(weights)
        cumulative += weights[k]
        cumulative > threshold && return k
    end
    lastindex(weights)  # protects against floating-point summation at the last state
end

function _coupled_transition(Pi, z_base, z_shock, u1, u2, u3)
    p_base = view(Pi,z_base,:)
    p_shock = view(Pi,z_shock,:)
    common = min.(p_base,p_shock)
    common_mass = sum(common)
    if u1 < common_mass
        z = _draw_state(common,common_mass,u2)
        return z,z
    end
    remainder = 1-common_mass
    z_base_next = _draw_state(p_base.-common,remainder,u2)
    z_shock_next = _draw_state(p_shock.-common,remainder,u3)
    z_base_next,z_shock_next
end

function _girf_pair(m, hh, forecast, F0, baseline_z, shocked_z, uniforms)
    horizon = size(uniforms,2) + 1
    base_levels = Matrix{Float64}(undef,6,horizon)
    shock_levels = similar(base_levels)
    base_states = Vector{Int}(undef,horizon)
    shock_states = similar(base_states)
    F_base = F0; F_shock = F0
    z_base = baseline_z; z_shock = shocked_z
    for h in 1:horizon
        base = clear_market(m,hh,forecast,F_base,z_base)
        shock = clear_market(m,hh,forecast,F_shock,z_shock)
        base_states[h] = z_base; shock_states[h] = z_shock
        base_levels[:,h] .= (m.yT[z_base],m.yN[z_base],base.price,base.cT,base.Bnext,base.binding)
        shock_levels[:,h] .= (m.yT[z_shock],m.yN[z_shock],shock.price,shock.cT,shock.Bnext,shock.binding)
        F_base = base.Fnext; F_shock = shock.Fnext
        if h < horizon
            z_base,z_shock = _coupled_transition(m.Pi,z_base,z_shock,
                                                    uniforms[1,h],uniforms[2,h],uniforms[3,h])
        end
    end
    (;base_levels,shock_levels,base_states,shock_states)
end

"""
    generalized_irf(result; F0=result.simulation.Faverage, baseline_z=13, shocked_z=12,
                    horizon=40, n_paths=2000, seed=20260929)

Compare paired paths from a common inherited distribution. Maximal coupling and
a fixed-distribution control variate reduce Monte Carlo noise. Response rows are
yT, yN, price, cT (percent changes), and next assets (level difference).
Constrained-household shares are returned separately as population percentages,
with separate standard errors for the baseline and shock paths.
"""
function generalized_irf(result;F0=result.simulation.Faverage,
                             baseline_z=13,shocked_z=12,horizon=40,
                             n_paths=2000,seed=20260929)
    result.converged || error("Solve the aggregate forecasting rules before computing the GIRF")
    m = result.m
    1 <= baseline_z <= m.nz && 1 <= shocked_z <= m.nz || error("Invalid aggregate state")
    horizon >= 2 && n_paths >= 2 || error("Use at least two dates and two paired paths")
    length(F0)==m.nf && all(F0 .>= -1e-12) && isapprox(sum(F0),1.0;atol=1e-9) ||
        error("F0 must be a normalized asset distribution on the model grid")

    # Outcomes with the inherited distribution fixed at F0: one market clearing per z.
    reference = Matrix{Float64}(undef, 6, m.nz)
    for z in 1:m.nz
        v = clear_market(m,result.hh,result.forecast,F0,z)
        reference[:,z] .= (m.yT[z],m.yN[z],v.price,v.cT,v.Bnext,v.binding)
    end

    # Exact conditional state probabilities, starting from each impact state.
    base_prob = zeros(m.nz,horizon); shock_prob = zeros(m.nz,horizon)
    base_prob[baseline_z,1] = 1.0; shock_prob[shocked_z,1] = 1.0
    for h in 2:horizon
        base_prob[:,h] = m.Pi' * view(base_prob,:,h-1)
        shock_prob[:,h] = m.Pi' * view(shock_prob,:,h-1)
    end
    exact_base = reference * base_prob
    exact_shock = reference * shock_prob

    # Maximal coupling preserves both Markov marginals and maximizes reconvergence.
    uniforms = rand(Random.Xoshiro(seed),3,horizon-1,n_paths)
    residual_base = Array{Float64}(undef,6,horizon,n_paths)
    residual_shock = similar(residual_base)
    matched = zeros(Int,horizon)
    for n in 1:n_paths
        draws = view(uniforms,:,:,n)
        pair = _girf_pair(m,result.hh,result.forecast,F0,baseline_z,shocked_z,draws)
        for h in 1:horizon, j in 1:6
            residual_base[j,h,n] = pair.base_levels[j,h]-reference[j,pair.base_states[h]]
            residual_shock[j,h,n] = pair.shock_levels[j,h]-reference[j,pair.shock_states[h]]
        end
        for h in 1:horizon
            matched[h] += pair.base_states[h]==pair.shock_states[h]
        end
    end

    mean_base = exact_base + dropdims(mean(residual_base,dims=3),dims=3)
    mean_shock = exact_shock + dropdims(mean(residual_shock,dims=3),dims=3)
    response = zeros(5,horizon); mc_se = zeros(5,horizon)
    share_base_pct = zeros(horizon); share_shock_pct = zeros(horizon)
    share_se_base = zeros(horizon); share_se_shock = zeros(horizon)
    for h in 1:horizon
        for j in 1:4
            ratio = mean_shock[j,h]/mean_base[j,h]
            response[j,h] = 100*(ratio-1)
            influence = view(residual_shock,j,h,:) .- ratio .* view(residual_base,j,h,:)
            mc_se[j,h] = 100*std(influence)/(sqrt(n_paths)*mean_base[j,h])
        end
        response[5,h] = mean_shock[5,h]-mean_base[5,h]
        mc_se[5,h] = std(view(residual_shock,5,h,:) .- view(residual_base,5,h,:))/sqrt(n_paths)
        share_base_pct[h] = 100*mean_base[6,h]
        share_shock_pct[h] = 100*mean_shock[6,h]
        share_se_base[h] = 100*std(view(residual_base,6,h,:))/sqrt(n_paths)
        share_se_shock[h] = 100*std(view(residual_shock,6,h,:))/sqrt(n_paths)
    end
    (;response,mc_se,share_base_pct,share_shock_pct,share_se_base,share_se_shock,
      mean_base,mean_shock,reference,match_rate=matched./n_paths,horizon,n_paths,seed,
      baseline_z,shocked_z,B0=dot(m.bf,F0))
end


"""Export the five responses, both constrained-share levels, and their Monte Carlo errors."""
function save_girf_csv(girf, path)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "h,yT_pct,yN_pct,price_pct,cT_pct,Bnext_level,share_base_pct,share_shock_pct,se_yT,se_yN,se_price,se_cT,se_Bnext,se_share_base,se_share_shock")
        for i in 1:girf.horizon
            println(io, join((i-1, girf.response[:,i]...,
                              girf.share_base_pct[i], girf.share_shock_pct[i],
                              girf.mc_se[:,i]...,
                              girf.share_se_base[i], girf.share_se_shock[i]), ","))
        end
    end
    path
end
