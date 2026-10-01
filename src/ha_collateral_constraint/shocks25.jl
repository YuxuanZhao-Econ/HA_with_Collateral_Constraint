"""Gaussian-Hermite extension of the supplied 4×4 joint endowment chain."""
function gh_joint_shocks(n::Int)
    n >= 2 || error("Need at least two nodes per aggregate shock dimension")

    # The supplied 16 states are an affine image of a 4×4 standard-normal
    # Gaussian-Hermite grid. Infer that affine map from the original nodes.
    original_yT, original_yN, _ = original_shocks()
    x4, _ = gh_normal_nodes_weights(4)
    X4 = reduce(vcat, ([x4[i] x4[j]] for i in 1:4 for j in 1:4))
    M = (X4 \ hcat(original_yT .- 1, original_yN .- 1))'
    S = M * M'
    Sinv = inv(Symmetric(S))

    # Column-vector convention. These coefficients reproduce the original
    # 16×16 transition matrix to floating-point precision with the formula
    # below; the three-decimal VAR printed in the paper is insufficient.
    A = [0.901009 -0.453210; 0.495189 0.225130]

    x, w = gh_normal_nodes_weights(n)
    X = reduce(vcat, ([x[i] x[j]] for i in 1:n for j in 1:n))
    W = [w[i] * w[j] for i in 1:n for j in 1:n]
    Y = ones(n * n, 2) + X * M'
    Pi = Matrix{Float64}(undef, n * n, n * n)
    scores = Vector{Float64}(undef, n * n)
    for i in 1:n*n
        conditional_mean = A * (view(Y, i, :) .- 1)
        for j in 1:n*n
            xj = view(Y, j, :) .- 1
            residual = xj - conditional_mean
            scores[j] = log(W[j]) - 0.5 * dot(residual, Sinv * residual) +
                        0.5 * dot(xj, Sinv * xj)
        end
        scores .-= maximum(scores)
        Pi[i, :] .= exp.(scores)
        Pi[i, :] ./= sum(view(Pi, i, :))
    end
    vec(Y[:, 1]), vec(Y[:, 2]), Pi
end

"""Nodes and probabilities for integration against a standard normal."""
function gh_normal_nodes_weights(n::Int)
    J = SymTridiagonal(zeros(n), sqrt.(Float64.(1:n-1)))
    E = eigen(J)
    E.values, vec(E.vectors[1, :].^2)
end
