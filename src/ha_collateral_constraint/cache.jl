"""Fingerprint numerical inputs and solver code; plotting changes do not invalidate a solve."""
function equilibrium_cache_signature(m, solver_options)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    files = [
        joinpath(@__DIR__, name) for name in
        ("model.jl", "forecasts.jl", "shocks25.jl", "households.jl", "equilibrium.jl")
    ]
    append!(files, [joinpath(root, "src", "bianchi2011", "shocks.jl"),
                    joinpath(root, "Project.toml"), joinpath(root, "Manifest.toml")])
    source_sha = bytes2hex(SHA.sha256(join(read.(files, String), "\n")))
    (;format_version=2, model=m, solver_options,
      julia_version=string(VERSION), source_sha)
end

"""
    load_or_solve(m, cache_path; force_compute=false, verbose=true, solver_options...)

Load a converged equilibrium with matching model, solver settings, Julia version,
dependencies, and numerical source. Otherwise solve and save a converged result.
The return value is the same as solve_equilibrium.
"""
function load_or_solve(m, cache_path; force_compute=false, verbose=true, solver_options...)
    signature = equilibrium_cache_signature(m, (;solver_options...))
    if !force_compute && isfile(cache_path)
        try
            saved = Serialization.deserialize(cache_path)
            if saved.signature == signature && saved.result.converged
                verbose && println("Loaded solved equilibrium from ", cache_path)
                return saved.result
            end
            verbose && println("Saved equilibrium is incompatible; recomputing.")
        catch err
            @warn "Could not load saved equilibrium; recomputing" exception=(err, catch_backtrace())
        end
    end

    result = solve_equilibrium(m; solver_options..., verbose)
    if result.converged
        mkpath(dirname(cache_path))
        Serialization.serialize(cache_path, (;signature, result))
        verbose && println("Saved solved equilibrium to ", cache_path)
    end
    result
end
