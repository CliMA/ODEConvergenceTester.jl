module ODEConvergenceTester

import LinearAlgebra
import Logging
import SciMLBase

compute_err_norm_default(u::FT, v::FT) where {FT <: Real} = LinearAlgebra.norm(u - v)
compute_err_norm_default(u, v) = LinearAlgebra.norm(u .- v)

Base.@kwdef struct ConvergenceTestResults{FT,T}
    computed_order::Vector{FT}
    err::Vector{FT}
    dts::Vector{FT}
    expected_order::T
end

"""
    refinement_study(
        problem,
        alg,
        dts::Vector{<:Real};
        compute_err_norm::Function = compute_err_norm_default,
    )

    refinement_study(
        problem,
        alg;
        compute_err_norm::Function = compute_err_norm_default,
        refinement_range::UnitRange = 1:3,
    )

Estimates and reports the convergence rates.

## Arguments
 - `problem` a SciMLBase problem
 - `alg` an algorithm
 - `dts` a vector of timesteps
 - `compute_err_norm = (x,y) -> norm(x,y)` a function for computing
    the norm of the solution
 - `refinement_range = 1:3` a `UnitRange` of refinements.
   `2:4` is finer than `1:3`

## Theory /  analysis

Consider the error for a time-stepping scheme

```
    err = C Δt^p+ H.O.T
    err_k ≈ C Δt_k^p
    err_k/err_m ≈ Δt_k^p/Δt_m^p
    log(err_k/err_m) ≈ log((Δt_k/Δt_m)^p)
    log(err_k/err_m) ≈ p*log(Δt_k/Δt_m)
    log(err_k/err_m)/log(Δt_k/Δt_m) ≈ p
```

Or, put differently, for `refinement_range = 1:3`:

```
            |f3 - f2|   /
    p = log ---------  / log (rf)
            |f2 - f1| /
```

Where
 - `f1` is the finest resolution
 - `f3` is the coarsest resolution
 - `rf` (where `rf` > 1) is the refinement factor

## References

    Roache, P. J. Quantification of Uncertainty in Computational
    Fluid Dynamics. Annu. Rev. Fluid Mech. 29, 123–160 (1997).
    De Vahl Davis, G. Natural convection of air in a square cavity: a
    benchmark solution. Int. J. Num. Methods Fluids 3, 249–264 (1983).

    http://www.grc.nasa.gov/WWW/wind/valid/tutorial/spatconv.html
"""
function refinement_study(
        problem,
        alg,
        dts::Vector{dtType};
        expected_order = SciMLBase.alg_order(alg),
        compute_err_norm::Function = compute_err_norm_default,
        verbose::Bool = false,
        integrator_kwargs...
    ) where {dtType <: Real}

    n_refinements = length(dts)
    # This can be pretty expensive, depending on the problem
    # so let's limit the number of refinements
    @assert 3 ≤ n_refinements "Need at least 3 runs to compute convergence order"

    # Create deepcopy of integrators to prevent accidental mutation
    integrators = map(1:n_refinements) do n
        SciMLBase.init(deepcopy(problem), deepcopy(alg); dt=dts[n], integrator_kwargs...)
    end

    # Refinement factors
    refinement_factor = map(1:(n_refinements - 1)) do i
        dts[i] / dts[i + 1]
    end

    # Compute new t_final such that n_steps[i]*dt[i] = t_final
    t_final₀ = last(problem.tspan)
    n_steps₀ = round(t_final₀ / dts[1])
    n_steps = map(i -> n_steps₀ * 2^(i - 1), 1:n_refinements)
    t_final = map(enumerate(dts)) do (i, dt)
        n_steps[i] * dt
    end

    if verbose
        @info "------ Convergence parameters ------"
        @info "tspan (original)          : $(problem.tspan)"
        @info "nsteps                    : $n_steps"
        @info "refinement factors        : $refinement_factor"
        @info "dts                       : $dts"
        @info "tfinal (rounded)          : $t_final"
        @info "--- Running convergence study... ---"
    end

    for i in 1:n_refinements
        dt = dts[i]
        if verbose
            print("@timing iteration $i:")
            @time Logging.with_logger(Logging.NullLogger()) do
                for n in 1:n_steps[i]
                    SciMLBase.step!(integrators[i], dt)
                end
            end
        else
            Logging.with_logger(Logging.NullLogger()) do
                for n in 1:n_steps[i]
                    SciMLBase.step!(integrators[i], dt)
                end
            end
        end
        u = integrators[i].u
    end
    err = map(1:(n_refinements - 1)) do i
        u = deepcopy(integrators[i].u)
        v = deepcopy(integrators[i + 1].u)
        compute_err_norm(u, v)
    end
    p = map(1:(length(err) - 1)) do i
        log(err[i] / err[i + 1]) / log(max(refinement_factor...))
    end

    FT = eltype(p)
    if verbose
        @info "----- Convergence study output -----"
        @info "Errors                    : $(err)"
        @info "Convergence orders        : $(p)"
        @info "Expected convergence order: $(expected_order)"
        @info "------------------------------------"
    end
    return ConvergenceTestResults(; computed_order=p, err, dts, expected_order)
end

function refinement_study(
        problem,
        alg;
        expected_order = SciMLBase.alg_order(alg),
        refinement_range::UnitRange = 1:3,
        compute_err_norm::Function = compute_err_norm_default,
        verbose::Bool = false,
        integrator_kwargs...,
    )
    t_final₀ = last(problem.tspan)
    # For this analysis, Δt must be ~ 2^n, so we must
    # find a scale that will result in a reasonable timestep:
    # 2^(scale) = t_final
    # i.e., dt_max should scale linearly with t_final, while
    # still being of the form 2^n
    scale = round(log(t_final₀) / log(2))
    buffer = 3 # use a buffer to avoid extremely coarse timesteps
    dt_max = 2^(scale - buffer)
    refinement_factor = 2
    dts = map(refinement_range) do i
        dt_max / (refinement_factor^(i - 1))
    end
    return refinement_study(
        problem,
        alg,
        dts;
        expected_order,
        verbose,
        compute_err_norm,
        integrator_kwargs...,
    )
end

end # module
