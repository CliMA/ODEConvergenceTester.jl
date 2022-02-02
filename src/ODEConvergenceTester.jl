module ODEConvergenceTester

import LinearAlgebra
import Logging
import OrdinaryDiffEq

compute_err_norm_default(u::FT, v::FT) where {FT <: Real} = LinearAlgebra.norm(u - v)
compute_err_norm_default(u, v) = LinearAlgebra.norm(u .- v)

Base.@kwdef struct ConvergenceTestResults{FT}
    p::Vector{FT}
    err::Vector{FT}
    dts::Vector{FT}
    p_expected::FT
end

"""
    refinement_study(
        integrator₀,
        dts::Vector{<:Real};
        compute_err_norm::Function = compute_err_norm_default,
    )

    refinement_study(
        integrator₀;
        compute_err_norm::Function = compute_err_norm_default,
        refinement_range::UnitRange = 1:3,
    )

Estimates and reports the convergence rates.

## Arguments
 - `integrator₀` a OrdinaryDiffEq integrator
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
        integrator₀,
        dts::Vector{dtType};
        compute_err_norm::Function = compute_err_norm_default,
        print_report::Bool = true,
    ) where {dtType <: Real}

    n_refinements = length(dts)
    # This can be pretty expensive, depending on the problem
    # so let's limit the number of refinements
    @assert 3 ≤ n_refinements "Need at least 3 runs to compute convergence order"

    # Create deepcopy of integrators to prevent accidental mutation
    integrators = map(1:n_refinements) do n
        deepcopy(integrator₀)
    end

    # Refinement factors
    refinement_factor = map(1:(n_refinements - 1)) do i
        dts[i] / dts[i + 1]
    end

    # Compute new t_final such that n_steps[i]*dt[i] = t_final
    t_final₀ = last(integrator₀.sol.prob.tspan)
    n_steps₀ = round(t_final₀ / dts[1])
    n_steps = map(i -> n_steps₀ * 2^(i - 1), 1:n_refinements)
    t_final = map(enumerate(dts)) do (i, dt)
        n_steps[i] * dt
    end

    if print_report
        @info "------ Convergence parameters ------"
        @info "nsteps                    : $n_steps"
        @info "refinement factors        : $refinement_factor"
        @info "dts                       : $dts"
        @info "tfinal (rounded)          : $t_final"
        @info "--- Running convergence study... ---"
    end

    for i in 1:n_refinements
        dt = dts[i]
        OrdinaryDiffEq.set_proposed_dt!(integrators[i], dt)
        print("@timing iteration $i:")
        @time Logging.with_logger(Logging.NullLogger()) do
            for n in 1:n_steps[i]
                OrdinaryDiffEq.step!(integrators[i], dt)
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
    p_expected = FT(OrdinaryDiffEq.alg_order(integrator₀.sol.alg))
    if print_report
        @info "----- Convergence study output -----"
        @info "Errors                    : $(err)"
        @info "Convergence orders        : $(p)"
        @info "Expected convergence order: $(p_expected)"
        @info "------------------------------------"
    end
    return ConvergenceTestResults(; p, err, dts, p_expected)
end

function refinement_study(
        integrator₀;
        refinement_range::UnitRange = 1:3,
        compute_err_norm::Function = compute_err_norm_default,
        print_report::Bool = true,
    )
    t_final₀ = last(integrator₀.sol.prob.tspan)
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
        integrator₀,
        dts;
        print_report=print_report,
        compute_err_norm=compute_err_norm,
    )
end

end # module
