using Test
import ODEConvergenceTester

# These problem configurations were borrowed, and slightly modified, from
# https://github.com/SciML/OrdinaryDiffEq.jl/blob/master/test/algconvergence/ode_ssprk_tests.jl

using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

Random.seed!(100)

f = (u,p,t)->cos(t)
prob_ode_sin = ODEProblem(ODEFunction(f; analytic=(u0,p,t)->sin(t)), 0.,(0.0,1.0))

f = (du,u,p,t)->du[1]=cos(t)
prob_ode_sin_inplace = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->[sin(t)]), [0.], (0.0,1.0))

f = (u,p,t)->sin(u)
prob_ode_nonlinear = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->2*acot(exp(-t)*cot(0.5))), 1.,(0.,0.5))

f = (du,u,p,t)->du[1]=sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(ODEFunction(f;analytic=(u0,p,t)->[2*acot(exp(-t)*cot(0.5))]),[1.],(0.,0.5))

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]

f_ssp = (u,p,t) -> begin
  sin(10t) * u * (1-u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0., 8.))

f_ssp_inplace = (du,u,p,t) -> begin
  @. du = sin(10t) * u * (1-u)
end
test_problem_ssp_inplace = ODEProblem(f_ssp_inplace, rand(3,3), (0., 8.))

all_problems = [
    test_problems_only_time...,
    test_problems_linear...,
    test_problems_nonlinear...,
    test_problem_ssp_inplace,
    test_problem_ssp,
]

all_algos = [Euler(), SSPRK22()]

@testset "ODEConvergenceTester" begin
    tc = ODEConvergenceTester.refinement_study
    for alg in all_algos
        for (i, prob) in enumerate(all_problems)
            integrator = OrdinaryDiffEq.init(prob, alg; dt = Float64(0.1))
            ctr = tc(integrator; refinement_range = 9:13, print_report = false)
            @test last(ctr.p) ≈ OrdinaryDiffEq.alg_order(alg) rtol = 5e-2
        end
    end

    # Make sure printing doesn't break:
    integrator = OrdinaryDiffEq.init(
        first(all_problems),
        first(all_algos);
        dt = Float64(0.1)
    )
    ctr = tc(integrator; refinement_range = 9:13)

    a_err = "Need at least 3 runs to compute convergence order"
    @test_throws AssertionError(a_err) tc(integrator; refinement_range = 1:2)

end

