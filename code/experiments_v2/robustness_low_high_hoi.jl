include("src_functions.jl")

using HomotopyContinuation, Serialization, Statistics, LinearAlgebra

const nsim      = 2
const n         = 2
const d         = 2            # polynomial order
const n_alphas  = 50           # number of alphas to test
const n_pert    = 25           # number of perturbation directions
const k_ratio   = 100.0
const pert_size = 10.0
const d_pert    = 0            # perturb growth‐rates
const param_dir = joinpath(@__DIR__, "parameter_sets", "n_$(n)")

# 0) Pre‐generate all parameter sets (one per seed)
isdir(param_dir) || mkpath(param_dir)
# U.generate_and_save_parameter_sets(nsim, n, n, n, param_dir;
#                                    lower = -1.0, upper = 1.0,
#                                    prefix = "seed")

# 1) HomotopyContinuation setup
@var x[1:n]
@var α[1:2]
pert_dirs   = U.points_hypersphere(n, 1.0, n_pert, false)
ref_eqs, coeffs_mat = U.get_ref_polynomials(x, d, n; coeff_name=:c)
init_sol    = fill(1.0 + 0im, n)

# Path to store pairwise distances
pairwise_path = joinpath("simulation_results",
                         "solution_pairwise_distances_n_$(n).csv")


# 2) Open CSV and loop
open("simulation_results/boundary_distances_n_$(n)_many_alphas_many_perts.csv","w") do io
  # write header
  println(io, "seed,pert_i,alpha,δ_p,δ_x,flag")

  for seed_i in 2:nsim
    seed_i = 3
    let filepath = joinpath(param_dir, "seed_seed$(seed_i).bin")
      # load parameters
      (r, A0, B0) = open(filepath, "r") do fio
        deserialize(fio)
      end

      # precompute for baseline
      invA0 = inv(transpose(A0))

      # compute the two HOI alphas
      σA, σB = std(vec(A0)), std(vec(B0))
      α0 = U.alpha_k_ratio(σA, σB, k_ratio)
      α1 = U.alpha_k_ratio(σA, σB, 1.0)

      # build the generic continuation system once
      num_eqs      = U.build_glvhoi((r,A0,B0), x)
      syst_generic = U.get_parametrized_system(num_eqs,
                                               ref_eqs,
                                               coeffs_mat,
                                               d_pert,
                                               x,
                                               α)

      # now loop over alpha = 0 (baseline), α0, and α1
      #alphavec = vcat(0.0, range(α0, α1; length=n_alphas))
      alphavec = range(0.0, 1.0; length=n_alphas)
      for αval in alphavec
        if αval == 0.0
          # --- analytic‐GLV baseline ---
          tmaxs = U.baseline_robustness_GLV(transpose(A0), r, pert_dirs, pert_size)
          for pert_i in 1:size(pert_dirs,1)
            v = pert_dirs[pert_i, :]
            δp = tmaxs[pert_i]
            # reconstruct the critical equilibrium
            y = invA0 * v
            xstar_crit = ones(n) .- δp .* y
            δx = norm(real.(xstar_crit) .- ones(n))
            flag = "baseline"  # mark analytic baseline
            println(io, "$(seed_i),$(pert_i),0.0,$(δp),$(δx),$(flag)")
          end
        else
          # --- continuation run ---
          syst = U.evaluate_pars(syst_generic, α, [1-αval, αval])
          for pert_i in 1:size(pert_dirs,1)
            println("Seed: $seed_i, Perturbation: $pert_i, α: $αval")
            v           = pert_dirs[pert_i, :]
            init_pars   = r
            target_pars = r .+ pert_size*v
            pars_crit, xstar_crit, flag =
              U.findparscrit(syst, init_sol, init_pars, target_pars)

            # Compute metrics for the main results file
            δp     = norm(pars_crit .- init_pars)
            δx     = norm(real.(xstar_crit) .- ones(n))
            δx_rel = δx / max(δp, eps())  # guard against division by zero

            # Append to the main results file
            println(io, "$(seed_i),$(pert_i),$(αval),$(δp),$(δx_rel),$(flag)")

            # Build the fully numeric system at the critical parameters, solve, and INCLUDE complex solutions
            syst_at_crit = U.plug_parameters_into_system(syst, pars_crit)
            all_sols     = solutions(solve(syst_at_crit))  # 1) include complex + real

            # 2–3) Distances use only real parts; only distances to xstar_crit are computed
            pairs = U.pairwise_distances(xstar_crit, all_sols; atol=1e-6)

            # 4) Append rows and flag types (i = real, j = real/complex)
            if isempty(pairs)
                U.append_no_pairs!(pairwise_path, seed_i, pert_i, αval; nvars = n, nsols = length(all_sols))
            else
                U.append_pairwise_distances!(pairwise_path, seed_i, pert_i, αval, pairs;
                                            nvars = n, nsols = length(all_sols))
            end
          end
        end
      end
    end
  end
end

println("✅ Test complete – see boundary_distances_n_$(n)_test.csv")
