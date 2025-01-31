include("functions.jl")

#Generates data to plot average feasibility critical radius 
#for a given n (number of species) and alpha (HOIs strength).

constrain_type = 2 #tractability and equilibrium constrains are on.
#number of species
nmax = 5
#order of interactions
#to be coded...
#values of alpha
alphavec = 0.01:0.1:0.99
#number of points to check is fix for all n
nperturbations = 1000
#number of simulations for each parameter set
nsim = 100
tol = 1e-9
npointsbeyond = 10
rhobeyond = 1.5
#random seed
rng = MersenneTwister(1)
mode = "follow" #either "follow" or "all"

@time begin
    for sim in 1:nsim #perform many replicates
        #load the vector of perturbations
        for n in 2:nmax #loop through system sizes
            #declare symboilic vector of species abundances
            @var x[1:n]
            #sample parameters
            r0, A, B = sampleparameters(n, rng, constrain_type)
            perturbations = points_hypersphere(n, 1.0, nperturbations)
            for alpha in alphavec #loop through HOI strength
                println("Simulation ", sim, " for n ", n, " alpha ", alpha)
                #form list of parameters for current simulation
                pars = (alpha, r0, A, B)
                target = get_first_target(pars, 1e-4, perturbations, nperturbations, x, n, tol)
                #find critical radius
                println("Find rho critical (mode is ", mode, ")")
                crit_radius = findmaxperturbation(tol, 10.0, perturbations, 
                                                    nperturbations, 
                                                    pars, n, x,
                                                    target,
                                                    mode, tol)
                #now extend from critical radius to 1.5*critical radius to see how fast does the 
                #feasibility breaking happens.
                rho_vec = range(crit_radius, rhobeyond*crit_radius, length = npointsbeyond)
                for rho in rho_vec
                    rpert = pars[2] .+ rho
                    parsnew = (pars[1], rpert, pars[3], pars[4])
                    #get proportion of feasible states
                    equilibria = perturbondisc(perturbations, rhob, parsnew, n, x, false, false, mode)
                    if rho == crit_radius
                        if mode == "follow"
                            all_equilibria = perturbondisc(perturbations, rhob, parsnew, n, x, false, false, "all")
                            #compute and store average number of positive equilibria at critical_radius
                            n_eq_av = average_number_positive_rows(all_equilibria)
                            tosave_n_eq = [sim n alpha crit_radius n_eq_av mode]
                            open("../data/results_critical_radius_n_eq_"*string(mode)*".csv", :"a") do io
                                writedlm(io, tosave_n_eq, ' ')
                            end
                        else
                            #compute and store average number of positive equilibria at critical_radius
                            n_eq_av = average_number_positive_rows(equilibria)
                            tosave_n_eq = [sim n alpha crit_radius n_eq_av mode]
                            open("../data/results_critical_radius_n_eq_"*string(mode)*".csv", :"a") do io
                                writedlm(io, tosave_n_eq, ' ')
                            end
                        end
                    end
                    positive_equilibria = select_feasible_equilibria(equilibria, target)
                    prop_feas = length(getfeasrowinds(positive_equilibria))/nperturbations
                    tosave_prop_feas_states = [sim n alpha rho prop_feas mode]
                    #store
                    open("../data/results_prop_feasible_states_"*string(mode)*".csv", "a") do io
                        writedlm(io, tosave_prop_feas_states, ' ')
                    end
                end
            end
        end
    end
end