#Code to see whats the maximum fully feasible region

include("functions.jl")

#baseline number of perturbations
npertbase = 10
nmax = 3
#set number of simulations
nsim = 2
#set seed for reproducibility
rng = MersenneTwister(2)
rhovec = 10 ./ collect(1:1:1000)
for n in 2:nmax
    #specific number of perturbations given system dimension need
    nperts = 10^n
    #define variables for polynomial construction
    @var x[1:n]
    #set type of constraints for simulations
    constrain_type = 1
    #loop over simulations
    for sim in 1:nsim
        #sample parameters with appropriate constraints
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        #loop through all perturbation magnitudes
        for alpha in 0:0.1:0.99
            #loop through all perturbations on this shell
            full_feasible = false
            for rho in rhovec
                #generate all perturbations on surface of hypershpere of radius rho
                perts_rho = points_hypersphere(n, rho, nperts) #ADD OPTION OF N=2, ORDERED PERTURBATIONS
                #for this shell of perturbations, loop through all values of alphas
                for pert in 1:nperts 
                    #get specific perturbation and build new vector
                    pert_rho_i = perts_rho[pert,:]
                    rpert = r0 + pert_rho_i
                    println("Perturbing for system ", sim, ", n  = ", n, ", radius = ", rho, 
                            ",  alpha = ", alpha, " pert: ", pert, " of ", nperts)
                    #form parameters of perturbed system with alpha value
                    pars = (alpha, rpert, A, B)
                    #solve system (get real solutions)
                    solmat = makeandsolve(x, pars)
                    #count number of solutions
                    nsols = size(solmat)[1]
                    if nsols == 0 #no real solutions
                        println("not real solutions, breaking")
                        break
                    elseif !anyrowfeasible(solmat) #no feasible solutions
                        println("not feasible solutions, breaking")
                        break
                    else
                        continue
                    end
                    if pert == nperts
                        println("found fully feasible, breaking")
                        full_feasible = true
                    end
                end
                if full_feasible == true
                    #save data for the fully feasible orbit
                    tosave = [sim, n, alpha, rho]
                    open("../data/feasibility_boundary.csv", "a") do io
                        writedlm(io, tosave, ' ')
                    end
                    break
                end
            end
        end
    end
end