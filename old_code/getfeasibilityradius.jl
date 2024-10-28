#Code to see whats the maximum fully feasible region

include("functions.jl")

#baseline number of perturbations
npertbase = 10
nmax = 3
#set number of simulations
nsim = 2
#set seed for reproducibility
rng = MersenneTwister(2)
rhovec = 10 ./ collect(1:1000:100000)
global full_feasible = false
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
        if sim != 43
            continue
        end
        #loop through all alphas for simulation i
        for alpha in 0:0.1:0.99
            rhoind = 1
            while !full_feasible
                #decrease radius until fully feasible
                rho = rhovec[rhoind]
                #generate all perturbations on surface of hypershpere of radius rho
                perts_rho = points_hypersphere(n, rho, nperts) #ADD OPTION OF N=2, ORDERED PERTURBATIONS
                #loop through all perturbations on this shell
                for pert in 1:nperts 
                    #get specific perturbation and build new vector
                    pert_rho_i = perts_rho[pert,:]
                    rpert = r0 + pert_rho_i
                    println("Perturbing for system ", sim, ", n  = ", n, ", rho = ", rho, 
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
                        if pert == nperts
                            println("found fully feasible")
                            #declare feasible if all perturbations traversed yielded feasible
                            global full_feasible = true
                        end
                    end
                end
                #see if all perturbations werer traversed
                if full_feasible == true
                    print("saving data")
                    #save data for the fully feasible orbit
                    tosave = [sim n alpha rho]
                    open("../data/feasibility_boundarysim43.csv", "a") do io
                        writedlm(io, tosave, ' ')
                    end
                    global full_feasible = false
                    rhoind = 1
                    break
                end
                rhoind += 1
            end
        end
    end
end