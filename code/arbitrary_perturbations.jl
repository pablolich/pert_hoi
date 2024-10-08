#Code to perform all simulations needed for the paper

#load needed packages and functions into script
include("functions.jl")

#baseline number of perturbations
npertbase = 10
nmax = 2
#set seed for reproducibility
rng = MersenneTwister(2)
#set number of simulations
nsim = 43
#set type of constraints for simulations
constrain_type = 1
rhovec = collect(0.11:0.2:1.5)
for n in 2:nmax
    #specific number of perturbations given system dimension need
    nperts = npertbase^n
    #define variables for polynomial construction
    @var x[1:n]
    #set solution without perturbations
    global xstar = repeat([1], n)
    #loop over simulations
    for sim in 1:nsim
        #sample parameters with appropriate constraints
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        if sim != 43
            continue
        end
        #loop through all perturbation magnitudes
        for alpha in 0.01:0.1:0.99 #WHEN ALPHA IS 0, THE SYSTEM CAN BE SOLVED VERY FAST WITH MATRIX INVERSION
            #get maximum radius for feasibility
            println("Searching maximum radius for sim ", sim, " alpha: ", alpha)
            pars0 = (alpha, r0, A, B)
            #loop through all perturbations on this shell
            for rho in rhovec
                #generate all perturbations on surface of hypershpere of radius rho
                perts_rho = points_hypersphere(n, rho, nperts)
                #for this shell of perturbations, loop through all values of alphas
                for pert in 1:nperts
                    #get specific perturbation
                    pert_rho_i = perts_rho[pert,:]
                    rpert = r0 + pert_rho_i
                    dtheta = angle(r0, rpert)
                    println("Perturbing for system ", sim, ", n  = ", n, ", radius = ", rho, 
                            ",  alpha = ", alpha, " pert: ", pert, " of ", nperts)
                    #form parameters of perturbed system with alpha value
                    pars = (alpha, rpert, A, B)
                    #get linear approximation of response to perturbation
                    xstarlinear = linearresponse(x, pars, ones(n), r0)
                    #solve system (get real solutions)
                    solmat = makeandsolve(x, pars)
                    #count number of solutions
                    nsols = size(solmat)[1]
                    feasible = anyrowfeasible(solmat)
                    if nsols == 0
                        println("No real solutions found")
                        continue
                    #when we hit a boundary, skip iteration
                    elseif feasible == false
                        println("No feasible solutions found")
                        continue
                    else
                        #pick the real solution closest to the one for previous alpha
                        ind_eq_new = identifyequilibrium(xstar, solmat)
                        #update xstar with new closest equilibrium
                        global xstar = solmat[ind_eq_new,:]
                        #prepare to save data in long format
                        simvec = repeat([sim], n)
                        nvec = repeat([n], n)
                        alphavec = repeat([alpha], n)
                        sppid = collect(1:n)
                        rhovec = repeat([rho], n)
                        pertvec = repeat([pert], n)
                        dthetavec = repeat([dtheta], n)
                        tosave = hcat(nvec, simvec, rhovec,
                                      alphavec, pertvec, sppid, r0, rpert,
                                      dthetavec, xstar, xstarlinear)
                        #save data
                        open("../data/simulations_arbitrary_perturbationssim43.csv", "a") do io
                            writedlm(io, tosave, ' ')
                        end
                    end
                end
            end
        end
    end
end
