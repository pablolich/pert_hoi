#Code to see whats the maximum fully feasible region

include("functions.jl")

#baseline number of perturbations
npertbase = 10
nmax = 4
#set number of simulations
nsim = 100
#set seed for reproducibility
rng = MersenneTwister(2)
for n in 2:nmax
    #specific number of perturbations given system dimension need
    nperts = npertbase^n
    #define variables for polynomial construction
    @var x[1:n]
    #set type of constraints for simulations
    constrain_type = 1
    #loop over simulations
    for sim in 1:nsim
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        #loop through all alphas for simulation i
        for alpha in 0.01:0.1:0.99
            println("Searching for n = ", n, " sim = ", sim, " alpha: ", alpha)
            pars = (alpha, r0, A, B)
            rmax = findmaxperturbation(0, 10, pars, n, nperts, x, 1e-9)
            tosave = [sim n alpha rmax]
            open("../data/feasibility_boundary_radius_2.csv", "a") do io
                writedlm(io, tosave, ' ')
            end
            # if n == 2 && sim == 43
            #     if alpha == 0.01 || alpha == 0.61 || alpha == 0.81
            #         eqmat = perturbondisc(rmax, pars, n, nperts, x)
            #         neqs = size(eqmat, 1)
            #         #add a column with the alphas
            #         alphavec = repeat([alpha], neqs)
            #         tosave = hcat(alphavec, eqmat)
            #         #draw the points and save them
            #         open("../data/domainstodraw.csv", "a") do io 
            #             writedlm(io, tosave, ' ')
            #         end
            #     end
            # end
        end
    end
end