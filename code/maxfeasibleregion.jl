#Code to see whats the maximum fully feasible region

include("functions.jl")

#baseline number of perturbations
npertbase = 10
nmax = 2
#set number of simulations
nsim = 1
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
        #sample parameters with appropriate constraints
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        #loop through all alphas for simulation i
        for alpha in 0:0.1:0.99
            println("Searching for sim ", sim, " alpha: ", alpha)
            pars = (alpha, r0, A, B)
            rmax = findmaxperturbation(0, 10, pars, n, nperts, x, 1e-9)
            tosave = [sim n alpha rmax]
            open("../data/feasibility_boundary.csv", "a") do io
                writedlm(io, tosave, ' ')
            end
        end
    end
end