#Generates data to plot average feasibility critical radius 
#for a given n (number of species) and alpha (HOIs strength).

constrain_type = 1 #tractability and equilibrium constrains are on.
#number of species
nmax = 2#4
#order of interactions
#to be coded...
#values of alpha
alphavec = [0.2, 0.4]#0.01:0.1:0.99
#number of points to check is fix for all n
nperts = 2#1000
#number of simulations for each parameter set
nsim = 2#1000
tol = 1e-9

for sim in 1:nsim #perform many replicates
    for n in 2:nmax #loop through system sizes
        #declare symboilic vector of species abundances
        @var x[1:n]
        #sample parameters
        r0, A, B = sampleparameters(n, rng, constrain_type)
        #load the vector of perturbations
        ################################
        #TODO
        ################################
        for alpha in alphavec #loop through HOI strength
            #form list of parameters for current simulation
            pars = (alpha, r0, A, B)
            #find first equilibria around the (1, 1, 1, ... , 1) point
            target_equilibira = get_first_target(pars, tol, perturbations, nperturbations, x, n)
            #find critical radius
            critical_radius = findmaxperturbation(rho1, rho2, perturbations, nperturbations, 
                                                  pars, x, target_equilibria, mode, tol)
            #save results
            tosave = [sim n alpha critical_radius]
            open("../data/results_critical_radius.csv", "a") do io
                writedlm(io, tosave, ' ')
            #now extend from critical radius to 1.5*critical radius to see how fast does the 
            #feasibility breaking happens.
            rho_vec = range(critical_radius, 1.5*critical_radius, length = 10)
            for rho in rho_vec
                #get proportion of feasible states
                all_equlibria = perturbondisc(perturbations, rhob, parameters, x, false)
                prop_feas = proportion_of_positive_rows(xperts)
            end
        end
    end
end