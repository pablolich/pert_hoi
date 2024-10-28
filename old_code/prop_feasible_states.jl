include("functions.jl")

#set baseline parameters
npertbase = 10
maxn = 5
maxsim = 100
maxrho = 2
alpha_vec = [0.3, 0.7]
rho_vec = [0.1]
tol = 1e-9
rng = MersenneTwister(1)
constrain_type = 1

for sim in 1:maxsim
    for n in 2:maxn
        #declare symboilic vector of species abundances
        @var x[1:n]
        #sample parameters
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        #set number of perturbations based on dimension n
        nperts = npertbase^n
        for alpha in alpha_vec#0.01:0.1:0.99
            #update parameters
            pars = (alpha, r0, A, B)
            for rho in rho_vec#tol:0.02:0.5
                println("Proportion of feasible states for sim = ", sim, 
                " n =  ", n, " alpha = ", alpha, " rho = ", rho)
                xperts = perturbondisc(rho, pars, n, nperts, x, false)
                prop_feas = proportion_of_positive_rows(xperts)
                tosave = [prop_feas rho alpha n sim]
                open("../data/prop_feasible_states_large_n.csv", "a") do io
                    writedlm(io, tosave, ' ')
                end
            end 
        end
    end
end

