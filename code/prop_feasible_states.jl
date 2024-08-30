include("functions.jl")

#set baseline parameters
npertbase = 10
maxn = 4
maxsim = 100
maxrho = 2
tol = 1e-9
rng = MersenneTwister(2)
constrain_type = 1

for sim in 1:maxsim
    for n in 2:maxn
        #declare symboilic vector of species abundances
        @var x[1:n]
        #sample parameters
        r0, A, B = sampleparameters(n, rng, constrain_type) #CAN I SAMPLE DIFFERENT RS PER SPP?
        #set number of perturbations based on dimension n
        nperts = npertbase^n
        for alpha in 0.01:0.1:0.99
            #update parameters
            pars = (alpha, r0, A, B)
            for rho in tol:0.02:0.5
                println("Proportion of feasible states for sim = ", sim, 
                " n =  ", n, " alpha = ", alpha, " rho = ", rho)
                xperts = perturbondisc(rho, pars, n, nperts, x)
                prop_feas = proportion_of_positive_rows(xperts)
                tosave = [prop_feas rho alpha n sim]
                open("../data/prop_feasible_states_2.csv", "a") do io
                    writedlm(io, tosave, ' ')
                end
            end 
        end
    end
end

