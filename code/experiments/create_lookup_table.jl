#this script creates a lookup table of what is the 
#minimum number of points needed to get an accurate representation of
#histogram of distances for each value of n, alpha, d=2, d_pert = 0,
#averaged over histograms for many parameter sets.

#for each value of n
    #for each parameter set
        #for each value of n_pert
            #generate directions
            #for each alpha
                #build histogram
            #caluclate the Kullback–Leibler divergence between the current 
            #and previous ones.
            #if the distance is less than a tolerance, flag that value of alpha
            #as finished, and skip in future iterations of n_pert
        #record the n_pert_min for all alphas and current parameter set, and n value 
#return a table with n, seed, alpha, n_pert_min 