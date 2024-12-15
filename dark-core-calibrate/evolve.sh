for lgM in -2 -3 -4 -5 -6 -7; do
    shmesa change inlist_calibrate new_core_mass "1d$lgM"
    ./star inlist_evolve 
    mv LOGS grid_evol/LOGS_$lgM
done
