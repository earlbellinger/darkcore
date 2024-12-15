mkdir grid 
shmesa change inlist_calibrate "x_ctrl(1)" 1

mkdir grid
for lgM in -8 -7 -6 -5 -4 -3 -2; do
    shmesa change inlist_calibrate new_core_mass "1d$lgM"
    python3 calibrate_diffusion.py
    cp inlist_calibrate LOGS
    mv LOGS grid/LOGS_$lgM
done

mkdir griddense
for lgM in -3.0 -2.9 -2.8 -2.7 -3.1 -3.2 -2.95; do
    shmesa change inlist_calibrate new_core_mass "$(echo "scale=10; e($lgM*l(10))" | bc -l)"
    python3 calibrate_diffusion.py
    cp inlist_calibrate LOGS
    mv LOGS griddense/LOGS_$lgM
done

exit 

mkdir grid_rb_-1
shmesa change inlist_calibrate "x_ctrl(1)" 0.1

for lgM in -7 -6 -5 -4 -3 -2; do
    shmesa change inlist_calibrate new_core_mass "1d$lgM"
    python3 calibrate_diffusion.py
    cp inlist_calibrate LOGS
    mv LOGS grid_rb_-1/LOGS_$lgM
done

exit 


mkdir grid_rb_-2
shmesa change inlist_calibrate "x_ctrl(1)" 0.01

for lgM in -7 -6 -5 -4; do
    shmesa change inlist_calibrate new_core_mass "1d$lgM"
    python3 calibrate_diffusion.py
    cp inlist_calibrate LOGS
    mv LOGS grid_rb_-2/LOGS_$lgM
done
