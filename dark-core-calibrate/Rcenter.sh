cd grid
for LOGS in $(ls); do 
    mkdir $LOGS/solar-freqs
    cd $LOGS/solar-freqs 
    cp ../../../gyre.in .
    $GYRE_DIR/bin/gyre gyre.in
    cp solar-freqs.dat ..
    cd -
done 

for LOGS in $(ls); do 
    mkdir $LOGS/solar-freqs
    cd $LOGS/solar-freqs 
    Rcenter=$(python -c "import pandas as pd; data = pd.read_table('../solar.data', nrows=1, skiprows=1, sep='\s+'); print((data['R_center']/data['rsun']).to_string(index=False))")
    cp ../../../gyre.in .
    shmesa change gyre.in inner_bound "'ZERO_R'"
    shmesa change gyre.in x_i $Rcenter
    $GYRE_DIR/bin/gyre gyre.in
    cp solar-freqs.dat ..
    cd -
done 

for LOGS in $(ls); do 
    mkdir $LOGS/solar-freqs
    cd $LOGS/solar-freqs 
    Rcenter=$(python -c "import pandas as pd; data = pd.read_table('../solar.data', nrows=1, skiprows=1, sep='\s+'); print((data['R_center']/data['rsun']).to_string(index=False))")
    cp ../../../gmode.in .
    shmesa change gmode.in inner_bound "'ZERO_R'"
    shmesa change gmode.in x_i $Rcenter
    $GYRE_DIR/bin/gyre gmode.in
    cp gmodes.dat ..
    cd -
done 

for LOGS in $(ls); do 
    mkdir $LOGS/solar-freqs
    cd $LOGS/solar-freqs 
    mv ../gmodes_xih.dat ../gmodes_xih2.dat
    #Rcenter=$(python -c "import pandas as pd; data = pd.read_table('../solar.data', nrows=1, skiprows=1, sep='\s+'); print((data['R_center']/data['rsun']).to_string(index=False))")
    Rcenter=$(python -c "import pandas as pd; data = pd.read_table('../history.data', skiprows=5, sep='\s+'); print((data['R_B'].values[-1]))")
    cp ../../../gmode_xih.in gmode_xih.in
    shmesa change gmode_xih.in inner_bound "'ZERO_H'"
    shmesa change gmode_xih.in x_i $Rcenter
    $GYRE_DIR/bin/gyre gmode_xih.in
    cp gmodes_xih.dat ..
    cd -
done 

mkdir LOGS/solar-freqs
cd LOGS/solar-freqs
cp ../../gmode.in .
$GYRE_DIR/bin/gyre gmode.in
cp gmodes.dat ..
cp ../../gyre.in .
$GYRE_DIR/bin/gyre gyre.in
cp solar-freqs.dat ..
cd - 
