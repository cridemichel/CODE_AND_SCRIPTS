#start_run.sh 0.20 500 500
cp Phi0.20/ellipsoid.cor Phi0.40/
start_run.sh 0.40 500 500 
cp Phi0.40/ellipsoid.cor Phi0.50/
start_run.sh 0.50 500 1000 
cp Phi0.50/ellipsoid.cor Phi0.53/
start_run.sh 0.53 2000 2000 
cp Phi0.53/ellipsoid.cor Phi0.56/
start_run.sh 0.56 10000 10000 
cp Phi0.56/ellipsoid.cor Phi0.58/
start_run.sh 0.58 20000 20000 
