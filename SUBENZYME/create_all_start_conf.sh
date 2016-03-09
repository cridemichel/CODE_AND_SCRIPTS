SIGS="1.0"
SIGE="4.5"
SIGC="9.4"
NE="1"
if [ "$1" == "" ]
then
NS="1000"
else
NS="$1"
fi
PHI="0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40"
for f in `echo $PHI`
do
#startCROW_phi0.40.cnf 15 15 15 0.40 4.5 1.0 9.4 4.5 1.0 9.4 1 1000
../create_SUBENZ_conf start_CROW_PHI_${f}_SIGC_${SIGC}.cnf 24 24 24 $f $SIGE $SIGS $SIGC $SIGE $SIGS $SIGC 5 5000
done
