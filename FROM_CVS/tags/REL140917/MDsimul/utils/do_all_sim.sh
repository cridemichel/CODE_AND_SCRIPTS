i=0
while [ $i -lt 16 ]
do
cat $HOME/simdat/pars/BMsoft3.templ | awk -v I=$i -F : '{if ($1=="inifile") printf ("inifile:/infm/cmpnapa3/simdat4/mdtmp/bimixT0.78_PP12PR.cor_R%d\n",I); else print $0}' > $HOME/simdat/pars/BMsoft3.par
$HOME/MDsimul/bin/bimix3 -f $HOME/simdat/pars/BMsoft3.par
mv $HOME/simdat3/mdtmp/CorT0.78-beoPR.1 $HOME/simdat3/mdtmp/CorT0.78PP12_R$i
i=$[$i+1]
done
