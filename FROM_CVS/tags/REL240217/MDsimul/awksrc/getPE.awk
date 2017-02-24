BEGIN { split(fn, parti, "_R");
 if (parti[2] != "")
   rango = parti[2];
 #print ("rango: " rango)
}

$1 == "ENmin:" {Emin = $2}
$1 == "ENmax:" {Emax = $2}
$1 == "T:" { parT = $2}

$1 == "lambda0:" { 
  for (i = 2 ; i <= NF; i++ ) {
    lambda0[i-2] = $i;
  }
}

$1 == "lambdat:" {
  for (i = 2 ; i <= NF; i++ ) {
    lambdat[i-2] = $i;
#print ($i)
  }
}
$1 == "PE:" {
  PE_PTS = NF - 1;
  for (i = 2 ; i <= NF; i++ ) {
    E = Emin + i*(Emax-Emin)/PE_PTS; 
    PE[i] = $i;
    EV[i] = E;
  }
}

END {
  #print ("rango: " rango " lt: " lambdat[rango] );
  lt = lambdat[rango];
  temp = parT / lambda0[lt];
  #temp = 0.485;
  #print ("lambda0[lt]: " lambda0[lt] "temp:" temp);;
  for (i = 2; i <= (PE_PTS+1); i++) {
    Esc = EV[i] - 1.5 * temp;
    print (Esc " " PE[i]);
  }
}
  

