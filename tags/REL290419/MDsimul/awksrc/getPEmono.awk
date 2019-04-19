BEGIN { split(fn, parti, "_R");
 if (parti[2] != "")
   rango = parti[2];
 #print ("rango: " rango)
}

$1 == "ENmin:" {Emin = $2}
$1 == "ENmax:" {Emax = $2}
$1 == "T:" { parT = $2}

$1 == "PE:" {
  PE_PTS = NF - 1;
  for (i = 2 ; i <= NF; i++ ) {
    E = Emin + i*(Emax-Emin)/PE_PTS; 
    PE[i] = $i;
    EV[i] = E;
  }
}

END {
  temp = parT;
  for (i = 2; i <= (PE_PTS+1); i++) {
    Esc = EV[i];# - 1.5 * temp;
    print (Esc " " PE[i]);
  }
}
  

