BEGIN {
  print "";
  #angstrom units (radius)
  Orad=1.52;
  Hrad=1.2;
  Crad=1.7;
  Prad=1.8;
  Nrad=1.55;
  msms=1;
  comx= 13.933;
  comy= 16.4632;
  comz= 67.8486;
}
{
 format="%7.5G    %7.5G    %7.5G  %G";
 got=0;

 rx = $7;
 ry = $8;
 rz = $9;
 #rx = rx - comx;
 #ry = ry - comy;
 #rz = rz - comz; 
 if ($3~"^H.*") {printf(format, rx,ry,rz,Hrad); got=1};
 if ($3~"^O.*") {printf(format, rx,ry,rz,Orad)  ; got=1};
 if ($3~"^P.*") {printf(format, rx,ry,rz,Prad)  ; got=1};
 if ($3~"^C.*") {printf(format, rx,ry,rz,Crad)  ; got=1};
 if ($3~"^N.*") {printf(format, rx,ry,rz,Nrad)  ; got=1};
 if ($3~"^1H.*"){printf(format, rx,ry,rz,Hrad)  ; got=1};
 if ($3~"^2H.*"){printf(format, rx,ry,rz,Hrad)  ; got=1};
 if (got==1)
  {
    if (msms==0)
      print $3;
    else 
      print "";
  }
}
