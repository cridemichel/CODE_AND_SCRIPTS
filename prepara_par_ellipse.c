#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>



int main(int argc, char **argv){
  
  int parnum=256;
  int parnumA=256;
  int stepnum = 10000;
  char endfile[1000];
  char inifile[1000];
  char tmpPath[1000];
  char misPath[1000];
  char localPath[1000];
  char *Phi;
  char *El; 
  double epsd = 1.0e-5; 
  double zbrentTol = 1.0e-7;
  double A0, B0=1.0, C0=1.0, rcut;
  
  double Dt=0.01;
  int    Nsave=1000; /* number of configurations to save */
  double storerate;

  double rNebrShell = 0.12;
  double diagNN;

#ifdef _NSAVE
  Nsave = _NSAVE;
#endif
#ifdef _PARNUM
  parnum=_PARNUM;
  parnumA=_PARNUM;
#endif
  
  if(argc<3){
    fprintf(stderr,"\nTo use the program, I need Phi and El(ongation)\n\n");
    fprintf(stderr,"I can use also stepnum\n\n");
    fprintf(stderr,
	    "Compile with -D_GROWTH=endphi to grow a configuration\n");
    fprintf(stderr,
	    "Compile with  -D_NEWCONF -D_GROWTH=endphi to grow a new conf\n");
    fprintf(stderr,"Compile with -D_STEEPDESC to use steepest descent\n");
    exit(0);
  }


  Phi=argv[1];
  El=argv[2];
  if(argc>3) { stepnum=atoi(argv[3]); if(stepnum<Nsave) Nsave=1+stepnum/10;}

  sprintf(localPath,"CNF");

#ifdef _NEWCONF
//  sprintf(inifile,"*\nL: %g",60.2*(double)parnum/(double)256);
  sprintf(inifile,"*\nL: %g",60.2*pow((double)parnum/(double)256,1./3.));
//  sprintf(inifile,"*\nL: %g",60.0*pow((double)parnum/(double)256,1./3.));
#else
  sprintf(inifile,"%s/ellipsoidPhi%s_%sPR.cor",localPath,Phi,El);
#endif

  sprintf(endfile,"%s/ellipsoidPhi%s_%sPR.cor",localPath,Phi,El);

  sprintf(tmpPath,"/home/scala/simdat/mdtmp/ellipsoid/Phi%s_%s/",Phi,El);
  sprintf(misPath,tmpPath);

  A0=atof(El); while(A0<1.0) {A0*=2.0; B0*=2.0; C0*=2.0;}

  {
    double a0=A0+rNebrShell,b0=B0+rNebrShell,c0=C0+rNebrShell;
    diagNN = 2.0*sqrt(a0*a0+b0*b0+c0*c0);
  }

  /*  if(A0>B0) rcut=2.01*A0; else rcut=2.01*B0; */
  rcut = 1.01*diagNN; // If I don't use neighbour lists, change this!!
  
  storerate=(stepnum*Dt)/(double)Nsave;

  printf("parnum: %d\n",parnum);
  printf("parnumA: %d\n",parnumA);
  printf("stepnum: %d\n",stepnum);

  printf("inifile: %s\n",inifile);
  
  printf("endfile: %s\n",endfile);
  printf("endFormat: 2\n");
  
  printf("guessDistOpt: 1\n");
  printf("forceguess: 0\n");

  printf("inistep: 1\n");
  printf("bakSteps: 0\n");
  printf("bakStepsAscii: 0\n");

  printf("tmpPath: %s\n",tmpPath);
  printf("misPath: %s\n",misPath);
  
  printf("temperat: 1.0\n");
  printf("Dt: %g\n",Dt);
  printf("nRun: ell%s\n",El);
  //printf("seed: %d\n",0); /* 0 fixed, -1 random */
  printf("seed: %d\n",-1); /* 0 fixed,-1 random */
  printf("rescaleTime: 0.1\n");
  printf("scalevel: 1\n");
  printf("CMreset: 0\n");
  printf("tapeTimes: 0\n");
  if((stepnum/Nsave)>100) printf("energyCalc: %d\n",stepnum/Nsave);
  else printf("energyCalc: 100\n");
  printf("energyName: energy-\n");
  printf("tempSteps: 10000\n");
  printf("tempName: temp-\n");

  printf("DtrSteps: 0\n");
  printf("DtrCalc: 1\n");
  printf("rotMSDCalc: 1\n");
  printf("rotMSDName: rotMSD-\n");
  printf("DtrName: D-\n");

  //neighbour lists
  printf("rNebrShell: %lf\n",rNebrShell);
  printf("useNNL: %d \n",1);
  printf("epsdNL: %lf\n",2.0e3*epsd);
  printf("epsdFastNL: %lf\n",2.1e3*epsd);
  printf("epsdFastRNL: %lf\n",2.1e3*epsd);
  printf("epsdMaxNL: %lf\n",2.0e3*epsd);
  printf("nebrTabFac: %d\n",500); //max number of n.n. in a n.n.l.

#ifdef _STEEPDESC
  printf("springkSD: 1\n");
  if(A0>B0) printf("stepSDA: %g\n",B0);
  else printf("stepSDA: %g\n",A0);
  printf("tolSD: 0.20\n");
  printf("tolSDlong: 0.05\n");
  printf("tolSDconstr: 0.05\n");
  printf("maxitsSD:  500\n");
#endif

  printf("h: %g\n",1.0e-7); // minimum step
  printf("epsd: %g\n",epsd);
  printf("epsdFast: %g\n",1.2*epsd);
  printf("epsdFastR: %g\n",1.2*epsd);
  printf("epsdMax: %g\n",1.1*epsd);
  printf("zbrakn: 100\n");
  if(epsd<100.0*zbrentTol) zbrentTol=1.0e-2*epsd;
  printf("zbrentTol: %g\n",zbrentTol);


  // sistema ridotto:
  printf("dist5NL: %d\n",1);
  { 
    double elong = atof(El);
    if((elong<=4.0)&&(elong>=0.25)) printf("dist5: %d\n",1);
  }
      
#ifdef _GROWTH
  printf("targetPhi: %g\n",_GROWTH);
#else
  printf("targetPhi: 0.0\n");
#endif
  printf("phitol: 1E-10\n");
  printf("axestol: 1E-8\n");
  printf("scalfact: %g\n",0.2);
  printf("reducefact: %g\n",0.4);
  printf("A0: %lf\n",A0);
  printf("B0: %lf\n",B0);
  printf("C0: %lf\n",C0);
  printf("Ia: 0.4\n");
  printf("Ib: 0.4\n");
  printf("mass0: 1.0\n");
  printf("rcut: %lf\n",rcut);
  printf("overlaptol: 0.01\n");
  printf("intervalSum: 0.2\n");
  printf("eventMult: 200\n");
  printf("bakSaveMode: 0\n");
  printf("storerate: %lf\n",storerate);
  printf("NN: 1\n");
  printf("base: 1\n");
  
 return 0;
}
