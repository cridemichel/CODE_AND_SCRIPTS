/** do_infile.cpp 
 * 
 * creates the input file for the ED-MD code. it implements a primitive 
 * parser that uses argv for changing the default values of the variables.
 * for example,
 * 
 * 		./do_infile.x --N 100  --maxpf 0.49 --growthrate 0.01 \ 
 *          --readfile new --writefile first.cnf > in.txt
 * 
 * creates the input file to generate the new configuration first.cnf of 100 
 * particles at phi=0.49, while
 * 
 * 		./do_infile.x --N 100 --readfile first.cnf --writefile \ 
 *          scnd.cnf > in.txt
 * 
 * creates the input file to evolve first.cnf for a time = 10.7 and saves the
 * configuration in scnd.cnf */

#include <iostream>
#include <string>    
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <typeinfo>

#define INF    100000000
#define dblINF 100000000.

using namespace std;

string line; ///< here goes the flattening of **argv
string::size_type idx; ///< used to indicate a position in a string
istringstream istr; ///< stream for the string/something converstion 

/// (overloaded): searchs the token search in the global string line; 
/// returns the first value after the token in the variable x 
bool search_token(string search,double &x){
  idx = line.find(search); if (idx!=string::npos) {
    istr.str( line.substr(idx+search.size()) );
    istr >> x; return true; 
  }
  return false;
}
bool search_token(string search,int &x){
  idx = line.find(search); if (idx!=string::npos) {
    istr.str( line.substr(idx+search.size()) );
    istr >> x; return true; 
  }
  return false;
}
bool search_token(string search,string &x){
  idx = line.find(search); if (idx!=string::npos) {
    istr.str( line.substr(idx+search.size()) );
    istr >> x; return true; 
  }
  return false;
}


/***********************************************************************/


int main(int argc, char **argv)
{
  // flattening of argv into line 
  line=""; 
  for(int i = 1; i<argc; i++) line += string(argv[i])+" ";



  //////////////////////////////////////////////////////////////////
  // parse the parameters
  //////////////////////////////////////////////////////////////////

  int parnum(128); search_token("--parnum",parnum);
  int parnumA(128); search_token("--parnumA",parnumA);
  int stepnum(10000000); search_token("--stepnum",stepnum);
  string endfile("ellipsoid-3.cor"); search_token("--endfile",endfile);
  int endFormat(2); search_token("--endFormat",endFormat);
  string inifile("ellipsoid-2.cor"); search_token("--inifile",inifile);
  int inistep(1); search_token("--inistep",inistep);
  int bakSteps(0); search_token("--bakSteps",bakSteps);
  int bakStepsAscii(0); search_token("--bakStepsAscii",bakStepsAscii);
  string tmpPath("./"); search_token("--tmpPath",tmpPath);
  string misPath("./"); search_token("--misPath",misPath);
  double temperat(1.0); search_token("--temperat",temperat);
  double ENmin(-6.5); search_token("--ENmin",ENmin);
  double ENmax(-5.0); search_token("--ENmax",ENmax);
  double Dt(0.05); search_token("--Dt",Dt);
  string nRun("ellips"); search_token("--nRun",nRun);
  int seed(-1); search_token("--seed",seed);
  double rescaleTime(0.001); search_token("--rescaleTime",rescaleTime);
  int scalevel(0); search_token("--scalevel",scalevel);
  int CMreset(0); search_token("--CMreset",CMreset);
  int tapeTimes(0); search_token("--tapeTimes",tapeTimes);
  int energyCalc(50); search_token("--energyCalc",energyCalc);
  string energyName("energy-"); search_token("--energyName",energyName);
  int tempCalc(1); search_token("--tempCalc",tempCalc);
  int tempSteps(1); search_token("--tempSteps",tempSteps);
  string tempName("temp-"); search_token("--tempName",tempName);
  int DtrSteps(10); search_token("--DtrSteps",DtrSteps);
  int DtrCalc(1); search_token("--DtrCalc",DtrCalc);
  string DtrName("D-"); search_token("--DtrName",DtrName);
  int rotMSDCalc(1); search_token("--rotMSDCalc",rotMSDCalc);
  string rotMSDName("rotMSD-"); search_token("--rotMSDName",rotMSDName);
  string VName("V-"); search_token("--VName",VName);
  int VSteps(0); search_token("--VSteps",VSteps);
  double partDiss(0.98); search_token("--partDiss",partDiss);
  double tc(1e-5); search_token("--tc",tc);
  double h(1e-08); search_token("--h",h);
  double epsd(1e-5); search_token("--epsd",epsd);
  double epsdFast(1.1e-5); search_token("--epsdFast",epsdFast);
  double epsdFastR(1.1e-5); search_token("--epsdFastR",epsdFastR);
  double epsdMax(1e-5); search_token("--epsdMax",epsdMax);
  double epsdNL(3.0e-2); search_token("--epsdNL",epsdNL);
  double epsdFastNL(3.1e-2); search_token("--epsdFastNL",epsdFastNL);
  double epsdFastRNL(3.1e-2); search_token("--epsdFastRNL",epsdFastRNL);
  double epsdMaxNL(3.0e-2); search_token("--epsdMaxNL",epsdMaxNL);
  int guessDistOpt(1); search_token("--guessDistOpt",guessDistOpt);
  int forceguess(0); search_token("--forceguess",forceguess);
  double springkSD(1.0); search_token("--springkSD",springkSD);
  double toldxNR(0.35); search_token("--toldxNR",toldxNR);
  int SDmethod(0); search_token("--SDmethod",SDmethod);
  double stepSDA(0.4); search_token("--stepSDA",stepSDA);
  double stepSDB(1.0); search_token("--stepSDB",stepSDB);
  double tolSDgrad(0.4); search_token("--tolSDgrad",tolSDgrad);
  int maxitsSD(1000); search_token("--maxitsSD",maxitsSD);
  int dist8stps(0); search_token("--dist8stps",dist8stps);
  double tolSD(0.01); search_token("--tolSD",tolSD);
  double tolSDlong(0.01); search_token("--tolSDlong",tolSDlong);
  double epsdSD(0.01); search_token("--epsdSD",epsdSD);
  double tolSDconstr(0.002); search_token("--tolSDconstr",tolSDconstr);
  int zbrakn(50); search_token("--zbrakn",zbrakn);
  double zbrentTol(1e-7); search_token("--zbrentTol",zbrentTol);
  double targetPhi(0.0); search_token("--targetPhi",targetPhi);
  double phitol(1e-7); search_token("--phitol",phitol);
  double axestol(1e-6); search_token("--axestol",axestol);
  double scalfact(0.1); search_token("--scalfact",scalfact);
  double reducefact(0.1); search_token("--reducefact",reducefact);
  double A0(3.0); search_token("--A0",A0);
  double B0(1.0); search_token("--B0",B0);
  double C0(1.0); search_token("--C0",C0);
  double Ia(0.4); search_token("--Ia",Ia);
  double mass0(1.0); search_token("--mass0",mass0);
  double rcut(7.1); search_token("--rcut",rcut);
  double overlaptol(0.01); search_token("--overlaptol",overlaptol);
  double intervalSum(0.2); search_token("--intervalSum",intervalSum);
  int eventMult(200); search_token("--eventMult",eventMult);
  int useNNL(1); search_token("--useNNL",useNNL);
  int paralNNL(1); search_token("--paralNNL",paralNNL);
  int dist5(1); search_token("--dist5",dist5);
  int dist5NL(1); search_token("--dist5NL",dist5NL);
  double rNebrShell(0.14); search_token("--rNebrShell",rNebrShell);
  int nebrTabFac(500); search_token("--nebrTabFac",nebrTabFac);
  int bakSaveMode(0); search_token("--bakSaveMode",bakSaveMode);
  double storerate(0.01); search_token("--storerate",storerate);
  int NN(43); search_token("--NN",NN);
  double base(1.3); search_token("--base",base);

  //////////////////////////////////////////////////////////////////
  //adjust parameters
  //////////////////////////////////////////////////////////////////
  string junk_s;

#ifdef _LOG_SAVING
  /* logarithmic saving of configurations */
  double Tau(1); search_token("--Tau",Tau);
  int NTau(10); search_token("--NTau",NTau);
  storerate = Dt;
  stepnum   = NTau*static_cast<int>( ceil(Tau/Dt) );
  NN        = static_cast<int>( ceil( log(Tau/storerate) / log(base) ) );
#else
  /* linear saving of configurations */
  int Nsave(1000); search_token("--Nsave",Nsave);
  if(stepnum<Nsave) Nsave=1+stepnum/10;
  storerate=(stepnum*Dt)/(double)Nsave;
#endif
  
  parnumA=parnum; // only 1 species

  if ( search_token("--elastic",junk_s) ) {partDiss=1.0; tc=0.0;}
  
  // if a newconf, should be very diluted 
  if(search_token("--NewConf",junk_s)){
    inifile="*"; 
    rescaleTime=0.1;
    scalevel=1;
    stepnum=100; storerate=stepnum*Dt/2.0;
    partDiss=1.0; tc=0.0;
    cout <<"L: "<<60.2*pow((double)parnum/(double)256,1./3.)<<endl;
  }  

  /*
  // A0 should always be ~1
  while(A0<1.0) {A0*=2.0; B0*=2.0; C0*=2.0;}
  */

  double Elongation = A0/B0;

  if(Elongation<1) rNebrShell = 0.2; // better for oblate
  double a0=A0+rNebrShell,b0=B0+rNebrShell,c0=C0+rNebrShell;
  double diagNN = 2.0*sqrt(a0*a0+b0*b0+c0*c0);
  rcut = 1.01*diagNN; // If I don't use neighbour lists, change this!!



  //////////////////////////////////////////////////////////////////
  //write the parameters
  //////////////////////////////////////////////////////////////////

  cout<<"parnum: "<<parnum<<endl; 
  cout<<"parnumA: "<<parnumA<<endl; 
  cout<<"endfile: "<<endfile<<endl; 
  cout<<"endFormat: "<<endFormat<<endl; 
  cout<<"inifile: "<<inifile<<endl; 
  cout<<"inistep: "<<inistep<<endl; 
  cout<<"bakSteps: "<<bakSteps<<endl; 
  cout<<"bakStepsAscii: "<<bakStepsAscii<<endl; 
  cout<<"tmpPath: "<<tmpPath<<endl; 
  cout<<"misPath: "<<misPath<<endl; 
  cout<<"temperat: "<<temperat<<endl; 
  //  cout<<"ENmin: "<<ENmin<<endl; 
  //  cout<<"ENmax: "<<ENmax<<endl; 
  cout<<"Dt: "<<Dt<<endl; 
  cout<<"nRun: "<<nRun<<endl; 
  cout<<"seed: "<<seed<<endl; 
  cout<<"rescaleTime: "<<rescaleTime<<endl; 
  cout<<"scalevel: "<<scalevel<<endl; 
  cout<<"CMreset: "<<CMreset<<endl; 
  cout<<"tapeTimes: "<<tapeTimes<<endl; 
  cout<<"energyCalc: "<<energyCalc<<endl; 
  cout<<"energyName: "<<energyName<<endl; 
  cout<<"tempCalc: "<<tempCalc<<endl; 
  cout<<"tempSteps: "<<tempSteps<<endl; 
  cout<<"tempName: "<<tempName<<endl; 
  cout<<"DtrSteps: "<<DtrSteps<<endl; 
  cout<<"DtrCalc: "<<DtrCalc<<endl; 
  cout<<"DtrName: "<<DtrName<<endl; 
  cout<<"rotMSDCalc: "<<rotMSDCalc<<endl; 
  cout<<"rotMSDName: "<<rotMSDName<<endl; 
  //  cout<<"VName: "<<VName<<endl; 
  //  cout<<"VSteps: "<<VSteps<<endl; 
  if(partDiss<1.0&&partDiss>0.0){ 
	cout<<"partDiss: "<<partDiss<<endl; 
  	cout<<"tc: "<<tc<<endl;
  }	
  cout<<"h: "<<h<<endl; 
  cout<<"epsd: "<<epsd<<endl; 
  cout<<"epsdFast: "<<epsdFast<<endl; 
  cout<<"epsdFastR: "<<epsdFastR<<endl; 
  cout<<"epsdMax: "<<epsdMax<<endl; 
  cout<<"epsdNL: "<<epsdNL<<endl; 
  cout<<"epsdFastNL: "<<epsdFastNL<<endl; 
  cout<<"epsdFastRNL: "<<epsdFastRNL<<endl; 
  cout<<"epsdMaxNL: "<<epsdMaxNL<<endl; 
  cout<<"guessDistOpt: "<<guessDistOpt<<endl; 
  cout<<"forceguess: "<<forceguess<<endl; 
  cout<<"springkSD: "<<springkSD<<endl; 
  cout<<"toldxNR: "<<toldxNR<<endl; 
  cout<<"SDmethod: "<<SDmethod<<endl; 
  cout<<"stepSDA: "<<stepSDA<<endl; 
  cout<<"stepSDB: "<<stepSDB<<endl; 
  cout<<"tolSDgrad: "<<tolSDgrad<<endl; 
  cout<<"maxitsSD: "<<maxitsSD<<endl; 
  cout<<"dist8stps: "<<dist8stps<<endl; 
  cout<<"tolSD: "<<tolSD<<endl; 
  cout<<"tolSDlong: "<<tolSDlong<<endl; 
  cout<<"epsdSD: "<<epsdSD<<endl; 
  cout<<"tolSDconstr: "<<tolSDconstr<<endl; 
  cout<<"zbrakn: "<<zbrakn<<endl; 
  cout<<"zbrentTol: "<<zbrentTol<<endl; 
  cout<<"targetPhi: "<<targetPhi<<endl; 
  cout<<"phitol: "<<phitol<<endl; 
  cout<<"axestol: "<<axestol<<endl; 
  cout<<"scalfact: "<<scalfact<<endl; 
  cout<<"reducefact: "<<reducefact<<endl; 
  cout<<"A0: "<<A0<<endl; 
  cout<<"B0: "<<B0<<endl; 
  cout<<"C0: "<<C0<<endl; 
  cout<<"Ia: "<<Ia<<endl; 
  cout<<"mass0: "<<mass0<<endl; 
  cout<<"rcut: "<<rcut<<endl; 
  cout<<"overlaptol: "<<overlaptol<<endl; 
  cout<<"intervalSum: "<<intervalSum<<endl; 
  cout<<"eventMult: "<<eventMult<<endl; 
  cout<<"useNNL: "<<useNNL<<endl; 
  cout<<"paralNNL: "<<paralNNL<<endl; 
  cout<<"dist5: "<<dist5<<endl; 
  cout<<"dist5NL: "<<dist5NL<<endl; 
  cout<<"rNebrShell: "<<rNebrShell<<endl; 
  cout<<"nebrTabFac: "<<nebrTabFac<<endl; 
  cout<<"bakSaveMode: "<<bakSaveMode<<endl;
 
  cout<<"stepnum: "<<stepnum<<endl; 
  cout<<"storerate: "<<storerate<<endl; 
  cout<<"NN: "<<NN<<endl; 
  cout<<"base: "<<base<<endl; 


  return 0;
}

