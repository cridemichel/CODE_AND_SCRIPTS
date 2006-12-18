#include "Riostream.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
gROOT->Reset();
TH2F *h2;
TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y");
TCanvas *c1;
TGraph *grafico;
Int_t logx=0, logy=0, hXpnts=100, hYpnts=100, autox=100, autoy=100;
TGraph* readascii(char* name, Float_t *minX, Float_t *maxX, Float_t *minY, Float_t *maxY)
{
  //   example of macro to read data from an ascii file and
  //   create a root file with an histogram and an ntuple.
  //ifstream in;
  // we assume a file basic.dat in the current directory
  // this file has 3 columns of float data
  //in.open(name);
  printf("reading: %s\n", name);
  Float_t *xarr, *yarr, erry, avg;
  Float_t x, y, foo;
  Int_t nlines = 0, i;

#if 0
  while (1) {
    in >> x >> y >> foo;
    if (!in.good()) break;
    if (nlines < 5) printf("x=%8f, y=%8f \n",x,y);
    ntuple->Fill(x,y);
    nlines++;
  }
#else
  FILE *f;
  char dummy[4096];
  f=fopen(name,"r");
  while (!feof(f))
    {
      fscanf(f, "%f %f %[^\n]\n", &x, &y, &dummy);
      ntuple->Fill(x,y);
      nlines++;
    }
  fclose(f);
#endif
  xarr = new Float_t[nlines];
  yarr = new Float_t[nlines];
  ntuple->SetBranchAddress("x", &x);
  for (i=0; i<nlines; i++) 
    {
      ntuple->GetEntry(i);
      xarr[i] = x;
      if (i==0 || x < *minX)
	*minX = x;
      if (i==0 || x > *maxX)
	*maxX = x;
    }
  ntuple->SetBranchAddress("y", &y);
  for (i=0; i<nlines; i++) 
    {
      ntuple->GetEntry(i);
      yarr[i] = y;
      if (i==0 || y < *minY)
	*minY = y;
      if (i==0 || y > *maxY)
	*maxY = y;
    }
  Double_t dx, dy;
  if (logx)
    {
      dx = 5.0*(*minX/100.0); 
    }
  else  
    {
      dx = 5.0*(*maxX-*minX)/100.0;
    }
  if (logy)
    {
      dy = 5.0*(*minY/100.0);
    }
  else
    {
      dy = 5.0*(*maxY-*minY)/100.0;
    }

#if 1
  if (autox > 0)
    hXpnts = autox;
  if (autoy > 0)
    hYpnts = autoy;
#endif
  h2 = new TH2F("h2","fit",hXpnts,(*minX-dx),(*maxX+dx), hYpnts, (*minY-dy), (*maxY+dy));
  grafico = new TGraph(nlines, xarr, yarr);
  printf(" found %d pointsn (min=%.15G max=%.15G miny=%.15G maxy=%.15G dx=%f dy=%f)\n",nlines, 
	 *minX, *maxX, *minY, *maxY, dx, dy);

  //in.close();
  return grafico;
}

Double_t fitFunctionStreched(Double_t *x, Double_t* par)
{
  return par[0]*exp(-pow(x[0]/par[1],par[2]));
}
Double_t fitFunctionExp(Double_t *x, Double_t* par)
{
  return par[0]*exp(-x[0]/par[1]);
}

TF1 *fitFcn; 
// type = 0 exponential
// type = 1 stretched
// type = 2 fa prima il fit esponenziale e poi usa i parametri ottenuti 
// per il fit stretched
void fitSE(char *fileName=NULL, Int_t type=0, Float_t beg=0.0, Float_t end=0.0, Int_t llogx=0, Int_t llogy=0,
	   Int_t xpnts=-1, Int_t ypnts=-1)
{
  TFile *f = new TFile("basic.root","RECREATE");
  c1 = new TCanvas("c1","fit with stretched exponential",10,10,700,500);
  TAxis *xax, *yax;
  Int_t npar;
  Double_t par[10], dx, dy;
  Float_t minXF, maxXF;
  Float_t minX, maxX, minY, maxY;
  if (llogx)
    gPad->SetLogx(1);
  if (llogy)
    gPad->SetLogy(1);

  c1->SetFillColor(33);
  c1->SetFrameFillColor(41);
  c1->SetGrid();
  logx = llogx;
  logy = llogy;
  ;autox = xpnts;
  autoy = ypnts;
  if (fileName==NULL)
    {
      printf("You have to supply the filename!\n");
      printf("Usage: root.exe 'fitSE(<filename>,<type>,<beg>,<end>,<logx>,<logy>,<xpnts>,<ypnts>)\n");
      printf("<filename>: file containing data to fit\n");
      printf("<type>: 0 = exponential fit, 1 = stretched exp fit, 2 = fit with exp and use result to fit with stretched exp\n");
      printf("<beg>: start of subrange to fit\n");
      printf("<end>: end of subrange to fit\n");
      printf("<logx>: 1 log scale for x-axis 0 = linear\n");
      printf("<logy>: 1 log scale for y-axis 0 = linear\n");
      printf("<xpnts>: number of points along x-axis for histogram mesh\n");
      printf("<ypnts>: number of points along y-axis for histogram mesh\n");
      printf("<xpnts> and <ypnts> affect the zoom\n");
      printf("if <beg> and <end> are omitted then use all points\n");
      printf("if <beg> = <end> then start fitting from <beg>\n");  
      exit(-1);
    }
  //printf("fileName=%s\n", fileName);
  grafico = readascii(fileName, &minX, &maxX, &minY, &maxY);
#if 1
  /* q=800 guess */

  if (beg==end && beg==0.0)
    {
       minXF = minX;
       maxXF = maxX;   
    }
  else if (beg==end && beg !=0.0)
    {
      minXF = beg;
      maxXF = maxX;
    }
  else if (beg != 0.0 && end != 0.0)
    {
      minXF = beg;
      maxXF = end;
    }
  else
    {
      minXF = minX;
      maxXF = maxX;
    }
  fitFcn = new TF1("fitFcn",fitFunctionStreched, minXF, maxXF, 3);
 // create a TF1 with the range from 0 to 3 and 6 parameters
  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(2);
  fitFcn->SetLineColor(kBlue);
  /* q=800 guess */
  fitFcn->SetParLimits(2,minXF, maxXF);
  fitFcn->SetParameters(1.0,1.0,1.0);
  if (type==0||type==2)
    fitFcn->FixParameter(2,1.0); 
  else 
    fitFcn->SetParLimits(2,0.01,2.2);
	
#endif
  //fitFcn->Draw();
  //h2->Fit("fitFcn");
  //h2->SetMarkerStyle(21);
  //h2->SetMarkerSize(0.8);
  // writes the fit results into the par array
  //fitFcn->GetParameters(par);
  gStyle->SetOptFit(1111);	
  grafico->SetLineWidth(3); 
  grafico->SetMarkerStyle(21); 
  grafico->SetLineColor(2); 
  Double_t dx = 5*(maxX-minX)/100;
  Double_t dy = 5*(maxY-minY)/100;
#if 0
  yax = h2->GetYaxis();
  yax->SetRangeUser(minY-dy,maxY+dy);
  xax = h2->GetXaxis();
  xax->SetRangeUser(minX-dx,maxX+dx);
  // h2 server solo per scegliere la viewport e stampare gli assi
#endif
#if 1
  h2->SetStats(kFALSE);
  h2->Draw();
#endif
  grafico->Fit("fitFcn", "R");

  if (type==2)
    {
      fitFcn->GetParameters(par);
      fitFcn->SetParameters(par[0],par[1],1.0);
      fitFcn->SetParLimits(2,0.01,2.2);
      grafico->Fit("fitFcn", "R");
    }
  grafico->Draw("PL");

  printf("CHISQUARE: %.15G PROB: %.15G\n", fitFcn->GetChisquare(), fitFcn->GetProb());
 //prof->Draw("same");
  //prof->SetStats(kFALSE);
  //fitFcn->Draw("same");
  // draw the legend
#if 1
  TLegend *legend = new TLegend(0.6,0.45,0.88,0.65);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(grafico,"Data","lp");
  legend->AddEntry(fitFcn,"Global Fit","l");
  legend->Draw();
#endif
  f->Write();

}
