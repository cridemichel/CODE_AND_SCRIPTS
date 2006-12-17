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
TGraph* readascii(char* name, Float_t *minX, Float_t *maxX, Float_t *minY, Float_t *maxY)
{
  //   example of macro to read data from an ascii file and
  //   create a root file with an histogram and an ntuple.
  ifstream in;
  // we assume a file basic.dat in the current directory
  // this file has 3 columns of float data
  in.open(name);
  printf("reading: %s\n", name);
  Float_t *xarr, *yarr, erry, avg;
  Float_t x, y, foo;
  Int_t nlines = 0, i;

  while (1) {
    in >> x >> y >> foo;
    if (!in.good()) break;
    if (nlines < 5) printf("x=%8f, y=%8f \n",x,y);
    ntuple->Fill(x,y);
    nlines++;
  }
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
  Double_t dx = 5*(*maxX-*minX)/100;
  Double_t dy = 5*(*maxY-*minY)/100;
  h2 = new TH2F("h2","fit demo ",1000,(*minX-dx),(*maxX+dx), 1000, (*minY-dy), (*maxY+dy));
  //for (i=0; i < nlines; i++)
  //  {
  //    h2->Fill(xarr[i],yarr[i]);
  //  }
  grafico = new TGraph(nlines, xarr, yarr);
  printf(" found %d pointsn (min=%.15G max=%.15G miny=%.15G maxy=%.15G dx=%f dy=%f)\n",nlines, 
	 *minX, *maxX, *minY, *maxY, dx, dy);

  in.close();
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
void fitFqt(char *fileName=NULL, Int_t type=0, Float_t beg=0.0, Float_t end=0.0)
{
  TFile *f = new TFile("basic.root","RECREATE");
  c1 = new TCanvas("c1","fit with stretched exponential",10,10,700,500);
  TAxis *xax, *yax;
  Int_t npar;
  Double_t par[10], dx, dy;
  Float_t minxXF, maxXF;
  Float_t minX, maxX, minY, maxY;
  c1->SetFillColor(33);
  c1->SetFrameFillColor(41);
  c1->SetGrid();
  if (fileName==NULL)
    {
      printf("You have to supply the filename!\n");
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
  else
    {
      minXF = minX;
      maxXF = maxX;
    }
  if (type==0)
    {
      npar = 2;
    }  
  else
    {
      npar = 3;
    }
  if (type==0)
    fitFcn = new TF1("fitFcn",fitFunctionExp, minXF, maxXF, npar);
  else
    fitFcn = new TF1("fitFcn",fitFunctionStreched, minXF, maxXF, npar);
  // create a TF1 with the range from 0 to 3 and 6 parameters
  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(2);
  fitFcn->SetLineColor(kBlue);
  /* q=800 guess */
  if (npar == 3)
    fitFcn->SetParameters(1.0,1.0,1.0);
  else
    fitFcn->SetParameters(1.0,1.0);
#endif
  //fitFcn->Draw();
  //h2->Fit("fitFcn");
  //h2->SetMarkerStyle(21);
  //h2->SetMarkerSize(0.8);
  // writes the fit results into the par array
  //fitFcn->GetParameters(par);
  grafico->SetLineWidth(3); 
  grafico->SetMarkerStyle(21); 
  grafico->SetLineColor(2); 
  //use a TProfile to convert the 2-d to 1-d problem
  //TProfile *prof = h2->ProfileX();
  //prof->Fit("fitFcn");
  Double_t dx = 5*(maxX-minX)/100;
  Double_t dy = 5*(maxY-minY)/100;
#if 1
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
  grafico->Draw("PL");
 //prof->Draw("same");
  //prof->SetStats(kFALSE);
  //fitFcn->Draw("same");
  // draw the legend
#if 1
  TLegend *legend = new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(grafico,"Data","lp");
  legend->AddEntry(fitFcn,"Global Fit","l");
  legend->Draw();
#endif
  f->Write();

}
