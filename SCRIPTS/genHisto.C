#include "Riostream.h"
#include "TH2.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
gROOT->Reset();
TH1F *h1, *h2;
TCanvas *c1;
char dummy[1024];
TNtuple* readascii(char* name, Int_t *nlines, Int_t ncol)
{
  //   example of macro to read data from an ascii file and
  //   create a root file with an histogram and an ntuple.
  //ifstream in;
  // we assume a file basic.dat in the current directory
  // this file has 3 columns of float data
  //in.open(name);
  if (ncol==2)
    TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y");
  else
    TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x");
  printf("reading: %s\n", name);

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
  Float_t x, y;
  char dummy[4096];
  *nlines=0;
  f=fopen(name,"r");
  while (!feof(f))
    {
      if (ncol==2)
	{
  	  fscanf(f, "%f %f %[^\n]\n", &x, &y, dummy);
  	  ntuple->Fill(x,y);
	}
      else 
	{
	  fscanf(f, "%f %[^\n]\n", &x, dummy);
  	  ntuple->Fill(x);
	}
      (*nlines)++;
    }
  fclose(f);
#endif
  printf(" found %d pointsn \n",*nlines);

  //in.close();
  return ntuple;
}

void genHisto(char *fileName=NULL, Int_t ncol=2, Int_t nbins=100, Int_t minx=-1, Int_t maxx=-1)
{
  TNtuple *ntuple; 
  TFile *f = new TFile("basic.root","RECREATE");
  c1 = new TCanvas("c1","fit with stretched exponential",10,10,700,500);
  TAxis *xax, *yax;
  Int_t npar, nlines, i;
  Float_t minX, maxX, x, y, *xarr, *yarr;
  
  c1->SetFillColor(33);
  c1->SetFrameFillColor(41);
  c1->SetGrid();
  if (ncol==2)
    c1->Divide(1,2);
  if (fileName==NULL)
    {
      printf("You have to supply the filename!\n");
      printf("Usage: root.exe 'getHisto(<filename>,<nbins>,<minx>,<maxy>)\n");
      printf("<filename>: file containing data to fit\n");
      printf("<nbins>: number of bins\n");
      printf("<minx>: minimum x for binning\n");
      printf("<maxx>: maximum x for binning\n");
      exit(-1);
    }
  //printf("fileName=%s\n", fileName);
  ntuple = readascii(fileName, &nlines, ncol);
  xarr = new Float_t[nlines];
  yarr = new Float_t[nlines];
  ntuple->SetBranchAddress("x", &x);
  for (i=0; i<nlines; i++) 
    {
      ntuple->GetEntry(i);
      xarr[i] = x;
      if (i==0 || x < minX)
	minX = x;
      if (i==0 || x > maxX)
	maxX = x;
    }
  if (ncol == 2)
    {
      ntuple->SetBranchAddress("y", &y);
      for (i=0; i<nlines; i++) 
	{
	  ntuple->GetEntry(i);
	  yarr[i] = y;
	}
    }
  if (minx == maxx  && maxx == -1)
    {
      minx = minX;
      maxx = maxX;
    }
  if (ncol ==2)
    c1->cd(1);	 
  h1 = new TH1F("Phi1","Histogram", nbins, minx, maxx);
  for (i = 0; i < nlines; i++)
    h1->Fill(xarr[i]);
 
  //h1->SetStats(kFALSE);
  h1->SetFillColor(kRed);
  //gStyle->SetHistLineColor(4);
  h1->Draw("bar3");
  if (ncol==2)
    {
      c1->cd(2);	 
      h2 = new TH1F("Phi2","Histogram", nbins, minx, maxx);
      for (i = 0; i < nlines; i++)
	h2->Fill(yarr[i]);
      h2->SetFillColor(kGreen);
      //gStyle->SetHistLineColor(6);
      h2->Draw("bar3"); 
    }
#if 1
  TLegend *legend = new TLegend(0.8,0.45,1.0,0.65);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(Phi1,"Phi1");
  if (ncol==2)
    legend->AddEntry(Phi2,"Phi2");
  legend->Draw();
#endif
  f->Write();

}
