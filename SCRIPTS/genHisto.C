#include "Riostream.h"
#include "TH2.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
gROOT->Reset();
TH1F *h1, *h2;
TCanvas *c1;
TNtuple* readascii(char* name, Int_t *nlines)
{
  //   example of macro to read data from an ascii file and
  //   create a root file with an histogram and an ntuple.
  //ifstream in;
  // we assume a file basic.dat in the current directory
  // this file has 3 columns of float data
  //in.open(name);
  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y");
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
      fscanf(f, "%f %f\n", &x, &y);
      ntuple->Fill(x,y);
      (*nlines)++;
    }
  fclose(f);
#endif
  printf(" found %d pointsn \n",*nlines);

  //in.close();
  return ntuple;
}

void genHisto(char *fileName=NULL, Int_t nbins=100, Int_t minx=-1, Int_t maxx=-1, Int_t col=0)
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
  c1->Divide(1,2);
  if (fileName==NULL)
    {
      printf("You have to supply the filename!\n");
      printf("Usage: root.exe 'getHisto(<filename>,<logx>,<logy>,<xpnts>,<ypnts>)\n");
      printf("<filename>: file containing data to fit\n");
      printf("<logx>: 1 log scale for x-axis 0 = linear\n");
      printf("<logy>: 1 log scale for y-axis 0 = linear\n");
      printf("<xpnts>: number of points along x-axis for histogram mesh\n");
      printf("<ypnts>: number of points along y-axis for histogram mesh\n");
      printf("<xpnts> and <ypnts> affect the zoom\n");
      exit(-1);
    }
  //printf("fileName=%s\n", fileName);
  ntuple = readascii(fileName, &nlines);
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
  ntuple->SetBranchAddress("y", &y);
  for (i=0; i<nlines; i++) 
    {
      ntuple->GetEntry(i);
      yarr[i] = y;
    }
  if (minx == maxx  && maxx == -1)
    {
      minx = minX;
      maxx = maxX;
    }
  c1->cd(1);	 
  h1 = new TH1F("Phi1","Histogram", nbins, minx, maxx);
  for (i = 0; i < nlines; i++)
    h1->Fill(xarr[i]);
 
  //h1->SetStats(kFALSE);
  h1->Draw();
  
  c1->cd(2);	 
  h2 = new TH1F("Phi2","Histogram", nbins, minx, maxx);
  for (i = 0; i < nlines; i++)
    h2->Fill(yarr[i]);
  h2->Draw(); 
#if 1
  TLegend *legend = new TLegend(0.8,0.45,1.0,0.65);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(Phi1,"Phi1","lp");
  legend->AddEntry(Phi2,"Phi2","l");
  legend->Draw();
#endif
  f->Write();

}
