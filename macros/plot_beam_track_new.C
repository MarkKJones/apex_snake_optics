#include <iostream>
#include <fstream>
#include <iomanip>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TGraph.h>
#include <TPolyLine.h>
#include <TLegend.h>
void plot_beam_track_new() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 //
  TString fstart="track_tz095";
  Int_t npt[4]={118,648,1086,1589};// last entry number for each magnet region
 //  TString fstart="track_tz0";
 //Int_t npt[4]={374,766,1130,1443};// last entry number for each magnet region
 // TString fstart="track_tz080";
 // Int_t npt[4]={344,872,1292,1895};// last entry number for each magnet region
 // TString fstart="track_tz110";
 //  Int_t npt[4]={342,866,1286,1887};// last entry number for each magnet region
TFile *ffield =  new TFile("rootfiles/"+fstart+".root");
TFile *hfield = new TFile("rootfiles/"+fstart+"_hist.root","recreate");
TTree *tfield = (TTree*)ffield->Get("ntuple");
//
 int nentries; //number of entries in file
 Float_t x,y,z,bx,by,bz; // position in mm and field in Tesla
 //
 tfield->SetBranchAddress("x",&x);
 tfield->SetBranchAddress("y",&y);
 tfield->SetBranchAddress("z",&z);
 tfield->SetBranchAddress("bx",&bx);
 tfield->SetBranchAddress("by",&by);
 tfield->SetBranchAddress("bz",&bz);
 //
 Double_t ang[4]={7.0,8.5,8.5,8.5}; // angles of rotation for 4 magnet regions
 Double_t ycen[4]={1760.,4136.7,2660.,2650.}; // y center loc with q2 and q3 relative to q1 
 Double_t zcen[4]={174.5,520.,0.,0.,};// z center loc with q2 and q3 relative to q1 
 // Int_t npt[4]={378,838,1362,1673};// last entry number for each magnet region
   Double_t yoff[4],zoff[4];
 Double_t bdl[4]={0.,0.,0.,0.};
 Double_t dis[4]={0.,0.,0.,0.};
            for (int n = 0; n < 4; n++)
      {
            if ( n <= 1 ) yoff[n]=ycen[n];
            if ( n > 1) yoff[n]=yoff[n-1]+ycen[n]*TMath::Cos(ang[n]/180.*3.14159);
            if ( n <= 1 ) zoff[n]=zcen[n];
            if ( n > 1) zoff[n]=zoff[n-1]+ycen[n]*TMath::Sin(ang[n]/180.*3.14159);
      }
      nentries = (int)tfield->GetEntries();
Double_t ybeam[10000],zbeam[10000],bxsave[10000];
      Int_t nreg=0;
            for (int ie = 0; ie < nentries; ie++)
      {
            tfield->GetEntry(ie);
	    if (ie == npt[nreg]) printf(" nreg = %i x = %f y = %f z = %f  yoff = %f \n",nreg,x,y,z,yoff[nreg]) ;
            if (ie == npt[nreg]) nreg++;
            ybeam[ie] = ((y*TMath::Cos(ang[nreg]/180.*3.14159)-z*TMath::Sin(ang[nreg]/180.*3.14159))+yoff[nreg])/10.;
            zbeam[ie] = ((y*TMath::Sin(ang[nreg]/180.*3.14159)+z*TMath::Cos(ang[nreg]/180.*3.14159))+zoff[nreg])/10.;
            bxsave[ie]=bx*10000.; // convert to gauss
            if ( ie >= 1) {
	      bdl[nreg]=bdl[nreg]+bxsave[ie]*(ybeam[ie]-ybeam[ie-1]);
            dis[nreg]=dis[nreg]+(ybeam[ie]-ybeam[ie-1]);
            }
            //cout  << y << " " << z << " " << ybeam[ie] << " " << zbeam[ie] << " " << bxsave[ie] << endl;
      }
            for (int n = 0; n < 4; n++) { 
              cout << n << " bdl = " << bdl[n] << " Gauss cm " << endl ;
            }
	    //
TLegend *myLegend=new TLegend(0.4,0.8,.98,.9,"");
TCanvas *cplot = new TCanvas("cplot"," ",800,800);
TGraph  *gr=new TGraph(nentries,ybeam,bxsave);
TGraph  *gr2=new TGraph(nentries,ybeam,zbeam);
 cplot->Divide(1,2);
 cplot->cd(1);

 gr->Draw("AC*");
 gr->SetTitle(" Vertical field vs beam distance");
   gr->GetYaxis()->SetTitle("Vertical field (Gauss) ");
   gr->GetXaxis()->SetTitle("Beam position (cm)"); 
 myLegend->AddEntry(gr,"SHMS =  5.5 deg","l");
 myLegend->SetTextSize(.025);
 myLegend->Draw();
 cplot->cd(2);

 gr2->SetTitle(" Horizontal beam position vs beam distance");
   gr2->GetYaxis()->SetTitle("Horizontal position (cm) ");
   gr2->GetXaxis()->SetTitle("Beam position (cm)"); 
 gr2->Draw("AC*");
 gr2->Fit("pol1","","",800.,850.);
 //
}
