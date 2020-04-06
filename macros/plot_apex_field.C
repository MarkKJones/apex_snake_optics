#include <iostream.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
void plot_dipole_field() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
  //
TFile *fdip =  new TFile("rootfiles/.root");
TTree *tdip = (TTree*)fdip->Get("ntuple");
//
 float_t x,y,z,bx,by,bz;
 tdip->SetBranchAddress("x",&x);
 tdip->SetBranchAddress("y",&y);
 tdip->SetBranchAddress("z",&z);
 tdip->SetBranchAddress("bx",&bx);
 tdip->SetBranchAddress("by",&by);
 tdip->SetBranchAddress("bz",&bz);
 //
 const int nx = 61;
 const int ny = 126;
 const int nz = 161;
 TGraph *bxvzpos[nx] ;
 TGraph *bxvzpos_yp30[nx] ;
 TGraph *bxvzpos_ym30[nx] ;
 TGraph *bzvxpos[nx] ;
 TGraph *bzvxpos_yp30[nx] ;
 TGraph *bzvxpos_ym30[nx] ;
 for ( int i = 0; i < nx ; i++) {
 bxvzpos[i]= new TGraph(nz);
 bxvzpos_yp30[i]= new TGraph(nz);
 bxvzpos_ym30[i]= new TGraph(nz);
 }
 for ( int i = 0; i < nz ; i++) {
 bzvxpos[i]= new TGraph(nx);
 bzvxpos_yp30[i]= new TGraph(nx);
 bzvxpos_ym30[i]= new TGraph(nx);
 }
 //
 int i=0;
      for (int ixx = 0; ixx < nx; ixx++) {
      for (int iyy = 0; iyy < ny; iyy++) {
      for (int izz = 0; izz < nz; izz++) {
        tdip->GetEntry(i);
	i++;
	if ( y == 0 ) { 
	  bxvzpos[ixx]->SetPoint(izz,z,bx/10000.);
	  bzvxpos[izz]->SetPoint(ixx,x,bz/10000.);
        }
	if ( y == -30 ) { 
	  bxvzpos_ym30[ixx]->SetPoint(izz,z,bx/10000.);
	  bzvxpos_ym30[izz]->SetPoint(ixx,x,bz/10000.);
        }
	if ( y == 30 ) { 
	  bxvzpos_yp30[ixx]->SetPoint(izz,z,bx/10000.);
	  bzvxpos_yp30[izz]->SetPoint(ixx,x,bz/10000.);
        }
      }
      }	
      }
      //
    TCanvas *cdip = new TCanvas("cdip"," Dipole field",800,800);
    cdip->Divide(1,1);
    cdip->cd(1);
    bxvzpos[0]->Draw("AL*");
      for (int ixx = 1; ixx < nx; ixx++) {
       bxvzpos[ixx]->Draw("L*");
      }
      for (int ixx = 0; ixx < nx; ixx++) {
       bxvzpos_yp30[ixx]->Draw("L*");
       bxvzpos_ym30[ixx]->Draw("L*");
       bxvzpos_ym30[ixx]->SetLineColor(2);
      }
      //
    TCanvas *cdip2 = new TCanvas("cdip2"," Dipole field",800,800);
    cdip2->Divide(1,1);
    cdip2->cd(1);
    bzvxpos[0]->Draw("AL*");
      for (int ixx = 1; ixx < nx; ixx++) {
       bzvxpos[ixx]->Draw("L*");
      }
      for (int ixx = 0; ixx < nx; ixx++) {
       bzvxpos_yp30[ixx]->Draw("L*");
       bzvxpos_ym30[ixx]->Draw("L*");
       bzvxpos_ym30[ixx]->SetLineColor(2);
      }
    //
  // end brace
}
