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
void plot_apex_track() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 //
  TString fstart="track";
 TFile *ffield =  new TFile("rootfiles/"+fstart+".root");
TFile *hfield = new TFile("rootfiles/"+fstart+"_hist.root","recreate");
TTree *tfield = (TTree*)ffield->Get("track");
 int nentries; //number of entries in file
 Float_t trnum,pmom,x,y,z,bx,by,bz; // position in mm and field in Tesla
 //
 tfield->SetBranchAddress("trnum",&trnum);
 tfield->SetBranchAddress("tr_mom",&pmom);
 tfield->SetBranchAddress("tr_x",&x);
 tfield->SetBranchAddress("tr_y",&y);
 tfield->SetBranchAddress("tr_z",&z);
 tfield->SetBranchAddress("tr_bx",&bx);
 tfield->SetBranchAddress("tr_by",&by);
 tfield->SetBranchAddress("tr_bz",&bz);
 //
 // tracking three momentum
 //
      nentries = (int)tfield->GetEntries();
      float_t prev_trnum=1;
      vector<double > vecX;
      vector<double > vecDZDY;
      vector<double > vecY;
      vector<double > vecZ;
      vector<double > vecBX;
      vector<double > vecBY;
      vector<double > vecBZ;
      vector<TGraph*> h_Y_X;
      vector<TGraph*> h_Y_DZDY;
      vector<TGraph*> h_Y_Z;
      vector<TGraph*> h_Y_BX;
      vector<TGraph*> h_Y_BY;
      vector<TGraph*> h_Y_BZ;
      Double_t max_Z=0;
      Double_t prev_z=0,prev_y=0;
      Double_t angz;
            for (int ie = 0; ie < nentries; ie++) {
            tfield->GetEntry(ie);
	    if (trnum != prev_trnum || ie == nentries-1) {
	      TGraph *temp = new TGraph(vecX.size(),&vecY[0],&vecX[0]) ;
	      h_Y_X.push_back(temp);
	      temp = new TGraph(vecZ.size(),&vecY[0],&vecZ[0]) ;
	      h_Y_Z.push_back(temp);
	      temp = new TGraph(vecY.size(),&vecY[0],&vecBX[0]) ;
	      h_Y_BX.push_back(temp);
	      auto max_elem=max_element(vecZ.begin(),vecZ.end());
	      if (*max_elem>max_Z) max_Z=*max_elem;
	      temp = new TGraph(vecY.size(),&vecY[0],&vecBY[0]) ;
	      h_Y_BY.push_back(temp);
	      temp = new TGraph(vecY.size(),&vecY[0],&vecBZ[0]) ;
	      h_Y_BZ.push_back(temp);
	      temp = new TGraph(vecY.size(),&vecY[0],&vecDZDY[0]) ;
	      h_Y_DZDY.push_back(temp);
	      vecX.clear();
	      vecY.clear();
	      vecZ.clear();
	      vecDZDY.clear();
	      vecBX.clear();
	      vecBY.clear();
	      vecBZ.clear();
	      prev_trnum=trnum;
	      vecX.push_back(x);
	      vecY.push_back(y);
	      vecZ.push_back(z);
	      vecDZDY.push_back(0.);
	      prev_z=z;
	      prev_y=y;
	      vecBX.push_back(bx);
	      vecBY.push_back(by);
	      vecBZ.push_back(bz);
	    } else {
	      vecX.push_back(x);
	      vecY.push_back(y);
	      vecZ.push_back(z);
	      if (ie==0) angz=0;
	      if (ie!=0) angz=(z-prev_z)/(y-prev_y);
	      vecDZDY.push_back(angz);
	      prev_z=z;
	      prev_y=y;
	      vecBX.push_back(bx);
	      vecBY.push_back(by);
	      vecBZ.push_back(bz);
	    }
	    }
	    //
	    TCanvas *canpos = new TCanvas("canpos","canpos",700,700);
	    canpos->Divide(1,3);
	    canpos->cd(1);
	    for (UInt_t ntr=0;ntr<h_Y_X.size();ntr++) {
	      h_Y_X[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_X[ntr]->Draw("AL");
	      if (ntr!=0) h_Y_X[ntr]->Draw("same");
	    }
	    canpos->cd(2);
	    for (UInt_t ntr=0;ntr<h_Y_Z.size();ntr++) {
	      h_Y_Z[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_Z[ntr]->Draw("AL");
	      if (ntr==0) h_Y_Z[ntr]->SetMaximum(max_Z*1.1);
	      if (ntr!=0) h_Y_Z[ntr]->Draw("same");
	    }
	    canpos->cd(3);
	    for (UInt_t ntr=0;ntr<h_Y_DZDY.size();ntr++) {
	      h_Y_DZDY[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_DZDY[ntr]->Draw("AL");
	      if (ntr!=0) h_Y_DZDY[ntr]->Draw("same");
	    }
	    //
	    TCanvas *canb = new TCanvas("canb","canb",700,700);
	    canb->Divide(1,3);
	    canb->cd(1);
	    for (UInt_t ntr=0;ntr<h_Y_BX.size();ntr++) {
	      h_Y_BX[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_BX[ntr]->Draw("AL");
	      if (ntr!=0) h_Y_BX[ntr]->Draw("same");
	    }
	    canb->cd(2);
	    for (UInt_t ntr=0;ntr<h_Y_BY.size();ntr++) {
	      h_Y_BY[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_BY[ntr]->Draw("AL");
	      if (ntr!=0) h_Y_BY[ntr]->Draw("same");
	    }
	    canb->cd(3);
	    for (UInt_t ntr=0;ntr<h_Y_BZ.size();ntr++) {
	      h_Y_BZ[ntr]->SetLineColor(ntr+1);
	      if (ntr==0) h_Y_BZ[ntr]->Draw("AL");
	      if (ntr!=0) h_Y_BZ[ntr]->Draw("same");
	    }
	    //
}
