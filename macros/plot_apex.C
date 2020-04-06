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
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TGraph.h>
void plot_apex_snake() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 TString fstart;
    Double_t fac;
    Double_t step_size,miny,maxy; Int_t nbiny;
    Double_t target_cent=0.0;
    Double_t zcut[3]={150.,150.,150.}; // mm z_snake is horizontal
    static const  Int_t nxcut=11; // x_snake is vertical
      Double_t step=5.5/nxcut;
      Double_t xcut[nycut];
      for (Int_t nx=0;nx<nxcut;ny++) {
	xcut[ny] = nx*step;
      }
     fstart="apex_snake";
     fac=1.0;
     miny=-75; //cm
     maxy=175;
     step_size=.02; // along y direction in meters
     nbiny=((maxy-miny)/(step_size*100));
    //
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
      nentries = (int)tfield->GetEntries();
      Double_t deriv,vfield;
      Double_t IntBdl=0.;
      Double_t MaxField=0.;
    //
      TH1F *hbz[nycut];
      TH1F *hby[nycut];
      TH1F *hbx[nycut];
      for (Int_t ny=0;ny<nycut;ny++) {
	hbz[ny] = new TH1F(Form("hbz_%d",ny),"  ;Ysnake (cm) ; Bhorz (Gauss) ",nbinz,minz,maxz);
      hbx[ny] = new TH1F(Form("hbx_%d",ny),"  ;Ysnake (cm) ; BVert (Gauss) ",nbinz,minz,maxz);
      hby[ny] = new TH1F(Form("hby_%d",ny),"  ;Ysnake (cm) ; BLong (Gauss) ",nbinz,minz,maxz);
      }
       for (int ie = 0; ie < nentries; ie++)
      {
            tfield->GetEntry(ie);
            if ( TMath::Abs(z-zcut[0]) <= .25*10.) {
                   for (Int_t nx=0;nx<nxcut;nx++) {
		     if ( TMath::Abs(x-xcut[nx]) <= .25*10.) {
             vfield=by*fac;
	      if (ny==0) {
		if (TMath::Abs(vfield) > MaxField) MaxField=TMath::Abs(vfield);
	      deriv= step_size*vfield;
              IntBdl = IntBdl + deriv;
	      }
              hby[ny]->Fill(z-target_cent,vfield);
              hbx[ny]->Fill(z-target_cent,bx*fac/1.);
              hbz[ny]->Fill(z-target_cent,bz*fac/1.);
		     }
		   }
	    }
	    //
       }
      //
      printf("For x= %4.2f and z= %4.2f \n  Int Bz = %5.4f Maxfieldz = %5.4f  L_eff = %5.4f \n",xcut[0]/10,zcut[0]/10,IntBdl/1.,MaxField/1.,IntBdl/MaxField);
      //
TLegend *myLegendz=new TLegend(0.4,0.7,.8,.9,"");
TLegend *myLegendx=new TLegend(0.4,0.8,.98,.8,"");
TLegend *myLegendy=new TLegend(0.4,0.8,.98,.8,"");
TCanvas *cplot = new TCanvas("cplot"," ",1000,1000);
 TMultiGraph *mgr ;
 mgr= new TMultiGraph("mgr_0","; Distance (cm); Field ");
 cplot->Divide(1,3);
 cplot->cd(1);
      for (Int_t nx=0;nx<nxcut;nx++) {
	if (nx==0) {
         hbz[nx]->SetMinimum(-.4);
         hbz[nx]->SetMaximum(.4);
         hbz[nx]->Draw("hist");	  
	  } else {
           hbz[nx]->Draw("hist same");
            hbz[nx]->SetLineColor(nx+1);
	  }
 myLegendz->AddEntry(hbz[nx],Form("horz = %4.2f cm vert= %4.2f cm",xcut[0],ycut[nx]),"l");
	  
      }
myLegendz->SetTextSize(.05);
 myLegendz->Draw();
 cplot->cd(2);
     for (Int_t nx=0;nx<nxcut;nx++) {
	if (nx==0) {
         hbx[nx]->SetMinimum(-.4);
         hbx[nx]->SetMaximum(.4);
         hbx[nx]->Draw("hist");	  
	  } else {
           hbx[nx]->Draw("hist same");
            hbx[nx]->SetLineColor(nx+1);
	  }
      }
cplot->cd(3);
      for (Int_t nx=0;nx<nxcut;nx++) {
	if (nx==0) {
         hby[nx]->SetMinimum(0.);
         hby[nx]->SetMaximum(1.0);
         hby[nx]->Draw("hist");	  
	  } else {
           hby[nx]->Draw("hist same");
            hby[nx]->SetLineColor(nx+1);
	  }
      }

    //
}
