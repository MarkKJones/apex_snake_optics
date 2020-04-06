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
      Double_t step=55./nxcut;
      Double_t xcut[nxcut];
      for (Int_t nx=0;nx<nxcut;nx++) {
	xcut[nx] = nx*step;
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
      TH1F *hbz[nxcut];
      TH1F *hby[nxcut];
      TH1F *hbx[nxcut];
      for (Int_t nx=0;nx<nxcut;nx++) {
	hbz[nx] = new TH1F(Form("hbz_%d",nx),"  ;Ysnake (cm) ; Bhorz (Gauss) ",nbiny,miny,maxy);
      hbx[nx] = new TH1F(Form("hbx_%d",nx),"  ;Ysnake (cm) ; BVert (Gauss) ",nbiny,miny,maxy);
      hby[nx] = new TH1F(Form("hby_%d",nx),"  ;Ysnake (cm) ; BLong (Gauss) ",nbiny,miny,maxy);
      }
       for (int ie = 0; ie < nentries; ie++)
      {
            tfield->GetEntry(ie);
            if ( TMath::Abs(z-zcut[0]) <= .25*10.) {
                   for (Int_t nx=0;nx<nxcut;nx++) {
		     if ( TMath::Abs(x-xcut[nx]) <= .25*10.) {
             vfield=bx*fac;
	      if (nx==0) {
		if (TMath::Abs(vfield) > MaxField) MaxField=TMath::Abs(vfield);
	      deriv= step_size*vfield;
              IntBdl = IntBdl + deriv;
	      cout << y << " " << x << " " << z << " " << by  << " " << bx << " " << bz << endl;
	      }
              hby[nx]->Fill((-y-target_cent)/10.,by*fac/1.);
              hbx[nx]->Fill((-y-target_cent)/10.,bx*fac/1.);
              hbz[nx]->Fill((-y-target_cent)/10.,bz*fac/1.);
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
 myLegendz->AddEntry(hbz[nx],Form("horz = %4.2f cm vert= %4.2f cm",zcut[0]/10.,xcut[nx]/10.),"l");
	  
      }
myLegendz->SetTextSize(.05);
 myLegendz->Draw();
 cplot->cd(2);
     for (Int_t nx=0;nx<nxcut;nx++) {
	if (nx==0) {
         hbx[nx]->SetMinimum(0.);
         hbx[nx]->SetMaximum(1.0);
         hbx[nx]->Draw("hist");	  
	  } else {
           hbx[nx]->Draw("hist same");
            hbx[nx]->SetLineColor(nx+1);
	  }
      }
cplot->cd(3);
      for (Int_t nx=0;nx<nxcut;nx++) {
	if (nx==0) {
         hby[nx]->SetMinimum(-.4);
         hby[nx]->SetMaximum(.4);
         hby[nx]->Draw("hist");	  
	  } else {
           hby[nx]->Draw("hist same");
            hby[nx]->SetLineColor(nx+1);
	  }
      }

    //
}
