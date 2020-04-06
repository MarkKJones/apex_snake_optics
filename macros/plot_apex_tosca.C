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
void plot_apex_tosca() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 TString fstart;
    Double_t fac;
    Double_t step_size,minz,maxz; Int_t nbinz;
    Double_t target_cent=0.0;
    //Double_t xcut[3]={15.0,15.0,15.0}; // cm x horizontal steps of 5mm, -15 to 15cm
      Double_t ycut[3]={0.0,4.0,8.0}; //cm y vertical steps of 5mm, 
      static const  Int_t nxcut=6;
      Double_t start=-71;
      Double_t end=nxcut*10+start;
      Double_t step = (end-start)/nxcut;
      Double_t xcut[nxcut];
      for (Int_t nx=0;nx<nxcut;nx++) {
	xcut[nx] = (nx)*step+start;
        cout << xcut[nx] << endl;
      }
     fstart="apex_tosca";
     fac=1.0;
     minz=0; //cm
     maxz=40;
     step_size=.005; // along y direction in meters
     nbinz=((maxz-minz)/(step_size*100));
    //
TFile *ffield =  new TFile("rootfiles/"+fstart+".root");
TFile *hfield = new TFile("rootfiles/"+fstart+"_hist.root","recreate");
TTree *tfield = (TTree*)ffield->Get("ntuple");
//
 int nentries; //number of entries in file
 Float_t x,y,z,bx,by,bz; // position in cm and field in Gauss
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
	hbz[nx] = new TH1F(Form("hbz_%d",nx),"  ;X (cm) ; Blong (Gauss) ",nbinz,minz,maxz);
      hbx[nx] = new TH1F(Form("hbx_%d",nx),"  ;X (cm) ; Bhorz (Gauss) ",nbinz,minz,maxz);
      hby[nx] = new TH1F(Form("hby_%d",nx),"  ;X (cm) ; Bvert (Gauss) ",nbinz,minz,maxz);
      }
       //
      for (int ie = 0; ie < nentries; ie++)
      {
            tfield->GetEntry(ie);
            if ( TMath::Abs(y-ycut[0]) <= .25) {
                   for (Int_t nx=0;nx<nxcut;nx++) {
		     if ( TMath::Abs(z-xcut[nx]) <= .25) {
             vfield=by*fac/10000.;
	     /*
	      if (ny==0) {
		if (TMath::Abs(vfield) > MaxField) MaxField=TMath::Abs(vfield);
	      deriv= step_size*vfield;
              IntBdl = IntBdl + deriv;
	      }
	     */
              hby[nx]->Fill(x,vfield);
              hbx[nx]->Fill(x,bx*fac/10000.);
              hbz[nx]->Fill(x,bz*fac/10000.);
		     }
		   }
	    }
       }
      //
      printf("For x= %4.2f and z= %4.2f \n  Int Bz = %5.4f Maxfieldz = %5.4f  L_eff = %5.4f \n",xcut[0],ycut[0],IntBdl/1.,MaxField/1.,IntBdl/MaxField);
      //
TLegend *myLegendz=new TLegend(0.4,0.7,.8,.9,"");
TLegend *myLegendx=new TLegend(0.4,0.8,.98,.8,"");
TLegend *myLegendy=new TLegend(0.4,0.8,.98,.8,"");
TCanvas *cplot = new TCanvas("cplot"," ",1000,1000);
 TMultiGraph *mgr ;
 mgr= new TMultiGraph("mgr_0","; Distance (cm); Field ");
 cplot->Divide(1,1);
 cplot->cd(1);
      for (Int_t ny=0;ny<nxcut;ny++) {
	if (ny==0) {
         hby[ny]->SetMinimum(0.);
         hby[ny]->SetMaximum(2.0);
         hby[ny]->Draw("hist");	  
	  } else {
           hby[ny]->Draw("hist same");
	   Int_t ncol=ny+1;
	   if (ncol>=5) ncol=ncol+1;
            hby[ny]->SetLineColor(ncol);
	  }
 myLegendz->AddEntry(hby[ny],Form("long = %4.2f cm",xcut[ny]),"l");
	  
      }
 myLegendz->SetTextSize(.05);
 myLegendz->Draw();
 //

    //
}
