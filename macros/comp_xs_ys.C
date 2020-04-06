#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void comp_xs_ys(Int_t nrun=4648) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 //
     outputpdf=Form("plots/run%d.pdf",nrun);
  TString inputroot;
   TFile *fhistroot;
     inputroot="rootfiles/apex_data_4648_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   //
 static const Int_t ndelcut=10;
 TH2F *hxs_ys_orig_delcut[ndelcut];
	for (Int_t i = 0; i < ndelcut; i++) {
	  hxs_ys_orig_delcut[i] = (TH2F*)fhistroot->Get(Form("hxs_ys_orig_delcut_%d",i));
	}
	//
	TCanvas *cxsys[ndelcut];
	for (Int_t i = 0; i < ndelcut; i++) {
	  cxsys[i] = new TCanvas(Form("can_%d",i),Form("can_%d",i),700,700);
	  cxsys[i]->Divide(1,1);
	  cxsys[i]->cd(1);
	  hxs_ys_orig_delcut[i]->Draw("colz");
	}	
	//
}
