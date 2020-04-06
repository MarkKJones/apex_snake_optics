#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
void plot_apex_data(Int_t nrun=4648) {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 //
   TString inputroot;
   inputroot="data_rootfiles/new_apex_4648.root";
   TString outputhist;
   outputhist= "rootfiles/apex_data_4648_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t  ntr;
   tsimc->SetBranchAddress("R.tr.n",&ntr);
 Double_t  xfp[100];
   tsimc->SetBranchAddress("R.tr.x",&xfp);
 Double_t  yfp[100];
   tsimc->SetBranchAddress("R.tr.y",&yfp);
 Double_t  xpfp[100];
   tsimc->SetBranchAddress("R.tr.th",&xpfp);
 Double_t  ypfp[100];
   tsimc->SetBranchAddress("R.tr.ph",&ypfp);
 Double_t  xptg[100];
   tsimc->SetBranchAddress("R.tr.tg_th",&xptg);
 Double_t  yptg[100];
   tsimc->SetBranchAddress("R.tr.tg_ph",&yptg);
 Double_t  delta[100];
   tsimc->SetBranchAddress("R.tr.tg_dp",&delta);
   //
  TH2F *hxyfp = new TH2F("hxyfp","xfp vs yfp (cm)",20,-5.,5.,100,-100.,100.);
  HList.Add(hxyfp);
 TH2F *hxpypfp = new TH2F("hxpypfp","xpfp vs ypfp (rad)",100,-.02,.02,100,-.2,.2);
  HList.Add(hxpypfp);
 TH2F *hxpxfp = new TH2F("hxpxfp","xfp vs xpfp ",90,-.15,.15,100,-100.,100.);
  HList.Add(hxpxfp);
 TH2F *hypyfp = new TH2F("hypyfp","yfp vs ypfp ",100,-.02,.02,100,-5.,5.);
 HList.Add(hypyfp);
 TH2F *hxptg_yptg_new = new TH2F("hxptg_yptg_new","New xptg vs yptg  ",100,-.05,.05,100,-.05,.05);
 HList.Add(hxptg_yptg_new);
 TH2F *hxptg_yptg_orig = new TH2F("hxptg_yptg_orig","Orig xptg vs yptg  ",100,-.05,.05,100,-.05,.05);
 HList.Add(hxptg_yptg_orig);
 static const Int_t ndelcut=10;
 Double_t delcutlo[ndelcut];
 Double_t delcuthi[ndelcut];
 Double_t delstart=-.5;
 Double_t delstop=5.;
 Double_t delstep=(delstop-delstart)/ndelcut;
 TH2F *hxptg_yptg_orig_delcut[ndelcut];
 TH2F *hxptg_yptg_new_delcut[ndelcut];
 TH2F *hxs_ys_orig_delcut[ndelcut];
	for (Int_t i = 0; i < ndelcut; i++) {
	  delcutlo[i]= delstart+ i*delstep;
	  delcuthi[i]=delcutlo[i]+ delstep;
	  hxptg_yptg_orig_delcut[i]   = new TH2F(Form("hxptg_yptg_orig_delcut_%d",i),"Orig xptg vs yptg delcut ",100,-.05,.05,100,-.05,.05);
 HList.Add(hxptg_yptg_orig_delcut[i]);
	  hxptg_yptg_new_delcut[i]   = new TH2F(Form("hxptg_yptg_new_delcut_%d",i),"New xptg vs yptg delcut ",100,-.05,.05,100,-.05,.05);
 HList.Add(hxptg_yptg_new_delcut[i]);
 hxs_ys_orig_delcut[i]   = new TH2F(Form("hxs_ys_orig_delcut_%d",i),Form("Orig xs vs ys delcut %d ",i),100,-3,3,100,-5,5);
 HList.Add(hxs_ys_orig_delcut[i]);
	}
 //
  string oldcoeffsfilename="matrix-files/hms_newfit5_Rhrs_save.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  int num_recon_terms_old;
  vector<Double_t> xptarcoeffs_old;
  vector<Double_t> yptarcoeffs_old;
  vector<Double_t> ytarcoeffs_old;
  vector<Double_t> deltacoeffs_old;
  vector<Int_t> xfpexpon_old;
  vector<Int_t> xpfpexpon_old;
  vector<Int_t> yfpexpon_old;
  vector<Int_t> ypfpexpon_old;
  vector<Int_t> xtarexpon_old;
  TString currentline;
  num_recon_terms_old = 0;
  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    //    cout << currentline.Data() << endl;
    //extract the coeffs and exponents from the line:

    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs_old.push_back(sc1.Atof());
    ytarcoeffs_old.push_back(sc2.Atof());
    yptarcoeffs_old.push_back(sc3.Atof());
    deltacoeffs_old.push_back(sc4.Atof());
    int expontemp[5];

    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

    xfpexpon_old.push_back(expontemp[0]);
    xpfpexpon_old.push_back(expontemp[1]);
    yfpexpon_old.push_back(expontemp[2]);
    ypfpexpon_old.push_back(expontemp[3]);
    xtarexpon_old.push_back(expontemp[4]);
    num_recon_terms_old++;
  }
 //
  Double_t sieveDis=31.23*2.54;
Long64_t nentries = tsimc->GetEntries();
	for (Int_t i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (ntr==1) {
		  hxyfp->Fill(yfp[0]*100,xfp[0]*100);
		  hxpxfp->Fill(xpfp[0],xfp[0]*100);
		  hxpypfp->Fill(ypfp[0],xpfp[0]);
		  hypyfp->Fill(ypfp[0],yfp[0]*100);
	         for (Int_t i = 0; i < ndelcut; i++) {
 		     hxptg_yptg_orig->Fill(yptg[0],xptg[0]);
		   if (delta[0]<=delcuthi[i]/100.&&delta[0]>delcutlo[i]/100.) {
                     hxptg_yptg_orig_delcut[i]->Fill(yptg[0],xptg[0]);
                    hxs_ys_orig_delcut[i]->Fill(yptg[0]*sieveDis,xptg[0]*sieveDis);
		}
		 }
          Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
	  Double_t xtar=0,etemp=0.0;
                for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
        	 etemp= 
	  TMath::Power( xfp[0], xfpexpon_old[icoeffold] ) * 
	  TMath::Power( yfp[0], yfpexpon_old[icoeffold] ) * 
	  TMath::Power( xpfp[0], xpfpexpon_old[icoeffold] ) * 
	  TMath::Power( ypfp[0], ypfpexpon_old[icoeffold] ) * 
	  TMath::Power( xtar, xtarexpon_old[icoeffold] );
        	deltatemp += deltacoeffs_old[icoeffold] * etemp;
        	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	        yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	        xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
	           }
		hxptg_yptg_new->Fill(yptartemp,xptartemp);
	         for (Int_t i = 0; i < ndelcut; i++) {
		   if (delta[0]<=delcuthi[i]/100.&&delta[0]>delcutlo[i]/100.) {
                     hxptg_yptg_new_delcut[i]->Fill(yptartemp,xptartemp);
		   }		 }
}
	}
	//
    TCanvas *c3 = new TCanvas("c3","Focal Plane",800,800);
    TPad *c3_1 = new TPad("c3_1", "c3_1",0.01,0.51,0.49,0.99);
    TPad *c3_3 = new TPad("c3_3", "c3_3",0.51,0.51,0.99,0.99);
    TPad *c3_2 = new TPad("c3_2", "c3_2",0.01,0.01,0.49,0.49);
    TPad *c3_4 = new TPad("c3_4", "c3_4",0.51,0.01,0.99,0.49);
    c3_1->Draw();
    c3_2->Draw();
    c3_3->Draw();
    c3_4->Draw();
    c3_1->cd();
    c3_1->SetGridx(1);
    c3_1->SetGridy(1);
    hxyfp->Draw("colz");
    c3_2->cd();
    c3_2->SetGridx(1);
    c3_2->SetGridy(1);
    hxpypfp->Draw("colz");
    c3_3->cd();
    c3_3->SetGridx(1);
    c3_3->SetGridy(1);
    hxpxfp->Draw("colz");
    c3_4->cd();
    c3_4->SetGridx(1);
    c3_4->SetGridy(1);
    hypyfp->Draw("colz");
  	//
    TCanvas *ctar = new TCanvas("ctar","Target",800,800);
    TPad *ctar_1 = new TPad("ctar_1", "ctar_1",0.01,0.51,0.49,0.99);
    TPad *ctar_3 = new TPad("ctar_3", "ctar_3",0.51,0.51,0.99,0.99);
    TPad *ctar_2 = new TPad("ctar_2", "ctar_2",0.01,0.01,0.49,0.49);
    TPad *ctar_4 = new TPad("ctar_4", "ctar_4",0.51,0.01,0.99,0.49);
    ctar_1->Draw();
    ctar_2->Draw();
    ctar_3->Draw();
    ctar_4->Draw();
    ctar_1->cd();
    ctar_1->SetGridx(1);
    ctar_1->SetGridy(1);
    ctar_1->SetLogz();
    hxptg_yptg_orig->Draw("colz");
    ctar_2->cd();
    ctar_2->SetGridx(1);
    ctar_2->SetGridy(1);
   ctar_2->SetLogz();
    hxptg_yptg_orig_delcut[0]->Draw("colz");
    ctar_3->cd();
    ctar_3->SetGridx(1);
    ctar_3->SetGridy(1);
   ctar_3->SetLogz();
    hxptg_yptg_new->Draw("colz");
    ctar_4->cd();
    ctar_4->SetGridx(1);
    ctar_4->SetGridy(1);
   ctar_4->SetLogz();
    hxptg_yptg_new_delcut[0]->Draw("colz");
    
 //
 TFile hsimc(outputhist,"recreate");
	HList.Write();
}
