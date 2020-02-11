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
void plot_apex_septum_q1(TString fstart) {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 Double_t cmom=0.837595,delta,ytar,ztemp,thcent=5./180.*3.14159;
 Double_t yptar_cent=0.08749; // tan(5) deg
 Double_t ypfp_cent=0.; // q1 is already rotated
   Double_t solang_fac=0.10*0.20*1000.;
TFile *fsnake =  new TFile("rootfiles/"+fstart+".root");
TFile *hsnake = new TFile("rootfiles/"+fstart+"_standard_hist.root","recreate");
//
//
 int nentries; //number of entries in file
Double_t evnum,epnum,xabs,yabs,zabs,cxabs,cyabs,czabs,mom,pathl,live,xrel,yrel,zrel,cxrel,cyrel,czrel;
 //
 TH1F *hEndPlane = new TH1F("hEndPlane"," Endplane Number",6,-.5,7.5);
 TH1F *hMom = new TH1F("hMom","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomp = new TH1F("hMomp","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomf = new TH1F("hMomf","Momentum (GeV)",50,.5,1.5);
 TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon ",40,-.12,.12);
 TH1F *hDelta = new TH1F("hDelta","Delta ",20,-.1,.1);
 TH1F *hDeltap = new TH1F("hDeltap","Delta ",20,-.1,.1);
 TH1F *hDeltaf = new TH1F("hDeltaf","Delta ",20,-.1,.1);
 TH1F *hDeltar = new TH1F("hDeltar",";Delta ;Solid angle (msr)  ",20,-.1,.1);
 TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
 TH1F *hytarp = new TH1F("hytarp","ytar (cm)",70,-7.,7.);
 TH1F *hytarf = new TH1F("hytarf","ytar (cm)",70,-7.,7.);
 TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
 TH1F *hxptarp = new TH1F("hxptarp","xptar ",100,-.12,.12);
 TH1F *hxptarf = new TH1F("hxptarf","xptar ",100,-.12,.12);
 TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
 TH1F *hyptarp = new TH1F("hyptarp","yptar ",100,-.1,.1);
 TH1F *hyptarf = new TH1F("hyptarf","yptar ",100,-.1,.1);
 TH2F *hxyfp = new TH2F("hxyfp","xfp vs yfp (cm)",100,-15.,15.,40,-20.,20.);
 TH2F *hxpypfp = new TH2F("hxpypfp","xpfp vs ypfp (rad)",100,-.1,.1,100,-.05,.05);
 TH2F *hxpyfp = new TH2F("hxpyfp","xpfp vs yfp ",100,0.,40.,100,-.1,.1);
 TH2F *hxpdel = new TH2F("hxpdel","xptar vs delta ",50,-.1,.1,55,-.1,.1);
 TH2F *hxpxfp = new TH2F("hxpxfp","xfp vs xpfp ",90,-.05,.55,100,-20.,20.);
 TH2F *hypyfp = new TH2F("hypyfp","ypfp;yfp ",100,-.1,.1,100,-10.,20.);
 TH2F *hypfpyptar = new TH2F("hypfpyfptar",";ypfp;yptar ",100,-.1,.1,100,-.1,.1);
 TH2F *hypfpxptar = new TH2F("hypfpxfptar",";ypfp;xptar ",100,-.1,.1,100,-.1,.1);
 TH2F *hyfpyptar = new TH2F("hyfpyfptar",";yptar;yfp ",100,-.1,.1,100,-10.,20.);
 TH2F *hypdel = new TH2F("hypdel","yptar vs delta ",100,-.1,.1,400,-.05,.05);
 TH2F *hypfpdel = new TH2F("hypfpdel",";ypfp;delta ",100,-.08,.08,400,-.06,.06);
 TH2F *hypytar = new TH2F("hypytar","yptar vs ytar;ytar (cm); yptar (rad) ",50,-6.,6.,40,-.04,.04);
 TH2F *hxpyp = new TH2F("hxpyp","xptar vs yptar;yptar;xptar ",40,-.04,.04,55,-.1,.1);
 TH2F *hxptarep = new TH2F("hxptarep","xptar vs endplane ",38,-.5,37.5,50,-.055,.055);
 TH2F *hyptarep = new TH2F("hyptarep","yptar vs endplane ",38,-.5,37.5,50,-.04,.04);
 //
 //
const int NumEndPl=7;
const int NumFocalPlane=6;
 Double_t xmin[NumEndPl]= {
-12.,-15.,0.,-0,-0,-0,-0
};
 Double_t xmax[NumEndPl]= {
+12.,+15.,+40.,40.,40.,40.,40.
};
 Double_t ymin[NumEndPl]= {
-12.,-15.,-15.,-15.,-15.,-15.,-15.
};
 Double_t ymax[NumEndPl]= {
+12.,+15.,+15.,+15.,+15.,+15.,+15.
};
 TString EndPlname[NumEndPl]={"target","10cm from target"," sieve front of septum"," mid septum"," 70 behind septum","q1 -307mm","q1 +470mm"};
 TH2F *hxyEndPl[NumEndPl],*hxyEndPlp[NumEndPl],*hxyEndPlp2[NumEndPl];
 TH1F *hxEndPl[NumEndPl],*hyEndPl[NumEndPl];
 TH1F *hxEndPlp[NumEndPl],*hyEndPlp[NumEndPl];
 for ( int i = 0; i < NumEndPl ; i++) {
   hxEndPl[i] = new TH1F(Form("hxEndPl%02d", i),Form("%s; X (cm) ; Counts",EndPlname[i].Data()),60,xmin[i],xmax[i]);
   hyEndPl[i] = new TH1F(Form("hyEndPl%02d", i),Form("%s; Y (cm) ; Counts",EndPlname[i].Data()),60,ymin[i],ymax[i]);
   hxEndPlp[i] = new TH1F(Form("hxEndPlp%02d", i),Form("%s; X (cm) ; Counts",EndPlname[i].Data()),60,xmin[i],xmax[i]);
   hyEndPlp[i] = new TH1F(Form("hyEndPlp%02d", i),Form("%s; Y (cm) ; Counts",EndPlname[i].Data()),60,ymin[i],ymax[i]);
   hxyEndPl[i] = new TH2F(Form("hxyEndPl%02d", i),EndPlname[i].Data(),80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
   hxyEndPlp[i] = new TH2F(Form("hxyEndPlp%02d", i),Form("%s; Y (cm) ; X (cm)  ",EndPlname[i].Data()),80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
   hxyEndPlp2[i] = new TH2F(Form("hxyEndPlp2%02d", i),Form("%s; Y (cm) ;  X (cm)",EndPlname[i].Data()),80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
 }
 //
          Double_t etemp;
 //
TTree *tsnake = (TTree*)fsnake->Get("ntuple");
 tsnake->SetBranchAddress("evnum",&evnum);
tsnake->SetBranchAddress("epnum",&epnum);
tsnake->SetBranchAddress("xabs",&xabs);
tsnake->SetBranchAddress("yabs",&yabs);
tsnake->SetBranchAddress("zabs",&zabs);
tsnake->SetBranchAddress("cxabs",&cxabs);
tsnake->SetBranchAddress("cyabs",&cyabs);
tsnake->SetBranchAddress("czabs",&czabs);
tsnake->SetBranchAddress("mom",&mom);
tsnake->SetBranchAddress("pathl",&pathl);
tsnake->SetBranchAddress("live",&live);
tsnake->SetBranchAddress("xrel",&xrel);
tsnake->SetBranchAddress("yrel",&yrel);
tsnake->SetBranchAddress("zrel",&zrel);
tsnake->SetBranchAddress("cxrel",&cxrel);
tsnake->SetBranchAddress("cyrel",&cyrel);
tsnake->SetBranchAddress("czrel",&czrel);
//

 Double_t xptar,yptar,xpfp,ypfp,xtar,xfp,yfp;
 int ntracks;
 Double_t xsave[NumEndPl],zsave[NumEndPl];
   // In Snake abs coordinate system +x is down,+y into SHMS,+z to large angle
//
 int EndPl;
 int last_plane= NumFocalPlane;
 int nfit=0;
      nentries = (int)tsnake->GetEntries();
           for (int ie = 0; ie < nentries; ie++)
            {
	   	      if (ie%10000==0)    	      cout << ie << " nfit = " << nfit << endl;
	         tsnake->GetEntry(ie);
   	  delta=(mom-cmom)/cmom;
          EndPl = epnum;
          if ( EndPl > NumEndPl) {
	    EndPl =NumEndPl;
            cout << " endplane number too large" ;
          }
	if ( epnum == 0  ) {
          ntracks++;
          yptar=czabs/cyabs-yptar_cent;
	  xptar=cxabs/cyabs;
          ytar=zabs/10.;
          hMom->Fill(mom);
          hDelta->Fill(delta);
          hytar->Fill(ytar);
	  hxptar->Fill(xptar);
	  hyptar->Fill(yptar);
	  if (last_plane != NumFocalPlane) {
          for (int ii = 0; ii < last_plane; ii++) {
          hxyEndPlp2[ii]->Fill(zsave[ii]/10.,xsave[ii]/10.);
          }
	  }
	  last_plane=0;
        }
	  xsave[EndPl]=xrel;
	  zsave[EndPl]=zrel;
	if ( live == 0  ) {
	  if (epnum <= NumFocalPlane) last_plane=epnum;
	  hxptarep->Fill(epnum,xptar);
	   hyptarep->Fill(epnum,yptar);
	  hxptarf->Fill(xptar);
	  hyptarf->Fill(yptar);
	  hMomf->Fill(mom);
          hDeltaf->Fill(delta);
          hytarf->Fill(ytar);
          hEndPlane->Fill(epnum);
          hxyEndPl[EndPl]->Fill(zrel/10.,xrel/10.);
	}
	if ( live == 1  ) {
          hxEndPl[EndPl]->Fill(xrel/10.);
          hyEndPl[EndPl]->Fill(zrel/10.);
        }
	if ( live == 1 && epnum==NumFocalPlane   ) {
	  last_plane=NumFocalPlane;
          for (int ii = 0; ii < NumEndPl; ii++) {
          hxyEndPlp[ii]->Fill(zsave[ii]/10.,xsave[ii]/10.);
          hyEndPlp[ii]->Fill(zsave[ii]/10.);
          hxEndPlp[ii]->Fill(xsave[ii]/10.);
          }
	  hxptarp->Fill(xptar);
	  hyptarp->Fill(yptar);
          hMomp->Fill(mom);
          hDeltap->Fill(delta);
          hytarp->Fill(ytar);
	  hxyfp->Fill(zrel/10.,xrel/10.);
	  yfp=zrel/10.; // focal plane in centimeters
	  xfp=xrel/10.; // focal plane in centimeters
          ztemp=-ytar*(cos(thcent)/tan(thcent+yptar)+sin(thcent));
	  xtar=-xptar*ztemp*cos(thcent)/100.;
          ypfp=czrel/cyrel-ypfp_cent;
	  xpfp=cxrel/cyrel;
	  hxpypfp->Fill(ypfp,xpfp);
	  hxpxfp->Fill(xpfp,xrel/10.);
	  hxpyfp->Fill(zrel/10.,xpfp);
	  hypyfp->Fill(ypfp,zrel/10.);
	  hypfpdel->Fill(ypfp,delta);
	  hypfpyptar->Fill(ypfp,yptar);
	  hypfpxptar->Fill(ypfp,xptar);
	  hyfpyptar->Fill(yptar,zrel/10.);
	  hxpdel->Fill(delta,xptar);
	  hxpyp->Fill(yptar,xptar);
	  hypdel->Fill(delta,yptar);
	  hypytar->Fill(ytar,yptar);
	//        
	}
// end loop over tree
	    }
  //
//
    TCanvas *cep = new TCanvas("cep"," Endplanes versus xp  and yp",800,800);
    cep->Divide(1,2);
    cep->cd(1);
    hxptarep->Draw("colz");
    cep->cd(2);
    hyptarep->Draw("colz");
    //
  //
    TCanvas *c7 = new TCanvas("c7","X vs Y (cm) for endplane 1-4",800,800);
    c7->Divide(2,2);
    for (int i = 1; i <5 ; i++){
      EndPl=i;
    c7->cd(i);
    hxyEndPlp[EndPl]->Draw("box");
    //    hxyEndPlp2[EndPl]->Draw("box,same");
    if (i !=2) hxyEndPlp2[EndPl]->SetLineColor(2);
    if (i !=2) hxyEndPl[EndPl]->Draw("colz,same");
    }
    TCanvas *c2 = new TCanvas("c2","Compare target",800,800);
    c2->Divide(2,2);
    c2->cd(1);
    hxptar->Draw();   
    hxptarp->SetLineColor(2);
    hxptarp->Draw("same");   
    hxptarf->SetLineColor(4);
    hxptarf->Draw("same");   
    c2->cd(2);
    hyptar->Draw();   
    hyptarp->SetLineColor(2);
    hyptarp->Draw("same");   
    hyptarf->SetLineColor(4);
    hyptarf->Draw("same");   
    c2->cd(3);
    hDelta->Draw();   
    hDeltap->SetLineColor(2);
     hDeltap->Draw("same");   
    hDeltaf->SetLineColor(4);
    hDeltaf->Draw("same");   
    c2->cd(4);
   hytar->Draw();   
    hytarp->SetLineColor(2);
     hytarp->Draw("same");   
    hytarf->SetLineColor(4);
    hytarf->Draw("same");   
    //
   TCanvas *c33 = new TCanvas("c33","FP v target",800,800);
   c33->Divide(2,2);
   c33->cd(1);
   hypfpyptar->Draw("colz");
   c33->cd(2);
   hyfpyptar->Draw("colz");
   c33->cd(3);
   hypfpdel->Draw("colz");
   c33->cd(4);
   hypfpxptar->Draw("colz");
    //
   TCanvas *c3 = new TCanvas("c3","Focal Plane",800,800);
    TPad *c3_1 = new TPad("c3_1", "c3_1",0.01,0.51,0.49,0.99);    TPad *c3_3 = new TPad("c3_3", "c3_3",0.51,0.51,0.99,0.99);
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
  hDelta->Write();
  hDeltap->Write();
  hytarp->Write();
  hyptarp->Write();
  hxptarp->Write();
 for ( int i = 0; i < NumEndPl ; i++) {
  hxEndPl[i]->Write();
  hyEndPl[i]->Write();
  hxyEndPl[i]->Write();
  hxyEndPlp[i]->Write();
  hxEndPlp[i]->Write();
  hyEndPlp[i]->Write();
  cout <<  " EndPl num = " << i << " " << EndPlname[i] << "               " << hxyEndPl[i]->Integral() << endl;
 }
   }
