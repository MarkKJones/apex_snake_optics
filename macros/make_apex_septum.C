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
void make_apex_septum(TString fstart) {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 Double_t cmom=1.08,delta,ytar,ztemp,thcent=40./180.*3.14159;
 Double_t yptar_cent=0.08749; // tan(5) deg
 Double_t ypfp_cent=0.22169; //tan(12.5) deg
   Double_t solang_fac=0.10*0.20*1000.;
TFile *fsnake =  new TFile("rootfiles/"+fstart+".root");
TFile *hsnake = new TFile("rootfiles/"+fstart+"_standard_hist.root","recreate");
 Bool_t fit_matrix=kFALSE;
   string newcoeffsfilename="hms_newfit_save.dat";
 TRandom3 err; 
//
  string oldcoeffsfilename="hrs_recon_cosy.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
  int num_recon_terms_old;
  int num_recon_terms_fit=0;

  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_fit;
  vector<double> yptarcoeffs_fit;
  vector<double> ytarcoeffs_fit;
  vector<double> deltacoeffs_fit;
  vector<int> xfpexpon_fit;
  vector<int> xpfpexpon_fit;
  vector<int> yfpexpon_fit;
  vector<int> ypfpexpon_fit;
  vector<int> xtarexpon_fit;
  vector<int> expon_fit_flag;

  vector<double> xptarcoeffs_new;
  vector<double> yptarcoeffs_new;
  vector<double> ytarcoeffs_new;
  vector<double> deltacoeffs_new;
  vector<int> xfpexpon_new;
  vector<int> xpfpexpon_new;
  vector<int> yfpexpon_new;
  vector<int> ypfpexpon_new;
  vector<int> xtarexpon_new;
  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  TString currentline;

  num_recon_terms_old = 0;
  cout << " cur = " << currentline.Data() << endl;

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
    if( 1==1) {
    //        if ( ((expontemp[4]+expontemp[0]+expontemp[1]+expontemp[2]+expontemp[3])==1) || ((expontemp[0]+expontemp[1]+expontemp[2]+expontemp[3]+expontemp[4])==2) ) {
	  //   if ( ((expontemp[4]+expontemp[0]+expontemp[1]+expontemp[2]+expontemp[3])==1) ) {
    expon_fit_flag.push_back(1);
      num_recon_terms_fit++;     
    xfpexpon_fit.push_back(expontemp[0]);
    xpfpexpon_fit.push_back(expontemp[1]);
    yfpexpon_fit.push_back(expontemp[2]);
    ypfpexpon_fit.push_back(expontemp[3]);
    xtarexpon_fit.push_back(expontemp[4]);
    cout << expontemp[0] << " " << expontemp[1]<< " " <<expontemp[2]<< " " << expontemp[3] << " " << expontemp[4] << endl;
    } else {
    expon_fit_flag.push_back(0);
    }
    // cout << "(C_theta, C_y, C_phi, C_delta) = (" 
    //	 << xptarcoeffs_old[num_recon_terms_old] << ", "
    //	 << ytarcoeffs_old[num_recon_terms_old] << ", " 
    //	 << yptarcoeffs_old[num_recon_terms_old] << ", "
    //	 << deltacoeffs_old[num_recon_terms_old] << "), expon = "
    //	 << xfpexpon_old[num_recon_terms_old] << xpfpexpon_old[num_recon_terms_old] 
    //	 << yfpexpon_old[num_recon_terms_old] << ypfpexpon_old[num_recon_terms_old] 
    //	 << xtarexpon_old[num_recon_terms_old] << endl;
    
    num_recon_terms_old++;
  }

  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  cout << "num recon terms in OLD matrix to fit = " << num_recon_terms_fit << endl;
  int npar,nfit_max=3000;
  int npar_old;
  npar=  num_recon_terms_fit;
  npar_old=  num_recon_terms_old;
  TVectorD b_ytar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  TVectorD b_delta(npar);
  TMatrixD lambdaold(npar_old,nfit_max);
  TMatrixD lambda(npar,nfit_max);

  //TMatrixD A_yptar(npar,npar);
  TMatrixD Ay(npar,npar);
 
//
 int nentries; //number of entries in file
Double_t evnum,epnum,xabs,yabs,zabs,cxabs,cyabs,czabs,mom,pathl,live,xrel,yrel,zrel,cxrel,cyrel,czrel;
 //
 TH1F *hEndPlane = new TH1F("hEndPlane"," Endplane Number",4,-.5,5.5);
 TH1F *hMom = new TH1F("hMom","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomp = new TH1F("hMomp","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomf = new TH1F("hMomf","Momentum (GeV)",50,.5,1.5);
 TH1F *hDeltarecon = new TH1F("hDeltarecon","Delta Recon ",40,-.12,.12);
 TH1F *hDeltadiff = new TH1F("hDeltadiff","Delta Diff % ",40,-.3,.3);
 TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon ",40,-.12,.12);
 TH1F *hDeltanewdiff = new TH1F("hDeltanewdiff","Delta New Diff % ",40,-.3,.3);
 TH2F *hDeltadiffvsxpfp = new TH2F("hDeltadiffvxptar",";Delta Diff %;xpfp ",40,-.1,.1,40,-.1,.1);
 TH1F *hDelta = new TH1F("hDelta","Delta ",20,-.1,.1);
 TH1F *hDeltap = new TH1F("hDeltap","Delta ",20,-.1,.1);
 TH1F *hDeltaf = new TH1F("hDeltaf","Delta ",20,-.1,.1);
 TH1F *hDeltar = new TH1F("hDeltar",";Delta ;Solid angle (msr)  ",20,-.1,.1);
 TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
 TH1F *hytarrecon = new TH1F("hytarrecon","ytar recon(cm)",70,-7.,7.);
 TH1F *hztarrecon = new TH1F("hztarrecon","ztar recon(cm)",200,-10.,10.);
 TH1F *hytarnew = new TH1F("hytarnew","ytar new(cm)",70,-7.,7.);
 TH1F *hytardiff = new TH1F("hytardiff","ytar diff(cm)",70,-10.,10.);
 TH1F *hytarnewdiff = new TH1F("hytarnewdiff","ytar new diff(cm)",70,-10.,10.);
 TH1F *hytarp = new TH1F("hytarp","ytar (cm)",70,-7.,7.);
 TH1F *hytarf = new TH1F("hytarf","ytar (cm)",70,-7.,7.);
 TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
 TH1F *hxptarrecon = new TH1F("hxptarrecon","xptar recon",100,-.12,.12);
 TH1F *hxptardiff = new TH1F("hxptardiff","xptar diff (mr)",100,-100,100);
 TH1F *hxptarnew = new TH1F("hxptarnew","xptar new recon",100,-.12,.12);
 TH1F *hxptarnewdiff = new TH1F("hxptarnewdiff","xptar new diff (mr)",100,-100,100);
 TH1F *hxptarp = new TH1F("hxptarp","xptar ",100,-.12,.12);
 TH1F *hxptarf = new TH1F("hxptarf","xptar ",100,-.12,.12);
 TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
 TH1F *hyptarrecon = new TH1F("hyptarrecon","yptar recon ",100,-.1,.1);
 TH1F *hyptardiff = new TH1F("hyptardiff","yptar diff (mr) ",100,-100,100);
 TH1F *hyptarnew = new TH1F("hyptarnew","yptar new recon ",100,-.1,.1);
 TH1F *hyptarnewdiff = new TH1F("hyptarnewdiff","yptar new diff (mr) ",100,-100,100);
 TH1F *hyptarp = new TH1F("hyptarp","yptar ",100,-.1,.1);
 TH1F *hyptarf = new TH1F("hyptarf","yptar ",100,-.1,.1);
 TH2F *hxyfp = new TH2F("hxyfp","xfp vs yfp (cm)",100,0.,50.,40,-20.,20.);
 TH2F *hxpypfp = new TH2F("hxpypfp","xpfp vs ypfp (rad)",100,-.1,.1,100,-.05,.05);
 TH2F *hxpyfp = new TH2F("hxpyfp","xpfp vs yfp ",100,0.,40.,100,-.1,.1);
 TH2F *hxpdel = new TH2F("hxpdel","xptar vs delta ",50,-.1,.1,55,-.1,.1);
 TH2F *hxpxfp = new TH2F("hxpxfp","xfp vs xpfp ",90,-.15,.15,100,-20.,20.);
 TH2F *hypyfp = new TH2F("hypyfp","ypfp;yfp ",100,-.1,.1,100,0.,50.);
 TH2F *hypfpyptar = new TH2F("hypfpyfptar",";ypfp;yptar ",100,-.1,.1,100,-.1,.1);
 TH2F *hypfpxptar = new TH2F("hypfpxfptar",";ypfp;xptar ",100,-.1,.1,100,-.1,.1);
 TH2F *hyfpyptar = new TH2F("hyfpyfptar",";yptar;yfp ",100,-.1,.1,100,0.,50.);
 TH2F *hypdel = new TH2F("hypdel","yptar vs delta ",100,-.1,.1,400,-.05,.05);
 TH2F *hypfpdel = new TH2F("hypfpdel",";ypfp;delta ",100,-.04,.04,400,-.05,.05);
 TH2F *hypytar = new TH2F("hypytar","yptar vs ytar;ytar (cm); yptar (rad) ",50,-6.,6.,40,-.04,.04);
 TH2F *hypytarrecon = new TH2F("hypytarrecon","Recon yptar vs ytar; ytar recon (cm); yptar recon (rad) ",50,-6.,6.,40,-.04,.04);
 TH2F *hyfpvsydiff = new TH2F("hyfpvsydiff","y_fp vs ydiff; y_fp (cm); ytar diff (cm) ",50,-10.,10.,40,-2,2.);
 TH2F *hxpyp = new TH2F("hxpyp","xptar vs yptar;yptar;xptar ",40,-.04,.04,55,-.1,.1);
 TH2F *hxpyprecon = new TH2F("hxpyprecon","xptar vs yptar recon; yptar recon;xptar recon ",40,-.04,.04,55,-.1,.1);
 TH2F *hxtarvdeldiff = new TH2F("hxtarvdeldiff ","xtar vs delta diff ",30,-2.,2.,50,-.5,.5);
 TH2F *hypvdeldiff = new TH2F("hypvdeldiff ","yptar vs delta diff ",30,-.06,.06,30,-.6,.6);
 TH2F *hxpvdeldiff = new TH2F("hxpvdeldiff ","xptar vs delta diff; xptar; delta-delta_recon ",30,-.06,.06,30,-.6,.6);
 TH2F *hxpreconvdeldiff = new TH2F("hxpreconvdeldiff ","xptar recon vs delta diff; xptar recon (rad); delta - delta_recon ",30,-.06,.06,30,-.6,.6);
 TH2F *hypreconvxprecon = new TH2F("hypreconvxprecon","ytar = 0;xptar recon;yptar recon; ",30,-.05,.05,30,-.06,.06);
 TH2F *hxptarep = new TH2F("hxptarep","xptar vs endplane ",38,-.5,37.5,50,-.055,.055);
 TH2F *hyptarep = new TH2F("hyptarep","yptar vs endplane ",38,-.5,37.5,50,-.04,.04);
 //
 //
const int NumEndPl=5;
const int NumFocalPlane=4;
 Double_t xmin[NumEndPl]= {
-12.,-15.,0.,-0,-0
};
 Double_t xmax[NumEndPl]= {
+12.,+15.,+40.,40.,40.
};
 Double_t ymin[NumEndPl]= {
-12.,-15.,-15.,-15.,-15.
};
 Double_t ymax[NumEndPl]= {
+12.,+15.,+15.,+15.,+15.
};
 TString EndPlname[NumEndPl]={"target","10cm from target"," -75cm in front of septum"," mid septum"," 70 behind septum"};
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
	  /*
	  xpfp=xpfp+err.Gaus(0.0,0.001);
	  ypfp=ypfp+err.Gaus(0.0,0.001);         
	  xfp=xfp+err.Gaus(0.0,.03);
	  yfp=yfp+err.Gaus(0.0,.03);         
	  */
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
	  //	  cout << yptar << " " << cmom << " " << delta << " " << ypfp << endl;
          Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
          Double_t ytar_fixed = 0.0,yptar_fixed=0.0,xptar_fixed=0.0,delta_fixed=0.0;
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar, xtarexpon_old[icoeffold] );
        	deltatemp += deltacoeffs_old[icoeffold] * etemp;
        	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	        yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	         xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
		 if (expon_fit_flag[icoeffold]==0) {
        	delta_fixed+= deltacoeffs_old[icoeffold] * etemp;
        	ytar_fixed += ytarcoeffs_old[icoeffold] * etemp;
	        yptar_fixed += yptarcoeffs_old[icoeffold] * etemp;
	         xptar_fixed += xptarcoeffs_old[icoeffold] *etemp; 
		 }
	      if (nfit < nfit_max) {
              lambdaold[icoeffold][nfit] = etemp;
              }
	  } // for loop
          for( int icoefffit=0; icoefffit<num_recon_terms_fit; icoefffit++ ){
	    //	    cout << xfpexpon_fit[icoefffit] << " " << yfpexpon_fit[icoefffit] << " " << xpfpexpon_fit[icoefffit] << " " << ypfpexpon_fit[icoefffit]  << " " << xtarexpon_fit[icoefffit] << endl;
         	etemp= 
	  pow( xfp / 100.0, xfpexpon_fit[icoefffit] ) * 
	  pow( yfp / 100.0, yfpexpon_fit[icoefffit] ) * 
	  pow( xpfp, xpfpexpon_fit[icoefffit] ) * 
	  pow( ypfp, ypfpexpon_fit[icoefffit] ) * 
	  pow( xtar, xtarexpon_fit[icoefffit] );
	      if (nfit < nfit_max) {
              lambda[icoefffit][nfit] = etemp;
	      b_xptar[icoefffit] += (xptar-xptar_fixed) * etemp;
	      b_yptar[icoefffit] += (yptar-yptar_fixed) * etemp;
	      b_ytar[icoefffit] += ((ytar) /100.0-ytar_fixed) * etemp;
	      b_delta[icoefffit] += (delta-delta_fixed) * etemp;
	      }
	  }
	  if (nfit < nfit_max) {
             nfit++;
          }
	  if (fit_matrix && nfit == nfit_max) ie=nentries;
  	    xfptrue.push_back( xfp );
	    yfptrue.push_back( yfp );
	    xpfptrue.push_back( xpfp );
	    ypfptrue.push_back( ypfp );
	    xtartrue.push_back( xtar );
	    xptartrue.push_back( xptar );
	    ytartrue.push_back( ytar  );
	    yptartrue.push_back( yptar  );
	    deltatrue.push_back( delta  );
          ztemp=ytartemp*100.*(cos(thcent)/tan(thcent+yptartemp)+sin(thcent));
	  hztarrecon->Fill(ztemp);
	  hytarrecon->Fill(ytartemp*100.);
	  hyptarrecon->Fill(yptartemp);
	  hxptarrecon->Fill(xptartemp);
	  hDeltarecon->Fill(deltatemp);
	  hytardiff->Fill(ytartemp*100.-ytar);
	  hyptardiff->Fill(1000*(yptartemp-yptar));
	  hxptardiff->Fill(1000*(xptartemp-xptar));
	  hDeltadiff->Fill(100*(deltatemp-delta));
          hyfpvsydiff->Fill(zrel/10.,ytartemp*100.-ytar);
	  hDeltadiffvsxpfp->Fill(100*(deltatemp-delta),xpfp);         
	  hxtarvdeldiff->Fill(100*xtar,100*(deltatemp-delta));         
	  hypvdeldiff->Fill(yptar,100*(deltatemp-delta));         
	  hxpvdeldiff->Fill(xptar,100*(deltatemp-delta));         
	  hxpreconvdeldiff->Fill(xptartemp,100*(deltatemp-delta));  
	  hxpyprecon->Fill(yptartemp,xptartemp);
	  hypytarrecon->Fill(ytartemp*100,yptartemp);
	}
	//        
	    }
// end loop over tree
  //
if ( fit_matrix) {
 for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
 }
  for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 100 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
      }
    }
  }
  TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_ytar );
  cout << "ytar solution ok = " << ok << endl;
  //b_ytar.Print();
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  //b_yptar.Print();
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  //b_xptar.Print();
  ok = Ay_svd.Solve( b_delta );
  cout << "Delta solution ok = " << ok << endl;
  //b_delta.Print();
  // calculate target quantities with new fit parameter
  for( int ifit=0; ifit<nfit; ifit++){
          Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
     for( int ipar=0; ipar<npar; ipar++){
       etemp=lambda[ipar][ifit];
        	deltanew += b_delta[ipar] * etemp;
        	ytarnew += b_ytar[ipar] * etemp;
	        yptarnew += b_yptar[ipar] * etemp;
	         xptarnew += b_xptar[ipar] *etemp;        
    }
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
	    etemp=lambdaold[icoeffold][ifit];
		 if (expon_fit_flag[icoeffold]==0) {
        	deltanew+= deltacoeffs_old[icoeffold] * etemp;
        	ytarnew += ytarcoeffs_old[icoeffold] * etemp;
	        yptarnew += yptarcoeffs_old[icoeffold] * etemp;
	         xptarnew += xptarcoeffs_old[icoeffold] *etemp; 
		 }
	  }
	  hytarnew->Fill(ytarnew*100.);
	  hyptarnew->Fill(yptarnew);
	  hxptarnew->Fill(xptarnew);
	  hDeltanew->Fill(deltanew);
	  hytarnewdiff->Fill(ytarnew*100.-ytartrue.at(ifit));
	  hyptarnewdiff->Fill(1000*(yptarnew-yptartrue.at(ifit)));
	  hxptarnewdiff->Fill(1000*(xptarnew-xptartrue.at(ifit)));
	  hDeltanewdiff->Fill(100*(deltanew-deltatrue.at(ifit)));    
  }
  // write out coeff
  char coeffstring[100];
  int nfitted=0;
  cout << "writing new coeffs file" << endl;
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
      if (expon_fit_flag[icoeffold]==1) {
      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", b_xptar[nfitted] );
      newcoeffsfile << coeffstring; 
      //      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", b_ytar[nfitted] );
      newcoeffsfile << coeffstring;
      sprintf( coeffstring, "%16.9g", b_yptar[nfitted] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      sprintf( coeffstring, "%16.9g", b_delta[nfitted] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      newcoeffsfile << " ";
 	newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_fit[nfitted]; 
	newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_fit[nfitted]; 
	newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_fit[nfitted]; 
	newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_fit[nfitted]; 
	newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_fit[nfitted]; 
      newcoeffsfile << endl;
      nfitted++;
     } else{
      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", xptarcoeffs_old[icoeffold] );
      newcoeffsfile << coeffstring; 
      //      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", ytarcoeffs_old[icoeffold] );
      newcoeffsfile << coeffstring;
      sprintf( coeffstring, "%16.9g", yptarcoeffs_old[icoeffold] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      sprintf( coeffstring, "%16.9g", deltacoeffs_old[icoeffold] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      newcoeffsfile << " ";
      
	newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_old[icoeffold]; 
      newcoeffsfile << endl;
      }
	  }
  newcoeffsfile << " ---------------------------------------------" << endl;

  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
  //
	   } // if fitting
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
    if (fit_matrix ) {
    TCanvas *c2c = new TCanvas("c2c","Xtar versus Delta diff",800,800);
    c2c->Divide(1,1);
    c2c->cd(1);
    hxtarvdeldiff->Draw("colz");
    TCanvas *c2b = new TCanvas("c2b","Ytar versus Ytardiff",800,800);
    c2b->Divide(1,2);
    c2b->cd(1);
    hyfpvsydiff->Draw("colz");
    c2b->cd(2);
    hztarrecon->Draw();
//
//
    TCanvas *c4 = new TCanvas("c4","2d target plots",800,800);
    c4->Divide(2,2);
    c4->cd(1);
    hxpdel->Draw("colz");
    c4->cd(2);
    hypdel->Draw("colz");
    c4->cd(3);
    hxpyp->Draw("colz");
    c4->cd(4);
    hypytar->Draw("colz");
//
    TCanvas *crecon = new TCanvas("crecon","Recon target",800,800);
    crecon->Divide(2,2);
    crecon->cd(1);
    hytarrecon->Draw();
    hytarnew->Draw("same");
    hytarnew->SetLineColor(2);
    crecon->cd(2);
    hyptarrecon->Draw();
    hyptarnew->Draw("same");
    hyptarnew->SetLineColor(2);
    crecon->cd(3);
    hxptarrecon->Draw();
    hxptarnew->Draw("same");
    hxptarnew->SetLineColor(2);
    crecon->cd(4);
    hDeltarecon->Draw();
    hDeltanew->Draw("same");
    hDeltanew->SetLineColor(2);
    //
//
    TCanvas *cdiff2d = new TCanvas("cdiff2d","Diff delta",800,800);
    cdiff2d->Divide(2,3);
    cdiff2d->cd(1);
    hypytar->Draw("colz");
    cdiff2d->cd(2);
    hypytarrecon->Draw("colz");
    cdiff2d->cd(3);
    hxpvdeldiff->Draw("colz");
    cdiff2d->cd(4);
    hxpreconvdeldiff->Draw("colz");
    cdiff2d->cd(5);
    hxpyp->Draw("colz");
        cdiff2d->cd(6);
    hxpyprecon->Draw("colz");
    //cdiff2d->cd(7);
    // hDeltadiffvsxpfp->Draw("colz");
//
    TCanvas *cdiff = new TCanvas("cdiff","Diff target",800,800);
    cdiff->Divide(2,2);
    cdiff->cd(1);
    hytardiff->Draw();
 hytardiff->Fit("gaus");
 TF1 *fitcydiff=hytardiff->GetFunction("gaus");
    cdiff->cd(2);
    hyptardiff->Draw();
 hyptardiff->Fit("gaus");
 TF1 *fitcypdiff=hyptardiff->GetFunction("gaus");
    cdiff->cd(3);
    hxptardiff->Draw();
 hxptardiff->Fit("gaus");
 TF1 *fitcxpdiff=hxptardiff->GetFunction("gaus");
    cdiff->cd(4);
    //    hDeltadiffvsxfp->Draw("colz");
    hDeltadiff->Draw();
 hDeltadiff->Fit("gaus");
 TF1 *fitcdeldiff=hDeltadiff->GetFunction("gaus");
 cout << " sig_y  "<< fitcydiff->GetParameter(2) << " sig_yp " << fitcypdiff->GetParameter(2) << " sig_xp "<< fitcxpdiff->GetParameter(2) << " sig_delta " << fitcdeldiff->GetParameter(2) << endl;
    //
//
    TCanvas *cdiffnew = new TCanvas("cdiffnew","New Diff target",800,800);
    cdiffnew->Divide(2,2);
    cdiffnew->cd(1);
    hytarnewdiff->Draw();
 hytarnewdiff->Fit("gaus");
 TF1 *fitcynewdiff=hytarnewdiff->GetFunction("gaus");
    cdiffnew->cd(2);
    hyptarnewdiff->Draw();
 hyptarnewdiff->Fit("gaus");
 TF1 *fitcypnewdiff=hyptarnewdiff->GetFunction("gaus");
    cdiffnew->cd(3);
    hxptarnewdiff->Draw();
 hxptarnewdiff->Fit("gaus");
 TF1 *fitcxpnewdiff=hxptarnewdiff->GetFunction("gaus");
    cdiffnew->cd(4);
    //    hDeltadiffvsxfp->Draw("colz");
    hDeltanewdiff->Draw();
 hDeltanewdiff->Fit("gaus");
 TF1 *fitcdelnewdiff=hDeltanewdiff->GetFunction("gaus");
 cout << "new  sig_y  "<< fitcynewdiff->GetParameter(2) << " sig_yp " << fitcypnewdiff->GetParameter(2) << " sig_xp "<< fitcxpnewdiff->GetParameter(2) << " sig_delta " << fitcdelnewdiff->GetParameter(2) << endl;
    }

    //
  hyfpvsydiff->Write();
  hztarrecon->Write();
  hDelta->Write();
  hDeltap->Write();
  hytarp->Write();
  hyptarp->Write();
  hxptarp->Write();
  hDeltadiff->Write();
  hytardiff->Write();
  hyptardiff->Write();
  hxptardiff->Write();
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
