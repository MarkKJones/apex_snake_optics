#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
void make_track_tree(Int_t ntracks=3,TString fname="track") {
   TFile *f = new TFile(Form("rootfiles/%s.root",fname.Data()),"RECREATE");
   TTree *T = new TTree("ntuple","APEX");
   float_t nreg,pmom,x,y,z,bx,by,bz;
   TTree *TTrack = new TTree("track","APEX");
   Int_t ntot=0;
   float_t trnum,tr_mom,tr_x,tr_y,tr_z,tr_bx,tr_by,tr_bz;
   TTrack->Branch("trnum",&trnum,"trnum/F");
   TTrack->Branch("tr_mom",&tr_mom,"tr_mom/F");
   TTrack->Branch("tr_x",&tr_x,"tr_x/F");
   TTrack->Branch("tr_y",&tr_y,"tr_y/F");
   TTrack->Branch("tr_z",&tr_z,"tr_z/F");
   TTrack->Branch("tr_bx",&tr_bx,"tr_bx/F");
   TTrack->Branch("tr_by",&tr_by,"tr_by/F");
   TTrack->Branch("tr_bz",&tr_bz,"tr_bz/F");
  
   for (Int_t nt=1;nt<ntracks+1;nt++) {
     Long64_t nlines = T->ReadFile(Form("%s_%d.dat",fname.Data(),nt),"nreg:pmom:x:y:z:bx:by:bz");
   printf(" found %lld points\n",nlines);
   T->SetBranchAddress("nreg",&nreg);
   T->SetBranchAddress("pmom",&pmom);
   T->SetBranchAddress("x",&x);
   T->SetBranchAddress("y",&y);
   T->SetBranchAddress("z",&z);
   T->SetBranchAddress("bx",&bx);
   T->SetBranchAddress("by",&by);
   T->SetBranchAddress("bz",&bz);
   for (Int_t nev=ntot;nev<nlines+ntot;nev++) {
     T->GetEntry(nev);
     trnum = nt;
     tr_mom=pmom;
     tr_x = x;
     tr_y = y;
     tr_z = z;
     tr_bx = bx;
     tr_by = by;
     tr_bz = bz;
     TTrack->Fill();
   }
   ntot=ntot+nlines;
   }
   TTrack->Write();
}
