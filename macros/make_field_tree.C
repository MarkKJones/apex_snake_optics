#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TCint.h>
#include <TNtuple.h>
void make_field_tree() {
   TFile *f = new TFile("rootfiles/hb_midplane_withpos.root","RECREATE");
   TTree *T = new TTree("ntuple","hb field data");
   Long64_t nlines = T->ReadFile("show_field/hb_midplane_withpos.dat","x:y:z:bx:by:bz");
      printf(" found %lld points\n",nlines);
   T->Write();
   //   T->Close();
}
