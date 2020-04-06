#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <iostream>
using namespace std;
void make_dipole_field_tree() {
   TFile *f = new TFile("rootfiles/apex_snake.root","RECREATE");
   TTree *T = new TTree("ntuple","apex field data");
   Long64_t nlines = T->ReadFile("show_field/apex_withpos.dat","x:y:z:bx:by:bz");
   printf(" found %lld points\n",nlines);
   T->Write();
   //   T->Close();
}
