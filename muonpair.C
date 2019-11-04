/*
  Simple macro showing how to access branches from the delphes output root file,
  loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
  mass.

  root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void muonpair(const char *inputFile) {
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  //  TH1F *m1 = new TH1F("m1","m1",200,60.0,140.0);
  //  TH1F *m2 = new TH1F("m2","m2",200,60.0,140.0);

  Float_t mass;

  // Loop over all events
  for (Int_t entry = 0; entry < numberOfEntries; ++entry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    Muon *muon1temp, *muon2temp;

    std::vector< std::vector<const Muon*> > Pairs;  // to store all muon pairs
    std::vector< std::vector<const Muon*> > PairsOpCharge;  // to store pairs with opposite charge
    
    // Select events with at least 4 muons
    if(branchMuon->GetEntries() > 3) {
      
      for(Int_t one = 0; one < branchMuon->GetEntries(); ++one) { // loop muons
	
	muon1temp = (Muon*) branchMuon->At(one);
	
	for (Int_t two = one; two < branchMuon->GetEntries(); ++two) { // loop muons


	  if(one != two){
	  
	    muon2temp = (Muon*) branchMuon->At(two);
	  
	    std::vector<const Muon*> pairOfMuons;

	    pairOfMuons.push_back(muon1temp);
	    pairOfMuons.push_back(muon2temp);
	    
	    Pairs.push_back(pairOfMuons);
	  
	  }
	} // close loop for two 
      } // close loop for one 
      

      cout<<" Unique number of pairs  "<<Pairs.size()<<endl;
      

      for(int j = 0; j< Pairs.size(); ++j){  // to check pairs of muons with opposite charge

	if(Pairs[j][0]->Charge != Pairs[j][1]->Charge){

	  PairsOpCharge.push_back(Pairs[j]);
	
	}
      }

      cout<<" Pairs with opposite charge  "<<PairsOpCharge.size()<<endl;
      
    }  // close conditional >= 4 muons
    Pairs.clear();
    PairsOpCharge.clear();
    
  } // close loop for events
} // close  funtion

