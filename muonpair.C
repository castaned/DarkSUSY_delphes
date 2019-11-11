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

  TH1F *m1 = new TH1F("m1","m1",200,0.0,240.0);
  TH1F *m2 = new TH1F("m2","m2",200,0.0,240.0);


  TH1F *ptmuon1 = new TH1F("ptmuon1","ptmuon1",100,0.0,300.0);
  TH1F *ptmuon2 = new TH1F("ptmuon2","ptmuon2",100,0.0,300.0);
  TH1F *ptmuon3 = new TH1F("ptmuon3","ptmuon3",100,0.0,300.0);
  TH1F *ptmuon4 = new TH1F("ptmuon4","ptmuon4",100,0.0,300.0);

  
  Float_t mass;

  // Loop over all events
  for (Int_t entry = 0; entry < numberOfEntries; ++entry) {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    Muon *muon1temp, *muon2temp;
    Muon *muon1,*muon2,*muon3,*muon4;
    
    std::vector< std::vector<const Muon*> > Pairs;  // to store all muon pairs
    std::vector< std::vector<const Muon*> > PairsOpCharge;  // to store pairs with opposite charge
    std::vector< std::vector<const Muon*> > FinalPairs;  // to store pairs with opposite charge
    
    // Select events with at least 4 muons
    if(branchMuon->GetEntries() > 3) {

      muon1 = (Muon*) branchMuon->At(0);
      muon2 = (Muon*) branchMuon->At(1);
      muon3 = (Muon*) branchMuon->At(2);
      muon4 = (Muon*)branchMuon->At(3);
      
      ptmuon1->Fill(muon1->PT);
      ptmuon2->Fill(muon2->PT);
      ptmuon3->Fill(muon3->PT);
      ptmuon4->Fill(muon4->PT);

      
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
      

      //      cout<<" Unique number of pairs  "<<Pairs.size()<<endl;
      

      for(int j = 0; j< Pairs.size(); ++j){  // to check pairs of muons with opposite charge

	if(Pairs[j][0]->Charge != Pairs[j][1]->Charge){

	  PairsOpCharge.push_back(Pairs[j]);
	
	}
      }

      // m1->Fill( ( PairsOpCharge[0][0]->P4() + PairsOpCharge[0][1]->P4()).M()   );
      // m2->Fill( ( PairsOpCharge[1][0]->P4() + PairsOpCharge[1][1]->P4()).M()  );

      
      double deltaMassTmp = 13000.0;
      
      // delta mass b/t each of two pairs
      std::vector<float> AbsdMass;
      std::vector<int> Alice;
      std::vector<int> Bob;
      
      
      if(PairsOpCharge.size()>1){

	for(unsigned int Ajet=0; Ajet < PairsOpCharge.size(); Ajet++){
	  for(unsigned int Bjet=0; Bjet< PairsOpCharge.size(); Bjet++){

	    deltaMassTmp = fabs(  (PairsOpCharge[Ajet][0]->P4() + PairsOpCharge[Ajet][1]->P4()).M() - (PairsOpCharge[Bjet][0]->P4()+PairsOpCharge[Bjet][1]->P4()).M());
	    
	    AbsdMass.push_back(deltaMassTmp);

	    cout<<"  Alice  "<<Ajet<<" Bob  "<<Bjet<<endl;

	    Alice.push_back(Ajet);
	    Bob.push_back(Bjet);
	    
	  }
	}
      
	std::vector<const Muon*> PairOne;
	std::vector<const Muon*> PairTwo;
      
      
	if ( ( AbsdMass.size() == Alice.size() ) &&
	     ( AbsdMass.size() == Bob.size() ) ) {
	
	
	  unsigned int Findcount = 0;
	  while( FinalPairs.size() != 2 && Findcount < AbsdMass.size() ){
	  
	    //find two pairs with min |dM|
	    double MinDeltaMass = 13000.;
	    unsigned int Index = -1;
	    Findcount++;

	  
	    for (unsigned int X = 0; X < AbsdMass.size(); X++) {
	      //std::cout << ">>>>>> AbsdMass.at("<<X<<") = "<< AbsdMass.at(X) << std::endl;
	      if( AbsdMass.at(X) < MinDeltaMass){
		MinDeltaMass = AbsdMass[X];

		PairOne = PairsOpCharge[Alice[X]];
		PairTwo = PairsOpCharge[Bob[X]];
		Index = X;
	      }//end if
	      //	      cout<<" Alice  "<<Alice[X]<<"  Bob   "<<Bob[X]<<endl;

	    }//end find pairs with min |dM|
	    

	    //	  std::cout << ">>> Found two pairs with min |dM|, [mass #"<<Alice.at(Index)<<", mass #"<<Bob.at(Index)<<"  "<<AbsdMass.at(Index)<<", " << Index <<"]"<< std::endl;
	  }



	  cout<<" m1   "<<  (PairOne[0]->P4() + PairOne[1]->P4()).M() <<endl;
	  cout<<" m2   "<<  (PairTwo[0]->P4() + PairTwo[1]->P4()).M() <<endl;

	  
	  FinalPairs.push_back(PairOne);
	  FinalPairs.push_back(PairTwo);

	}
      
      }
      ///      cout<<" Pairs with opposite charge  "<<PairsOpCharge.size()<<endl;
      //      cout<<" Pairs with min dM  "<<FinalPairs.size()<<endl;


      //      cout<<" m1   "<<(FinalPairs[0][0]->P4() + FinalPairs[0][1]->P4()).M()<<endl;
      //      cout<<" m2   "<<(FinalPairs[1][0]->P4() + FinalPairs[1][1]->P4()).M()<<endl;
      
      m1->Fill( ( FinalPairs[0][0]->P4() + FinalPairs[0][1]->P4()).M()   );
      m2->Fill( ( FinalPairs[1][0]->P4() + FinalPairs[1][1]->P4()).M()  );
       
    }  // close conditional >= 4 muons
    Pairs.clear();
    PairsOpCharge.clear();
    FinalPairs.clear();
    
  }


  TCanvas *c = new TCanvas("c","c");
  m1->Draw();

  TCanvas *c1 = new TCanvas("c1","c1");
  m2->Draw();


  TLegend *leg = new TLegend(0.5,0.7);
  leg->AddEntry(ptmuon1,"ptmuon1");
  leg->AddEntry(ptmuon2,"ptmuon2");
  leg->AddEntry(ptmuon3,"ptmuon3");
  leg->AddEntry(ptmuon4,"ptmuon4");
  
  TCanvas *c2 = new TCanvas("c2","c2");
  //  c2->SetLogy();
  ptmuon1->SetLineColor(1);
  ptmuon3->Draw();

  ptmuon2->SetLineColor(2);
  ptmuon2->Draw("same");
  ptmuon3->SetLineColor(3);
  ptmuon1->Draw("same"); 
  ptmuon4->SetLineColor(4);
  ptmuon4->Draw("same");
  leg->Draw("same");

  
}


