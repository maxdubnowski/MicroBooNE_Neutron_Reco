#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>


void test(){
   TH1D::SetDefaultSumw2();
   TH2D::SetDefaultSumw2();

   // Output file
   TString FileName = "Test.root";
   TFile* OutputFile = new TFile(FileName,"recreate");
   std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

   TH1D* Plot;
   TH1D* Plot1;

   Plot = new TH1D("Title", ";X;Y", 10, 1, 10);
   
   Plot->Fill(7);
   Plot->Fill(5);
   Plot->Fill(7);
   Plot->Fill(5);
   Plot->Fill(7);
   Plot->Fill(5);

   OutputFile->cd();
   OutputFile->Write();
   OutputFile->Close();
}
