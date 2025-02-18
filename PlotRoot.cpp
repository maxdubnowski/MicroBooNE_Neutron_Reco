#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include "Constants.h"

using namespace std;
using namespace Constants;

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

void PlotRoot() {


  
        TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			
	double LegendTextSize = 0.03;

	TString OutFilePath = "/uboone/app/users/maxd/MicroBooNE_Neutron_Reco/";
	
	bool stackedHist = false;


	std::vector<TString> PlotNames;
	///*
	//Choose which plots to display in the files
	
	// PlotNames.push_back("TrackMultiplicity");
	// PlotNames.push_back("ShowerMultiplicity");
	// PlotNames.push_back("BlipPDG");
	// PlotNames.push_back("BlipMultiplicity");
	// PlotNames.push_back("XAxisBlipLocation");
	// PlotNames.push_back("YAxisBlipLocation");
	// PlotNames.push_back("ZAxisBlipLocation");
	// PlotNames.push_back("BlipVertexDist");
	// PlotNames.push_back("ProtonBlipDist");
	// PlotNames.push_back("BlipCharge");
	// PlotNames.push_back("BlipProxTrkDist");
	// PlotNames.push_back("BlipTime");
	// PlotNames.push_back("BlipMaxWireSpan");
	// PlotNames.push_back("BlipEnergy");

	// PlotNames.push_back("ProtonBlipCharge");
	// PlotNames.push_back("ProtonBlipProxTrkDist");
	// PlotNames.push_back("ProtonBlipTime");
	// PlotNames.push_back("ProtonBlipMaxWireSpan");
	// PlotNames.push_back("ProtonBlipEnergy");

	PlotNames.push_back("RecoMuonCosThetaPlot");
	PlotNames.push_back("RecoDeltaPtPlot");
	PlotNames.push_back("RecoDeltaAlphaTPlot");
	// PlotNames.push_back("RecoVertexXPlot");
	// PlotNames.push_back("RecoVertexYPlot");
	// PlotNames.push_back("RecoVertexZPlot");
	
	// PlotNames.push_back("RecoPMissMagPlot");
	// PlotNames.push_back("RecoPMissDirPlot");
	
	// PlotNames.push_back("BacktrackedPMissMagPlot");
	// PlotNames.push_back("BacktrackedPMissDirPlot");

	// PlotNames.push_back("TruePMissMagPlot");
	// PlotNames.push_back("TruePMissDirPlot");
		
	// PlotNames.push_back("BacktrackedDeltaPtPlot");
	// PlotNames.push_back("BacktrackedDeltaAlphaTPlot");
	

	//PlotNames.push_back("TrueMuonEnergyPlot");
	//PlotNames.push_back("TrueNeutronMultiplicityPlot");

	const int NPlots = PlotNames.size();
	






	if(!stackedHist){

	  // -------------------- Method 1: Overlayed Line Histograms ---------------------- //


	  std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors;

	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("All");
	  Colors.push_back(kBlack);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("CC1p0#pi");
	  Colors.push_back(kCyan+2);	
	  
	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("nonCC1p0#pi");
	  Colors.push_back(kRed+2);		

	  const int NSamples = Names.size();
	  std::vector<TFile*> Files; Files.resize(NSamples);


	  // Loop over the samples to open the files and the TTree
	  for (int iSample = 0; iSample < NSamples; iSample++) {
	    Files[iSample] = new TFile(Names[iSample],"readonly");
	  } // End of the loop over the samples




	  vector<int> LineStyle;
	  LineStyle = {1, 1, 2};



	  // Loop over the plots to be compared
	  for (int iPlot = 0; iPlot < NPlots; iPlot++) {
	    
	    TString CanvasName = "Canvas_" + PlotNames[iPlot];
	    TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	    PlotCanvas->cd();
	    PlotCanvas->SetTopMargin(0.12);
	    PlotCanvas->SetLeftMargin(0.2);
	    PlotCanvas->SetBottomMargin(0.15);		
	    PlotCanvas->Draw();	
	    
	    TLegend* leg  = new TLegend(0.35,0.89, 0.85, 0.99 ); //On the right side of plot and above border
	    //TLegend* leg  = new TLegend(0.2,0.7,0.55,0.83); //Original from Afro
	    leg->SetBorderSize(0);
	    leg->SetNColumns(1);
	    leg->SetTextSize(LegendTextSize); 
	    leg->SetTextFont(FontStyle);						
	    
	    // Loop over the samples to open the files and to get the corresponding plot
	    double all_integral =0;
	    std::vector<double> integral; integral.resize(NSamples);


	    std::vector<TH1D*> Histos; Histos.resize(NSamples);
	    
	    for (int iSample = 0; iSample < NSamples; iSample++) {	
	      
	      Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotNames[iPlot]+CC1p0piLabels[iSample] ));
	      
	      
	      Histos[iSample]->SetLineWidth(4);
	      Histos[iSample]->SetLineStyle(LineStyle[iSample]);
	      Histos[iSample]->SetLineColor( Colors.at(iSample) );	
	      
	      Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
	      Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
	      Histos[iSample]->GetXaxis()->SetNdivisions(8);
	      Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
	      Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
	      Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
	      Histos[iSample]->GetXaxis()->CenterTitle();						
	      
	      Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
	      Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
	      Histos[iSample]->GetYaxis()->SetNdivisions(6);
	      Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
	      //Histos[iSample]->GetYaxis()->SetTitle("Weighted Events / 1.62E20 POT ");
	      Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
	      Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
	      Histos[iSample]->GetYaxis()->SetTickSize(0);
	      Histos[iSample]->GetYaxis()->CenterTitle();	
	      
	      double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());			
	      Histos[iSample]->GetYaxis()->SetRangeUser(0.,1.1*imax);
	      Histos[0]->GetYaxis()->SetRangeUser(0.,1.1*imax);			
	      
	      PlotCanvas->cd();
	      Histos[iSample]->Draw("hist same");// e");
	      Histos[0]->Draw("hist same");// e");	
	      
	      if (iSample == 0) all_integral = Histos[iSample]->Integral();
	      
	      leg->AddEntry(Histos[iSample],Labels[iSample] + "  (" + to_string(Histos[iSample]->Integral() / all_integral) +  ") "  ,"l");
	      
		
	    } // End of the loop over the samples grabing the plots	
		
	    PlotCanvas->cd();
	    leg->Draw();
	    
	    if (PMissingCut == false) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"_MicroBooNE_Reco.pdf");
	    else if (PMissingCut == true) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"NeutronCuts_MicroBooNE_Reco.pdf");
	    //PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"1MillionEvt_MicroBooNE_Reco.pdf"); 
	    
	    
	  } // End of the loop over the plots
	  
	  
	  
	  // ---------------------------- End of Method 1 ------------------------------ //
	}
	




	else{

	  // --------------------- Method 2: Stacked Histograms ----------------------- //


	  std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors; std::vector<TString> CC1p0piFindLabels; std::vector<TString> InterFindLabels;
	
	  // Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  // Labels.push_back("CC1p0#pi");
	  // Colors.push_back(kCyan+2);	
	  

	  Names.push_back(OutFilePath + "MicroBooNE_Reco.root");
	  Labels.push_back("QE-CC1p0#pi");
	  Colors.push_back(kBlue+2);	
	  CC1p0piFindLabels.push_back(CC1p0piLabels[1]);
	  InterFindLabels.push_back(InteractionLabels[1]);
	  
	  Names.push_back(OutFilePath + "MicroBooNE_Reco.root");
	  Labels.push_back("MEC-CC1p0#pi");
	  Colors.push_back(kOrange+2);	
	  InterFindLabels.push_back(InteractionLabels[1]);
	  CC1p0piFindLabels.push_back(CC1p0piLabels[2]);

	  Names.push_back(OutFilePath + "MicroBooNE_Reco.root");
	  Labels.push_back("RES-CC1p0#pi");
	  Colors.push_back(kGreen+2);	
	  CC1p0piFindLabels.push_back(CC1p0piLabels[1]);
	  InterFindLabels.push_back(InteractionLabels[3]);

	  Names.push_back(OutFilePath + "MicroBooNE_Reco.root");
	  Labels.push_back("DIS-CC1p0#pi");
	  Colors.push_back(kViolet+2);	
	  CC1p0piFindLabels.push_back(CC1p0piLabels[1]);
	  InterFindLabels.push_back(InteractionLabels[4]);


	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("QE-nonCC1p0#pi");
	  Colors.push_back(kPink+2);		
	  CC1p0piFindLabels.push_back(CC1p0piLabels[2]);
	  InterFindLabels.push_back(InteractionLabels[1]);


	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("MEC-nonCC1p0#pi");
	  Colors.push_back(kTeal+2);		
	  CC1p0piFindLabels.push_back(CC1p0piLabels[2]);
	  InterFindLabels.push_back(InteractionLabels[2]);

	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("RES-nonCC1p0#pi");
	  Colors.push_back(kMagenta+2);		
	  CC1p0piFindLabels.push_back(CC1p0piLabels[2]);
	  InterFindLabels.push_back(InteractionLabels[3]);

	  Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	  Labels.push_back("DIS-nonCC1p0#pi");
	  Colors.push_back(kAzure+2);		
	  CC1p0piFindLabels.push_back(CC1p0piLabels[2]);
	  InterFindLabels.push_back(InteractionLabels[4]);


	  const int NSamples = Names.size();
	  std::vector<TFile*> Files; Files.resize(NSamples);

	
      	  // Loop over the samples to open the files and the TTree
	  for (int iSample = 0; iSample < NSamples; iSample++) {
	    Files[iSample] = new TFile(Names[iSample],"readonly");
	  } // End of the loop over the samples



	  for (int iPlot = 0; iPlot < NPlots; iPlot++) {
	  
	    TString CanvasName = "Canvas_" + PlotNames[iPlot];
	    TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	    PlotCanvas->cd();
	    PlotCanvas->SetTopMargin(0.12);
	    PlotCanvas->SetLeftMargin(0.2);
	    PlotCanvas->SetBottomMargin(0.15);		
	    PlotCanvas->Draw();	
	  
	    TLegend* leg  = new TLegend(0.25,0.91, 0.85, 0.99 ); //On the right side of plot and above border
	    //TLegend* leg  = new TLegend(0.2,0.7,0.55,0.83); //Original from Afro
	    leg->SetBorderSize(0);
	    leg->SetNColumns(3);
	    leg->SetTextSize(LegendTextSize); 
	    leg->SetTextFont(FontStyle);						

	  
	    THStack* HistStack = new THStack("HistStack", " ");

	    // Loop over the samples to open the files and to get the corresponding plot 
	    std::vector<TH1D*> Histos; Histos.resize(NSamples);
  
	    for (int iSample = 0; iSample < NSamples; iSample++) {	
	    
	    
	      Histos[iSample] = (TH1D*)(Files[iSample]->Get(InterFindLabels[iSample]+PlotNames[iPlot]+CC1p0piFindLabels[iSample] ));
	    
	      //Histos[iSample]->SetLineWidth(4);
	      // Histos[iSample]->SetLineStyle(LineStyle[iSample]);
	      Histos[iSample]->SetFillColor(Colors.at(iSample) );	
	      Histos[iSample]->SetLineColor(Colors.at(iSample));
	      Histos[iSample]->SetMarkerStyle(21);
	      Histos[iSample]->SetMarkerColor(Colors.at(iSample));
	      // Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
	      // Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
	      // Histos[iSample]->GetXaxis()->SetNdivisions(8);
	      // Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
	      // Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
	      // Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
	      // Histos[iSample]->GetXaxis()->CenterTitle();						
	    
	      // Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
	      // Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
	      // Histos[iSample]->GetYaxis()->SetNdivisions(6);
	      // Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
	      // //Histos[iSample]->GetYaxis()->SetTitle("Weighted Events / 1.62E20 POT ");
	      // Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
	      // Histos[iSample]->GetYaxis()->SetTitleOffset(1.2);
	      // Histos[iSample]->GetYaxis()->SetTickSize(0);
	      // Histos[iSample]->GetYaxis()->CenterTitle();	
	    
	      double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());			
	      Histos[iSample]->GetYaxis()->SetRangeUser(0.,1.1*imax);
	      Histos[0]->GetYaxis()->SetRangeUser(0.,1.1*imax);			
	    
	      PlotCanvas->cd();
	      HistStack->Add(Histos[iSample]);

	      //Histos[iSample]->Draw("hist same");// e");
	      //Histos[0]->Draw("hist same");// e");	
	    
	      leg->AddEntry(Histos[iSample],Labels[iSample],"F");
	    
	    
	    } // End of the loop over the samples grabing the plots	
	  
	    PlotCanvas->cd();
	    HistStack->Draw("hist");
	    leg->Draw();
	  
	  
	    
	    if (PMissingCut == false) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"StackedHisto_MicroBooNE_Reco.pdf");
	    else if (PMissingCut == true) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"StackedHistoNeutronCuts_MicroBooNE_Reco.pdf"); 
	    //PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"1MillionEvtStackedHisto_MicroBooNE_Reco.pdf"); 
	  
	  
	  } // End of the loop over the plots
        

	  //*/
	  // ---------------------- End of Method 2 ------------------------- //
	}
	



	// TCanvas* test = new TCanvas("test", " ");
	// TH1D* h1 = new TH1D("h1", "Test hist" , 10, -4,4);
	// h1->FillRandom("gaus", 2000);
	// h1->SetFillColor(kRed);

	//test->Draw();
	//h1->Draw("hist");
	

} // End of the program
