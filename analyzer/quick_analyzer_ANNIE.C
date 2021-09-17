///////////  Reads a RAT-PAC output file and does a quick analysis ////////
/////////// Author: Vincent FISCHER  ///////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TVector3.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

/* #include "rootstart.h"*/

// Header file for the classes stored in the TTree if any.
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Calib.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DSReader.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

using namespace std;
using namespace TMath;

void quick_analyzer(const char* filename_ratpac) {
	
	// ------------------------------------------------------------------------------------------- //
	// Load RAT libraries (for dsReader)
	gSystem->Load("$(RATPAC_PATH)/lib/libRATEvent.so");
	
	// Initialization
	RAT::DSReader *dsReader;
	RAT::DS::Root *ds;
	RAT::TrackNav *nav;
	RAT::TrackCursor *cursor;
	RAT::TrackNode *node;
	TChain* tri;
	TChain* runtri;
	
	RAT::DS::Run* run;
	RAT::DS::PMTInfo* pmtInfo;
	
	std::clock_t start;
	double duration;
	
	ULong64_t NbEntries;
	
	TVector3 muon_momentum;
	TVector3 unit_vect(0.,0.,1.);
	
	TVector3 init_pos;
	TVector3 fin_pos;
	Double_t init_time, fin_time;
	TString nucl_cap_pdg_code;
	
	TH1::SetDefaultSumw2();
	//SetMyStyle();
	
	// ------------------------------------------------------------------------------------------- //
	// Starts the timer
	start = clock();
	
	gRandom = new TRandom3();
	
	// Output file
	TFile f_output("Analyzer_output.root","RECREATE");
	
	// Load the files
	dsReader = new RAT::DSReader(filename_ratpac);
	NbEntries = dsReader->GetTotal();
	
	// Load the trees 
	// load the ratpac trees and DS
	tri = new TChain("T");
	runtri = new TChain("runT");
	
	if (TString(filename_ratpac).MaybeWildcard()) {
		// Assume there is a runT in all files
		runtri->Add(filename_ratpac);
		RAT::DS::RunStore::SetReadTree(runtri);
	} else {
		// In single file case, we can check
		TFile *ftemp = TFile::Open(filename_ratpac);
		if (ftemp->Get("runT")) {
			runtri->Add(filename_ratpac);
			RAT::DS::RunStore::SetReadTree(runtri);
		} // else, no runT, so don't register runtri with RunStore
		
		delete ftemp;
	}
	
	RAT::DS::Root *branchDS = new RAT::DS::Root();
	tri->SetBranchAddress("ds", &branchDS);
	RAT::DS::RunStore::GetRun(branchDS);
	
	TH1::SetDefaultSumw2(kTRUE);
	
	// Create some histograms
	TH1F* h_HitsPrompt = new TH1F("h_HitsPrompt","Hits for prompt event",500,0, 500);
	TH1F* h_HitsDelayed = new TH1F("h_HitsDelayed","Hits for delayed event",500,0, 500);
	TH1F* h_HitsDelayed_withEdep = new TH1F("h_HitsDelayed_withEdep","Hits for delayed event with an Edep",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min2MeV = new TH1F("h_HitsDelayed_withEdep_min2MeV","Hits for delayed event with an Edep_min250keV",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min5MeV = new TH1F("h_HitsDelayed_withEdep_min5MeV","Hits for delayed event with an Edep_min500keV",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min8MeV = new TH1F("h_HitsDelayed_withEdep_min8MeV","Hits for delayed event with an Edep_min1MeV",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min10MeV = new TH1F("h_HitsDelayed_withEdep_min10MeV","Hits for delayed event with an Edep_min2MeV",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min15MeV = new TH1F("h_HitsDelayed_withEdep_min15MeV","Hits for delayed event with an Edep_min3MeV",50,0, 500);
	TH1F* h_HitsDelayed_withEdep_min20MeV = new TH1F("h_HitsDelayed_withEdep_min20MeV","Hits for delayed event with an Edep_min4MeV",50,0, 500);
	
	TH1F* h_Pair_DeltaT = new TH1F("h_Pair_DeltaT","Pair DeltaT",10000,0, 200000);
	TH1F* h_Pair_DeltaR = new TH1F("h_Pair_DeltaR","Pair DeltaR",2000,0, 2000);
	
	TH1F* h_DeltaE_1MeV_e = new TH1F("h_DeltaE_1MeV_e","1MeV e DeltaE",2000,0, 2000);
	TH1F* h_DeltaR_1MeV_e = new TH1F("h_DeltaR_1MeV_e","1MeV e DeltaR",150,0, 1500);
	TH1F* h_DeltaR_1MeV_e_cumul = new TH1F("h_DeltaR_1MeV_e_cumul","1MeV e DeltaR cumulative",150,0, 1500);
	
	TH2F* h_Edep_to_PE = new TH2F("h_Edep_to_PE","Edep (in BGO) to PE",2000, 0, 20, 500, 0, 500);
	
	
	// Some global variables
	Double_t prompt_window_time_low = 0, prompt_window_time_high = 200;
	Double_t delayed_window_time_low = 1000, delayed_window_time_high = 1000000;
	vector<Double_t> v_prompt_hits_times, v_delayed_hits_times;
	
	UInt_t cluster_threshold = 0;
	UInt_t captures_not_in_water = 0;
	UInt_t captures_after_cuts = 0;
	UInt_t captures_on_AmBe = 0;
	UInt_t double_neutrons = 0.;
	bool b_all_cap_on_water = true;
	bool b_cap_on_AmBe = false;
	bool b_neutron_Edep_BGO = false;
	
	// ------------------------------------------------------------------------------------------- //
	// Analysis loop over all the events
	for (ULong64_t entry=0; entry<NbEntries; ++entry) {
		
		cout << "Entry: " << entry << endl;
		ds = dsReader->GetEvent(entry);    
		run = RAT::DS::RunStore::Get()->GetRun(ds);
		
		// Some initializations
		v_prompt_hits_times.clear(); v_delayed_hits_times.clear();
		b_all_cap_on_water = true;
		b_cap_on_AmBe = false;
		b_neutron_Edep_BGO = false;
		
		// ---------- PMT loop ---------------- //
		for(long iPMT = 0; iPMT < ds->GetMC()->GetMCPMTCount(); iPMT++ ){
			for(long iPhot = 0; iPhot < ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount(); iPhot++){
				
				// Fill prompt and delayed times vector
				if ( ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime() > prompt_window_time_low && 	ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime() < prompt_window_time_high) {
					v_prompt_hits_times.push_back(ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime());
				}
				if ( ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime() > delayed_window_time_low && ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime() < delayed_window_time_high) {
					v_delayed_hits_times.push_back(ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime());
				}
			}
		}
		// -------------------------------------- //
		
		// Skips events with less than XX tubes hits (delayed)
		
		if (v_delayed_hits_times.size() < cluster_threshold) {
			//cout << "Less than " << cluster_threshold << " tubes hit in delayed window... Skipping...\n";
			continue;
		}
		
		captures_after_cuts++;
		
		
		//if (TMath::Sqrt( TMath::Power(ds->GetMC()->GetMCParticle(0)->GetPosition().X(),2) +  TMath::Power(ds->GetMC()->GetMCParticle(0)->GetPosition().Y(),2)) < 1000 && TMath::Abs(ds->GetMC()->GetMCParticle(0)->GetPosition().Z()) < 1650 && !v_prompt_hits_times.empty() && !v_delayed_hits_times.empty()) {
		
		//}
		
		// Calculate the mean times for prompt and delayed
		Double_t mean_prompt_hits_times = 0;
		for(unsigned int i=0; i < v_prompt_hits_times.size(); i++){
			mean_prompt_hits_times += v_prompt_hits_times.at(i);
		}
		mean_prompt_hits_times /= v_prompt_hits_times.size();
		
		Double_t mean_delayed_hits_times = 0;
		for(unsigned int i=0; i < v_delayed_hits_times.size(); i++){
			mean_delayed_hits_times += v_delayed_hits_times.at(i);
		}
		mean_delayed_hits_times /= v_delayed_hits_times.size();
		
		std::sort(v_delayed_hits_times.begin(), v_delayed_hits_times.end());
		if ( v_delayed_hits_times.size() != 0 ){
			if ( (v_delayed_hits_times.at(v_delayed_hits_times.size()-1) - v_delayed_hits_times.at(0)) > 1000) {
				double_neutrons++;
				continue;
			}
		}
		
		

		
		// Find the initial position and final positions of the pair (vertex and n-Gd)
		TVector3 Pair_vertex_start(0.,0.,0.);
		TVector3 Pair_vertex_stop(0.,0.,0.);
		
		
		if (v_delayed_hits_times.size() > 60 && v_delayed_hits_times.size() < 100) {
			cout << "\n\nDelayed PE: " << v_delayed_hits_times.size() << ", time is " << mean_delayed_hits_times << endl;	
			
			cout << ds->GetMC()->GetMCParticleCount() << endl;
			for(int ipart = 0; ipart < ds->GetMC()->GetMCParticleCount(); ipart++){
				cout << ds->GetMC()->GetMCParticle(ipart)->GetPDGCode() <<endl;
			}
		}
		
		for(int ipart = 0; ipart < ds->GetMC()->GetMCParticleCount(); ipart++){
			
			if(ds->GetMC()->GetMCParticle(ipart)->GetPDGCode() == 11 && ds->GetMC()->GetMCParticle(ipart)->GetKE() == 1.){
				h_DeltaE_1MeV_e->Fill(v_prompt_hits_times.size());
				
				std::cout << "PDG ID: "<< ds->GetMC()->GetMCParticle(ipart)->GetPDGCode() << v_prompt_hits_times.size() << std::endl;

				h_DeltaR_1MeV_e->Fill((ds->GetEV(0)->GetCentroid()->GetPosition() - ds->GetMC()->GetMCParticle(ipart)->GetPosition()).Mag());
				h_DeltaR_1MeV_e->Fill(ds->GetEV(0)->GetCentroid()->GetPosition().Mag());
			}
			
			if(ds->GetMC()->GetMCParticle(ipart)->GetPDGCode() == -11){
				Pair_vertex_start = ds->GetMC()->GetMCParticle(ipart)->GetPosition();
			}
			
		}
		
		if (v_prompt_hits_times.size() > 0){
			for(int itrack = 0; itrack < ds->GetMC()->GetMCTrackCount(); itrack++){
				if(ds->GetMC()->GetMCTrack(itrack)->GetParticleName() != "opticalphoton")
					cout << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetProcess() << endl;
			}
		}
			
			for(int itrack = 0; itrack < ds->GetMC()->GetMCTrackCount(); itrack++){
				
				if(ds->GetMC()->GetMCTrack(itrack)->GetPDGCode() == 2112 && 
					ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetProcess() != "Transportation" && 
					ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_BGO" ) {
					cout << "Neutron Edep in BGO! Edep before =" << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetKE() << ", Edep after = " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(1)->GetKE() << ", process = " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetProcess() << endl;
					b_neutron_Edep_BGO = true;
				}
				
				if (ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetProcess() == "nCapture"){
					// 					cout << "Capture in volume -> " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() << " " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetKE() << endl;
					if (ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() != "detector") {
						b_all_cap_on_water = false;
						cout << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() << endl;
					}
					if (v_delayed_hits_times.size() > 60 && v_delayed_hits_times.size() < 100) {
						cout << "Capture location is " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() << endl;
						cout << "Gamma energy: " << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetKE() << ", ID is " << ds->GetMC()->GetMCTrack(itrack)->GetID() << ", parent ID is " << ds->GetMC()->GetMCTrack(itrack)->GetParentID() << endl;
						for(int ipart = 0; ipart < ds->GetMC()->GetMCParticleCount(); ipart++){
							cout << "Initial neutron energy: " << ds->GetMC()->GetMCParticle(ipart)->GetKE() << endl;
						}
					}
					if (ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_housing_cap" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_housing_body" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_housing_air" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_BGO" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_placer" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_case" ||
						ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() == "AmBe_source" ) {
						b_cap_on_AmBe = true;
						cout << ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetVolume() << endl;
					}
					Pair_vertex_stop = ds->GetMC()->GetMCTrack(itrack)->GetMCTrackStep(0)->GetEndpoint();
				}
			}
		
		if (!b_all_cap_on_water) {
			captures_not_in_water++;
			cout << "Not a capture on water" << endl;
		}
		if (b_cap_on_AmBe) {
			captures_on_AmBe++;
			cout << "Capture on AmBe source assembly" << endl;
		}
		if (b_neutron_Edep_BGO){
			continue;
		}
		
		//if (TMath::Sqrt( TMath::Power(ds->GetMC()->GetMCParticle(0)->GetPosition().X(),2) +  TMath::Power(ds->GetMC()->GetMCParticle(0)->GetPosition().Y(),2)) < 1000 && TMath::Abs(ds->GetMC()->GetMCParticle(0)->GetPosition().Z()) < 1650 && !v_prompt_hits_times.empty() && !v_delayed_hits_times.empty()) {
		h_Pair_DeltaT->Fill(mean_delayed_hits_times - 0);
		h_Pair_DeltaR->Fill((Pair_vertex_stop - Pair_vertex_start).Mag()*gRandom->Gaus(0.0, 450.0/TMath::Sqrt(v_prompt_hits_times.size()/140.)));
		//}
		
		h_HitsPrompt->Fill(v_prompt_hits_times.size());
		h_HitsDelayed->Fill(v_delayed_hits_times.size());
		
		h_Edep_to_PE->Fill(ds->GetMC()->GetMCSummary()->GetTotalScintEdep(), v_prompt_hits_times.size());
		
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 0){
			h_HitsDelayed_withEdep->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 2.0){
			h_HitsDelayed_withEdep_min2MeV->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 5.0){
			h_HitsDelayed_withEdep_min5MeV->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 8.0){
			h_HitsDelayed_withEdep_min8MeV->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 10.0){
			h_HitsDelayed_withEdep_min10MeV->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 15.0){
			h_HitsDelayed_withEdep_min15MeV->Fill(v_delayed_hits_times.size());
		}
		if (ds->GetMC()->GetMCSummary()->GetTotalScintEdep() > 20.0){
			h_HitsDelayed_withEdep_min20MeV->Fill(v_delayed_hits_times.size());
		}
		
		// kind of progress bar...
		if (NbEntries > 10) {
			if ( entry%(NbEntries/10) == 0 ) { 
				cout << "Evt " << entry << " out of " << NbEntries << " events ===> " << Floor(Double_t(entry)/Double_t(NbEntries)*100.) << " %\n";
			}
		}
		// ---------- //
		
		// 		cout << "Prompt hits: " << v_prompt_hits_times.size() << endl;
		// 		cout << "Delayed hits: " << v_delayed_hits_times.size() << endl;
		// 		cout << "Prompt time: " << mean_prompt_hits_times << endl;
		// 		cout << "Delayed time: " << mean_delayed_hits_times << endl;
		// 		cout << "Prompt pos: " << Pair_vertex_start.X() << " " << Pair_vertex_start.Y() << " " << Pair_vertex_start.Z() << endl;
		// 		cout << "Delayed pos: " << Pair_vertex_stop.X() << " " << Pair_vertex_stop.Y() << " " << Pair_vertex_stop.Z() << endl;
		// 		cout << "-----------------------\n";
		
	}
	// ------------------------------------------------------------------------------------------- //
	
	delete run;
	delete tri, runtri, branchDS;
	delete dsReader;
	
	//Fill cumulative DeltaR (1 MeV e-) hist
	Double_t cumul = 0.;
	for(int bin = 0; bin < h_DeltaR_1MeV_e->GetNbinsX(); bin++){
		cumul += h_DeltaR_1MeV_e->GetBinContent(bin);
		h_DeltaR_1MeV_e_cumul->Fill(h_DeltaR_1MeV_e->GetBinCenter(bin), cumul);
	}
	h_DeltaR_1MeV_e_cumul->Scale(100./h_DeltaR_1MeV_e_cumul->GetBinContent(h_DeltaR_1MeV_e->GetNbinsX()-1));
	
	f_output.cd();
	h_HitsPrompt->SetXTitle("Photoelectrons");
	h_HitsPrompt->SetYTitle("Events");
	h_HitsPrompt->Write();
	
	h_HitsDelayed->SetXTitle("Photoelectrons");
	h_HitsDelayed->SetYTitle("Events");
	h_HitsDelayed->Write();
	
	h_HitsDelayed_withEdep->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep->SetYTitle("Events");
        h_HitsDelayed_withEdep->SetOption("hist");
	h_HitsDelayed_withEdep->Write();
	
	h_HitsDelayed_withEdep_min2MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min2MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min2MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min2MeV->Write();
	
	h_HitsDelayed_withEdep_min5MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min5MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min5MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min5MeV->Write();
	
	h_HitsDelayed_withEdep_min8MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min8MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min8MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min8MeV->Write();
	
	h_HitsDelayed_withEdep_min10MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min10MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min10MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min10MeV->Write();
	
	h_HitsDelayed_withEdep_min15MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min15MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min15MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min15MeV->Write();
	
	h_HitsDelayed_withEdep_min20MeV->SetXTitle("Photoelectrons");
	h_HitsDelayed_withEdep_min20MeV->SetYTitle("Events");
	h_HitsDelayed_withEdep_min20MeV->SetOption("hist");
	h_HitsDelayed_withEdep_min20MeV->Write();
	
	h_Pair_DeltaT->SetXTitle("#DeltaT [ns]");
	h_Pair_DeltaT->SetYTitle("Events");
	h_Pair_DeltaT->Write();
	
	h_Pair_DeltaR->SetXTitle("#DeltaR [mm]");
	h_Pair_DeltaR->SetYTitle("Events");
	h_Pair_DeltaR->Write();
	
	h_DeltaR_1MeV_e->SetXTitle("Reconstructed #DeltaR [mm]");
	h_DeltaR_1MeV_e->SetYTitle("Events");
	h_DeltaR_1MeV_e->Write();
	
	h_DeltaE_1MeV_e->SetXTitle("Photoelectrons");
	h_DeltaE_1MeV_e->SetYTitle("Events");
	h_DeltaE_1MeV_e->Write();
	
	h_Edep_to_PE->SetXTitle("Energy deposited (scint) [MeV]");
	h_Edep_to_PE->SetYTitle("Photoelectrons");
	h_Edep_to_PE->Write();
	
	h_DeltaR_1MeV_e_cumul->SetXTitle("Reconstructed #DeltaR [mm]");
	h_DeltaR_1MeV_e_cumul->SetYTitle("Probability [%/10 mm]");
	h_DeltaR_1MeV_e_cumul->Write();
	
	f_output.Close();
	
	cout << "There were " << captures_not_in_water << " captures not on water (and " << captures_on_AmBe << " on the AmBe assembly) out of " << captures_after_cuts << " total captures after cuts\n";
	cout << "Event with what appears to be a double neutron capture : " << double_neutrons << endl;
	
	// Ends the timer
	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Execution time: " << duration << " seconds\n";
	
}
