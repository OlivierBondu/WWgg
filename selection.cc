// Basic WWgg selection
// O. Bondu (June 2013)
// C++ headers
#include <iostream>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

int main()
{
//	TFile *infile = TFile::Open("root://eoscms//eos/cms/store/group/phys_higgs/Resonant_HH/trees/radion_tree_v01/Data_full2012.root");
	TFile *infile = TFile::Open("Data_full2012.root");
	TFile *outfile = new TFile("data.root", "RECREATE");
	TTree *intree = (TTree*)infile->Get("Data");
	TTree *outtree = new TTree("data", "data preselected");

	// setup tree inputs
	float ph1_eta, ph2_eta, ph1_pt, ph2_pt, PhotonsMass, ph1_phi, ph2_phi, ph1_e, ph2_e;
	int ph1_ciclevel, ph2_ciclevel;
	float pu_n, nvtx, rho;
	float weight, evweight, pu_weight;
	float j1_e, j1_pt, j1_phi, j1_eta, j1_betaStarClassic, j1_dR2Mean;
	float j2_e, j2_pt, j2_phi, j2_eta, j2_betaStarClassic, j2_dR2Mean;
	float j3_e, j3_pt, j3_phi, j3_eta, j3_betaStarClassic, j3_dR2Mean;
	float j4_e, j4_pt, j4_phi, j4_eta, j4_betaStarClassic, j4_dR2Mean;
	// intermediate storage
	float jet_e, jet_pt, jet_phi, jet_eta;
	float jet_betaStarClassic, jet_dR2Mean;
	// setup tree outputs
	float pho1_pt, pho1_e, pho1_phi, pho1_eta, pho1_mass;
	float pho2_pt, pho2_e, pho2_phi, pho2_eta, pho2_mass;
	float jet1_pt, jet1_e, jet1_phi, jet1_eta, jet1_mass;
	float jet2_pt, jet2_e, jet2_phi, jet2_eta, jet2_mass;
	float jj_pt, jj_e, jj_phi, jj_eta, jj_mass;
	float gg_pt, gg_e, gg_phi, gg_eta, gg_mass;

	intree->SetBranchAddress("ph1_eta", &ph1_eta);
	intree->SetBranchAddress("ph2_eta", &ph2_eta);
	intree->SetBranchAddress("ph1_pt", &ph1_pt);
	intree->SetBranchAddress("ph2_pt", &ph2_pt);
	intree->SetBranchAddress("ph1_phi", &ph1_phi);
	intree->SetBranchAddress("ph2_phi", &ph2_phi);
	intree->SetBranchAddress("ph1_e", &ph1_e);
	intree->SetBranchAddress("ph2_e", &ph2_e);
	intree->SetBranchAddress("PhotonsMass", &PhotonsMass);
	intree->SetBranchAddress("ph1_ciclevel", &ph1_ciclevel);
	intree->SetBranchAddress("ph2_ciclevel", &ph2_ciclevel);
	intree->SetBranchAddress("pu_n", &pu_n);
	intree->SetBranchAddress("nvtx", &nvtx);
	intree->SetBranchAddress("rho", &rho);
	intree->SetBranchAddress("weight", &weight);
	intree->SetBranchAddress("evweight", &evweight);
	intree->SetBranchAddress("pu_weight", &pu_weight);

	intree->SetBranchAddress("j1_e", &j1_e);
	intree->SetBranchAddress("j1_pt", &j1_pt);
	intree->SetBranchAddress("j1_phi", &j1_phi);
	intree->SetBranchAddress("j1_eta", &j1_eta);
	intree->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic);
	intree->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean);
	intree->SetBranchAddress("j2_e", &j2_e);
	intree->SetBranchAddress("j2_pt", &j2_pt);
	intree->SetBranchAddress("j2_phi", &j2_phi);
	intree->SetBranchAddress("j2_eta", &j2_eta);
	intree->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic);
	intree->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean);
	intree->SetBranchAddress("j3_e", &j3_e);
	intree->SetBranchAddress("j3_pt", &j3_pt);
	intree->SetBranchAddress("j3_phi", &j3_phi);
	intree->SetBranchAddress("j3_eta", &j3_eta);
	intree->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic);
	intree->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean);
	intree->SetBranchAddress("j4_e", &j4_e);
	intree->SetBranchAddress("j4_pt", &j4_pt);
	intree->SetBranchAddress("j4_phi", &j4_phi);
	intree->SetBranchAddress("j4_eta", &j4_eta);
	intree->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic);
	intree->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean);

	outtree->Branch("weight", &weight, "weight/F");
	outtree->Branch("evweight", &evweight, "evweight/F");
	outtree->Branch("pu_weight", &pu_weight, "pu_weight/F");
	outtree->Branch("pho1_pt", &pho1_pt, "pho1_pt/F");
	outtree->Branch("pho1_e", &pho1_e, "pho1_e/F");
	outtree->Branch("pho1_phi", &pho1_phi, "pho1_phi/F");
	outtree->Branch("pho1_eta", &pho1_eta, "pho1_eta/F");
	outtree->Branch("pho1_mass", &pho1_mass, "pho1_mass/F");
	outtree->Branch("pho2_pt", &pho2_pt, "pho2_pt/F");
	outtree->Branch("pho2_e", &pho2_e, "pho2_e/F");
	outtree->Branch("pho2_phi", &pho2_phi, "pho2_phi/F");
	outtree->Branch("pho2_eta", &pho2_eta, "pho2_eta/F");
	outtree->Branch("pho2_mass", &pho2_mass, "pho2_mass/F");
	outtree->Branch("jet1_pt", &jet1_pt, "jet1_pt/F");
	outtree->Branch("jet1_e", &jet1_e, "jet1_e/F");
	outtree->Branch("jet1_phi", &jet1_phi, "jet1_phi/F");
	outtree->Branch("jet1_eta", &jet1_eta, "jet1_eta/F");
	outtree->Branch("jet1_mass", &jet1_mass, "jet1_mass/F");
	outtree->Branch("jet2_pt", &jet2_pt, "jet2_pt/F");
	outtree->Branch("jet2_e", &jet2_e, "jet2_e/F");
	outtree->Branch("jet2_phi", &jet2_phi, "jet2_phi/F");
	outtree->Branch("jet2_eta", &jet2_eta, "jet2_eta/F");
	outtree->Branch("jet2_mass", &jet2_mass, "jet2_mass/F");
	outtree->Branch("jj_pt", &jj_pt, "jj_pt/F");
	outtree->Branch("jj_e", &jj_e, "jj_e/F");
	outtree->Branch("jj_phi", &jj_phi, "jj_phi/F");
	outtree->Branch("jj_eta", &jj_eta, "jj_eta/F");
	outtree->Branch("jj_mass", &jj_mass, "jj_mass/F");
	outtree->Branch("gg_pt", &gg_pt, "gg_pt/F");
	outtree->Branch("gg_e", &gg_e, "gg_e/F");
	outtree->Branch("gg_phi", &gg_phi, "gg_phi/F");
	outtree->Branch("gg_eta", &gg_eta, "gg_eta/F");
	outtree->Branch("gg_mass", &gg_mass, "gg_mass/F");


	int nevents[30] = {0};
	string eventcut[30];
	int njets[20] = {0};
  int decade = 0;
  int totevents = intree->GetEntries();
  if(DEBUG) totevents = 500;
  cout << "#entries= " << totevents << endl;
  // loop over events
  for(int ievt=0 ; ievt < totevents ; ievt++)
  {
    double progress = 10.0*ievt/(1.0*totevents);
    int k = TMath::FloorNint(progress);
    if (k > decade) cout<<10*k<<" %"<<endl;
    decade = k;

    intree->GetEntry(ievt);
	
		// Apply photon ID cuts
		nevents[0]++; eventcut[0] = "Before photon ID";
		if( (fabs(ph1_eta) > 2.5) || (fabs(ph2_eta) > 2.5) ) continue;
		nevents[1]++; eventcut[1] = "After eta < 2.5";
		if( (fabs(ph1_eta) < 1.566) && (fabs(ph1_eta) >1.4442) ) continue;
		nevents[2]++; eventcut[2] = "After eta gap for photon 1";
		if( (fabs(ph2_eta) < 1.566) && (fabs(ph2_eta) >1.4442) ) continue;
		nevents[3]++; eventcut[3] = "After eta gap for photon 2";
		if( ph1_pt < (float)(40.*PhotonsMass)/(float)120. ) continue;
		nevents[4]++; eventcut[4] = "After floating pt cut for photon 1";
		if( ph2_pt < 25. ) continue;
		nevents[5]++; eventcut[5] = "After fixed pt cut for photon 2";
		if(DEBUG) cout << "ph1_ciclevel= " << ph1_ciclevel << "\tph2_ciclevel= " << ph2_ciclevel << endl;
		if( (ph1_ciclevel < 4) || (ph2_ciclevel < 4) ) continue;
		nevents[6]++; eventcut[6] = "After cic cut on both photons";
		if(DEBUG) cout << "PhotonsMass= " << PhotonsMass << endl;
		if( (PhotonsMass < 130.) ) continue;
		nevents[7]++; eventcut[7] = "After mgg > 130";

		int njets_passing_kLooseID = 0;
		if( j1_pt > 0.) njets_passing_kLooseID++;
		if( j2_pt > 0.) njets_passing_kLooseID++;
		if( j3_pt > 0.) njets_passing_kLooseID++;
		if( j4_pt > 0.) njets_passing_kLooseID++;
		if(DEBUG) cout << "njets_passing_kLooseID= " << njets_passing_kLooseID << endl;
		// take only the subset of events where at least two jets remains
		if( njets_passing_kLooseID < 2 ) continue;
		nevents[8]++; eventcut[8] = "After njet > 2";

		TLorentzVector jet;
		vector<float> jetPt;
		vector<float> jetE;
		vector<float> jetEta;
		vector<float> jetPhi;
		jetPt.clear();
		jetE.clear();
		jetEta.clear();
		jetPhi.clear();

		// loop over jets, store jet info + info on closest genjet / parton (no selection applied)
		for( int ijet = 0 ; ijet < min(njets_passing_kLooseID, 4); ijet ++ )
		{
			njets[0]++;
			if( ijet == 0 )
			{
				jet_e = j1_e;
				jet_pt = j1_pt;
				jet_phi = j1_phi;
				jet_eta = j1_eta;
				jet_betaStarClassic = j1_betaStarClassic;
				jet_dR2Mean = j1_dR2Mean;
				jet.SetPtEtaPhiE(j1_pt, j1_eta, j1_phi, j1_e);
			} // end if jet == 0

			if( ijet == 1 )
			{
				jet_e = j2_e;
				jet_pt = j2_pt;
				jet_phi = j2_phi;
				jet_eta = j2_eta;
				jet_betaStarClassic = j2_betaStarClassic;
				jet_dR2Mean = j2_dR2Mean;
				jet.SetPtEtaPhiE(j2_pt, j2_eta, j2_phi, j2_e);
			} // end if jet == 1

			if( ijet == 2 )
			{
				jet_e = j3_e;
				jet_pt = j3_pt;
				jet_phi = j3_phi;
				jet_eta = j3_eta;
				jet_betaStarClassic = j3_betaStarClassic;
				jet_dR2Mean = j3_dR2Mean;
				jet.SetPtEtaPhiE(j3_pt, j3_eta, j3_phi, j3_e);
			} // end if jet == 2

			if( ijet == 3 )
			{
				jet_e = j4_e;
				jet_pt = j4_pt;
				jet_phi = j4_phi;
				jet_eta = j4_eta;
				jet_betaStarClassic = j4_betaStarClassic;
				jet_dR2Mean = j4_dR2Mean;
				jet.SetPtEtaPhiE(j4_pt, j4_eta, j4_phi, j4_e);
			} // end if jet == 3

			if(DEBUG) cout << "\tjet_pt= " << jet_pt << "\tjet_eta= " << jet_eta << "\tjet_betaStarClassic/log(nvtx-0.64)= " << jet_betaStarClassic / log(nvtx-0.64) << "\tjet_dR2Mean= " << jet_dR2Mean << endl;
			// jet selection
			// ** acceptance + pu id **
			if( jet_pt < 25. ) continue;
			if(DEBUG) cout << "jet[" << ijet << "] survives pt cut" << endl;
			njets[1]++;
			if( fabs(jet_eta) > 4.7 ) continue;
			if(DEBUG) cout << "jet[" << ijet << "] survives eta cut" << endl;
			njets[2]++;
			if( (fabs(jet_eta)<2.5) && ((jet_betaStarClassic > 0.2*log(nvtx-0.64)) || (jet_dR2Mean>0.06)) ) continue;
			else if( (fabs(jet_eta)>2.5) && (fabs(jet_eta)<2.75) && ((jet_betaStarClassic > 0.3*log(nvtx-0.64)) || (jet_dR2Mean>0.05)) ) continue;
			else if( (fabs(jet_eta)>2.75) && (fabs(jet_eta)<3.) && (jet_dR2Mean>0.05) ) continue;
			else if( (fabs(jet_eta)>3.) && (fabs(jet_eta)<4.7) && (jet_dR2Mean>0.055) ) continue;
			if(DEBUG) cout << "jet[" << ijet << "] survives pu rejection" << endl;
			njets[3]++;
			// ** store 4-momentum + csv output for combinatorics **
			jetPt.push_back(jet_pt);
			jetE.push_back(jet_e);
			jetEta.push_back(jet_eta);
			jetPhi.push_back(jet_phi);
		} // end of loop over jets
	
		if(DEBUG) cout << "Njets left after jet cuts: jetPt.size()= " << jetPt.size() << endl;	
		// jet combinatorics
//		if( jetPt.size() < 2 ) continue;
//		nevents[9]++; eventcut[9] = "njet> 2 After jet pt,eta cut and pu id";
		if( jetPt.size() != 2 ) continue;
		nevents[9]++; eventcut[9] = "njet == 2 After jet pt,eta cut and pu id";

		int ij1 = 0;
		int ij2 = 1;
		TLorentzVector pho1;
		TLorentzVector pho2;
		TLorentzVector jet1;
		TLorentzVector jet2;
		TLorentzVector regjet1;
		TLorentzVector regjet2;
		pho1.SetPtEtaPhiE(ph1_pt, ph1_eta, ph1_phi, ph1_e);
		pho2.SetPtEtaPhiE(ph2_pt, ph2_eta, ph2_phi, ph2_e);
		jet1.SetPtEtaPhiE(jetPt[ij1], jetEta[ij1], jetPhi[ij1], jetE[ij1]);
		jet2.SetPtEtaPhiE(jetPt[ij2], jetEta[ij2], jetPhi[ij2], jetE[ij2]);
		TLorentzVector jj = jet1 + jet2;
		TLorentzVector gg = pho1 + pho2;

		pho1_pt = pho1.Pt();
		pho1_e = pho1.E();
		pho1_phi = pho1.Phi();
		pho1_eta = pho1.Eta();
		pho1_mass = pho1.M();
		pho2_pt = pho2.Pt();
		pho2_e = pho2.E();
		pho2_phi = pho2.Phi();
		pho2_eta = pho2.Eta();
		pho2_mass = pho2.M();
		jet1_pt = jet1.Pt();
		jet1_e = jet1.E();
		jet1_phi = jet1.Phi();
		jet1_eta = jet1.Eta();
		jet1_mass = jet1.M();
		jet2_pt = jet2.Pt();
		jet2_e = jet2.E();
		jet2_phi = jet2.Phi();
		jet2_eta = jet2.Eta();
		jet2_mass = jet2.M();
		jj_pt = jj.Pt();
		jj_e = jj.E();
		jj_phi = jj.Phi();
		jj_eta = jj.Eta();
		jj_mass = jj.M();
		gg_pt = gg.Pt();
		gg_e = gg.E();
		gg_phi = gg.Phi();
		gg_eta = gg.Eta();
		gg_mass = gg.M();

		if(DEBUG) cout << "jj_mass= " << jj_mass << endl;
		if( jj_mass < 300 ) continue;
		nevents[10]++; eventcut[10] = "m_jj > 300";
		if( jj_mass < 400 ) continue;
		nevents[11]++; eventcut[11] = "m_jj > 400";
		if( jj_mass < 500 ) continue;
		nevents[12]++; eventcut[12] = "m_jj > 500";
		if( jj_mass < 600 ) continue;
		nevents[13]++; eventcut[13] = "m_jj > 600";
		if( jj_mass < 700 ) continue;
		nevents[14]++; eventcut[14] = "m_jj > 700";
		if( jj_mass < 800 ) continue;
		nevents[15]++; eventcut[15] = "m_jj > 800";
		if( jj_mass < 900 ) continue;
		nevents[16]++; eventcut[16] = "m_jj > 900";
		if( jj_mass < 1000 ) continue;
		nevents[17]++; eventcut[17] = "m_jj > 1000";
		if( jj_mass < 1100 ) continue;
		nevents[18]++; eventcut[18] = "m_jj > 1100";
		if( jj_mass < 1200 ) continue;
		nevents[19]++; eventcut[19] = "m_jj > 1200";
		if( jj_mass < 1300 ) continue;
		nevents[20]++; eventcut[20] = "m_jj > 1300";
		if( jj_mass < 1400 ) continue;
		nevents[21]++; eventcut[21] = "m_jj > 1400";
		if( jj_mass < 1500 ) continue;
		nevents[22]++; eventcut[22] = "m_jj > 1500";
		if( jj_mass < 1600 ) continue;
		nevents[23]++; eventcut[23] = "m_jj > 1600";
		if( jj_mass < 1700 ) continue;
		nevents[24]++; eventcut[24] = "m_jj > 1700";
		if( jj_mass < 1800 ) continue;
		nevents[25]++; eventcut[25] = "m_jj > 1800";
		if( jj_mass < 1900 ) continue;
		nevents[26]++; eventcut[26] = "m_jj > 1900";
		if( jj_mass < 2000 ) continue;
		nevents[27]++; eventcut[27] = "m_jj > 2000";


		outtree->Fill();

	} // end of loop over events


	for(int i=0 ; i < 28 ; i++)
    cout << "#nevents[" << i << "]= " << nevents[i] << "\teventcut[" << i << "]= " << eventcut[i] << endl;
	cout << endl;
	for(int i=0 ; i < 3 ; i++)
    cout << "#njets[" << i << "]= " << njets[i] << endl;

  outfile->cd();
  outtree->Write();
  outfile->Close();
  infile->Close();



	return 0;
}
