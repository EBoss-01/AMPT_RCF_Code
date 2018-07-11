//-------------------------------------------------------------------------------------
// This code does the participant plane calculations using partons and hadrons.
// The parton file used here is different than the npart-xy.dat file used other times.
// The files used can be seen below
//
// File 1: parton-initial-afterPropagation.dat
// File 2: ampt.dat
//
// Created: 07-10-18
// Last Updated: 07-11-18
//-------------------------------------------------------------------------------------

#include "TLatex.h"
#include "TChain.h"
#include "TH1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TRegexp.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TPaletteAxis.h"
#include "TTree.h"
#include "TF1.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//-----------------------------------------
// Structures for partons
//-----------------------------------------
// This stucture is used for File 1.
struct parton {
//Event Number
	int evtN;
// Parton ID
	int PID;
// Momenta at birth
	float px;
	float py;
	float pz;
// Mass of parton
	float m;
// Spacetime at birth
	float x;
	float y;
	float z;
	float t;
//other
	float eta;
	float pT;
};
// Possible structure for File 2. If using hadrons.
struct finalhadrons {
// Event Number
	int evtN;
// Parton ID
	int PID;
// Momenta
	float px;
	float py;
	float pz;
// Mass
	float m;
// Other
	float eta;
	float phi;
	float pT;
};

//-----------------------------------------------
// Variables
//-----------------------------------------------

// Event Counter
int eventnumber = 0;
// Vectors
vector<parton> initialpartons;
vector<finalhadrons> finalstatehadrons;
vector<float> psi2values;
vector<float> psi3values;
vector<float> epsilon2values;
vector<float> epsilon3values;
// Global variables
float psi2;
float psi3;
float AverageE2;
float AverageE3;
float TotalE2;
float TotalE3;

//-------------------------------------------------
// Functions
//-------------------------------------------------

void calculateParticipantPlanes() {

	// Variables for calculations
	float cmx = 0;
	float cmy = 0;
	float q2x = 0;
	float q2y = 0;
	float q3x = 0;
	float q3y = 0;
	float epsilon2 = 0;
	float epsilon3 = 0;
	float rsquare = 0;

	// vector for use.
	vector<parton> v = initialpartons;

	// Loop through partons.
	for (unsigned int k = 0; k < v.size(); k++) {
		cmx = cmx + v[k].x;
		cmy = cmy + v[k].y;
	}

	int partoncounter = v.size();

	cmx = cmx/(float)partoncounter;
	cmy = cmy/(float)partoncounter;

	for (unsigned int i = 0; i < v.size(); i++) {

		v[i].x = v[i].x - cmx;
		v[i].y = v[i].y - cmy;

		float phivalue = TMath::ATan2(v[i].y, v[i].x);
		float rvalue = TMath::Sqrt(pow(v[i].x,2) + pow(v[i].y,2));
		// Component calculations
		q2x = q2x + pow(rvalue,2)*TMath::Cos(2*phivalue);
		q2y = q2y + pow(rvalue,2)*TMath::Sin(2*phivalue);
		q3x = q3x + pow(rvalue,2)*TMath::Cos(3*phivalue);
		q3y = q3y + pow(rvalue,2)*TMath::Sin(3*phivalue);
		rsquare = rsquare + pow(rvalue,2);

	}

	// Some averages
	q2x = q2x/(float)partoncounter;
	q2y = q2y/(float)partoncounter;
	q3x = q3x/(float)partoncounter;
	q3y = q3y/(float)partoncounter;
	rsquare = rsquare/(float)partoncounter;

	// Participant Plane calculations
	psi2 = (TMath::ATan2(q2y, q2x)/2.0) + (TMath::Pi()/2.0);
	psi3 = (TMath::ATan2(q3y, q3x)/3.0) + (TMath::Pi()/3.0);
	// Eccentricity calculations
	epsilon2 = (TMath::Sqrt(pow(q2x,2) + pow(q2y,2)))/rsquare;
	epsilon3 = (TMath::Sqrt(pow(q3x,2) + pow(q3y,2)))/rsquare;
	// Pushing the values to more vectors
	psi2values.push_back(psi2);
	psi3values.push_back(psi3);
	epsilon2values.push_back(epsilon2);
	epsilon3values.push_back(epsilon3);
}

void calculateFlow(TProfile *v2plot, TProfile *v3plot) {

	float v2 = 0;
	float v3 = 0;

	for (unsigned int i = 0; i < psi2values.size(); i++) {

		for (unsigned int j = 0; j < finalstatehadrons.size(); j++) {

			if (finalstatehadrons[j].evtN == (i+1)) {

				if (fabs(finalstatehadrons[j].eta) < 0.35) {

					v2 = TMath::Cos(2*(finalstatehadrons[j].phi - psi2values[i]));
					v3 = TMath::Cos(3*(finalstatehadrons[j].phi - psi3values[i]));
					// Filling the TProfiles.
					v2plot->Fill(finalstatehadrons[j].pT,v2);
					v3plot->Fill(finalstatehadrons[j].pT,v3);
				}
			}
		}
	}
}

//-------------------------------------------------
// Main Code
//-------------------------------------------------

void partonParticipantPlaneCalculation(void) {

	cout << "Running partonParticipantPlaneCalculation()..." << endl;
	// Varialbes needed for reading in file 1.
	string line;
	int iterationN;
	int Npartons;
	int Nbaryons;
	int Nmesons;
	int particlesN;
	int nonpartparticles;
	// Needed for lines after header.
	int pid;
	float momenta[3];
	float mass;
	float spacetime[4];
	// File initialization
	ifstream myInitialFileInfo;
	ifstream myAmptFile;

	//Opening
	myInitialFileInfo.open("ana/parton-initial-afterPropagation.dat");

	if (!myInitialFileInfo) {
		// Informs user if file doesn't open
		cout << "Unable to open File 1." << endl;
		return;
	}
	else {
		// Tells when file opens.
		cout << "Successfully opened File 1." << endl;
	}

	// Begin reading the file
	while (std::getline(myInitialFileInfo, line)) {

		// Prevents double reading of last line
		if (!myInitialFileInfo) break;

		stringstream heading(line);

		if (heading >> eventnumber >> iterationN >> Npartons >> Nbaryons >> Nmesons >> particlesN >> nonpartparticles) {

			// Prevents duplicate events
			if (iterationN > 0) {
				continue;
			}

			// Loop that goes through all partons.
			for (int i = 0; i < Npartons; i++) {

				myInitialFileInfo >> pid >> momenta[0] >> momenta[1] >> momenta[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

				float energy1 = TMath::Sqrt(pow(momenta[0],2) + pow(momenta[1],2) + pow(momenta[2],2) + pow(mass,2));
				TLorentzVector cut(momenta[0], momenta[1], momenta[2], energy1);

				parton partinitial;
				partinitial.evtN = eventnumber;
				partinitial.PID = pid;
				partinitial.px = momenta[0];
				partinitial.py = momenta[1];
				partinitial.pz = momenta[2];
				partinitial.m = mass;
				partinitial.x = spacetime[0];
				partinitial.y = spacetime[1];
				partinitial.z = spacetime[2];
				partinitial.t = spacetime[3];
				partinitial.eta = cut.Eta();
				partinitial.pT = cut.Pt();

				if (fabs(partinitial.eta) < 3 && partinitial.t < 3) {
					initialpartons.push_back(partinitial);
				}
			}

			// Call to calculate planes
			calculateParticipantPlanes();

			initialpartons.clear();
		}
	}

	myInitialFileInfo.close();

	// Moving to the second file
	//Variables for reading File 2.
	string line1;
	int evt;
	int test;
	int npartevt;
	float imparm[2];
	int totalpart;
	int participantproj;
	int participanttarg;
	int collisions[4];
	// Variables for reading in everything after header.
	int particleid;
	float momentas[3];
	float pmass;
	float locationt[4];

	myAmptFile.open("ana/ampt.dat");

	if (!myAmptFile) {
		cout << "Unable to open File 2." << endl;
		return;
	}
	else {
		cout << "Successfully open File 2." << endl;
	}

	TCanvas *c = new TCanvas("c","v_{2} calculations",700,700);
	gStyle->SetOptStat(0);
	TProfile *v2plot = new TProfile("v2plot","v_{2} vs. p_{T}",20,0,2.5,-1,1);

	TCanvas *c1 = new TCanvas("c1","v_{3} calculation",700,700);
	gStyle->SetOptStat(0);
	TProfile *v3plot = new TProfile("v3plot","v_{3} vs. p_{T}",20,0,2.5,-1,1);

	while (std::getline(myAmptFile, line1)) {

		// to avoid reading last line twice
		if (!myAmptFile) break;

		stringstream header(line1);

		if (header >> evt >> test >> npartevt >> imparm[0] >> totalpart >> participantproj >> participanttarg >> collisions[0] >> collisions[1] >> collisions[2] >> collisions[3] >> imparm[1]) {

			// Looping through hadrons
			for (int j = 0; j < npartevt; j++) {

				myAmptFile >> particleid >> momentas[0] >> momentas[1] >> momentas[2] >> pmass >> locationt[0] >> locationt[1] >> locationt[2] >> locationt[3];

				if (TMath::Sqrt(pow(momentas[0],2) + pow(momentas[1],2)) < 0.0001) {
					continue;
				}

				if (abs(particleid) != 211 && abs(particleid) != 321 && abs(particleid) != 2212) continue;

				finalhadrons hadronpart;
				hadronpart.evtN = evt;
				hadronpart.PID = particleid;
				hadronpart.px = momentas[0];
				hadronpart.py = momentas[1];
				hadronpart.pz = momentas[2];
				hadronpart.m = pmass;

				float energy = TMath::Sqrt(pow(momentas[0],2) + pow(momentas[1],2) + pow(momentas[2],2) + pow(pmass,2));
				TLorentzVector ev(momentas[0], momentas[1], momentas[2], energy);

				hadronpart.eta = ev.Eta();
				hadronpart.phi = ev.Phi();
				hadronpart.pT = ev.Pt();

				finalstatehadrons.push_back(hadronpart);
			}

			calculateFlow(v2plot,v3plot);

			finalstatehadrons.clear();
		}
	}

	c->cd();
	v2plot->SetLineWidth(2);
	v2plot->SetLineColor(kBlack);
	v2plot->Draw();

	c1->cd();
	v3plot->SetLineWidth(2);
	v3plot->SetLineColor(kBlack);
	v3plot->Draw();
	TFile *f = new TFile("out_partonParticipantPlaneCalculation.root","RECREATE");
	v2plot->Write();
	v3plot->Write();
	f->Close();

	for (unsigned int i = 0; i < psi2values.size(); i++) {

		cout << "Psi 2 for this event is: " << psi2values[i] << endl;
		cout << "Psi 3 for this event is: " << psi3values[i] << endl;
	}

	for (unsigned int i = 0; i < epsilon2values.size(); i++) {

		TotalE2 = TotalE2 + epsilon2values[i];
		TotalE3 = TotalE3 + epsilon3values[i];
	}

	int totalevents = epsilon2values.size();

	AverageE2 = TotalE2/float(totalevents);
	AverageE3 = TotalE3/float(totalevents);

	cout << "Average Epsilon 2: " << AverageE2 << endl;
	cout << "Average Epsilon 3: " << AverageE3 << endl;

	myAmptFile.close();

	cout << "Finished running analysis code..." << endl;
	return;
}
