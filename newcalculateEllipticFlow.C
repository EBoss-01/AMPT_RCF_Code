//-----------------------------------------------------------------------------------------
//This code is designed to read in files from AMPT to calculate the elliptic flow.
//This version uses partons for the calculation of v2 instead of final state particles.
//
//File 1: parton-initial-afterPropagation.dat
//File 2: parton-after-coalescence.dat
//
//07-05-2018
// Updated: 07-12-18
//------------------------------------------------------------------------------------------

#include "TLatex.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
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

//---------------------------------
//Structure for the Parton
//---------------------------------

struct parton {
// Event Number
	int evtN;
// Parton ID
	int PID;
// Mometa at birth
	float px;
	float py;
	float pz;
// Mass of parton
	float m;
// Parton initial location
	float x;
	float y;
	float z;
	float t;
// Other
	int eta;
	int pT;

};

// This is for the final state partons.
struct fpartons {
// Event Number
	int evtnumber;
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

//----------------------------------
// Variables
//----------------------------------

// Event Counter
int eventnumber = 0;
// Vectors
vector<parton> initialpartons;
vector<float> psi2values;
vector<float> epsilon2values;
vector<float> psi3values;
vector<float> epsilon3values;
// This vector holds the information needed for the v2 calculations.
vector<fpartons> finalstate;

float psi2;
float psi3;
float AverageE2;
float AverageE3;
float TotalE2;
float TotalE3;

//----------------------------------
//Functions
//----------------------------------

// This function does the calculation for the Participant plane to be used in the calculation of v2.
void calculateParticipantPlane() {

	// Variables for calculations.
	float cmx = 0;
	float cmy = 0;
	float q2x = 0;
	float q2y = 0;
	float q3x = 0;
	float q3y = 0;
	float epsilon2 = 0;
	float epsilon3 = 0;
	float rsquare = 0;

	// Vector used for CM adjustment
	vector<parton> v = initialpartons;

	for (unsigned int k = 0; k < v.size(); k++) {
		cmx = cmx + v[k].x;
		cmy = cmy + v[k].y;
	}

	int nparton = v.size();

	cmx = cmx/(float)nparton;
	cmy = cmy/(float)nparton;

	for (unsigned int i = 0; i < v.size(); i++) {

		v[i].x = v[i].x - cmx;
		v[i].y = v[i].y - cmy;

	 	float phivalue = TMath::ATan2(v[i].y, v[i].x);
	 	float rvalue = TMath::Sqrt(pow(v[i].x,2) + pow(v[i].y,2));
	 	//Compnent calculations for psi.
	 	q2x = q2x + pow(rvalue,2)*TMath::Cos(2*phivalue);
	 	q2y = q2y + pow(rvalue,2)*TMath::Sin(2*phivalue);
		q3x = q3x + pow(rvalue,2)*TMath::Cos(3*phivalue);
		q3y = q3y + pow(rvalue,2)*TMath::Sin(3*phivalue);
		rsquare = rsquare + pow(rvalue,2);
	}

	// Averaging components
	q2x = q2x/(float)nparton;
	q2y = q2y/(float)nparton;
	q3x = q3x/(float)nparton;
	q3y = q3y/(float)nparton;
    rsquare = rsquare/(float)nparton;

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
// This function does the actual v2 calculation.
void calculateFlow(int &eventcounter, TProfile *v2plot, TProfile *v3plot) {

	float v2 = 0;
	float v3 = 0;
	// This value used to set the participant plane by event.
	int increment = eventcounter;

	// This loop goes through all the partons for an event.
	for (unsigned int j = 0; j < finalstate.size(); j++) {

		if (finalstate[j].pT < 0.0001) {
			continue;
		}

		// This is the rapidity cut.
		if (fabs(finalstate[j].eta) < 2) {

			v2 = TMath::Cos(2*(finalstate[j].phi - psi2values[increment-1]));
			// Filling the TProfile.
			v2plot->Fill(finalstate[j].pT,v2);

			v3 = TMath::Cos(3*(finalstate[j].phi - psi3values[increment-1]));
			// Filling v3 TProfile.
			v3plot->Fill(finalstate[j].pT,v3);

		}
	}
}
// This section of the code reads in the files and fills the first needed vectors. 
void newcalculateEllipticFlow(void) {

  cout << "Running calculateEllipticFlow()..." << endl;
	// Information needed to read the event header.
	string line;
	int iterationN;
	int partonsN;
	int baryonsN;
	int mesonsN;
	int particlesN;
	int nonpartparticles;
	// Needed for lines after header
	int pid;
	float momentas[3];
	float pmass;
	float spacetime[4];
	// File initialization
	ifstream myInitialFileInfo;
	ifstream myPacFile;

	myInitialFileInfo.open("ana/parton-initial-afterPropagation.dat");

	if (!myInitialFileInfo) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file 1." << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened file 1." << endl;
	}


	while (std::getline(myInitialFileInfo, line)) {

		// This line prevents the file from trying to read the last line multiple times.
		if (!myInitialFileInfo) break;

		stringstream heading(line);

		if (heading >> eventnumber >> iterationN >> partonsN >> baryonsN >> mesonsN >> particlesN >> nonpartparticles) {

			if (iterationN > 0) {
				continue;
			}

			// This loop fills in the initial information for each parton.
			for (int i = 0; i < partonsN; i++) {

				myInitialFileInfo >> pid >> momentas[0] >> momentas[1] >> momentas[2] >> pmass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

				float energy1 = TMath::Sqrt(pow(momentas[0],2) + pow(momentas[1],2) + pow(momentas[2],2) + pow(pmass,2));
				TLorentzVector cut(momentas[0], momentas[1], momentas[2], energy1);

				parton partinitial;
				partinitial.evtN = eventnumber;
				partinitial.PID = pid;
				partinitial.px = momentas[0];
				partinitial.py = momentas[1];
				partinitial.pz = momentas[2];
				partinitial.m = pmass;
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

			// This is where I will call the function to calculate the participant plane.
			calculateParticipantPlane();

			// Clearing vectors at the end of an event.
			initialpartons.clear();
		}

	}

	myInitialFileInfo.close();

	string line1;
	int evt;
	int npartons;
	int nbaryons;
	int nmesons;
	float imparm;
	int stuff[5];

	// Needed for the fpartons structure.
	int partonid;
	float momenta[3];
	float mass;
	int seqnumber;
	int hadronN;
	float extra;

	myPacFile.open("ana/parton-after-coalescence.dat");

	if (!myPacFile) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file 2." << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened file 2." << endl;
	}

	// Counters used for events and the amount of partons.
	int eventcounter = 0;
	int participantparton = 0;
	int counterparton = 0;

	// Creation of the TCanvas and TProfile.
	TCanvas *c = new TCanvas("c","v2 calculations",700,700);
	gStyle->SetOptStat(0);
	TProfile *v2plot = new TProfile("v2plot","v_{2} vs p_{T}",20,0,2.5,-1,1);

	TCanvas *c1 = new TCanvas("c1","v3 calculations",700,700);
	gStyle->SetOptStat(0);
	TProfile *v3plot = new TProfile("v3plot","v_{3} vs p_{T}",20,0,2.5,-1,1);

	// Reading in file 2.
	while (std::getline(myPacFile, line1)) {

		if (!myPacFile) break;

		// These strings are used to tokenize each line read in.
		stringstream header(line1);
		string intermediate;

		// Vector that will give the size of each line.
		vector<float> token;

		// This piece finds each value and puts it into the tokenized vector.
		while (getline(header, intermediate, ' ')) {

			if (intermediate != ""){
				token.push_back(atof(intermediate.c_str()));
			}
		}

		// Getting ready to fill fparton structure.
		fpartons finalpart;
		float energy = 0;

		// Used to recognize the headers per event.
		if (token.size() == 10) {
			eventcounter++;
			participantparton = token[1];

			//cout << eventcounter << "  " << participantparton << endl;
		}

		// Used for all lines with 8 or 7 elements.
		else if ((token.size() == 8) || (token.size() == 7)) {

			// Filling the fparton structure.
			finalpart.evtnumber = eventcounter;
			finalpart.PID = token[0];
			finalpart.px = token[1];
			finalpart.py = token[2];
			finalpart.pz = token[3];
			finalpart.m = token[4];

			// Energy and Lorentz vectors.
			energy = TMath::Sqrt(pow(token[1],2) + pow(token[2],2) + pow(token[3],2) + pow(token[4],2));
			TLorentzVector ev(token[1], token[2], token[3], energy);

			// Filling final pieces using TLorentz operations.
			finalpart.eta = ev.Eta();
			finalpart.phi = ev.Phi();
			finalpart.pT = ev.Pt();

			// Pushing structure inforamtion to global vector.
			finalstate.push_back(finalpart);
			counterparton++;

			//cout << counterparton << endl;

			if (counterparton == participantparton) {

			  calculateFlow(eventcounter,v2plot,v3plot);

				finalstate.clear();
				counterparton = 0;
			}

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
	TFile *f = new TFile("out_newcalculateEllipticFlow.root","RECREATE");
	v2plot->Write();
	v3plot->Write();
	f->Close();

	for (unsigned int i = 0; i < epsilon2values.size(); i++){
	  
	  TotalE2 = TotalE2 + epsilon2values[i];
	  TotalE3 = TotalE3 + epsilon3values[i];
	}

	int totalevents = epsilon2values.size();

	AverageE2 = TotalE2/float(totalevents);
	AverageE3 = TotalE3/float(totalevents);

	cout << "Average Epsilon 2: " << AverageE2 << endl;
	cout << "Average Epsilon 3: " << AverageE3 << endl;

	myPacFile.close();

	cout << "Finished running analysis code..." << endl;
	return;
}

