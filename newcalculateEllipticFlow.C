//-----------------------------------------------------------------------------------------
//This code is designed to read in files from AMPT to calculate the elliptic flow.
//This version uses partons for the calculation of v2 instead of final state particles.
//
//File 1: npart-xy.dat
//File 2: parton-after-coalescence.dat
//
//07-05-2018
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
// Parton initial location (from npart-xy.dat)
	float x;
	float y;
//Parton sequence and status
	int seq;
	int status;

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

// This is the counter to keep track of which event is being read.
int eventnumber = 0;

// Vectors for the projectile and the target.
vector<parton> particleproj;
vector<parton> particletarg;
// This vector will contain both the target and projectile in one vector.
vector<parton> allparticles;
vector<float> participantplane;
vector<float> eccentricity;

// This vector holds the information needed for the v2 calculations.
vector<fpartons> finalstate;

float psi2;
float phi;
float pT;
float smearpsi2;
float AverageE2;
float TotalE2;

//----------------------------------
//Functions
//----------------------------------

// This function does the calculation for the Participant plane to be used in the calculation of v2.
void calculateParticipantPlane() {

	// This is a counter to enable getting an average.
	int nparton = 0;
	// Variables used in the calculation of psi2.
	float q2x = 0;
	float q2y = 0;
	float e21 = 0;
	float e22 = 0;
	float e2 = 0;
	float rsqure = 0;

	TF1 *fradial = new TF1("fradial","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradial->SetParameter(0,0.4);

	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());

	vector<parton> v = allparticles;

	 for (unsigned int i = 0; i < v.size(); i++) {

	 	if (v[i].status > 0) {

	 		float phi1 = TMath::ATan2(v[i].y, v[i].x);
	 		float r = TMath::Sqrt(pow(v[i].x,2) + pow(v[i].y,2));

	 		float rtemp = fradial->GetRandom();
	 		float phitemp = fphi->GetRandom();

	 		float qxtemp = v[i].x + rtemp*TMath::Sin(phitemp);
	 		float qytemp = v[i].y + rtemp*TMath::Cos(phitemp);

	 		float sr = TMath::Sqrt(pow(qxtemp,2) + pow(qytemp,2));
	 		float sphi = TMath::ATan2(qytemp, qxtemp);

	 		q2x = q2x + pow(sr,2)*TMath::Cos(2*sphi);
	 		q2y = q2y + pow(sr,2)*TMath::Sin(2*sphi);

	 		// For eccentricity calculations

	 		e21 = e21 + sr*sr*TMath::Cos(2*sphi);
	 		e22 = e22 + sr*sr*TMath::Sin(2*sphi);
	 		rsqure = rsqure + pow(sr,2);

	 		nparton++;
	 	}
	 }

	 //cout << nparton << endl;

	 q2x = q2x/nparton;
	 q2y = q2y/nparton;
	 e21 = e21/nparton;
	 e22 = e22/nparton;
     rsqure = rsqure/nparton;

	 psi2 = (TMath::ATan2(e22, e21)/2.0) + (TMath::Pi()/2.0);

	 //cout << ":::::::" << psi2 << endl;
	 e2 = (TMath::Sqrt(e21*e21 + e22*e22)/rsqure); 

	 participantplane.push_back(psi2);
	 eccentricity.push_back(e2);

	//cout << participantplane.size() << endl;

}
// This function does the actual v2 calculation.
void calculateFlow(int &eventcounter, TProfile *v2plot) {

	float v2;
	// This value used to set the participant plane by event.
	int increment = eventcounter;

	// This loop goes through all the partons for an event.
	for (unsigned int j = 0; j < finalstate.size(); j++) {

		// This is the rapidity cut.
		if (finalstate[j].eta < 1 && finalstate[j].eta > -1) {

		  //cout << increment << "     " << participantplane[increment-1] << endl;

			v2 = TMath::Cos(2*(finalstate[j].phi - participantplane[increment-1]));
			// Filling the TProfile.
			v2plot->Fill(finalstate[j].pT,v2);

		}
	}
}
// This section of the code reads in the files and fills the first needed vectors. 
void newcalculateEllipticFlow(void) {

  cout << "Running calculateEllipticFlow()..." << endl;
	// Information needed to read the event header.
	string line;
	int iterationN;
	int atomicMproj;
	int atomicMtarg;
	int impactparm;
	float space[2];
	int sequence;
	int stat;
	float junk[3];
	ifstream myInitialFileInfo;
	ifstream myPacFile;

	myInitialFileInfo.open("ana/npart-xy.dat");

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

		if (heading >> eventnumber >> iterationN >> atomicMproj >> atomicMtarg >> impactparm) {

			if (iterationN > 0) {
				continue;
			}

			int partoncounter = 0;

			// This loop fills in the initial information for each parton.
			for (int i = 0; i < (atomicMproj + atomicMtarg); i++) {

				myInitialFileInfo >> space[0] >> space[1] >> sequence >> stat >> junk[0] >> junk[1] >> junk[2];

				parton partinfo;
				partinfo.evtN = eventnumber;
				partinfo.x = space[0];
				partinfo.y = space[1];
				partinfo.seq = sequence;
				partinfo.status = stat;

				allparticles.push_back(partinfo);

				if (sequence > 0) {
					particleproj.push_back(partinfo);
				}
				else if (sequence < 0) {
					particletarg.push_back(partinfo);
				}

				if (stat > 0 && sequence > 0) {
					partoncounter = partoncounter;
					partoncounter++;
				}
			}

			// This is where I will call the function to calculate the participant plane.
			calculateParticipantPlane();

			// Clearing vectors at the end of an event.
			allparticles.clear();
			particleproj.clear();
			particletarg.clear();
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

				calculateFlow(eventcounter,v2plot);

				finalstate.clear();
				counterparton = 0;
			}

		}
	}

	v2plot->SetLineWidth(2);
	v2plot->SetLineColor(kBlack);
	v2plot->Draw();
	TFile *f = new TFile("out_newcalculateEllipticFlow.root","RECREATE");
	v2plot->Write();
	f->Close();

	for (unsigned int i = 0; i < eccentricity.size(); i++){
	  
	  TotalE2 = TotalE2 + eccentricity[i];
	}

	int totalevents = eccentricity.size();

	AverageE2 = TotalE2/float(totalevents);
	cout << "Average Epsilon 2: " << AverageE2 << endl;

	myPacFile.close();

	cout << "Finished running analysis code..." << endl;
	return;
}

