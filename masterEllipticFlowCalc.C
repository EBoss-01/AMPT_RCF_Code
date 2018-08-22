//--------------------------------------------------------------------------------------------
// This code uses multiple AMPT files to calulate the elliptic flow.
// There is a flag that is used to determine how the participant plane is calculatied.
//
// File 1: parton-initial-afterPropagation.dat (uses partons for Psi)
// File 2: npart-xy.dat (uses nucleons for Psi for set tform = const)
// File 3: parton-after-coalescence.dat (uses partons for v2 calculation)
// File 4: ampt.dat (uses hadrons for v2 calculation)
// File 5: parton-collisionsHistory.dat (used to trace partons to formation time)
//
// Flagpsi is used to determine the files used for the plane calculation.
// Flagpsi = 0, uses partons and tform != const.
// Flagpsi = 1, uses nucleons and tform = const.
//
// 07-16-2018
// Updated 08-22-18 (created a new formation time distribution plot)
//--------------------------------------------------------------------------------------------

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
#include <map>

using namespace std;

//-------------------
// Structures
//-------------------

// This is for File 1.
struct parton {
// Event Number
	int evtN;
// Parton ID
	int PID;
// Parton Line Position
	int position;
// Momenta at birth
	float px;
	float py;
	float pz;
// Mass of parton
	float m;
// Initial location
	float x;
	float y;
	float z;
	float t;
// Additional
	float eta;
	float pT;
};

// This is for File 2.
struct nucleon {
// Event number
	int evtN;
// Nucleon location
	float x;
	float y;
// Sequence and status
	int seq;
	int status;
};

// This is for File 3.
struct fpartons {
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
// Additional
	float eta;
	float phi;
	float pT;
};

// This is for File 4.
struct hadrons {
// Event Number
	int evtN;
// Particle ID
	int PID;
// Momenta
	float px;
	float py;
	float pz;
// Mass
	float m;
// location
	float x;
	float y;
	float z;
	float t;
// Additional
	float eta;
	float phi;
	float pT;
};

// For parseCollisions
struct collision {
// Event number
	int evtN;
// Parton ID
	int PID;
// Parton number
	int nline;
// Momenta
	float px;
	float py;
	float pz;
// Mass
	float m;
// Spacetime
	float x;
	float y;
	float z;
	float t;
};

// For findFormationTime
struct newparton {
// Event number
	int evtN;
// Parton ID
	int PID;
// New Parton Number
	int Nparton;
// Momenta
	float px;
	float py;
	float pz;
// Mass
	float m;
// Additional
	float eta;
	float phi;
	float pT;
// Formation Time
	float t;
};

//-------------------------------------
// Variables
//-------------------------------------
int setFlag = 0;
// Event Counter
int eventnumber = 0;
// Vectors
vector<parton> initialpartons;
vector<nucleon> initialnucleon;
// New vectors for extracting formation times.
vector<collision> partoncollisions;
vector<newparton> partonWTform;
vector<int> partoncounts;
vector<parton> allinitialpartons;

vector<float> psi2valuespartons;
vector<float> psi2valuesnucleon;
vector<float> epsilon2valuespartons;
vector<float> epsilon2valuesnucleon;
vector<float> psi3valuespartons;
vector<float> psi3valuesnucleon;
vector<float> epsilon3valuespartons;
vector<float> epsilon3valuesnucleon;
// This vector holds the information needed for the v2 calculations.
vector<fpartons> finalstate;
vector<hadrons> finalhadrons;

float psi2partons;
float psi2nucleon;
float psi3partons;
float psi3nucleon;
float AverageE2partons;
float AverageE2nucleon;
float AverageE3partons;
float AverageE3nucleon;
float TotalE2partons;
float TotalE2nucleon;
float TotalE3partons;
float TotalE3nucleon;
// Variables needed for multimap
int Pevent = 0;
float Ppx = 0;
float Ppy = 0;
float Ppz = 0;
string pkey;
typedef multimap<string, int> PartonMap;

//-------------------------------------
// Functions
//-------------------------------------

void calculateParticipantPlanePartons() {

	int cmcounter = 0;
	float cmx = 0;
	float cmy = 0;
	float q2x = 0;
	float q2y = 0;
	float q3x = 0;
	float q3y = 0;
	float epsilon2 = 0;
	float epsilon3 = 0;
	float rsquare = 0;

	vector<parton> v = initialpartons;

	for (unsigned int k = 0; k < v.size(); k++) {
		cmx = cmx + v[k].x;
		cmy = cmy + v[k].y;
	}

	int ncounter = v.size();

	cmx = cmx/(float)ncounter;
	cmy = cmy/(float)ncounter;

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
	q2x = q2x/(float)ncounter;
	q2y = q2y/(float)ncounter;
	q3x = q3x/(float)ncounter;
	q3y = q3y/(float)ncounter;
    rsquare = rsquare/(float)ncounter;

    // Participant Plane calculations
	psi2partons = (TMath::ATan2(q2y, q2x)/2.0) + (TMath::Pi()/2.0);
	psi3partons = (TMath::ATan2(q3y, q3x)/3.0) + (TMath::Pi()/3.0);
	// Eccentricity calculations
	epsilon2 = (TMath::Sqrt(pow(q2x,2) + pow(q2y,2)))/rsquare;
	epsilon3 = (TMath::Sqrt(pow(q3x,2) + pow(q3y,2)))/rsquare;
	// Pushing the values to more vectors
	psi2valuespartons.push_back(psi2partons);
	psi3valuespartons.push_back(psi3partons);
	epsilon2valuespartons.push_back(epsilon2);
	epsilon3valuespartons.push_back(epsilon3);
}

void calculateParticipantPlaneNucleons() {

	int nucleoncounter = 0;
	int Ncmcounter = 0;
	float Ncmx = 0;
	float Ncmy = 0;
	float Nq2x = 0;
	float Nq2y = 0;
	float Nq3x = 0;
	float Nq3y = 0;
	float Nepsilon2 = 0;
	float Nepsilon3 = 0;
	float Nrsquare = 0;

	TF1 *fradial = new TF1("fradial","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradial->SetParameter(0,0.4);

	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());

	vector<nucleon> vN = initialnucleon;

	for (unsigned int k = 0; k < vN.size(); k++) {
		if (vN[k].status > 0) {
			Ncmx = Ncmx + vN[k].x;
			Ncmy = Ncmy + vN[k].y;

			Ncmcounter++;
		}
	}

	Ncmx = Ncmx/(float)Ncmcounter;
	Ncmy = Ncmy/(float)Ncmcounter;

	for (unsigned int i = 0; i < vN.size(); i++) {

		if (vN[i].status > 0) {

	 		vN[i].x = vN[i].x - Ncmx;
	 		vN[i].y = vN[i].y - Ncmy;

	 		float rtemp = fradial->GetRandom();
	 		float phitemp = fphi->GetRandom();

	 		float qxtemp = vN[i].x + rtemp*TMath::Sin(phitemp);
	 		float qytemp = vN[i].y + rtemp*TMath::Cos(phitemp);

	 		float sr = TMath::Sqrt(pow(qxtemp,2) + pow(qytemp,2));
	 		float sphi = TMath::ATan2(qytemp, qxtemp);

	 		Nq2x = Nq2x + pow(sr,2)*TMath::Cos(2*sphi);
	 		Nq2y = Nq2y + pow(sr,2)*TMath::Sin(2*sphi);
	 		// psi3 components with smear.
	 		Nq3x = Nq3x + pow(sr,2)*TMath::Cos(3*sphi);
	 		Nq3y = Nq3y + pow(sr,2)*TMath::Sin(3*sphi);
	 		Nrsquare = Nrsquare + pow(sr,2);

	 		nucleoncounter++;
	 	}
	}

	Nq2x = Nq2x/(float)nucleoncounter;
	Nq2y = Nq2y/(float)nucleoncounter;
	Nq3x = Nq3x/(float)nucleoncounter;
	Nq3y = Nq3y/(float)nucleoncounter;
	Nrsquare = Nrsquare/(float)nucleoncounter;

	psi2nucleon = (TMath::ATan2(Nq2y, Nq2x)/2.0) + (TMath::Pi()/2.0);
	psi3nucleon = (TMath::ATan2(Nq3y, Nq3x)/3.0) + (TMath::Pi()/3.0);

	Nepsilon2 = (TMath::Sqrt(pow(Nq2x,2) + pow(Nq2y,2)))/Nrsquare;
	Nepsilon3 = (TMath::Sqrt(pow(Nq3x,2) + pow(Nq3y,2)))/Nrsquare;

	psi2valuesnucleon.push_back(psi2nucleon);
	psi3valuesnucleon.push_back(psi3nucleon);
	epsilon2valuesnucleon.push_back(Nepsilon2);
	epsilon3valuesnucleon.push_back(Nepsilon3);

}

void findFormationTime(PartonMap &keymomentums, string &pkey, TLorentzVector ev) {

	newparton Update;

	TLorentzVector v = ev;

	PartonMap::iterator it = keymomentums.find(pkey);

	int initialpartonvalue = keymomentums.find(pkey)->second;

	for (unsigned int i = 0; i < allinitialpartons.size(); i++) {

		if (allinitialpartons[i].evtN == Pevent) {

			if (allinitialpartons[i].position == initialpartonvalue) {

				Update.evtN = allinitialpartons[i].evtN;
				Update.PID = allinitialpartons[i].PID;
				Update.Nparton = allinitialpartons[i].position;
				Update.px = Ppx;
				Update.py = Ppy;
				Update.pz = Ppz;
				Update.m = allinitialpartons[i].m;
				Update.eta = v.Eta();
				Update.phi = v.Phi();
				Update.pT = v.Pt();
				Update.t = allinitialpartons[i].t;

				partonWTform.push_back(Update);
			}
		}
	}
}

void calculateFlowPartons(TProfile *v2plotPpartons, TProfile *v2plotNpartons, TProfile *v3plotPpartons, TProfile *v3plotNpartons, TProfile *v2plotPpartonseta1, TProfile *v2plotPpartonsetap5, TProfile *v2plotNpartonseta1, TProfile *v2plotNpartonsetap5) {

	float v2Ppartons = 0;
	float v2Npartons = 0;
	float v3Ppartons = 0;
	float v3Npartons = 0;
	int increment = 0;

	// Also don't use the flags in this section to ensure both calculations happen.

	for (unsigned int j = 0; j < finalstate.size(); j++) {

		if (finalstate[j].pT < 0.0001) {
			continue;
		}

		if (fabs(finalstate[j].eta) < 2) {

			increment = finalstate[j].evtN;
			// psi using partons
			v2Ppartons = TMath::Cos(2*(finalstate[j].phi - psi2valuespartons[increment-1]));
			v2plotPpartons->Fill(finalstate[j].pT,v2Ppartons);

			v3Ppartons = TMath::Cos(3*(finalstate[j].phi - psi3valuespartons[increment-1]));
			v3plotPpartons->Fill(finalstate[j].pT,v3Ppartons);

			// psi using nucleons
			v2Npartons = TMath::Cos(2*(finalstate[j].phi - psi2valuesnucleon[increment-1]));
			v2plotNpartons->Fill(finalstate[j].pT,v2Npartons);

			v3Npartons = TMath::Cos(3*(finalstate[j].phi - psi3valuesnucleon[increment-1]));
			v3plotNpartons->Fill(finalstate[j].pT,v3Npartons);
		}

		if (fabs(finalstate[j].eta) < 1) {

			increment = finalstate[j].evtN;

			// psi using partons
			v2Ppartons = TMath::Cos(2*(finalstate[j].phi - psi2valuespartons[increment-1]));
			v2plotPpartonseta1->Fill(finalstate[j].pT,v2Ppartons);

			// psi using nucleons
			v2Npartons = TMath::Cos(2*(finalstate[j].phi - psi2valuesnucleon[increment-1]));
			v2plotNpartonseta1->Fill(finalstate[j].pT,v2Npartons);

		}

		if (fabs(finalstate[j].eta) < 0.5) {

			increment = finalstate[j].evtN;

			// psi using partons
			v2Ppartons = TMath::Cos(2*(finalstate[j].phi - psi2valuespartons[increment-1]));
			v2plotPpartonsetap5->Fill(finalstate[j].pT,v2Ppartons);

			// psi using nucleons
			v2Npartons = TMath::Cos(2*(finalstate[j].phi - psi2valuesnucleon[increment-1]));
			v2plotNpartonsetap5->Fill(finalstate[j].pT,v2Npartons);

		}

	}

}

void calculateFlowHadrons(TProfile *v2plotPhadrons, TProfile *v2plotNhadrons, TProfile *v3plotPhadrons, TProfile *v3plotNhadrons) {

	float v2Phadrons = 0;
	float v2Nhadrons = 0;
	float v3Phadrons = 0;
	float v3Nhadrons = 0;

	// Partons
	for (unsigned int i = 0; i < psi2valuespartons.size(); i++) {

		for (unsigned int j = 0; j < finalhadrons.size(); j++) {

			if (finalhadrons[j].evtN == (i+1)) {

				if (fabs(finalhadrons[j].eta) < 0.35) {

					v2Phadrons = TMath::Cos(2*(finalhadrons[j].phi - psi2valuespartons[i]));
					v2plotPhadrons->Fill(finalhadrons[j].pT,v2Phadrons);

					v3Phadrons = TMath::Cos(3*(finalhadrons[j].phi - psi3valuespartons[i]));
					v3plotPhadrons->Fill(finalhadrons[j].pT,v3Phadrons);
				}
			}
		}
	}
	// nucleons
	for (unsigned int i = 0; i < psi2valuesnucleon.size(); i++) {

		for (unsigned int j = 0; j < finalhadrons.size(); j++) {

			if (finalhadrons[j].evtN == (i+1)) {

				if (fabs(finalhadrons[j].eta) < 0.35) {

					v2Nhadrons = TMath::Cos(2*(finalhadrons[j].phi - psi2valuesnucleon[i]));
					v2plotNhadrons->Fill(finalhadrons[j].pT,v2Nhadrons);

					v3Nhadrons = TMath::Cos(3*(finalhadrons[j].phi - psi3valuesnucleon[i]));
					v3plotNhadrons->Fill(finalhadrons[j].pT,v3Nhadrons);
				}
			}
		}
	}
}

void parsePartons(PartonMap &keymomentums) {

	// Needed to read header
	string line;
	int iterationN;
	int partonsN;
	int baryonsN;
	int mesonsN;
	int particlesN;
	int nonpartparticles;
	// Needed for additional lines
	int pid;
	float momentas[3];
	float pmass;
	float spacetime[4];
	// File initialization
	ifstream myFileOne;
	myFileOne.open("ana/parton-initial-afterPropagation.dat");

	if (!myFileOne) {
		cout << "Unable to open parton-initial-afterPropagation." << endl;
		return;
	}
	else {
		cout << "Successfully opened parton-initial-afterPropagation." << endl;
	}

	while (std::getline(myFileOne, line)) {
		// Prevents duplication of last line
		if (!myFileOne) break;

		stringstream heading(line);

		if (heading >> eventnumber >> iterationN >> partonsN >> baryonsN >> mesonsN >> particlesN >> nonpartparticles) {

			if (iterationN > 0) {
				continue;
			}

			// This counter is used to determine actual parton number.
			int linecounter = 0;

			// This loop files the initial information
			for (int i = 0; i < partonsN; i++) {

				myFileOne >> pid >> momentas[0] >> momentas[1] >> momentas[2] >> pmass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

				linecounter++;

				float energy1 = TMath::Sqrt(pow(momentas[0],2) + pow(momentas[1],2) + pow(momentas[2],2) + pow(pmass,2));
				TLorentzVector cut(momentas[0], momentas[1], momentas[2], energy1);

				parton partinitial;
				partinitial.evtN = eventnumber;
				partinitial.PID = pid;
				partinitial.position = linecounter;
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

				// Setting the multimap variables.
				Pevent = eventnumber;
				Ppx = momentas[0];
				Ppy = momentas[1];
				Ppz = momentas[2];

				pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);
				keymomentums.insert(make_pair(pkey, linecounter));

				allinitialpartons.push_back(partinitial);

				if (fabs(partinitial.eta) < 3 && partinitial.t < 3) {
					initialpartons.push_back(partinitial);
				}
			}

			calculateParticipantPlanePartons();

			initialpartons.clear();
		}
	}

	myFileOne.close();
}

void parseCollisions (PartonMap &keymomentums) {
	string line4;
	int Pevt;
	int extrajunk1;
	int partonindex1;
	int partonindex2;
	// Initial parton information
	int parton1_id_initial;
	float parton1_momenta_initial[3];
	float parton1_mass_initial;
	double parton1_spacetime_initial[4];
	int parton2_id_initial;
	float parton2_momenta_initial[3];
	float parton2_mass_initial;
	double parton2_spacetime_initial[4];
	// Final parton information
	int parton1_final_id;
	float parton1_final_momenta[3];
	float parton1_final_mass;
	double parton1_final_spacetime[4];
	int parton2_final_id;
	float parton2_final_momenta[3];
	float parton2_final_mass;
	double parton2_final_spacetime[4];

	ifstream myFileFive;

	myFileFive.open("ana/parton-collisionsHistory.dat");
	if (!myFileFive) {
		cout << "Unable to open parton-collisionsHistory." << endl;
		return;
	}
	else {
		cout << "Successfully opened parton-collisionsHistory." << endl;
	}

	while (std::getline(myFileFive, line4)) {

		if (!myFileFive) break;

		stringstream headline(line4);
		string description;

		if (headline >> description >> Pevt >> extrajunk1 >> partonindex1 >> partonindex2) {

			myFileFive >> parton1_id_initial >> parton1_momenta_initial[0] >> parton1_momenta_initial[1] >> parton1_momenta_initial[2] >> parton1_mass_initial >> parton1_spacetime_initial[0] >> parton1_spacetime_initial[1] >> parton1_spacetime_initial[2] >> parton1_spacetime_initial[3];

			myFileFive >> parton2_id_initial >> parton2_momenta_initial[0] >> parton2_momenta_initial[1] >> parton2_momenta_initial[2] >> parton2_mass_initial >> parton2_spacetime_initial[0] >> parton2_spacetime_initial[1] >> parton2_spacetime_initial[2] >> parton2_spacetime_initial[3];

			myFileFive >> parton1_final_id >> parton1_final_momenta[0] >> parton1_final_momenta[1] >> parton1_final_momenta[2] >> parton1_final_mass >> parton1_final_spacetime[0] >> parton1_final_spacetime[1] >> parton1_final_spacetime[2] >> parton1_final_spacetime[3];

			myFileFive >> parton2_final_id >> parton2_final_momenta[0] >> parton2_final_momenta[1] >> parton2_final_momenta[2] >> parton2_final_mass >> parton2_final_spacetime[0] >> parton2_final_spacetime[1] >> parton2_final_spacetime[2] >> parton2_final_spacetime[3];

			collision part1;
			part1.evtN = Pevt;
			part1.PID = parton1_final_id;
			part1.nline = partonindex1;
			part1.px = parton1_final_momenta[0];
			part1.py = parton1_final_momenta[1];
			part1.pz = parton1_final_momenta[2];
			part1.m = parton1_final_mass;
			part1.x = parton1_final_spacetime[0];
			part1.y = parton1_final_spacetime[1];
			part1.z = parton1_final_spacetime[2];
			part1.t = parton1_final_spacetime[3];

			// Setting the multimap variables again.
			Pevent = Pevt;
			Ppx = parton1_final_momenta[0];
			Ppy = parton1_final_momenta[1];
			Ppz = parton1_final_momenta[2];

			pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);
			keymomentums.insert(make_pair(pkey, partonindex1));

			collision part2;
			part2.evtN = Pevt;
			part2.PID = parton2_final_id;
			part2.nline = partonindex2;
			part2.px = parton2_final_momenta[0];
			part2.py = parton2_final_momenta[1];
			part2.pz = parton2_final_momenta[2];
			part2.m = parton2_final_mass;
			part2.x = parton2_final_spacetime[0];
			part2.y = parton2_final_spacetime[1];
			part2.z = parton2_final_spacetime[2];
			part2.t = parton2_final_spacetime[3];

			// Setting the multimap variables once more.
			Pevent = Pevt;
			Ppx = parton2_final_momenta[0];
			Ppy = parton2_final_momenta[1];
			Ppz = parton2_final_momenta[2];

			pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);
			keymomentums.insert(make_pair(pkey, partonindex2));

			partoncollisions.push_back(part1);
			partoncollisions.push_back(part2);
		}
	}

	cout << "Closing parton-collisionsHistory." << endl;

	myFileFive.close();
}

void parseNucleons() {

	// Variables needed to read header.
	string line1;
	int Niteration;
	int atomicMproj;
	int atomicMtarg;
	int impactparm;
	float space[2];
	int sequence;
	int stat;
	float junk[3];
	ifstream myFileTwo;

	myFileTwo.open("ana/npart-xy.dat");

	if (!myFileTwo) {
		cout << "Unable to open npart-xy." << endl;
		return;
	}
	else {
		cout << "Successfully opened npart-xy." << endl;
	}

	while (std::getline(myFileTwo, line1)) {

		// Prevents duplicating last line.
		if (!myFileTwo) break;

		stringstream header(line1);

		if (header >> eventnumber >> Niteration >> atomicMproj >> atomicMtarg >> impactparm) {

			if (Niteration > 0) {
				continue;
			}

			int partoncounter = 0;

			// Loop for filling
			for (int i = 0; i < (atomicMproj + atomicMtarg); i++) {

				myFileTwo >> space[0] >> space[1] >> sequence >> stat >> junk[0] >> junk[1] >> junk[2];

				nucleon nucinitial;
				nucinitial.evtN = eventnumber;
				nucinitial.x = space[0];
				nucinitial.y = space[1];
				nucinitial.seq = sequence;
				nucinitial.status = stat;

				initialnucleon.push_back(nucinitial);

				if (stat > 0 && sequence >> 0) {
					partoncounter = partoncounter;
					partoncounter++;
				}
			}

			calculateParticipantPlaneNucleons();

			initialnucleon.clear();
		}
	}

	myFileTwo.close();

}

void parseFinalPartons(TProfile *v2plotPpartons,TProfile *v2plotNpartons,TProfile *v3plotPpartons,TProfile *v3plotNpartons, TProfile *v2plotPpartonseta1, TProfile *v2plotPpartonsetap5, TProfile *v2plotNpartonseta1, TProfile *v2plotNpartonsetap5, PartonMap &keymomentums) {

	string line2;
	int evt;
	int npartons;
	int nbaryons;
	int nmesons;
	float imparm;
	int stuff[5];
	// Needed for fpartons struct.
	int partonid;
	float momenta[3];
	float mass;
	int seqnumber;
	int hadronN;
	float extra;

	ifstream myFileThree;

	myFileThree.open("ana/parton-after-coalescence.dat");

	if (!myFileThree) {
		cout << "Unable to open parton-after-coalescence." << endl;
		return;
	}
	else {
		cout << "Successfully opened parton-after-coalescence." << endl;
	}

	int eventcounter = 0;
	int participantparton = 0;
	int counterparton = 0;

	while (std::getline(myFileThree, line2)) {

		if (!myFileThree) break;

		stringstream eventheader(line2);
		string intermediate;

		vector<float> token;

		while (getline(eventheader, intermediate,' ')) {

			if (intermediate != "") {
				token.push_back(atof(intermediate.c_str()));
			}
		}

		fpartons finalpartons;
		float energy2 = 0;

		if (token.size() == 10) {
			eventcounter++;
			participantparton = token[1];
		}

		else if ((token.size() == 8) || (token.size() == 7)) {

			finalpartons.evtN = eventcounter;
			finalpartons.PID = token[0];
			finalpartons.px = token[1];
			finalpartons.py = token[2];
			finalpartons.pz = token[3];
			finalpartons.m = token[4];

			energy2 = TMath::Sqrt(pow(token[1],2) + pow(token[2],2) + pow(token[3],2) + pow(token[4],2));
			TLorentzVector ev(token[1], token[2], token[3], energy2);

			finalpartons.eta = ev.Eta();
			finalpartons.phi = ev.Phi();
			finalpartons.pT = ev.Pt();

			if (ev.Pt() < 0.0001) {
				continue;
			}

			finalstate.push_back(finalpartons);
			counterparton++;

			// Adjusting multimap variables.
			Pevent = eventcounter;
			Ppx = token[1];
			Ppy = token[2];
			Ppz = token[3];

			pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);
			// Call to new function
			findFormationTime(keymomentums, pkey, ev);

			if (counterparton == participantparton) {

				calculateFlowPartons(v2plotPpartons,v2plotNpartons,v3plotPpartons,v3plotNpartons,v2plotPpartonseta1,v2plotPpartonsetap5,v2plotNpartonseta1,v2plotNpartonsetap5);

				finalstate.clear();
				counterparton = 0;
			}
		}
	}

	cout << "Closing parton-after-coalescence." << endl;

	myFileThree.close();

}

void parseFinalHadrons(TProfile *v2plotPhadrons, TProfile *v2plotNhadrons, TProfile *v3plotPhadrons, TProfile *v3plotNhadrons) {

	string line3;
	int evtnumber;
	int test;
	int npartevt;
	float imparms[2];
	int totalpart;
	int participantproj;
	int participanttarg;
	int collisions[4];
	// Needed for struct
	int particleid;
	float momentum[3];
	float Hmass;
	float spaces[4];

	ifstream myFileFour;

	myFileFour.open("ana/ampt.dat");

	if (!myFileFour) {
		cout << "Unable to open ampt.dat." << endl;
		return;
	}
	else {
		cout << "Successfully opened ampt.dat." << endl;
	}

	while (std::getline(myFileFour, line3)) {

		if (!myFileFour) break;

		stringstream hadronheader(line3);

		if (hadronheader >> evtnumber >> test >> npartevt >> imparms[0] >> totalpart >> participantproj >> participanttarg >> collisions[0] >> collisions[1] >> collisions[2] >> collisions[3] >> imparms[1]) {

			for (int j = 0; j < npartevt; j++) {

				myFileFour >> particleid >> momentum[0] >> momentum[1] >> momentum[2] >> Hmass >> spaces[0] >> spaces[1] >> spaces[2] >> spaces[3];

				if (momentum[0] == 0.000 && momentum[1] == 0.000) {
					continue;
				}

				if ((fabs(particleid) != 211) && (fabs(particleid) != 2212) && (fabs(particleid) != 321)) continue;
					
				hadrons finalH;
				finalH.evtN = evtnumber;
				finalH.PID = particleid;
				finalH.px = momentum[0];
				finalH.py = momentum[1];
				finalH.pz = momentum[2];
				finalH.m = Hmass;
				finalH.x = spaces[0];
				finalH.y = spaces[1];
				finalH.z = spaces[2];
				finalH.t = spaces[3];

				float energy3 = TMath::Sqrt(pow(momentum[0],2) + pow(momentum[1],2) + pow(momentum[2],2) + pow(Hmass,2));
				TLorentzVector hv(momentum[0], momentum[1], momentum[2], energy3);

				finalH.eta = hv.Eta();
				finalH.phi = hv.Phi();
				finalH.pT = hv.Pt();

				finalhadrons.push_back(finalH);
			}

			calculateFlowHadrons(v2plotPhadrons,v2plotNhadrons,v3plotPhadrons,v3plotNhadrons);

			finalhadrons.clear();
		}
	}

	//calculateFlowHadrons(v2plotPhadrons,v2plotNhadrons,v3plotPhadrons,v3plotNhadrons);

	cout << "Closing ampt.dat." << endl;

	myFileFour.close();

}

//-------------------------------------
// Main Code
//-------------------------------------
void masterEllipticFlowCalc(int flagpsi = 0) {

	setFlag = flagpsi;

	// Empty container for use with multimap.
	PartonMap keymomentums;

	cout << "Running masterEllipticFlowCalc()..." << endl;

	// Call to File 1.
	parsePartons(keymomentums);

	// Call to File 5.
	parseCollisions(keymomentums);

	// Call to File 2.
	parseNucleons();

	// Creation of the TCanvas and TProfile.
	TCanvas *c = new TCanvas("c","v2 calculations",700,700);
	gStyle->SetOptStat(0);
	TProfile *v2plotPpartons = new TProfile("v2plotPpartons","v_{2} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v2plotPhadrons = new TProfile("v2plotPhadrons","v_{2} vs. p_{T}",20,0,2.5,-1,1);
	TProfile *v2plotNpartons = new TProfile("v2plotNpartons","v_{2} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v2plotNhadrons = new TProfile("v2plotNhadrons","v_{2} vs. p_{T}",20,0,2.5,-1,1);

	TProfile *v2plotPpartonseta1 = new TProfile("v2plotPpartonseta1","v_{2} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v2plotNpartonseta1 = new TProfile("v2plotNpartonseta1","v_{2} vs p_{T}",20,0,2.5,-1,1);

	TProfile *v2plotPpartonsetap5 = new TProfile("v2plotPpartonsetap5","v_{2} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v2plotNpartonsetap5 = new TProfile("v2plotNpartonsetap5","v_{2} vs p_{T}",20,0,2.5,-1,1);

	TCanvas *c1 = new TCanvas("c1","v3 calculations",700,700);
	gStyle->SetOptStat(0);
	TProfile *v3plotPpartons = new TProfile("v3plotPpartons","v_{3} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v3plotPhadrons = new TProfile("v3plotPhadrons","v_{3} vs. p_{T}",20,0,2.5,-1,1);
	TProfile *v3plotNpartons = new TProfile("v3plotNpartons","v_{3} vs p_{T}",20,0,2.5,-1,1);
	TProfile *v3plotNhadrons = new TProfile("v3plotNhadrons","v_{3} vs. p_{T}",20,0,2.5,-1,1);

	TCanvas *c2 = new TCanvas("c2","Plot of nucleon plane",700,700);
	gStyle->SetOptStat(0);
	TH1F* hist1 = new TH1F("hist1",";#delta#psi_{2};Counts",200,-4,4);
		hist1->Sumw2();
	TH1F* formationdist1 = new TH1F("formationdist1",";t_{form};Counts",4000,0,4);
		formationdist1->Sumw2();
	TH1F* formationdist2 = new TH1F("formationdist2",";t_{form};Counts",4000,0,4);
		formationdist2->Sumw2();
	TH1F* formationdist3 = new TH1F("formationdist3",";t_{form};Counts",4000,0,4);
		formationdist3->Sumw2();
	TH1F* formationdist4 = new TH1F("formationdist4",";t_{form};Counts",4000,0,4);
		formationdist4->Sumw2();
	TH1F* fullformationdist = new TH1F("fullformationdist",";t_{form};Counts",4000,0,4);
		fullformationdist->Sumw2();

	// Call to file 3.
	parseFinalPartons(v2plotPpartons,v2plotNpartons,v3plotPpartons,v3plotNpartons,v2plotPpartonseta1,v2plotPpartonsetap5,v2plotNpartonseta1,v2plotNpartonsetap5,keymomentums);

	for (unsigned int i = 0; i < partonWTform.size(); i++) {

		if ((partonWTform[i].pT > 0.2) && (partonWTform[i].pT <0.4)) {
			formationdist1->Fill(partonWTform[i].t);
		}
		if ((partonWTform[i].pT > 0.4) && (partonWTform[i].pT < 0.6)) {
			formationdist2->Fill(partonWTform[i].t);
		}
		if ((partonWTform[i].pT > 0.6) && (partonWTform[i].pT < 0.8)) {
			formationdist3->Fill(partonWTform[i].t);
		}
		if ((partonWTform[i].pT > 0.8) && (partonWTform[i].pT < 1.0)) {
			formationdist4->Fill(partonWTform[i].t);
		}
	}

	for (unsigned int j = 0; j < partonWTform.size(); j++) {

		fullformationdist->Fill(partonWTform[j].t);
	}

	// Call to File 4.
	parseFinalHadrons(v2plotPhadrons,v2plotNhadrons,v3plotPhadrons,v3plotNhadrons);

	c->cd();
	v2plotPpartons->SetLineWidth(2);
	v2plotPpartons->SetLineColor(kBlack);
	v2plotPpartons->Draw();

	v2plotNpartons->SetLineWidth(2);
	v2plotNpartons->SetLineColor(kBlack);
	v2plotNpartons->Draw();

	v2plotPhadrons->SetLineWidth(2);
	v2plotPhadrons->SetLineColor(kBlack);
	v2plotPhadrons->Draw();

	v2plotNhadrons->SetLineWidth(2);
	v2plotNhadrons->SetLineColor(kBlack);
	v2plotNhadrons->Draw();

	v2plotPpartonseta1->SetLineWidth(2);
	v2plotPpartonseta1->SetLineColor(kBlack);
	v2plotPpartonseta1->Draw();

	v2plotNpartonseta1->SetLineWidth(2);
	v2plotNpartonseta1->SetLineColor(kBlack);
	v2plotNpartonseta1->Draw();

	v2plotPpartonsetap5->SetLineWidth(2);
	v2plotPpartonsetap5->SetLineColor(kBlack);
	v2plotPpartonsetap5->Draw();

	v2plotNpartonsetap5->SetLineWidth(2);
	v2plotNpartonsetap5->SetLineColor(kBlack);
	v2plotNpartonsetap5->Draw();

	c1->cd();
	v3plotPpartons->SetLineWidth(2);
	v3plotPpartons->SetLineColor(kBlack);
	v3plotPpartons->Draw();

	v3plotNpartons->SetLineWidth(2);
	v3plotNpartons->SetLineColor(kBlack);
	v3plotNpartons->Draw();

	v3plotPhadrons->SetLineWidth(2);
	v3plotPhadrons->SetLineColor(kBlack);
	v3plotPhadrons->Draw();

	v3plotNhadrons->SetLineWidth(2);
	v3plotNhadrons->SetLineColor(kBlack);
	v3plotNhadrons->Draw();
	TFile *f = new TFile("out_masterEllipticFlowCalc.root","RECREATE");
	//partons
	v2plotPpartons->Write();
	v2plotPhadrons->Write();
	v2plotPpartonseta1->Write();
	v2plotPpartonsetap5->Write();
	v3plotPpartons->Write();
	v3plotPhadrons->Write();
	//nucleons
	v2plotNpartons->Write();
	v2plotNhadrons->Write();
	v2plotNpartonseta1->Write();
	v2plotNpartonsetap5->Write();
	v3plotNpartons->Write();
	v3plotNhadrons->Write();

	if (setFlag == 0) {
		for (unsigned int i = 0; i < psi2valuespartons.size(); i++) {
			cout << "Psi 2 for partons is: " << psi2valuespartons[i] << endl;
		}
	}
	else if (setFlag == 1) {
		for (unsigned int i = 0; i < psi2valuesnucleon.size(); i++) {
			cout << "Psi 2 for nucleons is: " << psi2valuesnucleon[i] << endl;
		}
	}

	for (unsigned int j = 0; j < psi2valuesnucleon.size(); j++) {
			float difference = psi2valuesnucleon[j] - psi2valuespartons[j];
			hist1->Fill(difference);
	}

	c2->cd();
	hist1->Draw();
	formationdist1->Draw();
	formationdist2->Draw();
	formationdist3->Draw();
	formationdist4->Draw();
	fullformationdist->Draw();

	hist1->Write();
	formationdist1->Write();
	formationdist2->Write();
	formationdist3->Write();
	formationdist4->Write();
	fullformationdist->Write();

	f->Close();

	cout << "Finished running analysis code..." << endl;
	return;
	
}