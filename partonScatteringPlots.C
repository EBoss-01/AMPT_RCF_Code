//--------------------------------------------------------------------------------------------
// This code matches the initial partons to their scatterings using two AMPT files.
// The code then uses a series of pT cuts to plot the scattering probabilities.
//
// File 1: parton-initial-afterPropagation.dat
// File 2: parton-collisionsHistory.dat
//
// Created: 08-08-2018
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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//-------------------------------------------
// Structures
//-------------------------------------------

struct parton {
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
// Spacetime
	float x;
	float y;
	float z;
	float t;
// Additional
	float pT;
};

//---------------------------------------------
// Variables
//---------------------------------------------

// Counters
int eventnumber = 0;
int totalcounter = 0;
int noptbin = 0;
int fillcounter = 0;

//------------------
// Vectors
//------------------
vector<parton> ScatteringTime;
vector<vector<parton> > EventPartons;

// Energy
float energy1 = 0;
float energy2 = 0;
float energy3 = 0;

// Plots
TH1F* scattering_pt1 = new TH1F("scattering_pt1",";n scatterings;Probability",14,0,14);
TH1F* scattering_pt2 = new TH1F("scattering_pt2",";n scatterings;Probability",14,0,14);
TH1F* scattering_pt3 = new TH1F("scattering_pt3",";n scatterings;Probability",14,0,14);
TH1F* scattering_pt4 = new TH1F("scattering_pt4",";n scatterings;Probability",14,0,14);

//-------------------------------------------------------
// Functions
//-------------------------------------------------------

void myText(Double_t x,Double_t y,Color_t color,const char *text,Double_t tsize = 0.05,double angle=-1) {
	TLatex l;
	l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextColor(color);
	if (angle > 0) l.SetTextAngle(angle);
	l.DrawLatex(x,y,text);
}

void fillPlots(vector<parton> v) {

	int scatteringcounterpt1 = 0;
	int scatteringcounterpt2 = 0;
	int scatteringcounterpt3 = 0;
	int scatteringcounterpt4 = 0;

	for (unsigned int i = 0; i < v.size(); i++) {
		totalcounter++;

		if ((v[i].pT >= 0.2) && (v[i].pT < 0.4)) {
			scatteringcounterpt1++;
		}
		if ((v[i].pT >= 0.4) && (v[i].pT < 0.6)) {
			scatteringcounterpt2++;
		}
		if ((v[i].pT >= 0.6) && (v[i].pT < 0.8)) {
			scatteringcounterpt3++;
		}
		if ((v[i].pT >= 0.8) && (v[i].pT < 1.0)) {
			scatteringcounterpt4++;
		}

		if (i == (v.size()-1)) {
			fillcounter++;
			if ((scatteringcounterpt1-1) >= 0) {
				scattering_pt1->Fill(scatteringcounterpt1-1);
			}
			if ((scatteringcounterpt2-1) >= 0) {
				scattering_pt2->Fill(scatteringcounterpt2-1);
			}
			if ((scatteringcounterpt3-1) >= 0) {
				scattering_pt3->Fill(scatteringcounterpt3-1);
			}
			if ((scatteringcounterpt4-1) >= 0) {
				scattering_pt4->Fill(scatteringcounterpt4-1);
			}
		}
	}
}

void processEvent() {

	for (unsigned int p = 0; p < EventPartons.size(); p++) {

		std::vector<parton> v = EventPartons[p];

		fillPlots(v);
	}

	EventPartons.clear();
}

//-------------------------------------------------------
// Main Code
//-------------------------------------------------------

void partonScatteringPlots(void) {

	ifstream myFileOne;
	ifstream myFileTwo;

	myFileOne.open("ana/parton-initial-afterPropagation.dat");

	if (!myFileOne) {
		cout << "Unable to open parton-initial-afterPropagation." << endl;
		return;
	}
	else {
		cout << "Successfully opened parton-initial-afterPropagation." << endl;
	}

	while (myFileOne) {

		// Variables needed to read the event header.
		int iterationN;
		int nPartons;
		int nBaryons;
		int nMesons;
		int particleC;
		int particleNC;

		myFileOne >> eventnumber >> iterationN >> nPartons >> nBaryons >> nMesons >> particleC >> particleNC;

		if (!myFileOne) break;

		for (int i = 0; i < nPartons; i++) {

			int partID;
			float momenta[3];
			float mass;
			float Spacetime[4];

			myFileOne >> partID >> momenta[0] >> momenta[1] >> momenta[2] >> mass >> Spacetime[0] >> Spacetime[1] >> Spacetime[2] >> Spacetime[3];

			energy1 = TMath::Sqrt(pow(momenta[0],2) + pow(momenta[1],2) + pow(momenta[2],2) + pow(mass,2));
			TLorentzVector ev1(momenta[0],momenta[1],momenta[2],energy1);

			parton partinfo;
			partinfo.evtN = eventnumber;
			partinfo.PID = partID;
			partinfo.px = momenta[0];
			partinfo.py = momenta[1];
			partinfo.pz = momenta[2];
			partinfo.m = mass;
			partinfo.x = Spacetime[0];
			partinfo.y = Spacetime[1];
			partinfo.z = Spacetime[2];
			partinfo.t = Spacetime[3];
			partinfo.pT = ev1.Pt();

			vector<parton> auxillary;
			auxillary.push_back(partinfo);
			EventPartons.push_back(auxillary);

		}

		//-----------------------------------------
		// Variables for FileTwo
		//-----------------------------------------

		string line;
		int evt;
		int junk1;
		int partonindex1;
		int partonindex2;
		// Initial parton information
		int parton1_id_initial;
		float parton1_momenta_initial[3];
		float parton1_mass_initial;
		float parton1_spacetime_initial[4];

		int parton2_id_initial;
		float parton2_momenta_initial[3];
		float parton2_mass_initial;
		float parton2_spacetime_initial[4];
		// Final parton inforamtion
		int parton1_final_id;
		float parton1_final_momenta[3];
		float parton1_final_mass;
		float parton1_final_spacetime[4];

		int parton2_final_id;
		float parton2_final_momenta[3];
		float parton2_final_mass;
		float parton2_final_spacetime[4];

		myFileTwo.open("ana/parton-collisionsHistory.dat");

		if (!myFileTwo) {
			cout << "Unable to open parton-collisionsHistory." << endl;
			return;
		}

		while (std::getline(myFileTwo, line)) {

			if (!myFileTwo) break;

			stringstream heading(line);
			string description;

			if (heading >> description >> evt >> junk1 >> partonindex1 >> partonindex2) {

				if (evt == eventnumber) {

					myFileTwo >> parton1_id_initial >> parton1_momenta_initial[0] >> parton1_momenta_initial[1] >> parton1_momenta_initial[2] >> parton1_mass_initial >> parton1_spacetime_initial[0] >> parton1_spacetime_initial[1] >> parton1_spacetime_initial[2] >> parton1_spacetime_initial[3];

					myFileTwo >> parton2_id_initial >> parton2_momenta_initial[0] >> parton2_momenta_initial[1] >> parton2_momenta_initial[2] >> parton2_mass_initial >> parton2_spacetime_initial[0] >> parton2_spacetime_initial[1] >> parton2_spacetime_initial[2] >> parton2_spacetime_initial[3];

					myFileTwo >> parton1_final_id >> parton1_final_momenta[0] >> parton1_final_momenta[1] >> parton1_final_momenta[2] >> parton1_final_mass >> parton1_final_spacetime[0] >> parton1_final_spacetime[1] >> parton1_final_spacetime[2] >> parton1_final_spacetime[3];

					myFileTwo >> parton2_final_id >> parton2_final_momenta[0] >> parton2_final_momenta[1] >> parton2_final_momenta[2] >> parton2_final_mass >> parton2_final_spacetime[0] >> parton2_final_spacetime[1] >> parton2_final_spacetime[2] >> parton2_final_spacetime[3];

					energy2 = TMath::Sqrt(pow(parton1_final_momenta[0],2) + pow(parton1_final_momenta[1],2) + pow(parton1_final_momenta[2],2) + pow(parton1_final_mass,2));
					TLorentzVector ev2(parton1_final_momenta[0],parton1_final_momenta[1],parton1_final_momenta[2],energy2);

					parton part1;
					part1.evtN = evt;
					part1.PID = parton1_final_id;
					part1.px = parton1_final_momenta[0];
					part1.py = parton1_final_momenta[1];
					part1.pz = parton1_final_momenta[2];
					part1.m = parton1_final_mass;
					part1.x = parton1_final_spacetime[0];
					part1.y = parton1_final_spacetime[1];
					part1.z = parton1_final_spacetime[2];
					part1.t = parton1_final_spacetime[3];
					part1.pT =  ev2.Pt();

					energy3 = TMath::Sqrt(pow(parton2_final_momenta[0],2) + pow(parton2_final_momenta[1],2) + pow(parton2_final_momenta[2],2) + pow(parton2_final_mass,2));
					TLorentzVector ev3(parton2_final_momenta[0],parton2_final_momenta[1],parton2_final_momenta[2],energy3);

					parton part2;
					part2.evtN = evt;
					part2.PID = parton2_final_id;
					part2.px = parton2_final_momenta[0];
					part2.py = parton2_final_momenta[1];
					part2.pz = parton2_final_momenta[2];
					part2.m = parton2_final_mass;
					part2.x = parton2_final_spacetime[0];
					part2.y = parton2_final_spacetime[1];
					part2.z = parton2_final_spacetime[2];
					part2.t = parton2_final_spacetime[3];
					part2.pT = ev3.Pt();

					EventPartons[partonindex1 - 1].push_back(part1);
					EventPartons[partonindex2 - 1].push_back(part2);

				}
			}
		}

		myFileTwo.close();

		processEvent();

	}

	TCanvas* c = new TCanvas("c","Scattering Plots",700,700);
	gStyle->SetOptStat(0);

	scattering_pt1->Scale(1.0/(scattering_pt1->Integral()));
	scattering_pt1->Draw();
	scattering_pt2->Scale(1.0/(scattering_pt2->Integral()));
	scattering_pt2->Draw();
	scattering_pt3->Scale(1.0/(scattering_pt3->Integral()));
	scattering_pt3->Draw();
	scattering_pt4->Scale(1.0/(scattering_pt4->Integral()));
	scattering_pt4->Draw();

	TFile *f = new TFile("out_partonScatteringPlots.root","RECREATE");
	scattering_pt1->Write();
	scattering_pt2->Write();
	scattering_pt3->Write();
	scattering_pt4->Write();

	f->Close();

	myFileOne.close();
}
