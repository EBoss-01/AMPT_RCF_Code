//-----------------------------------------------------------------------------------------
//This code allows me to calculate the v2 and v3 from AMPT.
//This particular code utilizes initial parton information and hadrons.
//
//File 1: npart-xy.dat
//File 2: ampt.dat
//
//05-21-18
//Updated 07-10-18
//This update added CM adjustment and v3 calculations to this file.
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

// This is for the Final State particles.
struct particles {
// Event Number
	int evtnumber;
// Particle ID
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
};

struct phistuff {
// Event Number
	int evtN;
// Phi values
	float phi;
// pT values
	float pT;
  //eta values
  float eta;
};


//----------------------------------
// Variables
// Updated 7-10-18
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
vector<float> psi3values;
vector<float> epsilon3values;

vector<particles> finalstate;
vector<phistuff> phivalues;

float psi2;
float psi3;
float phi;
float pT;
float smearpsi2;
float AverageE2;
float AverageE3;
float TotalE2;
float TotalE3;

//----------------------------------
//Functions
// Updated 7-10-18
//----------------------------------

// This function does the calculation for the Participant plane to be used in the calculation of v2.
void calculateParticipantPlane() {

	// This is a counter to enable getting an average.
	int nparton = 0;
	int cmcounter = 0;
	// Variables used in the calculation of psi2.
	float q2x = 0;
	float q2y = 0;
	float q3x = 0;
	float q3y = 0;
	float e21 = 0;
	float e22 = 0;
	float e2 = 0;
	float epsilon3 = 0;
	float rsqure = 0;
	// Variables added for Center of Mass adjustment.
	float cmx = 0;
	float cmy = 0;

	TF1 *fradial = new TF1("fradial","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradial->SetParameter(0,0.4);

	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());

	vector<parton> v = allparticles;

	//This loop is used for the center of mass adjustment.
	for (unsigned int i = 0; i < v.size(); i++) {
		if (v[i].status > 0) {
			cmx = cmx + v[i].x;
			cmy = cmy + v[i].y;

			cmcounter++;
		}
	}

	cmx = cmx/(float)cmcounter;
	cmy = cmy/(float)cmcounter;

	 for (unsigned int i = 0; i < v.size(); i++) {

	 	if (v[i].status > 0) {

	 		v[i].x = v[i].x - cmx;
	 		v[i].y = v[i].y - cmy;

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
	 		// psi3 components with smear.
	 		q3x = q3x + pow(sr,2)*TMath::Cos(3*sphi);
	 		q3y = q3y + pow(sr,2)*TMath::Sin(3*sphi);

	 		// For eccentricity calculations

	 		e21 = e21 + pow(sr,2)*TMath::Cos(2*sphi);
	 		e22 = e22 + pow(sr,2)*TMath::Sin(2*sphi);
	 		rsqure = rsqure + pow(sr,2);

	 		nparton++;
	 	}
	 }

	 //cout << nparton << endl;

	 q2x = q2x/nparton;
	 q2y = q2y/nparton;
	 q3x = q3x/nparton;
	 q3y = q3y/nparton;
	 e21 = e21/nparton;
	 e22 = e22/nparton;
     rsqure = rsqure/nparton;

	 psi2 = (TMath::ATan2(e22, e21)/2.0) + (TMath::Pi()/2.0);
	 psi3 = (TMath::ATan2(q3y, q3x)/3.0) + (TMath::Pi()/3.0);

	 //cout << ":::::::" << psi2 << endl;
	 e2 = (TMath::Sqrt(pow(e21,2) + pow(e22,2))/rsqure);
	 epsilon3 = (TMath::Sqrt(pow(q3x,2) + pow(q3y,2)))/rsqure;

	 participantplane.push_back(psi2);
	 eccentricity.push_back(e2);
	 psi3values.push_back(psi3);
	 epsilon3values.push_back(epsilon3);

}

void calculatePhiValues() {

	int amountparticle = 0;

	vector<particles> p = finalstate;

	for (unsigned int i = 0; i < p.size(); i++) {

	  if ((fabs(p[i].PID) == 211) || (fabs(p[i].PID) == 2212) || (fabs(p[i].PID) == 321)) {
		  
		  if ((p[i].px == 0) && (p[i].py == 0)) {
		    continue;
		  }
		  float phienergy = TMath::Sqrt(pow(p[i].px,2) + pow(p[i].py,2) + pow(p[i].pz,2) + pow(p[i].m,2));

		  TLorentzVector ev(p[i].px, p[i].py, p[i].pz, phienergy);

		  //phi = TMath::ATan2(p[i].py, p[i].px);
		  //pT = TMath::Sqrt(pow(p[i].px,2) + pow(p[i].py,2));

			/*if (p[i].evtnumber == 1) {
				cout << p[i].PID << ",   " << phi << endl;
			}*/

			phistuff phitotal;
			phitotal.evtN = p[i].evtnumber;
			phitotal.phi = ev.Phi();
			phitotal.pT = ev.Pt();
			phitotal.eta = ev.Eta();

			phivalues.push_back(phitotal);

			amountparticle++;
		}
	}

	//cout << amountparticle << endl;

	//cout << phivalues[0].evtN << ", " << phivalues.size() << endl;
}

void calculateFlow() {

	TCanvas *c = new TCanvas("c","v2 calculations",700,700);
	gStyle->SetOptStat(0);

	TProfile *v2plot = new TProfile("v2plot","v_{2} vs p_{T}",20,0,2.5,-1,1);
	//TH2D *v2plot = new TH2D("v2plot","v_{2} vs. p_{T}",200,0,2.5,100,-1,1);

	TCanvas *c1 = new TCanvas("c1","v3 calculations",700,700);
	gStyle->SetOptStat(0);

	TProfile *v3plot = new TProfile("v3plot","v_{3} vs. p_{T}",20,0,2.5,-1,1);

	float v2;
	float v3;

	for (unsigned int i = 0; i < participantplane.size(); i++) {

		for (unsigned int j = 0; j < phivalues.size(); j++) {

			if (phivalues[j].evtN == (i+1)) {
			  
			  if (fabs(phivalues[j].eta) < 0.35) {

				v2 = TMath::Cos(2*(phivalues[j].phi - participantplane[i]));
				v2plot->Fill(phivalues[j].pT,v2);

				v3 = TMath::Cos(3*(phivalues[j].phi - psi3values[i]));
				v3plot->Fill(phivalues[j].pT,v3);
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
		TFile *f = new TFile("out_calculateEllipticFlow.root","RECREATE");
		v2plot->Write();
		v3plot->Write();
		f->Close();
		//c->Print(Form("v2plot_firsttry_%i.pdf",i));
	}
}
// This will be my attempt to read in the files that I will be using. 
void calculateEllipticFlow(void) {

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
	ifstream myAmptFile;

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
		  //cout << eventnumber << " " << iterationN << endl;

			int partoncounter = 0;

			// This loop fills in the initial information for each parton.
			for (int i = 0; i < (atomicMproj + atomicMtarg); i++) {


				myInitialFileInfo >> space[0] >> space[1] >> sequence >> stat >> junk[0] >> junk[1] >> junk[2];

				//cout << eventnumber << ", " << space[0] << ", " << space[1] << ", " << sequence << ", " << stat << endl;

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


			//cout << eventnumber << ", " << partoncounter << endl;

			calculateParticipantPlane();



			//cout << psi2 << endl;


			allparticles.clear();
			particleproj.clear();
			particletarg.clear();
		}

	}

	myInitialFileInfo.close();

	string line1;
	int evt;
	int test;
	int npartevt;
	float imparm[2];
	int totalpart;
	int participantproj;
	int participanttarg;
	int collisions[4];

	// Needed for the particle structure.
	int particleid;
	float momenta[3];
	float mass;
	float spacetime[4];

	myAmptFile.open("ana/ampt.dat");

	if (!myAmptFile) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file 2." << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened file 2." << endl;
	}

	while (std::getline(myAmptFile, line1)) {

		if (!myAmptFile) break;

		stringstream header(line1);

		if (header >> evt >> test >> npartevt >> imparm[0] >> totalpart >> participantproj >> participanttarg >> collisions[0] >> collisions[1] >> collisions[2] >> collisions[3] >> imparm[1]) {

		  //cout << "---------" << evt << " " << test << endl;

			for (int j = 0; j < npartevt; j++) {

				myAmptFile >> particleid >> momenta[0] >> momenta[1] >> momenta[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

				particles finalpart;
				finalpart.evtnumber = evt;
				finalpart.PID = particleid;
				finalpart.px = momenta[0];
				finalpart.py = momenta[1];
				finalpart.pz = momenta[2];
				finalpart.m = mass;
				finalpart.x = spacetime[0];
				finalpart.y = spacetime[1];
				finalpart.z = spacetime[2];
				finalpart.t = spacetime[3];

				finalstate.push_back(finalpart);

				//cout << finalstate[j].evtnumber << endl;

			}

			calculatePhiValues();

			finalstate.clear();
			//phivalues.clear();
		}
	}

	calculateFlow();

	for (unsigned int i = 0; i < eccentricity.size(); i++){
	  
	  TotalE2 = TotalE2 + eccentricity[i];
	  TotalE3 = TotalE3 + epsilon3values[i];
	}

	int totalevents = eccentricity.size();

	AverageE2 = TotalE2/float(totalevents);
	cout << "Average Epsilon 2: " << AverageE2 << endl;
	AverageE3 = TotalE3/float(totalevents);
	cout << "Average Epsilon 3: " << AverageE3 << endl;

	myAmptFile.close();

	cout << "Finished running analysis code..." << endl;
	return;
}

