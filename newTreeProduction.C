//-------------------------------------------------------------------------------------------------
// This code is used to create TTrees of the parton Evolution and calculate the participant planes.
//
// File 1: parton-initial-afterPropagation.dat (used for parton plane and initial values)
// File 2: parton-collisionsHistory.dat (used to follow parton evolution)
// File 3: npart-xy.dat (uses nucleons for Psi for set tform = const)
// File 4: parton-after-coalescence.dat (will be used for final momentum and hadron ID)
//
// Multimap just uses parton-after-coalescence.dat and the key is the last scattering.
//
// 09-26-18
// Updated 10-16-18 (added code to ignore additional iterations where applicable and comments)
//-------------------------------------------------------------------------------------------------

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
// TTree Creation
//---------------------------------

TFile* f1 = NULL;
TTree *tree = NULL;

//---------------------------------
//Structure for the Parton
//---------------------------------

struct parton {
// Event Number
	int evtN;
// Parton ID
	int pID;
// Nucleon Plane
	float Npsi;
// Parton Plane
	float Ppsi;
// Parton initial momentum
	float px;
	float py;
	float pz;
// Parton initial location
	float x;
	float y;
	float z;
	float t;
// Parton mass in GeV
	float m;
/*// Final State
	float Finalpx;
	float Finalpy;
	float Finalpz;
// Final hadron
	int hadron;*/

};

// This is for File 3.
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


//----------------------------------
// Variables
//----------------------------------

// This is the counter to keep track of which event is being read.
int eventnumber = 0;

// Vectors: 
// This vector is used to store all the parton information.
vector<vector<parton> > EventPartons;
// This vector is used to store nucleon information for psi calculation.
vector<nucleon> initialnucleon;
vector<parton> initialpartons;
vector<float> psi2valuespartons;
vector<float> psi2valuesnucleon;

float psi2partons;
float psi2nucleon;

// Vaiables needed for multimap
int Pevent = 0;
float Ppx = 0;
float Ppy = 0;
float Ppz = 0;
string pkey;
typedef multimap<string, int> PartonMap;

//----------------------------------
//Functions
//----------------------------------

void myText(Double_t x,Double_t y,Color_t color,const char *text,Double_t tsize = 0.05,double angle=-1) {
	TLatex l;
	l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextColor(color);
	if (angle > 0) l.SetTextAngle(angle);
	l.DrawLatex(x,y,text);
}

// This function searches the multimap with the determined key.
void findHadronValue(PartonMap &keymomentums, string &pkey, float &hadronnumber) {

	PartonMap::iterator it = keymomentums.find(pkey);
	hadronnumber = keymomentums.find(pkey)->second;
}

// This function fills the mutltimap in order to find hadron formation information. 
void parseCoalescenceFile(PartonMap &keymomentums) {

	string line2;
	int evt;
	int npartons;
	int nbaryons;
	int nmesons;
	float imparm;
	int stuff[5];

	ifstream myFileFour;

	myFileFour.open("ana/parton-after-coalescence.dat");

	if (!myFileFour) {
		cout << "Unable to open parton-after-coalescence file." << endl;
		return;
	}

	else {
		cout << "Successfully opened parton-after-coalescence file." << endl;
	}

	int eventcounter = 0;
	int participantparton = 0;
	int counterparton = 0;

	while (std::getline(myFileFour, line2)) {

		if (!myFileFour) break;

		stringstream eventheader(line2);
		string intermediate;

		vector<float> token;

		while (getline(eventheader, intermediate, ' ')) {

			if (intermediate != "") {
				token.push_back(atof(intermediate.c_str()));
			}
		}

		if (token.size() == 10) {
			eventcounter++;
			participantparton = token[1];
		}

		else if ((token.size() == 8) || (token.size() == 7)) {

			Pevent = eventcounter;
			Ppx = token[1];
			Ppy = token[2];
			Ppz = token[3];

			pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);
			keymomentums.insert(make_pair(pkey, token[6]));
		}
	}

	cout << "Closing parton-after-coalescence file." << endl;

	myFileFour.close();
}

// Function used to calculate the participant plane from Nucleons. Needed for set formation times. 
void calculateParticipantPlaneNucleons() {

	int nucleoncounter = 0;
	int Ncmcounter = 0;
	float Ncmx = 0;
	float Ncmy = 0;
	float Nq2x = 0;
	float Nq2y = 0;
	float Nrsquare = 0;

	// Needed for smearing.
	TF1 *fradial = new TF1("fradial","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradial->SetParameter(0,0.4);

	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());

	vector<nucleon> vN = initialnucleon;

	// All center of mass frame adjustments.
	for (unsigned int k = 0; k < vN.size(); k++) {
		if (vN[k].status > 0) {
			Ncmx = Ncmx + vN[k].x;
			Ncmy = Ncmy + vN[k].y;

			Ncmcounter++;
		}
	}

	// Averaging the center of mass.
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

	 		Nrsquare = Nrsquare + pow(sr,2);

	 		nucleoncounter++;
	 	}
	}

	Nq2x = Nq2x/(float)nucleoncounter;
	Nq2y = Nq2y/(float)nucleoncounter;

	Nrsquare = Nrsquare/(float)nucleoncounter;

	// Actual participant plane calculation. 
	psi2nucleon = (TMath::ATan2(Nq2y, Nq2x)/2.0) + (TMath::Pi()/2.0);

	psi2valuesnucleon.push_back(psi2nucleon);

}

// Funcation used to calculate the participant plane using Partons.
void calculateParticipantPlanePartons() {

	int cmcounter = 0;
	float cmx = 0;
	float cmy = 0;
	float q2x = 0;
	float q2y = 0;
	float rsquare = 0;

	vector<parton> v = initialpartons;

	// Used to begin adjusting into center of mass frame.
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
		//Component calculations for psi.
		q2x = q2x + pow(rvalue,2)*TMath::Cos(2*phivalue);
		q2y = q2y + pow(rvalue,2)*TMath::Sin(2*phivalue);
		rsquare = rsquare + pow(rvalue,2);
	}

	// Averaging components
	q2x = q2x/(float)ncounter;
	q2y = q2y/(float)ncounter;

    rsquare = rsquare/(float)ncounter;

    // Participant Plane calculations
	psi2partons = (TMath::ATan2(q2y, q2x)/2.0) + (TMath::Pi()/2.0);

	// Pushing the values to more vectors
	psi2valuespartons.push_back(psi2partons);
}

// This function is used to calculate the Nucleon Participant Plane.
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
	ifstream myFileThree;

	// File needed for nucleon geometry. 
	myFileThree.open("ana/npart-xy.dat");

	if (!myFileThree) {
		cout << "Unable to open npart-xy." << endl;
		return;
	}
	else {
		cout << "Successfully opened npart-xy." << endl;
	}

	while (std::getline(myFileThree, line1)) {

		// Prevents duplicating last line.
		if (!myFileThree) break;

		stringstream header(line1);

		// This reads the headers for each event.
		if (header >> eventnumber >> Niteration >> atomicMproj >> atomicMtarg >> impactparm) {

			// This prevents duplication of events.
			if (Niteration > 0) {
				continue;
			}

			// Loop for filling vector.
			for (int i = 0; i < (atomicMproj + atomicMtarg); i++) {

				myFileThree >> space[0] >> space[1] >> sequence >> stat >> junk[0] >> junk[1] >> junk[2];

				nucleon nucinitial;
				nucinitial.evtN = eventnumber;
				nucinitial.x = space[0];
				nucinitial.y = space[1];
				nucinitial.seq = sequence;
				nucinitial.status = stat;

				initialnucleon.push_back(nucinitial);

			}

			// Call to function used to calculate the psi value for nucleons.
			calculateParticipantPlaneNucleons();

			initialnucleon.clear();
		}
	}

	myFileThree.close();

}

// This function is used to pull in the positions for every parton and scattering. The banches of TTree are filled here.
void calculatePosition (vector<parton> v, int &actualevent, int &partonid, float &NPsi, float &PPsi, float &initialpx, float &initialpy, float &initialpz, float &formationt, float &initialx, float &initialy, float &initialz, float &initialmass, int &scatteringn, float scatteringpx[], float scatteringpy[], float scatteringpz[], float scatteringt[], float scatteringx[], float scatteringy[], float scatteringz[], float &finalpx, float &finalpy, float &finalpz, int &formedhadron, PartonMap &keymomentums) {

	float x0, y0;

	// The v vector in this loop is looking at the scattering's for each individual parton so the size varies in length.
	for (unsigned int i = 0; i < v.size(); i++) {

		scatteringn = v.size()-1;

		if (i == v.size()-1) {

			// Seting up the multimap key elements. 
			Pevent = eventnumber;
			Ppx = v[i].px;
			Ppy = v[i].py;
			Ppz = v[i].pz;

			pkey = Form("%d*%g*%g*%g", Pevent, Ppx, Ppy, Ppz);

			float hadronnumber;
			// Call to find the hadron value.
			findHadronValue(keymomentums, pkey, hadronnumber);

			// Fills the branches for the final parton position and hadron value.
			finalpx = v[i].px;
			finalpy = v[i].py;
			finalpz = v[i].pz;
			formedhadron = hadronnumber;

		}


		// Determine the initial position of the parton.
		if (i == 0) {

			x0 = v[i].x;
			y0 = v[i].y;

			// This is filling in the branch for initial conditions.
			actualevent = eventnumber;
			partonid = v[i].pID;
			NPsi = psi2valuesnucleon[eventnumber-1];
			PPsi = psi2valuespartons[0];
			initialpx = v[i].px;
			initialpy = v[i].py;
			initialpz = v[i].pz;
			formationt = v[i].t;
			initialx = x0;
			initialy = y0;
			initialz = v[i].z;
			initialmass = v[i].m;

		}
		else {

			// This is filling in the branches for scattering.
			scatteringpx[i-1] = v[i].px;
			scatteringpy[i-1] = v[i].py;
			scatteringpz[i-1] = v[i].pz;
			scatteringt[i-1] = v[i].t;
			scatteringx[i-1] = v[i].x; 
			scatteringy[i-1] = v[i].y; 
			scatteringz[i-1] = v[i].z;
		}
	}

}

// This function will loop over each parton and then fill the branches of the TTree. 
void processEvent(int &actualevent, int &partonid, float &NPsi, float &PPsi, float &initialpx, float &initialpy, float &initialpz, float &formationt, float &initialx, float &initialy, float &initialz, float &initialmass, int &scatteringn, float scatteringpx[], float scatteringpy[], float scatteringpz[], float scatteringt[], float scatteringx[], float scatteringy[], float scatteringz[], float &finalpx, float &finalpy, float &finalpz, int &formedhadron, PartonMap &keymomentums) {
	
	// This is a call for the calculation of parton plane.
	calculateParticipantPlanePartons();

	// This counter used to prevent the tree from being filled more than needed.
	unsigned int counter1 = 0;

		// Put Parton Loop here. Each parton is looped over here to calculate the needed information. 
		
		for (unsigned int p = 0; p < EventPartons.size(); p++) {

			std::vector<parton> v = EventPartons[p];

			// Call to function to fill the tree branches.
			calculatePosition(v,actualevent,partonid,NPsi,PPsi,initialpx,initialpy,initialpz,formationt,initialx,initialy,initialz,initialmass,scatteringn,scatteringpx,scatteringpy,scatteringpz,scatteringt,scatteringx,scatteringy,scatteringz,finalpx,finalpy,finalpz,formedhadron,keymomentums);

			counter1 = counter1;
			counter1++;

			// This statement is to ensure that the tree is only filled once for all partons.
			if (counter1 <= EventPartons.size()) {
				tree->Fill();
			}
			else {
				continue;
			}

		}

	// Clearing vectors to prevent problem information.
	EventPartons.clear();
	psi2valuespartons.clear();
}

//-----------------------------------------------------
// Main Code
//-----------------------------------------------------
void newTreeProduction(void) {

	// This sets up the multimap for use.
	PartonMap keymomentums;

	// Call to function that fills the multimap.
	parseCoalescenceFile(keymomentums);

	f1 = new TFile("New_AMPT_TTree_out.root", "RECREATE");

	const Int_t maxArrayLength = 10000;
	Int_t actualevent;
	Int_t partonid;
	Float_t NPsi;
	Float_t PPsi;
	Float_t initialpx;
	Float_t initialpy;
	Float_t initialpz;
	Float_t formationt;
	Float_t initialx;
	Float_t initialy;
	Float_t initialz;
	Float_t initialmass;
	Int_t scatteringn = 10000;
	Float_t scatteringpx[maxArrayLength];
	Float_t scatteringpy[maxArrayLength];
	Float_t scatteringpz[maxArrayLength];
	Float_t scatteringt[maxArrayLength];
	Float_t scatteringx[maxArrayLength];
	Float_t scatteringy[maxArrayLength];
	Float_t scatteringz[maxArrayLength];
	Float_t finalpx;
	Float_t finalpy;
	Float_t finalpz;
	Int_t formedhadron;

	tree = new TTree("tree", "A juicy apple tree");
	tree->Branch("event_number",&actualevent,"event_number/I");
	tree->Branch("parton_id",&partonid,"parton_id/I");
	tree->Branch("nucleon_plane",&NPsi,"nucleon_plane/F");
	tree->Branch("parton_plane",&PPsi,"parton_plane/F");
	tree->Branch("initial_px",&initialpx,"initial_px/F");
	tree->Branch("initial_py",&initialpy,"initial_py/F");
	tree->Branch("initial_pz",&initialpz,"initial_pz/F");
	tree->Branch("formation_t",&formationt,"formation_t/F");
	tree->Branch("initial_x",&initialx,"initial_x/F");
	tree->Branch("initial_y",&initialy,"initial_y/F");
	tree->Branch("initial_z",&initialz,"initial_z/F");
	tree->Branch("initial_m",&initialmass,"initial_m/F");
	tree->Branch("scattering_n",&scatteringn);
	tree->Branch("scattering_px",scatteringpx,"scattering_px[scattering_n]/F");
	tree->Branch("scattering_py",scatteringpy,"scattering_py[scattering_n]/F");
	tree->Branch("scattering_pz",scatteringpz,"scattering_pz[scattering_n]/F");
	tree->Branch("scattering_t",scatteringt,"scattering_t[scattering_n]/F");
	tree->Branch("scattering_x",scatteringx,"scattering_x[scattering_n]/F");
	tree->Branch("scattering_y",scatteringy,"scattering_y[scattering_n]/F");
	tree->Branch("scattering_z",scatteringz,"scattering_z[scattering_n]/F");
	tree->Branch("final_px",&finalpx,"final_px/F");
	tree->Branch("final_py",&finalpy,"final_py/F");
	tree->Branch("final_pz",&finalpz,"final_pz/F");
	tree->Branch("hadron_formed",&formedhadron,"hadron_formed/I");
	ifstream myInitialFileInfo;
	ifstream myEvolutionFile;

	// Call to function that begins Nucleon Plane calculation. 
	parseNucleons();

	myInitialFileInfo.open("ana/parton-initial-afterPropagation.dat");

	if (!myInitialFileInfo) {
	// This will let me know if the file fails to open and specifically which file.
		cout << "Unable to open file parton-initial-afterPropagation.dat" << endl;
		return;
	}
	else {
	// I've added this piece simply to confirm that the file did indeed open.
		cout << "Successfully opened parton-initial-afterPropagation.dat" << endl;
	}

	while (myInitialFileInfo) {

		if (eventnumber % 100 == 0) {

			cout << "Reading Event Number " << eventnumber << endl;
		}

		// Information needed to read the event header.
		int iterationN;
		int nPartons;
		int nBaryons;
		int nMesons;
		int particleC;
		int particleNC;

		myInitialFileInfo >> eventnumber >> iterationN >> nPartons >> nBaryons >> nMesons >> particleC >> particleNC;

		if (iterationN > 0) {
			continue;
		}

		// This line prevents the file from trying to read the last line multiple times.
		if (!myInitialFileInfo) break;

		// This loop fills in the initial information for each parton.
		for (int i = 0; i < nPartons; i++) {

			int partID;
			float momenta[3];
			float mass;
			double spacetime[4];

			myInitialFileInfo >> partID >> momenta[0] >> momenta[1] >> momenta[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

			parton partinfo;
			partinfo.evtN = eventnumber;
			partinfo.pID = partID;
			partinfo.px = momenta[0];
			partinfo.py = momenta[1];
			partinfo.pz = momenta[2];
			partinfo.m = mass;
			partinfo.x = spacetime[0];
			partinfo.y = spacetime[1];
			partinfo.z = spacetime[2];
			partinfo.t = spacetime[3];

			vector<parton> auxillary;
			auxillary.push_back(partinfo);
			initialpartons.push_back(partinfo);
			EventPartons.push_back(auxillary);

		}

		//-------------------------------------------
		// Additional variables for evolution
		//-------------------------------------------
		string line;
		int evt;
		int junk1;
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

		float benergy1;
		float benergy2;
		float benergy3;
		float benergy4;

		TLorentzVector com;
		TLorentzVector vb1;
		TLorentzVector vb2;
		TLorentzVector vb3;
		TLorentzVector vb4;

		// This is where the evolution file is opened each time.
		myEvolutionFile.open("ana/parton-collisionsHistory.dat");

		if (!myEvolutionFile) {
		// As before this will let me know if this specific file fails to open.
			cout << "Unable to open file parton-collisionsHistory.dat" << endl;
			return;
		}
		/*else {
		// This is here simply to confirm that this file opened successfully as well.
			cout << "Successfully opened parton-collisionsHistory.dat" << endl;
		}*/

		while (std::getline(myEvolutionFile,line)) {

			if (!myEvolutionFile) break;
			
			stringstream heading(line);
			string description;

			// Reading the file in.
			if (heading >> description >> evt >> junk1 >> partonindex1 >> partonindex2) {

				if (junk1 > 0) {
					continue;
				}

				if (evt == eventnumber) {

					myEvolutionFile >> parton1_id_initial >> parton1_momenta_initial[0] >> parton1_momenta_initial[1] >> parton1_momenta_initial[2] >> parton1_mass_initial >> parton1_spacetime_initial[0] >> parton1_spacetime_initial[1] >> parton1_spacetime_initial[2] >> parton1_spacetime_initial[3];

					myEvolutionFile >> parton2_id_initial >> parton2_momenta_initial[0] >> parton2_momenta_initial[1] >> parton2_momenta_initial[2] >> parton2_mass_initial >> parton2_spacetime_initial[0] >> parton2_spacetime_initial[1] >> parton2_spacetime_initial[2] >> parton2_spacetime_initial[3];

					myEvolutionFile >> parton1_final_id >> parton1_final_momenta[0] >> parton1_final_momenta[1] >> parton1_final_momenta[2] >> parton1_final_mass >> parton1_final_spacetime[0] >> parton1_final_spacetime[1] >> parton1_final_spacetime[2] >> parton1_final_spacetime[3];

					myEvolutionFile >> parton2_final_id >> parton2_final_momenta[0] >> parton2_final_momenta[1] >> parton2_final_momenta[2] >> parton2_final_mass >> parton2_final_spacetime[0] >> parton2_final_spacetime[1] >> parton2_final_spacetime[2] >> parton2_final_spacetime[3];


					parton part1;
					part1.evtN = evt;
					part1.pID = parton1_final_id;
					part1.px = parton1_final_momenta[0];
					part1.py = parton1_final_momenta[1];
					part1.pz = parton1_final_momenta[2];
					part1.x = parton1_final_spacetime[0];
					part1.y = parton1_final_spacetime[1];
					part1.z = parton1_final_spacetime[2];
					part1.t = parton1_final_spacetime[3];
					part1.m = parton1_final_mass;

					parton part2;
					part2.evtN = evt;
					part2.pID = parton2_final_id;
					part2.px = parton2_final_momenta[0];
					part2.py = parton2_final_momenta[1];
					part2.pz = parton2_final_momenta[2];
					part2.x = parton2_final_spacetime[0];
					part2.y = parton2_final_spacetime[1];
					part2.z = parton2_final_spacetime[2];
					part2.t = parton2_final_spacetime[3];
					part2.m = parton2_final_mass;

					EventPartons[partonindex1 - 1].push_back(part1);
					EventPartons[partonindex2 - 1].push_back(part2);
				}
			}
		}

		// This is where the evolution file is closed each time.
		myEvolutionFile.close();

		// Call function that Processes the Event.
		processEvent(actualevent,partonid,NPsi,PPsi,initialpx,initialpy,initialpz,formationt,initialx,initialy,initialz,initialmass,scatteringn,scatteringpx,scatteringpy,scatteringpz,scatteringt,scatteringx,scatteringy,scatteringz,finalpx,finalpy,finalpz,formedhadron,keymomentums);

	}

	f1->Write();
	f1->Close();

	myInitialFileInfo.close();
	return;
}