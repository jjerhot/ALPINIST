 /// \file DecayMCGlobal.h
 /// \Brief
 /// Initialization of global constants, production and decay modes and paths
 /// Definition of AxionParameters structure
 /// \EndBrief

#ifndef DECAYMCGLOBAL_H
#define DECAYMCGLOBAL_H

#include <iostream>
#include <TROOT.h>
#include "TRandom3.h"

extern TRandom3 rndmGen;

const std::vector<TString> exoLabels{"alp","hnl","dp","ds"}; ///< available exotics

const std::vector<TString> yVars {"Width","Tau","cTau"}; ///< available y-axis variables: 0-Width (GeV), 1-Tau (fs), 2-cTau (m)

const std::vector<std::vector<TString>> prodmodeNames {
						{"Bmeson", "Dmeson","Primakoff","PhotonFromMeson","MixingPi0","MixingEta","MixingEtaPrim"},
						{"Bmeson", "Dmeson"},
						{"Brems", "MesonDecay","MixingRho","MixingOmega","MixingPhi"},
						{"Brems","Bmeson","Bmeson2S"}}; ///< associated production modes

const std::vector<std::vector<TString>>  decaymodeNames {
							{"2Gamma","2El","2Mu","3Pi0","3Pi","2PiGamma","2Pi0Eta","2PiEta","2Pi0EtaPrim","2PiEtaPrim","2Pi","2K","2KPi0","2Pi2Pi0"},
							{"PiEl", "PiMu", "PiPiEl", "PiPiMu", "NuElEl", "NuMuMu", "NuElMu", "PiNu", "EtaNu", "PiPiNu"  /*"PiTau", "RhoTau", "EtaPrimNu", "KMu", "KMu", "KTau", "KEl"*/},
							{"2El","2Mu","2Pi","3Pi","4Pi","2Pi2Pi0","2K","2KPi0"},
							{"2El","2Mu","2Pi","2K","4Pi","2Pi2Pi0"}}; ///< associated decay modes

const std::vector<std::vector<TString>> activeCouplings {
							{""},
							{"-ElMixing", "-MuMixing", "-TauMixing"},
							{""},
							{""}}; ///< associated to-SM-mixing the production was run in (relevant only for hnls)

const std::vector<TString> expLabels{"CHARM","BEBC","NuCal","NuTeV","NA62","DarkQuest",
									"DarkQuestPhase2","DUNE","SHiP","KOTOdump","KOTOpnn",
									"KOTO2dump","KOTO2pnn","KOTOexclPnn","KLEVER","KLEVERext",
									"SHADOWS","SHiPecn4","ORCA","BEBCcuboid"}; ///< experiments

/// Paths
const TString sourcePath = "../../tab_prod";
const TString outPath = "../../tab_decay";
const TString widthPath = "../../widths/Dalitz";

/// Constants
const Double_t c = 299792458.; // m/s
const Double_t hc = 0.197E-15; // GeV m
const Double_t MKCh = 0.4937; // GeV
const Double_t MPiCh = 0.1397; // GeV
const Double_t MMu = 0.105658; // GeV
const Double_t MEl = 0.000511; // GeV
const Double_t MPi0 = 0.13498; // GeV
const Double_t MEta = 0.54786; // GeV
const Double_t MEtaPrim = 0.95778; //GeV
const Double_t MB0 = 5.27963 ; //GeV
const Double_t MDs = 1.96834; //Ds^+- mass in GeV 2020PDG
const Double_t MDCh = 1869.66; //D^+- mass in GeV 2020PDG
const Double_t MRho = 0.77526; //GeV
const Double_t MP = 0.93827208816; //GeV PdG

/// Energy-dependent constants
const std::map<int,Double_t> sigma_pp = {{30,38.38*1E9},{70,38.38*1E9},{120,38.54*1E9},{400,39.85*1E9},{800,41.14*1E9}};
const std::map<int,Double_t> sigma_cc = {{30,0.},
										{70,(2.2/* ± 1.2*/)*1E6},
										{120,(4.6/* ± 1.8*/)*1E6},
										{400,/*2.3**/(20./*±4.5*/)*1E6},// SHiP cascade factor for cc@400GeV is 2.3
										{800,(39.5/*±7.7*/)*1E6}};  //
const std::map<int,Double_t> sigma_bb = {{30,0.},
										{70,0.},
										{120,20./*± 20*/},
										{400,/*1.7**/(2.63/* ± 0.76*/)*1E3}, // SHiP cascade factor for bb@400GeV is 1.7
										{800,(18.6/* ± 4.6*/)*1E3}}; 


const std::map<std::string,std::array<Double_t,3>> targetMat
	= {	{"Cu",{63.546,29,6.*1e12}},{"Fe",{56,26,5.1*1e12}},{"Mo",{95,42,12.66*1e12}},
		{"C",{12,6,0.35*1e12}},{"Au",{197,79,38.1*1e12}},{"Be",{9.01218,4,1.67*1e11}},
		{"BeO",{12.5056,6.55872,1.45*1e11}},{"W",{183.84,74,3.31*1e13}}}; ///< material parameters (mass, Z, sigma_abs); sigma_abs from https://www.nist.gov/pml/xcom-photon-cross-sections-database (nb. using 4GeV for KOTO)

 /// \struct AxionParameters
 /// \Brief
 /// Exotic container
 /// \EndBrief
struct AxionParameters{
	Double_t massA;
	Double_t energyA;
	Double_t pA; ///< momentum
	Double_t betaA; ///< Lorentz
	Double_t thetaA; ///< initial angle
	Double_t phiA;
	Double_t widthExpA; ///< decay width exponential
	Double_t crossSecA;
	Double_t decayLengthA;
};

#endif
