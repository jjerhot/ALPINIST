#ifndef DECAYMCGLOBAL_H
#define DECAYMCGLOBAL_H

#include <iostream>
#include <TROOT.h>

// production modes
const std::vector<TString> prodmodeNames{"primakoff","photonfrommeson","mixingPi0","mixingEta","mixingEtaPrim","BmesonK","BmesonKstar","DmesonPi"};

// decay modes
const std::vector<TString> decaymodeNames{"2Gamma","2El","2Mu","3Pi0","3Pi","2PiGamma","2Pi0Eta","2PiEta","2Pi0EtaPrim","2PiEtaPrim","2Pi","2K"};

// experiments
const std::vector<TString> expLabels{"NA62","CHARM","NuCal","SHiP","DarkQuest","DUNE","SHADOWS","KOTOdump","KOTOpnn","KOTO2dump","KOTO2pnn","KOTOexclPnn","SHiPecn3"};

// paths
const TString sourcePath = "../../tab_prod";
const TString outPath = "../../tab_decay";
const TString widthPath = "../../widths/Dalitz";

// constants
const Double_t hc = 0.197E-15; // GeV m
const Double_t MKCh = 0.4937; // GeV
const Double_t MPiCh = 0.1397; // GeV
const Double_t MMu = 0.105; // GeV
const Double_t MEl = 0.000511; // GeV
const Double_t MPi0 = 0.13498; //GeV
const Double_t MEta = 0.54786; //GeV
const Double_t MEtaPrim = 0.95778; //GeV

// exotic container
struct AxionParameters{
	Double_t massA;
	Double_t energyA;
	Double_t pA; //momentum
	Double_t betaA; //Lorentz
	Double_t thetaA; //initial angle
	Double_t phiA;
	Double_t widthExpA;
	Double_t crossSecA;
	Double_t decayLengthA;
};

#endif
