 /// \class ExpParameters
 /// \Brief
 /// Contains geometry of experiments and selection cuts.
 /// \EndBrief
 /// \Detailed
 /// Defines experiment geometry and material in the constructor.
 /// Contains methods with analysis cuts which can then be used.
 /// \EndDetailed

#include "ExpParameters.h"

//ClassImp(ExpParameters)

 /// \fn ExpParameters
 /// \Brief
 /// Base constructor with ranges
 /// \EndBrief
ExpParameters::ExpParameters(Int_t nBinsX, Int_t nBinsY, Bool_t xIsLin, Bool_t yIsLin, Int_t yVar){
	//Intitializing common values:
	nMassFiles = 2;
	minMassFile = {"01","10"};
	maxMassFile = {"10","5310"};

	nY = nBinsY; //101;
	// GeV
	Double_t min = pow(10,-25);
	Double_t max = pow(10,-10);
	if(yVar == 1){ //fs
		min = 3*pow(10,4);
		max = 3*pow(10,5);
	} else if(yVar == 2){ //m
		min = 10;
		max = 100;
	}
	if(yIsLin){
		minY = min; //pow(10,-25)
		maxY = max; //pow(10,-10)
	} else{
		minY = log10(min); //log10(1.*pow(10,-25));
		maxY = log10(max); //log10(1.*pow(10,-10));
	}
	wdY = (maxY-minY)/(nY-1);

    // ncTau = 10;
	// cTauMin = 10;
	// cTauMax = 100;
	// cTauWidth = (cTauMax-cTauMin)/(ncTau-1);

	// nWidths = 1;
	// nWidths = 101;
	// widthMin = log10(1.*pow(10,-25));
	// widthMax = log10(1.*pow(10,-10));
	// // wdWidth = (widthMax-widthMin);
	// wdWidth = (widthMax-widthMin)/(nWidths-1);

	nMX = {nBinsX, nBinsX}; //number of bins for output
	if(xIsLin){
		massXMin = {0.0001,0.5}; //min mass in massfile[i] //log scale output
		massXMax = {0.01,2}; //max mass in massfile[i] //log scale output
	} else{
		massXMin = {log10(0.0001),log10(0.01)}; //min mass in massfile[i] //log scale output
		massXMax = {log10(0.01),log10(5.31)}; //max mass in massfile[i] //log scale output
	}
	for (Int_t ifile=0; ifile < nMassFiles; ifile++)
		massXWid[ifile] = (massXMax[ifile]-massXMin[ifile])/(nMX[ifile]-1); //log scale output

	PhiBeamEuler = 0;
	ThetaBeamEuler = 0;
	PsiBeamEuler = 0;

	nRegions = 1;
}

 /// \fn ExpParameters
 /// \Brief
 /// Full constructor with experiments, decay modes and selection cuts
 /// \EndBrief
 /// \Detailed
 /// Initializes final states as particle objects given the decay mode.
 /// Defines experiment geometry and material.
 /// Calculates normalization cross section for given production mode
 /// - target variable allows choosing a different target material than for which the input table was generated (but ALP photoproduction has also non-trivial kinematic dependency).
 /// \EndDetailed
ExpParameters::ExpParameters(Int_t ExoticType_, Int_t Experiment_, TString ProductionMode_, Int_t DecayMode_, Int_t nBinsX = 101, Int_t nBinsY = 101, Bool_t xIsLin = false, Bool_t yIsLin = false, Int_t yVar = 0) : ExpParameters::ExpParameters(nBinsX, nBinsY, xIsLin, yIsLin, yVar){
	//Intializing final states
	expNum = Experiment_;
	exoNum = ExoticType_;
	decayMode = DecayMode_;
	decayModeName = decaymodeNames[ExoticType_][DecayMode_];
	switch(ExoticType_){
		case 0: //alp
			switch(DecayMode_){ // so far "2Gamma", "2El","2Mu", "3Pi0", "3Pi", "2PiGamma", "2Pi0Eta", "2PiEta", "2Pi0EtaPrim", "2PiEtaPrim",
				case 0:
					finState.push_back(new SMParticleProperty("gamma"));
					finState.push_back(new SMParticleProperty("gamma"));
					break;
				case 1:
					finState.push_back(new SMParticleProperty("el+"));
					finState.push_back(new SMParticleProperty("el-"));
					break;
				case 2:
					finState.push_back(new SMParticleProperty("mu+"));
					finState.push_back(new SMParticleProperty("mu-"));
					break;
				case 3:
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
				case 4:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
				case 5:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("gamma"));
					break;
				case 6:
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("eta"));
					break;
				case 7:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("eta"));
					break;
				case 8:
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("eta_prime"));
					break;
				case 9:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("eta_prime"));
					break;
				case 10:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					break;
				case 11:
					finState.push_back(new SMParticleProperty("kaon+"));
					finState.push_back(new SMParticleProperty("kaon-"));
					break;
				case 12:
					finState.push_back(new SMParticleProperty("kaon+"));
					finState.push_back(new SMParticleProperty("kaon-"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
				case 13:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;

			} break;
		case 1: //hnl
			switch(DecayMode_){ // so far "PiEl", "PiMu", "PiPiEl", "PiMu", "NuElEl", "NuMuMu", "NuElMu", "PiNu", "EtaNu", "PiPiNu",
				case 0:
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("el+"));
					break;
				case 1:
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("mu+"));
					break;
				case 2:
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("el+"));
					break;
				case 3:
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("mu+"));
					break;
					break;
				case 4:
					finState.push_back(new SMParticleProperty("el+"));
					finState.push_back(new SMParticleProperty("el-"));
					finState.push_back(new SMParticleProperty("nu"));
					decayModeOpen = 1;
					break;
				case 5:
					finState.push_back(new SMParticleProperty("mu+"));
					finState.push_back(new SMParticleProperty("mu-"));
					finState.push_back(new SMParticleProperty("nu"));
					decayModeOpen = 1;
					break;
				case 6:
					finState.push_back(new SMParticleProperty("el+"));
					finState.push_back(new SMParticleProperty("mu-"));
					finState.push_back(new SMParticleProperty("nu"));
					decayModeOpen = 1;
					break;
				case 7:
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("nu"));
					decayModeOpen = 1;
					break;
				case 8:
					finState.push_back(new SMParticleProperty("eta"));
					finState.push_back(new SMParticleProperty("nu"));
					decayModeOpen = 1;
					break;
				case 9:
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("nu"));
					break;
					decayModeOpen = 1;
					break;
			} break;
		case 2: //dp
			switch(DecayMode_){ // so far "2Gamma", "2El","2Mu", "3Pi0", "3Pi", "2PiGamma", "2Pi0Eta", "2PiEta", "2Pi0EtaPrim", "2PiEtaPrim",
				case 0:
					finState.push_back(new SMParticleProperty("el+"));
					finState.push_back(new SMParticleProperty("el-"));
					break;
				case 1:
					finState.push_back(new SMParticleProperty("mu+"));
					finState.push_back(new SMParticleProperty("mu-"));
					break;
				case 2:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					break;
				case 3:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
				case 4:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					break;
				case 5:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
				case 6:
					finState.push_back(new SMParticleProperty("kaon+"));
					finState.push_back(new SMParticleProperty("kaon-"));
					break;
				case 7:
					finState.push_back(new SMParticleProperty("kaon+"));
					finState.push_back(new SMParticleProperty("kaon-"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
			} break;
		case 3: //ds
			switch(DecayMode_){ // so far "2Gamma", "2El","2Mu", "3Pi0", "3Pi", "2PiGamma", "2Pi0Eta", "2PiEta", "2Pi0EtaPrim", "2PiEtaPrim",
				case 0:
					finState.push_back(new SMParticleProperty("el+"));
					finState.push_back(new SMParticleProperty("el-"));
					break;
				case 1:
					finState.push_back(new SMParticleProperty("mu+"));
					finState.push_back(new SMParticleProperty("mu-"));
					break;
				case 2:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					break;
				case 3:
					finState.push_back(new SMParticleProperty("kaon+"));
					finState.push_back(new SMParticleProperty("kaon-"));
					break;
				case 4:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					break;
				case 5:
					finState.push_back(new SMParticleProperty("pi+"));
					finState.push_back(new SMParticleProperty("pi-"));
					finState.push_back(new SMParticleProperty("pi0"));
					finState.push_back(new SMParticleProperty("pi0"));
					break;
			} break;
	}

	// defining experiments
	switch(Experiment_){
		case 0 : //CHARM
			expLabel = expLabels[Experiment_];
			expName = "CHARM";
			beamEnergy = 400;
			target = targetMat.at("Cu");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.;  // m
			Y0 = -4.8;// m since the angle is ~10 mrad, but the detector is parallel to the beam line
			Z0 = 0.; //was 20 m
			ZFVIn = 480.; //m
			ZECal = 515.; //m
			ZFVEnd = ZECal;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal;
			ZMuonDetector = 515.; //  dummy, same as "ECal"

			Sigacceptance = 0.51; // for photons
			Sigacceptancemumu  = 0.85; // fof muons TO CHECK !!!

			POT = 2.4*1.E18;

			thetaMinInFile = 0.0025; 
			thetaMaxInFile = 0.0175; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 1 ://"BEBC" WA66
			expLabel = expLabels[Experiment_];
			expName = "BEBC";
			beamEnergy = 400;
			target = targetMat.at("Cu");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.;  // m
			Y0 = 0;// 
			Z0 = 0.; //
			ZECal = 405.9; //m 
			ZFVIn = ZECal-1.25; //m 404.
			ZFVEnd = ZFVIn+1.25;
			ZMuonDetector = ZECal+7.665; // EMI is 26x6 m semicircle at 7.665 m from bubble chamber center  
			Sigacceptance = 0.88; 
			Sigacceptancemumu  = 0.89;

			Magnets.push_back(Magnet(ZECal, 1.85, 3.5));
			POT = 2.72*1.E18;

			thetaMinInFile = 0.00017; 
			thetaMaxInFile = 0.00972; 
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); 
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 2://NuCal
			expLabel = expLabels[Experiment_];
			expName = "NuCal";
			beamEnergy = 70;
			target = targetMat.at("Fe");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.;  // m
			Y0 = 0.;  // m
			Z0 = ZTarget; //was 20 m
			ZFVIn = 64.; //m
			ZECal = 23.+64.;
			ZFVEnd = ZECal;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal;
			ZMuonDetector = 87.;

			Sigacceptance =0.7; // see Bluemlein for photons https://arxiv.org/pdf/1311.3870.pdf
			Sigacceptancemumu =0.8; // see Bluemlein for muons https://arxiv.org/pdf/1311.3870.pdf

			POT = 1.7*1.E18;

			thetaMinInFile = 0.00025; 
			thetaMaxInFile = 0.01525; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 1.1;
			energyMaxInFile = 64.9;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 3://NuTeV
			expLabel = expLabels[Experiment_];
			expName = "NuTeV";
			beamEnergy = 800;
			target = targetMat.at("BeO");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; 
			Y0 = 0.; 
			Z0 = 0.; 
			ZFVIn = 1400.; //m
			ZFVEnd = ZFVIn + 34.;
			ZTracker1 = ZFVIn + 35.5; // m, position of 1st tracker chamber   8
			ZTracker4 = ZTracker1 + 38.5;  // m, position of 4th tracker chamber
			ZECal = ZFVIn + 39.;
			// if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; // neutral decay (we allow decay downstream of spectrometer), comparing with NA62MC
			ZMuonDetector = ZECal + 17.6946 + 1.5*2.8973; // 6*(28*5.15+14*6.48+7*8.37)/100

			Sigacceptance = 0.23; // for photons
			Sigacceptancemumu = 0.23; // for muons 
			//these are the muon veto magnets: Magnets =  {Magnet(ZTracker1 + 17.6946 + 0.5*2.8973, 1.4, 1.9, (Bool_t)1, 4.86 ), Magnet(ZTracker1 + 17.6946 + 1.5*2.8973, 1.4, -1.9, (Bool_t)1,4.86 ), Magnet(ZTracker1 + 17.6946 + 2.5*2.8973, 1.4, 1.9, (Bool_t)1, 4.86 )}; // eff field radius = 4.86=  1.7/(1.9-1.55),  initial tracking: 6*(28*52.+ 14*64.8+7*83.7) = 17.6946 magnet unit thicknes = 9*83.7 + 8*64.8 + 8*203.2 = 2.8973

			POT = 2.54*1.E18; 

			thetaMinInFile =  0.0000295; 
			thetaMaxInFile =  0.00179; 
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); 
			energyMinInFile =  12.;
			energyMaxInFile = 650.;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 4 : //NA62-bd
			expLabel = expLabels[Experiment_];
			expName = "NA62";
			beamEnergy = 400;
			target = targetMat.at("Cu");
			refTarget = target;
			ZTarget = 23.070; //XTAX101024
			X0 = 0.; // m
			Y0 = -0.022; // m
			Z0 = ZTarget; //was 20 m
			ZFVIn = 102; //102.425; // 105, comparing with NA62MC
			ZTracker1 = 183. ; // m, position of 1st tracker chamber
			ZTracker4 = 218.7; // m, position of 4th tracker chamber
			// ZMagnet = 0.5*(196.350+197.650);
			// MagnetZLength       = 1.3; // m
			// MagnetFieldStrength = 0.6928; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV

			Magnets.push_back(Magnet(0.5*(196.350+197.650), 1.3, 0.6928));
			ZTrigger = 238.898;
			ZECal = 241.495;
			// ZFVEnd = ZTracker1;
			ZFVEnd = 182.;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; // neutral decay (we allow decay downstream of spectrometer), comparing with NA62MC
			ZMuonDetector = 246.850;

			Sigacceptance = 1; // better call that "efficiency?"
			Sigacceptancemumu = 1;

			acceptanceHole = 0.1 ; // m, half the hole length
			acceptanceSide = 1.; // m
			//POT = 1.3*1.E16; // 1 day
			POT = 1.*1.E18; //
			// POT = 1.*1.E19; //2026-2030 period (possible x4 update beyond)

			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 5: //DarkQuestPhase I
			expLabel = expLabels[Experiment_]; // DarkQuest as in 1804.00661v1
			expName = "DarkQuest";
			beamEnergy = 120;
			target = targetMat.at("Fe");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; // m
			Y0 = 0.; // m straight for mathematica comparison
			Z0 = 0.; //was 20 m
			ZFVIn = 5.; //m it used to be 7 according to text of 1804.00661v1, but now it is 5 accoring to 2008.08108
			ZTracker1 = 6.16 ; // m, position of 1st tracker of DarkQuest (x-view), https://www.sciencedirect.com/science/article/pii/S016890021930347X
			ZTracker4 = 18.79; // position of 3nd tracker x-view,
			// ZMagnet = 9;// // follow fig 1 in 2008.08108, check
			// MagnetZLength       = 3.; // m //
			// MagnetFieldStrength = 0.4; // Tesla, see https://www.sciencedirect.com/science/article/pii/S016890021930347X
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; //

			Magnets.push_back(Magnet(9, 3., 0.4));
			
			// assume that the second tracker is also hit when the third is hit, so 2nd is not included explcitye
			// 13.47; // m, position of 2nd tracker x-view, https://www.sciencedirect.com/science/article/pii/S016890021930347X
			ZECal = 18.5;// 25.; //m (total length is 25 m but Calo needs to go before wall)
			// THIS IS PHASE 1, in phase 2 the fiducial region would be larger
			ZFVEnd = ZTracker1; // it used to be 8, see fig 1 in 1804.00661v1 [hep-ph], but now it is 6 accoring to 2008.08108
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; /// gamma gamma
			ZMuonDetector = 21.; // according to 2008.08108

			Sigacceptance = 1;
			Sigacceptancemumu = 0.6; // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"

			POT = 1.44*1.E18; // 1.E20; ; // phase1, phase 2, cf caption of fig 10 in 2008.08108
			//POT = 2.144E18; // phase1

			thetaMinInFile = 0.0009; 
			thetaMaxInFile = 0.0549; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 1.9;
			energyMaxInFile = 112.1;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 6 : // DarkQuest PhaseII as proposed in [2008.08108]
			expLabel = expLabels[Experiment_];
			expLabel = expLabels[Experiment_]; // DarkQuest as in 1804.00661v1
			expName = "DarkQuest";
			beamEnergy = 120;
			target = targetMat.at("Fe");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; // m
			Y0 = 0.; // m 
			Z0 = 0.; // m
			ZFVIn = 7.; // m according to 2008.08108
			ZFVEnd = 12.; // m according to 2008.08108
			ZTracker1 = 13.5 ; // m, position of 2nd tracker of DarkQuest (x-view), https://www.sciencedirect.com/science/article/pii/S016890021930347X
			ZTracker4 = 18.79; // position of 3nd tracker x-view,
			Magnets.push_back(Magnet(9, 3., 0.4)); // MagnetZLength = 3 m  MagnetFieldStrength = 0.4 Tsee https://www.sciencedirect.com/science/article/pii/S016890021930347X
			// assume that the second tracker is also hit when the third is hit, so 2nd is not included explcitye
			// 13.47; // m, position of 2nd tracker x-view, https://www.sciencedirect.com/science/article/pii/S016890021930347X
			ZECal = 19;// 25.; //m (total length is 25 m but Calo needs to go before wall)
			// THIS IS PHASE 1, in phase 2 the fiducial region would be larger
			ZFVEnd = ZTracker1; // it used to be 8, see fig 1 in 1804.00661v1 [hep-ph], but now it is 6 accoring to 2008.08108
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; /// gamma gamma
			ZMuonDetector = 21.; // according to 2008.08108

			Sigacceptance = 1;
			Sigacceptancemumu = 0.6; // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"

			POT = 1.E20; // 1.E20; ; // phase1, phase 2, cf caption of fig 10 in 2008.08108
			//POT = 2.144E18; // phase1

			thetaMinInFile = 0.0009; 
			thetaMaxInFile = 0.0549; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 1.9;
			energyMaxInFile = 112.1;
			break;
		case 7: // based on https://arxiv.org/pdf/2103.13910.pdf and https://arxiv.org/pdf/2011.05995.pdf
			expLabel = expLabels[Experiment_];
			expName = "DUNE";
			beamEnergy = 120;
			target = targetMat.at("C");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.;
			Y0 = 0.; //by default but detector can move off-axis up to 30.5m (see DUNE-PRISM)
			Z0 = 0.;
			ZFVIn = 574.5; // front of LAr
			ZECal = 580;// front of GAr
			ZFVEnd = 584.; // end of GAr
			ZMuonDetector = ZFVEnd;

			Sigacceptance = 1;
			Sigacceptancemumu =1; //??

			POT = 1.1*1.E21; //1 year of data taking
//				POT = 1.1*1.E22; //10 years of data taking?
			break;
		case 8: // SHiP@ecn3 as in CERN-SPSC-2023-033
			expLabel = expLabels[Experiment_];
			expName = "SHiP";
			beamEnergy = 400;
			target = targetMat.at("Mo");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; // m
			Y0 = 0.; // m straight for mathematica comparison
			Z0 = 0.; //was 20 m
			ZFVIn = 33.5; //m
			ZFVEnd = ZFVIn+50;
			ZTracker1 = ZFVEnd+0.55; // m, position of 1st tracker chamber
			ZTracker4 = ZTracker1+9.55; // m, position of 4th tracker chamber
			// ZMagnet = ZFVEnd + 5.6;//
			// MagnetZLength       = 4;
			// MagnetFieldStrength = 0.1625; // Tesla			
			// MagnetFieldPhi = 0;
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; //
			Magnets.push_back(Magnet(ZFVEnd + 5.6, 4., 0.1625, 0.));
			ZECal = ZTracker4 + 0.55; //m 
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			// if(finState.front()->GetCharge() == 0)	ZFVEnd = ZTracker1;
			ZMuonDetector = ZFVEnd + ZECal + 3.5; // not condition

			acceptanceHole = 0. ; // m, half the hole length
			acceptanceSide = 6; // too big. In reality it's egg-shaped but the inCaloAcceptance is correct for ship

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 6.E20;

			thetaMinInFile = 0.00078; 
			thetaMaxInFile = 0.04758; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		
		case 9:
			// based on presentation https://indico.cern.ch/event/1055867/contributions/4437601/attachments/2280165/3874009/darkSectorKOTOv1.pdf
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			target = targetMat.at("Au");
			refTarget = target;
			ZTarget = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56; //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTarget;
			ZFVIn = ZTarget + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet = 0.5*(ZTracker1+ZTracker4); // no tracker
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			// Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 0., 0.));
			ZECal = 27;
			ZFVEnd = ZECal; //20(+3?) m decay volume for charged particles
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal + 1.; //summy

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 2.2*1.E17; //slide 13 in https://indico.cern.ch/event/1055867/contributions/4437601/attachments/2280165/3874009/darkSectorKOTOv1.pdf

			thetaMinInFile = 0.06; // factor 5 w.r.t. shadows
			thetaMaxInFile = 0.6; // factor 5 w.r.t. shadows
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 0.5;
			energyMaxInFile = 29.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 10:
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			target = targetMat.at("Au");
			refTarget = target;
			ZTarget = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56;  //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTarget;
			ZFVIn = ZTarget + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet = 0.5*(ZTracker1+ZTracker4); // no tracker
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			// Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 0., 0.));
			ZECal = 27;
			ZFVEnd = ZECal; //20(+3?) m decay volume for charged particles
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal + 1.; //summy

			Sigacceptance = 0.52;
			Sigacceptancemumu = 1;

			POT = 2.2*1.E19;

			thetaMinInFile = 0.06; // factor 5 w.r.t. shadows
			thetaMaxInFile = 0.6; // factor 5 w.r.t. shadows
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 0.5;
			energyMaxInFile = 29.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

			nRegions = 8;
			break;
		case 11:
			expLabel = expLabels[Experiment_];
			expName = "KOTO2";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.087266; //5∘ = 87.266 mrad
			PsiBeamEuler = 0;
			target = targetMat.at("Au");
			refTarget = target;
			ZTarget = 0; //presume that they dumped on T1 target
			X0 = 0.;//
			Y0 = 0.;
			Z0 = ZTarget;
			ZFVIn = ZTarget + 43 + 1 + 1.75; // 3-m long beam line and 1-m long space; The Front Barrel Counter is 1.75-m
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet = 0.5*(ZTracker1+ZTracker4); // no tracker
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			// Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 0., 0.));
			ZECal = 44 + 20;
			ZFVEnd = ZECal; 
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal; //summy

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 2.2*1.E17; //assuming they take the same statistics as in KOTOdump
			thetaMinInFile = 0.06;
			thetaMaxInFile = 0.12; 
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1);
			energyMinInFile = 0.5;
			energyMaxInFile = 29.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 12:
			expLabel = expLabels[Experiment_];
			expName = "KOTO2";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.087266; //5∘ = 87.266 mrad
			PsiBeamEuler = 0;
			target = targetMat.at("Au");
			refTarget = target;
			ZTarget = 0; //presume that they dumped on T1 target
			X0 = 0.;//
			Y0 = 0.;
			Z0 = ZTarget;
			ZFVIn = ZTarget + 43 + 1 + 1.75; // 3-m long beam line and 1-m long space; The Front Barrel Counter is 1.75-m
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet = 0.5*(ZTracker1+ZTracker4); // no tracker
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			// Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 0., 0.));
			ZECal = 44 + 20;
			ZFVEnd = ZECal; 
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal; //summy

			Sigacceptance = 0.73; // for photon clusters to reduce neutron clusters by applying cluster shape, pulse shape, and depth information of the hits in the calorimeter: p234 https://arxiv.org/pdf/2110.04462.pdf
			Sigacceptancemumu = 1;

			POT = 6.*1E20; //table 27: 3*1E7 s * 2.*1E13 POT in https://arxiv.org/pdf/2110.04462.pdf
			thetaMinInFile = 0.06;
			thetaMaxInFile = 0.12; 
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1);
			energyMinInFile = 0.5;
			energyMaxInFile = 29.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

			nRegions = 1; //we will work with one signal region for now (pT>400MeV/c and zVtx < 15m, with some expected 1.5 background events; possible also pT>350 MeV/c with 3 background events)
			break;
		case 13:
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			target = targetMat.at("Au");
			refTarget = target;
			ZTarget = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56;  //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTarget;
			ZFVIn = ZTarget + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet = 0.5*(ZTracker1+ZTracker4); // no tracker
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			// Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 0., 0.));
			ZECal = 27;
			ZFVEnd = ZECal; //20(+3?) m decay volume for charged particles
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal + 1.; //summy

			Sigacceptance = 0.9; // p248 KL halo selection: we also achieved to reduce it to be ×(2.1 × 10−2) with 90% signal efficiency with small correlation to the cluster and pulse shape cuts.
			Sigacceptancemumu = 1;

			POT = 14*1.E19; //assuming 

			thetaMinInFile = 0.06; // factor 5 w.r.t. shadows
			thetaMaxInFile = 0.6; // factor 5 w.r.t. shadows
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 0.5;
			energyMaxInFile = 29.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

			nRegions = 1;
			break;
		case 14 :
			expLabel = expLabels[Experiment_];
			expName = "NA62";
			beamEnergy = 400;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.008; 
			PsiBeamEuler = 0;
			target = targetMat.at("Be"); // possibility of Be, Cu or W
			refTarget = targetMat.at("Cu");// production yields are calculated assuming Cu
			ZTarget = 0.;
			X0 = 0.;  // m
			Y0 = 0.;  // m
			Z0 = ZTarget; //was 20 m
			ZFVIn = 120; //given by the position of upstream veto
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet     = 0;
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			ZECal = 241.495;
			ZFVEnd = ZTracker1;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; // neutral decay (we allow decay downstream of spectrometer), comparing with NA62MC
			ZMuonDetector = 246.850;

			Sigacceptance = 1; // better call that "efficiency?"
			Sigacceptancemumu = 1;

			acceptanceHole = 0.12 ; // can be extended to 15cm
			acceptanceSide = 1.25; // m
			POT = 6.*1.E19; 

			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 15 :
			expLabel = expLabels[Experiment_];
			expName = "NA62";
			beamEnergy = 400;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.008; 
			PsiBeamEuler = 0;
			target = targetMat.at("Be"); // possibility of Be, Cu or W
			refTarget = targetMat.at("Cu"); // production yields are calculated assuming Cu
			ZTarget = 0.;
			X0 = 0.;  // m
			Y0 = 0.;  // m
			Z0 = ZTarget; //was 20 m
			ZFVIn = 120+150.; //given by the position of upstream veto (extended by 150m)
			ZTracker1 = 0; // no tracker
			ZTracker4 = 0; // no tracker
			// ZMagnet     = 0;
			// MagnetZLength       = 0; // m
			// MagnetFieldStrength = 0; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV
			ZECal = 241.495+150.;
			ZFVEnd = ZTracker1;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; // neutral decay (we allow decay downstream of spectrometer), comparing with NA62MC
			ZMuonDetector = 246.850+150.;

			Sigacceptance = 1; // better call that "efficiency?"
			Sigacceptancemumu = 1;

			acceptanceHole = 0.12 ; // can be extended to 15cm
			acceptanceSide = 1.25; // m
			POT = 6.*1.E19; 

			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		
		case 16: //based on PBC proposal https://indico.cern.ch/event/1002356/contributions/4229628/attachments/2201710/3781990/shadows_PBC.pdf
			expLabel = expLabels[Experiment_];
			expName = "SHADOWS";
			beamEnergy = 400;
			PhiBeamEuler = TMath::Pi()/2; //SHADOWS is offaxis in X axis 
			// ThetaBeamEuler = 0.0624188; //(~62mrad) for calorimeter located 36m downstream the target and 2.25m shifted center w.r.t. beam axis
			PsiBeamEuler = 0;
			target = targetMat.at("Cu");
			refTarget = target;
			ZTarget = 23.070; //NA62
			X0 = 2.5;//(1+1.5);  //shifted 1m wrt NA62, 3m calo size (also option with 3m size in consideration)
			Y0 = -0.022; //as for NA62
			Z0 = ZTarget; //NA62
			ZFVIn = ZTarget + 10.; //10m from target
			ZTracker1 = ZFVIn + 20. + 6. - 3. ; // m, position of 1st tracker chamber (decay volume end + 6m Calo location - 3m spectrometer length)
			ZTracker4 = ZTracker1 + 2.5; // m, position of 4th tracker chamber (0.5m between 1,2 and 3,4 station and 1.5m magnet)
			// ZMagnet = 0.5*(ZTracker1+ZTracker4);
			// MagnetZLength       = 1.5; // m
			// MagnetFieldStrength = 1.; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; // GeV

			Magnets.push_back(Magnet(0.5*(ZTracker1+ZTracker4), 1.5, 1.));

			ZECal = ZTracker4 + 0.5; //~50cm behind 4th tracker station
			ZFVEnd = ZTracker1; //20(+3?) m decay volume for charged particles
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal; //for photons
			ZMuonDetector = ZECal + 1.; //1m distance between calo and muon detector as measured from the image in the SHADOWS paper

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 1.*1.E19; //2026-2030 period (possible x4 update beyond)

			thetaMinInFile = 0.012;
			thetaMaxInFile = 0.12;
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); // used to be log scale
			break;
		case 17: // SHiP@ecn4
			expLabel = expLabels[Experiment_];
			expName = "SHiP";
			beamEnergy = 400;
			target = targetMat.at("Mo");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; // m
			Y0 = 0.; // m straight for mathematica comparison
			Z0 = 0.; //was 20 m
			ZFVIn = 45.; //m
			ZTracker1 = ZFVIn+50.76 ; // m, position of 1st tracker chamber
			ZTracker4 = ZFVIn+60.16; // m, position of 4th tracker chamber
			// ZMagnet = ZFVIn + 53.45;//
			// MagnetZLength       = 5.; // m // in figure 1b https://cds.cern.ch/record/2005715/files/main.pdf they also state 0.75 Tm
			// MagnetFieldStrength = 0.15; // Tesla
			// MagKick = MagnetFieldStrength*MagnetZLength*c*1E-9; //

			Magnets.push_back(Magnet(ZFVIn + 53.45, 5., 0.15));
			
			ZECal = ZTracker4 + 0.5; //m
			ZFVEnd = ZTracker1;
			if(finState.front()->GetCharge() == 0)	ZFVEnd = ZECal;
			ZMuonDetector = 1E10; // not condition

			acceptanceHole = 0. ; // m, half the hole length
			acceptanceSide = 5; // too big. In reality it's egg-shaped but the inCaloAcceptance is correct for ship

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 1.E20;

			thetaMinInFile = 0.00078; 
			thetaMaxInFile = 0.04758; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 18:
			expLabel = expLabels[Experiment_];
			expName = "ORCA";
			beamEnergy = 400;
			target = targetMat.at("C");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.; 
			Y0 = 0.; 
			Z0 = 0.; 
			ZFVIn = 1700E3; //m
			ZFVEnd = ZFVIn + 210.;
			ZTracker1 = ZFVIn + 112.5; // m, point of maximum horizontal extend
			ZTracker4 = ZTracker1;       // m, point of maximum horizontal extend
			ZECal = ZTracker4; //triggering only on CHERENKOV 
			ZMuonDetector = 0;  //no distinction for muons

			Sigacceptance = 0.25; // for photons (maximum quantum effiency of photomultipliers)
			Sigacceptancemumu = 0.25; // for muons 
			POT = 8.E19; 

			thetaMinInFile = 1.5E-6; 
			thetaMaxInFile = 90E-6; 
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); 
			energyMinInFile = 6.3;
			energyMaxInFile = 371.7;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;
		case 19 : // BEBC with a asimplified geometry as proposed in [2208.00416]
			expLabel = expLabels[Experiment_];
			expName = "BEBC";
			beamEnergy = 400;
			target = targetMat.at("Cu");
			refTarget = target;
			ZTarget = 0;
			X0 = 0.;  // m
			Y0 = 0;// 
			Z0 = 0.; //
			ZFVIn = 404.; //m 
			ZECal = 405.92; //m /Chamber is 3.7 m tall and contains 19 m^3l -> could be approximating it as  2.27x2.27x3.7 cuboid,
			ZFVEnd = ZFVIn+1.85;
			ZMuonDetector = ZECal+7.665; // EMI is 26x6 m semicircle at 7.665 m from bubble chamber center  
			Sigacceptance = 0.88; 
			Sigacceptancemumu  = 0.89;

			POT = 2.72*1.E18;

			thetaMinInFile = 0.00017; 
			thetaMaxInFile = 0.00972; 
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); 
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale
			break;

		default:
			std::cout << "Invalid experiment label: " <<  Experiment_ << ". Use values ";
			for(Int_t iExp=0;iExp<expLabels.size();iExp++) std::cout << iExp << " (" << expLabels[iExp] << "), ";
			std::cout << std::endl;
			exit(1);

	}

	// calculating normalization cross section
	normCrossSec = 1.;
	
	switch(ExoticType_){
		case 0: //alp
			if(refTarget!=target){
				if(ProductionMode_ == "Bmeson") normCrossSec = TMath::Power(target[0]/refTarget[0], 1./3); 
				else if(ProductionMode_ == "Dmeson") normCrossSec = TMath::Power(target[0]/refTarget[0], 1./3); 
				else if(ProductionMode_ == "Primakoff"){
					normCrossSec = TMath::Power(target[0]/refTarget[0], 0.77);
					std::cout<<"[Warning] production mode " << ProductionMode_<< " cannot be trivially rescaled from reference please ensure to use an input distribution calculated for the desired target material." << std::endl;
				}
				else if(ProductionMode_ == "PhotonFromMeson"){
					normCrossSec = target[2]/refTarget[2];
					std::cout<<"[Warning] production mode " << ProductionMode_<< " cannot be trivially rescaled from reference please ensure to use an input distribution calculated for the desired target material." << std::endl;
				}
				else if (ProductionMode_ != "MixingPi0" && ProductionMode_ != "MixingEta" && ProductionMode_ != "MixingEtaPrim") 
					std::cout<<"[Warning] production mode " << ProductionMode_<< " not recognized for exo type alp. Defaulting normCrossSec to 1." << std::endl;
			} break;
		case 1: //hnl
			if(refTarget!=target){
				if(ProductionMode_ == "Bmeson") normCrossSec = TMath::Power(target[0]/refTarget[0], 1./3);
				else if(ProductionMode_ == "Dmeson") normCrossSec = TMath::Power(target[0]/refTarget[0], 1./3);
				else std::cout<<"[Warning] production mode " << ProductionMode_ << " not recognized for exo type hnl. Defaulting normCrossSec to 1." << std::endl;
			} break;
		case 2: //dp
			if (ProductionMode_ != "Brems" && ProductionMode_ != "MesonDecay" && ProductionMode_ != "MixingRho" && ProductionMode_ != "MixingOmega" && ProductionMode_ != "MixingPhi") 
				std::cout<<"[Warning] production mode " << ProductionMode_<< " not recognized for exo type dp. Defaulting normCrossSec to 1." << std::endl;
			break;
		case 3: //ds (B meson + Brems production)
			if (ProductionMode_ == "Bmeson" || ProductionMode_ == "Bmeson2S") 
				if(refTarget!=target) normCrossSec = TMath::Power(target[0]/refTarget[0], 1./3);
			if ( ProductionMode_ != "Brems" &&  ProductionMode_ != "Bmeson" && ProductionMode_ != "Bmeson2S") 
				std::cout<<"[Warning] production mode " << ProductionMode_<< " not recognized for exo type ds. Defaulting normCrossSec to 1." << std::endl;
			break;
	}
	ZDist = ZFVIn-ZTarget;
	BeamDecayLength = ZFVEnd - ZFVIn; // m
}

 /// \fn ~ExpParameters
 /// \Brief
 /// Destroy particle final states
 /// \EndBrief
ExpParameters::~ExpParameters(){
    while(!finState.empty()) {
        delete finState.back();
        finState.pop_back();
    }
	finState.clear();
}

 /// \fn inCaloAcceptance
 /// \Brief
 /// Checks if particle is in ECal acceptance
 /// \EndBrief
Bool_t ExpParameters::inCaloAcceptance(TVector3 posFVEnd, TVector3 posECal, Double_t zDecay){
	Double_t xa = fabs(posECal.X()*1000.); //in mm
	Double_t ya = fabs(posECal.Y()*1000.); //in mm

	if (expNum == 0) { // this works since you offset the proton beam
		if (xa > 1500.) return kFALSE; //for CHARM 3 m by 3 m
		if (ya > 1500.) return kFALSE; //for CHARM 3 m by 3 m
	}
	else if (expNum == 1){ // BEBC
		Double_t delta_z = TMath::Abs(ZECal - zDecay) *1000, r = 1250, h = 3500;

		if(r < delta_z || (h/2 < ya) || (r*r - delta_z*delta_z < xa*xa) ) return kFALSE;
	}
	else if (expNum == 2) { //nucal
		if(pow(xa*xa+ya*ya,0.5)>1300) return kFALSE;
	}
	else if (expNum == 3){// NuTeV
		if(1500 < xa || 1500 < ya) return kFALSE; // 3x3m scintilator planes
	}
	else if (expNum == 4) {
		if (xa > 1130)return kFALSE;
		if (ya > 1130)return kFALSE;
		if (xa+ya > 1598) return kFALSE;
		if (xa > 632 && ya > 847) return kFALSE; // !83.7?
		if (ya > 947 && xa > 522) return kFALSE;
		if(pow(xa*xa+ya*ya,0.5)<80) return kFALSE;
		// fECalDetectorOuterRadius = 1.825 * m
		// fRadiusSteelFoilColdWindow = 1.465 * m;
	}
	else if (expNum == 5 || expNum == 6) {
		if (xa > 1000.) return kFALSE; // 2 m by
		if (ya > 1000.) return kFALSE; // 2 m
	}
	else if (expNum == 7) { //DUNE
		if((zDecay > ZFVIn) && (zDecay < ZECal )){ //LAr
			if (abs(posECal.X()) > 3.) return kFALSE;
			if (abs(posECal.Y()) > 1.) return kFALSE;
		} else if((zDecay > ZECal) && (zDecay < ZFVEnd)) { //DUNE GAr: 2x4.8 m
			if (abs(posFVEnd.X()) > 2.4) return kFALSE;
			if (abs(posFVEnd.Y()) > 1.) return kFALSE;
		} else return kFALSE;
	}
	else if (expNum == 8) { // ship ecn3
		if (2590. < xa || 3310. < ya) return kFALSE; //ECN3 callorimeter is rectangle with semi axes x 2590 mm, y 3310 mm
	}
	else if (expNum == 9 || expNum == 10 || expNum == 13) { //KOTO
		if(pow(xa*xa+ya*ya,0.5)>1000) return kFALSE; // The CsI calorimeter consists of 2716 undoped CsI crystals stacked in a cylindrical shape of 2-m diameter
		if (abs(xa) < 75. && abs(ya) < 75.) return kFALSE; //A 15 x 15 cm2 center region of the stacked CSI permitted the beam particles to pass through the detector.
	}
	else if (expNum == 11 || expNum == 12) { //KOTO2
		if(pow(xa*xa+ya*ya,0.5)>1500) return kFALSE; //3m diameter
		if (abs(xa) < 100. && abs(ya) < 100.) return kFALSE; //20x20cm central hole
	}
	else if (expNum == 14 || expNum == 15) { //KLEVER
		if(pow(xa*xa+ya*ya,0.5)>1250) return kFALSE; //1.25m outer radius
		if(pow(xa*xa+ya*ya,0.5)<120) return kFALSE; //20x20cm central hole
	}
	else if (expNum == 16) { //SHADOWS
		if (xa > 1250.) return kFALSE; // 2.5 m by
		if (ya > 1250.) return kFALSE; // 2.5 m
	}
	else if (expNum == 17) { // shipecn4 // new setup has no ellipse anymore
		// Double_t xellip = 2500;
		// Double_t yellip = 5000;
		// Double_t rad= pow(xa/xellip,2) +pow(ya/yellip,2);
		// if (rad > 1.) return kFALSE; //for ship
		if (xa > 5000.) return kFALSE; // 10 m by
		if (ya > 2500.) return kFALSE; // 5 m

	}
	else if (expNum == 18){ // ORCA
		Double_t delta_z = TMath::Abs(ZECal - zDecay)*1000, r = 112500, h = 200000;
		if(r < delta_z || (h/2 < ya) || (r*r - delta_z*delta_z < xa*xa) ) return kFALSE;
	}
	else if (expNum == 19){ // BEBCcuboid
		if(1260. < xa || 1785. < ya  ) return kFALSE;
	}
	return kTRUE;
}

 /// \fn inTriggerAcceptance
 /// \Brief
 /// Trivial acceptance condition on NA62
 /// \EndBrief
Bool_t ExpParameters::inTriggerAcceptance(TVector3 posTrigger){

	Double_t xa = fabs(posTrigger.X()*1000.); //in mm
	Double_t ya = fabs(posTrigger.Y()*1000.); //in mm

	if (expNum == 4) {
		if(pow(xa*xa+ya*ya,0.5)>1210.) return kFALSE;
		// if(pow(xa*xa+ya*ya,0.5)<128.) return kFALSE;
	}
	return kTRUE;
}

 /// \fn inCaloOuterEdgeAcceptance
 /// \Brief
 /// Check for outer edge of calorimeter for more efficient simulation
 /// \EndBrief
Bool_t ExpParameters::inCaloOuterEdgeAcceptance(TVector3 posECal){
	Double_t xa = fabs(posECal.X()*1000.); //in mm
	Double_t ya = fabs(posECal.Y()*1000.); //in mm

	if (expNum == 0) { // this works since you offset the proton beam
		if (xa > 1500.) return kFALSE; //for CHARM 3 m by 3 m
		if (ya > 1500.) return kFALSE; //for CHARM 3 m by 3 m
	}
	else if (expNum == 1){ // BEBC r = 1850, h = 4000;
		if((1850<xa) || (2000<ya))	return kFALSE;
	}
	else if (expNum == 2) { //nucal
		if(pow(xa*xa+ya*ya,0.5)>1300) return kFALSE;
	}
	else if (expNum == 3) { //NuTeV
		if(1525<xa || 1525<ya) return kFALSE; //3.05m quadratic steel plates
	}
	else if (expNum == 4) {
		// if (xa > 1130)return kFALSE;
		// if (ya > 1130)return kFALSE;
		// if (xa+ya > 1598) return kFALSE;
		// if (xa > 632 && ya > 847) return kFALSE; // !83.7?
		// if (ya > 947 && xa > 522) return kFALSE;
		if(pow(xa*xa+ya*ya,0.5)>1465) return kFALSE; //fRadiusSteelFoilColdWindow
	}
	else if (expNum == 5 || expNum == 6) {
		if (xa > 1000.) return kFALSE; // 2 m by
		if (ya > 1000.) return kFALSE; // 2 m
	}
	else if (expNum == 7) { //DUNE
		if (abs(posECal.X()) > 3.) return kFALSE;
		if (abs(posECal.Y()) > 1.) return kFALSE;
	}
	else if (expNum == 9 || expNum == 10 || expNum == 13) { //KOTO
		if(pow(xa*xa+ya*ya,0.5)>1000) return kFALSE; // The CsI calorimeter consists of 2716 undoped CsI crystals stacked in a cylindrical shape of 2-m diameter
	}
	else if (expNum == 11 || expNum == 12) { //KOTO2
		if(pow(xa*xa+ya*ya,0.5)>1500) return kFALSE; //3m diameter
	}
	else if (expNum == 8) { // ship ecn3
		if (2650. < xa || 3350. < ya) return kFALSE; //ECN3 callorimeter is rectangle with semi axes x 2590 mm, y 3310 mm
	}
	else if (expNum == 14 || expNum == 15) { //KLEVER
		if(pow(xa*xa+ya*ya,0.5)>1250) return kFALSE; //1.25m outer radius
	}
	else if (expNum == 16) { //SHADOWS
		if (xa > 1250.) return kFALSE; // 2.5 m by
		if (ya > 1250.) return kFALSE; // 2.5 m
	}
	else if (expNum == 17) { // shipecn4 // new setup has no ellipse anymore
		// Double_t xellip = 2500;
		// Double_t yellip = 5000;
		// Double_t rad= pow(xa/xellip,2) +pow(ya/yellip,2);
		// if (rad > 1.) return kFALSE; //for ship
		if (xa > 5000.) return kFALSE; // 10 m by
		if (ya > 2500.) return kFALSE; // 5 m

	}
	else if (expNum == 18) {// ORCA
		if((112500.<xa) || (100000.<ya))	return kFALSE;
	}
	else if (expNum == 19){ // BEBCcuboid
		if(1260 < xa || 1785 < ya ) return kFALSE;
	}
	return kTRUE;
}

 /// \fn inSpectrometerAcceptance
 /// \Brief
 /// Check for spectrometer stations acceptance
 /// \EndBrief
Bool_t ExpParameters::inSpectrometerAcceptance(TVector3 posTracker){
	Double_t xa = fabs(posTracker.X()), ya = fabs(posTracker.Y());
	Double_t trackerRadius =pow(xa*xa + ya*ya, 0.5);
	
	if(expNum == 0) {return kFALSE; }//no trackers
	else if (expNum == 1){ // BEBC r = 1250, h = 3500
		if((xa<1.25) && (ya<1.75))	return kTRUE;}
	else if(expNum == 2) {return kFALSE;} //no trackers
	else if (expNum == 3){ // NuTeV
		if((xa<1.5) && (ya<1.5))	return kTRUE;}
	else if(expNum == 4) {// NA62 is round
		if((acceptanceHole < trackerRadius) && (trackerRadius<acceptanceSide)) return kTRUE;}
	else if(expNum == 5 || expNum == 6) {// DarkQuest is a rectangle, 2 by 2
		if((xa<1.) && (ya<1.))		return kTRUE;}
	else if(expNum == 7){ return kFALSE; }//no trackers
	else if (expNum == 8){ // ship ecn3
		if((xa<2.) && (ya<3.))		return kTRUE;}
	else if(expNum == 16){ // SHADOWS is a rectangle, 2.5 by 2.5
		if((xa<1.25) && (ya<1.25))	return kTRUE;}
	else if(expNum == 17) {// shipecn4 is a rectangle 5 by 10
		if((xa<2.5) && (ya<5.))		return kTRUE;}
	else if (expNum == 18){ // ORCA
		if((xa<112.5) && (ya<100.))	return kTRUE;}
	else if (expNum == 19){ // BEBCcuboid
		if(xa<1260.|| (ya<1785)  ) return kTRUE;}
	return kFALSE;
}

Bool_t ExpParameters::inMuonDetectorAcceptance(TVector3 posFVEnd, TVector3 posMuDet, TVector3 posDecay){
	Double_t xa = fabs(posMuDet.X()), ya = fabs(posMuDet.Y());
	
	if(expNum == 0) { // CHARM has no trackers. Sinmply ask both muons within 3x3 m at end of decay volume?
		if((xa<1.5) && (ya<1.5)) return kTRUE; // both Muons within 3x3m. Distance requirement?
	}
	else if (expNum == 1){// BEBC EMI angular coverage of r = 7.21 m 25x6m semicircle around (0,0,404.925) covers entire forward region 
		// return kTRUE;
		Double_t r = 7.665, zD = posDecay.Z() - ZECal, xM = fabs(posMuDet.X()), xD = fabs(posDecay.X()); 
		Double_t a = 2*(xM*xD +r*zD) - (xM*xM + xD*xD + zD*zD + r*r), b = 2*(xM*(xM-xD) + r *(r- zD));
		Double_t alpha_ =  (TMath::Sqrt(b*b + 4.*a*xM*xM) - b) / 2 / a;
		if( fabs(posMuDet.Y()*(1-alpha_) + alpha_ * posDecay.Y())<3 ) return kTRUE; // 
	}
	else if(expNum == 2) { // nucal, see https://arxiv.org/abs/1311.3870
		if(pow(xa*xa + ya*ya,0.5)<1.3) return kTRUE; // follow bluemlein
	}
	else if (expNum == 3){// NuTeV
		if(xa < 1.5 &&  ya<1.5) return kTRUE; // 3x3m scintilator planes
	}
	else if(expNum == 4){ // NA62 MuonDetector: 2.64 x 2.64 m with 0.22 diameter (central tile) hole
		if((xa<1.32) && (ya<1.32) && (0.0121<(xa*xa+ya*ya))) return kTRUE;
	}
	if(expNum == 5 || expNum == 6){ // DarkQuest, assume the same size as calo: 2 by 2
		if((xa<1) && (ya<1))	return kTRUE;
	}
	else if(expNum == 7) { // DUNE
		Double_t zDecay = posDecay.Z();
		if((zDecay > ZFVIn) && (zDecay < ZECal )){ //decays in DUNE LAr
			if((xa<3.) && (ya<1.))	return kTRUE;
		} else if((zDecay > ZECal) && (zDecay < ZFVEnd)){ //decays in DUNE GAr
			if((abs(posFVEnd.X())<2.4) && (abs(posFVEnd.Y())<1.))	return kTRUE;
		}
	}
	if (expNum == 8) { // ship ecn3
		if (xa < 2590. && ya < 3310. ) return kTRUE; //ECN3 muon/hadron callorimeter is rectangle with semi axes x 2590 mm, y 3310 mm

	}
	if(expNum == 16){ // SHADOWS, assume the same size as calo: 2.5 by 2.5
		if((xa<1.25) && (ya<1.25))	return kTRUE;
	}
	if(expNum == 17){ // ship, assume the same size as calo: 5 by 10
		if((xa<2.5) && (ya<5))	return kTRUE;
	}
	else if (expNum == 19){// BEBC EMI angular coverage of r = 7.21 m 25x6m semicircle around (0,0,404.925) covers entire forward region 
		Double_t r = 7.665, zD = posDecay.Z() - ZECal, xM = fabs(posMuDet.X()), xD = fabs(posDecay.X());
		Double_t a = 2*(xM*xD +r*zD) - (xM*xM + xD*xD + zD*zD + r*r), b = 2*(xM*(xM-xD) + r *(r- zD));
		Double_t alpha_ =  (TMath::Sqrt(b*b + 4.*a*xM*xM) - b) / 2 / a;
		if( fabs(posMuDet.Y()*(1-alpha_) + alpha_ * posDecay.Y())<3 ) return kTRUE; // 
	}
	return kFALSE;
}

 /// \fn twoPhotonCondition
 /// \Brief
 /// Experimental cut to pass two photon selection
 /// \EndBrief
Int_t ExpParameters::twoPhotonCondition(Int_t iOnECal, Double_t totalEnergyAcceptance, TVector3 posDecay, TVector3 posFVEnd[2], TVector3 posECal[2], TLorentzVector pgA[2]){
	TVector3 diffECal = posECal[1]-posECal[0];
	Double_t distance = diffECal.Mag();
	Double_t radiuspos0 = posECal[0].XYvector().Mod();
	Double_t radiuspos1 = posECal[1].XYvector().Mod();
	Double_t energy0 = pgA[0].E();
	Double_t energy1 = pgA[1].E();

	if(expNum == 0){ //CHARM
		if(iOnECal==1 && energy0>5. && energy0<50.) return 1;
		else if(iOnECal==2 && energy0>5. && energy1 > 5. && energy0<50. && energy1 < 50.) return 1; //iOnECal>=1 && ; improved CHARM
		return 0;
	}
	else if(expNum == 3){ //NuTeV
		return exoNum == 1; // NuTeV only cares about the charged tracks (true for hnls only (pipi0 decays))
	}
	else if(expNum == 2){ //NuCal
		if((iOnECal==1 && (totalEnergyAcceptance>10)) || (iOnECal ==2 && ((energy0>10.) || (energy1>10.)))) return 1; // NuCal
		return 0;
	}
	else if(expNum == 4){ //NA62
		if(iOnECal==2 && distance > 0.1 && energy0>1. && energy1 > 1. && radiuspos1> 0.15 && radiuspos0> 0.15) return 1; 
		return 0;
	}
	else if(expNum == 5 || expNum == 6){ //DarkQuest
		if(iOnECal==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; // DarkQuest
		return 0;
	}
	else if(expNum == 7){ //DUNE
		Double_t zDecay = posDecay.Z();
		if((zDecay > ZFVIn) && (zDecay < ZECal )){ //LAr
			if(iOnECal==2) return 1;
			return 0;
		} else if((zDecay > ZECal) && (zDecay < ZFVEnd)){ //GAr
			if(iOnECal==2 && energy0>0.02 && energy1 > 0.02 && pgA[0].Angle(pgA[1].Vect()) > 0.14) return 1; //at least 20MeV and 0.8deg separation angle
			return 0;
		}
		return 0;
	}
	else if(expNum == 8){ //SHiP in ecn3 same as original ship
		if(iOnECal==2  && distance > 0.1  && energy0>1. && energy1 > 1.) {
			if (!(this->GetDecayModeOpen() )) return kTRUE;
			TLorentzVector p_vis = (pgA[0]+pgA[1]);
			Double_t delta_z_FV = (posDecay.Z() - ZFVIn) / (ZFVEnd - ZFVIn);
			return ((posDecay.Z() - ZFVIn) > 1 )  &&  ( fabs(posDecay.X()) < (0.645 + delta_z_FV *  1.5) ) && ( fabs(posDecay.Y()) < (1.495 + delta_z_FV *  1.65 ) ) &&  fabs((posDecay.Perp() - p_vis.Z() / p_vis.Pt() * (posDecay.Z() - ZTarget )) < 2.5); //
		}
		return 0;
	}
	else if(expNum == 9){ //KOTOdump must resolve lower gamma energies
		// `dream cuts` similar to first attempt, see #34 on GITHUB.
		// for the photon cluster distance cut we follow the thesis by Melissa Hutcheson (30cm)
		if(iOnECal==2 && distance > 0.3 && energy0>0.05 && energy1 > 0.05) return kTRUE; // 50MeV

		return 0;
	}
	else if(expNum == 10){ //KOTO
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnECal==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posECal[0].X()),abs(posECal[0].Y())) > 0.15 && TMath::Min(abs(posECal[1].X()),abs(posECal[1].Y())) > 0.15)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posECal[0].X()*energy0+posECal[1].X()*energy1)/(energy0+energy1),
			(posECal[0].Y()*energy0+posECal[1].Y()*energy1)/(energy0+energy1),
			posECal[0].Z()
		);
		Double_t rCOE = posCOE.XYvector().Mod();
		if(rCOE < 0.2) return 0;
		// values assuming on-axis Pnn decay with pi0->2gamma
		Double_t angleSepPi0 = TMath::Pi();
		if(MPi0*MPi0/(energy0*energy1)<=4) angleSepPi0 = TMath::ACos(1-MPi0*MPi0/(2*energy0*energy1)); //photon separation angle pi0 hypothesis
		else return 0;
		Double_t bracket1 = distance*distance-TMath::Power(TMath::Sin(angleSepPi0),2)*(radiuspos0*radiuspos0+radiuspos1*radiuspos1);
		Double_t det1 = TMath::Power(bracket1,2)+TMath::Power(TMath::Sin(angleSepPi0),2)*(4*radiuspos0*radiuspos0*radiuspos1*radiuspos1*TMath::Power(TMath::Cos(angleSepPi0),2)-TMath::Power(radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance,2));
		if (det1 < 0) return 0;
		Double_t det2 =  bracket1 + TMath::Sqrt(det1);
		if (det2 < 0) return 0;
		Double_t zVtxDistPnn = TMath::Sqrt(det2)/(TMath::Sqrt2()*TMath::Sin(angleSepPi0));
		Double_t zDecayPnnInDV = ZECal-zVtxDistPnn-21; //21m - beginning of DV in pnn analysis

		Double_t ptmax = 0.25; //250MeV
		Double_t ptmin = 0.13;
		if(zDecayPnnInDV > 4) ptmin = 0.13 + (zDecayPnnInDV-4)*0.2/0.7; //for zDecay > 4m linear dependence of ptmin on zDecay

		Double_t angleDegPnnPhot0 = TMath::ATan(radiuspos0/zVtxDistPnn)*180./TMath::Pi();
		Double_t angleDegPnnPhot1 = TMath::ATan(radiuspos1/zVtxDistPnn)*180./TMath::Pi();

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		//Selection by regions
		if( zDecayPnnInDV > 2.9 && zDecayPnnInDV < 6 && angleDegPnnPhot0*energy0 > 2.5 && angleDegPnnPhot1*energy1 > 2.5){ // general selection: zVtx; photon angle (Pnn hypothesis) * cluster energy > 2500 MeV*deg
			if(zDecayPnnInDV > 3 && zDecayPnnInDV < 4.7 && pCOE.Pt() > ptmin && pCOE.Pt() < ptmax) return 1; //region 1 same as Pnn
			else if(zDecayPnnInDV < 5.1 && pCOE.Pt() < 0.12) return 2;
			else if(zDecayPnnInDV >= 5.1 && pCOE.Pt() < 0.12) return 3;
			else if(zDecayPnnInDV < 5.1 && pCOE.Pt() < 0.26) return 4;
			else if(zDecayPnnInDV >= 5.1 && pCOE.Pt() < 0.26) return 5;
			else if(zDecayPnnInDV < 5.1 && pCOE.Pt() < 0.5) return 6;
			else if(zDecayPnnInDV >= 5.1 && pCOE.Pt() < 0.5) return 7;
			else return 8; //region with pT > 500MeV
		} else return 0;

	}
	else if(expNum == 11){ //KOTO2dump
		if(iOnECal==2 && distance > 0.3 && energy0>0.1 && energy1 > 0.1 && (energy0 + energy1) > 0.5) return 1; // cuts 1-5 (p234) from Pnn ensuring photon cluster quality
		return 0;
	}
	else if(expNum == 12){ //KOTO2
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnECal==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0; //stays the same as for KOTO1
		if(!(radiuspos0 < 1.35 && radiuspos1 < 1.35 && distance > 0.3 && TMath::Min(abs(posECal[0].X()),abs(posECal[0].Y())) > 0.175 && TMath::Min(abs(posECal[1].X()),abs(posECal[1].Y())) > 0.175)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.5 && energy0 > 0.1 && energy1 > 0.1)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posECal[0].X()*energy0+posECal[1].X()*energy1)/(energy0+energy1),
			(posECal[0].Y()*energy0+posECal[1].Y()*energy1)/(energy0+energy1),
			posECal[0].Z()
		);
		Double_t rCOE = posCOE.XYvector().Mod();
		// if(rCOE < 0.2) return 0; no rCOE cut in KOTO2
		// values assuming on-axis Pnn decay with pi0->2gamma
		Double_t angleSepPi0 = TMath::Pi();
		if(MPi0*MPi0/(energy0*energy1)<=4) angleSepPi0 = TMath::ACos(1-MPi0*MPi0/(2*energy0*energy1)); //photon separation angle pi0 hypothesis
		else return 0;
		Double_t bracket1 = distance*distance-TMath::Power(TMath::Sin(angleSepPi0),2)*(radiuspos0*radiuspos0+radiuspos1*radiuspos1);
		Double_t det1 = TMath::Power(bracket1,2)+TMath::Power(TMath::Sin(angleSepPi0),2)*(4*radiuspos0*radiuspos0*radiuspos1*radiuspos1*TMath::Power(TMath::Cos(angleSepPi0),2)-TMath::Power(radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance,2));
		if (det1 < 0) return 0;
		Double_t det2 =  bracket1 + TMath::Sqrt(det1);
		if (det2 < 0) return 0;
		Double_t zVtxDistPnn = TMath::Sqrt(det2)/(TMath::Sqrt2()*TMath::Sin(angleSepPi0));

		Double_t zDecayPnnInDV = ZECal-zVtxDistPnn-44; //the Pnn decay volume starts after 44m beam line

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		if( zDecayPnnInDV < 15 && zDecayPnnInDV > 1.75 && pCOE.Pt() > 0.4) //exclude the 1.75m long front barrel counter
			return 1;
		else
			return 0;

	}
	else if(expNum == 13){ //KOTO Pnn region excluded
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnECal==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posECal[0].X()),abs(posECal[0].Y())) > 0.15 && TMath::Min(abs(posECal[1].X()),abs(posECal[1].Y())) > 0.15)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posECal[0].X()*energy0+posECal[1].X()*energy1)/(energy0+energy1),
			(posECal[0].Y()*energy0+posECal[1].Y()*energy1)/(energy0+energy1),
			posECal[0].Z()
		);
		Double_t rCOE = posCOE.XYvector().Mod();
		if(rCOE < 0.2) return 0;
		// values assuming on-axis Pnn decay with pi0->2gamma
		Double_t angleSepPi0 = TMath::Pi();
		if(MPi0*MPi0/(energy0*energy1)<=4) angleSepPi0 = TMath::ACos(1-MPi0*MPi0/(2*energy0*energy1)); //photon separation angle pi0 hypothesis
		else return 0;
		Double_t bracket1 = distance*distance-TMath::Power(TMath::Sin(angleSepPi0),2)*(radiuspos0*radiuspos0+radiuspos1*radiuspos1);
		Double_t det1 = TMath::Power(bracket1,2)+TMath::Power(TMath::Sin(angleSepPi0),2)*(4*radiuspos0*radiuspos0*radiuspos1*radiuspos1*TMath::Power(TMath::Cos(angleSepPi0),2)-TMath::Power(radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance,2));
		if (det1 < 0) return 0;
		Double_t det2 =  bracket1 + TMath::Sqrt(det1);
		if (det2 < 0) return 0;
		Double_t zVtxDistPnn = TMath::Sqrt(det2)/(TMath::Sqrt2()*TMath::Sin(angleSepPi0));
		Double_t zDecayPnnInDV = ZECal-zVtxDistPnn-21; //21m - beginning of DV in pnn analysis

		Double_t ptmax = 0.25; //250MeV
		Double_t ptmin = 0.13;
		if(zDecayPnnInDV > 4) ptmin = 0.13 + (zDecayPnnInDV-4)*0.2/0.7; //for zDecay > 4m linear dependence of ptmin on zDecay

		Double_t angleDegPnnPhot0 = TMath::ATan(radiuspos0/zVtxDistPnn)*180./TMath::Pi();
		Double_t angleDegPnnPhot1 = TMath::ATan(radiuspos1/zVtxDistPnn)*180./TMath::Pi();

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		//Selection by regions
		if( zDecayPnnInDV > 2.9 && zDecayPnnInDV < 6 && angleDegPnnPhot0*energy0 > 2.5 && angleDegPnnPhot1*energy1 > 2.5){ // general selection: zVtx; photon angle (Pnn hypothesis) * cluster energy > 2500 MeV*deg
			if(zDecayPnnInDV < 5.1 && pCOE.Pt() < 0.26) return 0; //exclude Pnn signal region
			else return 1;
		} else return 0;

	}
	
	else if(expNum == 14 || expNum == 15){ //KLEVER Pnn region excluded
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnECal==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(radiuspos0 < 0.35 || radiuspos1 < 0.35) return 0; //r_min requirement
		if(energyL < 2) return 0; //min energy requirement

		///// reconstruction of Pnn-like event
		TVector3 posCOE(
			(posECal[0].X()*energy0+posECal[1].X()*energy1)/(energy0+energy1),
			(posECal[0].Y()*energy0+posECal[1].Y()*energy1)/(energy0+energy1),
			posECal[0].Z()
		);
		// values assuming on-axis Pnn decay with pi0->2gamma
		Double_t angleSepPi0 = TMath::Pi();
		if(MPi0*MPi0/(energy0*energy1)<=4) angleSepPi0 = TMath::ACos(1-MPi0*MPi0/(2*energy0*energy1)); //photon separation angle pi0 hypothesis
		else return 0;
		Double_t bracket1 = distance*distance-TMath::Power(TMath::Sin(angleSepPi0),2)*(radiuspos0*radiuspos0+radiuspos1*radiuspos1);
		Double_t det1 = TMath::Power(bracket1,2)+TMath::Power(TMath::Sin(angleSepPi0),2)*(4*radiuspos0*radiuspos0*radiuspos1*radiuspos1*TMath::Power(TMath::Cos(angleSepPi0),2)-TMath::Power(radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance,2));
		if (det1 < 0) return 0;
		Double_t det2 =  bracket1 + TMath::Sqrt(det1);
		if (det2 < 0) return 0;
		Double_t zVtxDistPnn = TMath::Sqrt(det2)/(TMath::Sqrt2()*TMath::Sin(angleSepPi0));
		Double_t zDecayPnn = ZECal-zVtxDistPnn;

		Double_t ptmin = 0.3;
		if(expNum == 15 && zDecayPnn < 260) ptmin = 0.1; //can afford less for ext before pnn SR

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		if(pCOE.Pt() > ptmin && zDecayPnn > 120){// general selection: zVtx; pT (Pnn hypothesis)
			if(expNum == 14 && zDecayPnn < 220) //at least 20m before ECal
				return 1;
			if(expNum == 15 && zDecayPnn < 370)
				return 1;
		}
	}
	else if(expNum == 16){ //SHADOWS
		if(iOnECal==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; // SHADOWS
		return 0;
	}
	else if(expNum == 17){ //SHiPecn4
		if(iOnECal==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; //ship
		return 0;
	}
	else if(expNum == 18){ //ORCA
		if(iOnECal == 2 &&( 2E-9 < energy0 && energy0 < 4.5E-9 ) && ( 2E-9 < energy1 && energy1 < 5E-9 )  ) return 1 ; // Minimum distance between 2 stations should be considered and energies near visible spectrum ~275nm to 600nm
	}
	return 0;
}

 /// \fn multiplePhotonCondition
 /// \Brief
 /// Experimental cut to pass up to six photon selection
 /// \EndBrief
Bool_t ExpParameters::multiplePhotonCondition(Int_t nGammas, Double_t totalEnergyAcceptance, TVector3 posDecay, TVector3 posFVEnd[6], TVector3 posECal[6], TLorentzVector pgA[6]){

// evaluate number of distinct gammas seen
	Double_t clusterMerginDistance = 0.1, zDecay = posDecay.Z();
	if(expNum == 9 || expNum == 11)
		clusterMerginDistance = 0.3; //different threshold for KOTO
	int flagMerged[6] = {-1,-1,-1,-1,-1,-1};
	for (Int_t k1=0; k1<nGammas; k1++){
		if (!inCaloAcceptance(posFVEnd[k1],posECal[k1],zDecay)) continue;
		if (flagMerged[k1] >= 0) continue;

		for (Int_t k2=k1+1; k2<nGammas; k2++){
			if (!inCaloAcceptance(posFVEnd[k2],posECal[k2],zDecay)) continue;
			if (flagMerged[k2] >= 0) continue;

			if ((posECal[k1]-posECal[k2]).Mag() < clusterMerginDistance) {
				int idom = k1;
				int isub = k2;
				if (pgA[k2].E() > pgA[k1].E()) {
					idom = k2;
					isub = k1;
				}
				flagMerged[isub] = idom;
			}
		}
	}
	int nMerged = 0;
	for (Int_t i=0; i<6; i++) if(flagMerged[i] != -1) nMerged++;
	int nClusters = nGammas - nMerged;

	int iOnECal = 0;
	int i3OnECal = 0;
	int i5to50OnECal = 0;
	int indexAbove3[6] = {0,0,0,0,0,0};
	int index5to50[6] = {0,0,0,0,0,0};
	for (Int_t k1=0; k1<nGammas; k1++){
		if (!inCaloAcceptance(posFVEnd[k1],posECal[k1],zDecay)) continue;
		if (pgA[k1].E() > 3){
			indexAbove3[i3OnECal] = k1;
			i3OnECal++; // for nuCal
		}
		if (pgA[k1].E() > 5 && pgA[k1].E()<50){
			index5to50[i5to50OnECal] = k1;
			i5to50OnECal++; // for charm
		}
		if (flagMerged[k1] >= 0) continue;
  		if (expNum == 7 && zDecay > ZECal && zDecay < ZFVEnd && pgA[k1].E() < 0.02) continue; //DUNE GAr energy threshold
		if ((expNum == 9 && pgA[k1].E() < 0.05) || (expNum == 11 && pgA[k1].E() < 0.1)) continue; //KOTO(2)dump single cluster energy threshold
		iOnECal++;
	}

	//exp. conditions
	if(expNum == 4 || expNum == 17 || expNum == 5 || expNum == 16 || expNum == 8 || expNum == 6){ //NA62, SHiP, SHADOWS, DarkQuest
		if(nGammas==iOnECal && totalEnergyAcceptance>3) return kTRUE; // a minimum requirement on the individual gamma energy should be probably added here ?
		return kFALSE;
	}
	else if(expNum == 0){
		if(i5to50OnECal==1 || i5to50OnECal==2) return kTRUE;//CHARM, as in digamma, but not more than two photons allowed
		return kFALSE;
	}
	else if(expNum == 2){ //NuCal
		Int_t i3OnECaldetected=0;
		Double_t totalEnergyDetected=0;
		for(Int_t n3onECal=0; n3onECal< i3OnECal; n3onECal++){
			// take into account efficiency of photon detection
			if(rndmGen.Rndm()>0.7){
				i3OnECaldetected++;
				totalEnergyDetected+=pgA[indexAbove3[n3onECal]].E();
			}
		}
		// mimic the following statement:
		//Requiring events with an electromagnetic shower of Eel m > 3 GeV, 106 events remain. A large fraction of them has an energetic hadron shower. This part can be removed by excluding all events with hadron energies above 1.5 GeV. For the remaining 21 events the electromagnetic shower energy is shown in Fig. 3a.
		if(i3OnECaldetected>=1 &&totalEnergyDetected >10) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 7){ //DUNE
		if(nGammas==iOnECal) return kTRUE; //we drop the angle conditions
		return kFALSE;
	}
	else if(expNum == 9){ //KOTOdump
		if(nGammas==iOnECal && nMerged == 0 && totalEnergyAcceptance>0.5) return kTRUE; // to be checked
		return kFALSE;
	}
	else if(expNum == 10){ //KOTO - same conditions as for two photons but require 2 clusters (or more photons in acceptance if clusters overlap)
		if(nClusters != 2) return kFALSE;
		std::vector<Int_t> realIndex;
		for (Int_t k1=0; k1<nGammas; k1++) if(flagMerged[k1]==-1) realIndex.push_back(k1);
		if(realIndex.size()!=nClusters){
			std::cout << "[KOTO multiplePhotonCondition] Error: unexpected number of calorimeter clusters" << std::endl;
			exit(1);
		}
		TVector3 diffECal = posECal[realIndex[1]]-posECal[realIndex[0]];
		Double_t distance = diffECal.Mag();
		Double_t radiuspos0 = posECal[realIndex[0]].XYvector().Mod();
		Double_t radiuspos1 = posECal[realIndex[1]].XYvector().Mod();
		Double_t energy0 = pgA[realIndex[0]].E();
		Double_t energy1 = pgA[realIndex[1]].E();
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(angleSepXY*180/TMath::Pi() > 150) return kFALSE;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posECal[realIndex[0]].X()),abs(posECal[realIndex[0]].Y())) > 0.15 && TMath::Min(abs(posECal[realIndex[1]].X()),abs(posECal[realIndex[1]].Y())) > 0.15)) return kFALSE; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return kFALSE; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posECal[realIndex[0]].X()*energy0+posECal[realIndex[1]].X()*energy1)/(energy0+energy1),
			(posECal[realIndex[0]].Y()*energy0+posECal[realIndex[1]].Y()*energy1)/(energy0+energy1),
			posECal[realIndex[0]].Z()
		);
		Double_t rCOE = posCOE.XYvector().Mod();
		if(rCOE < 0.2) return kFALSE;
		// values assuming on-axis Pnn decay with pi0->2gamma
		Double_t angleSepPi0 = TMath::Pi();
		if(MPi0*MPi0/(energy0*energy1)<=4) angleSepPi0 = TMath::ACos(1-MPi0*MPi0/(2*energy0*energy1)); //photon separation angle pi0 hypothesis
		else return kFALSE;
		Double_t bracket1 = distance*distance-TMath::Power(TMath::Sin(angleSepPi0),2)*(radiuspos0*radiuspos0+radiuspos1*radiuspos1);
		Double_t det1 = TMath::Power(bracket1,2)+TMath::Power(TMath::Sin(angleSepPi0),2)*(4*radiuspos0*radiuspos0*radiuspos1*radiuspos1*TMath::Power(TMath::Cos(angleSepPi0),2)-TMath::Power(radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance,2));
		if (det1 < 0) return kFALSE;
		Double_t det2 =  bracket1 + TMath::Sqrt(det1);
		if (det2 < 0) return kFALSE;
		Double_t zVtxDistPnn = TMath::Sqrt(det2)/(TMath::Sqrt2()*TMath::Sin(angleSepPi0));
		Double_t zDecayPnnInDV = ZECal-zVtxDistPnn-21;

		Double_t ptmax = 0.25; //250MeV
		Double_t ptmin = 0.13;
		if(zDecayPnnInDV > 4) ptmin = 0.13 + (zDecayPnnInDV-4)*0.2/0.7; //for zDecay > 4m linear dependence of ptmin on zDecay

		Double_t angleDegPnnPhot0 = TMath::ATan(radiuspos0/zVtxDistPnn)*180./TMath::Pi();
		Double_t angleDegPnnPhot1 = TMath::ATan(radiuspos1/zVtxDistPnn)*180./TMath::Pi();

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		if(zDecayPnnInDV > 3 && zDecayPnnInDV < 4.7 //KOTO Pnn additional geometry conditions
		   && angleDegPnnPhot0*energy0 > 2.5 && angleDegPnnPhot1*energy1 > 2.5 //photon angle (Pnn hypothesis) * cluster energy > 2500 MeV*deg
		   && pCOE.Pt() > ptmin && pCOE.Pt() //photon pT (Pnn hypothesis) cuts
		  ) return kTRUE;

		return kFALSE;
	}

	else if(expNum == 11){ //KOTO2dump
		if(nGammas==iOnECal && nMerged == 0 && totalEnergyAcceptance>0.5) return kTRUE; // to be checked
		return kFALSE;
	}
	return kFALSE;
}

 /// \fn twoMuonCondition
 /// \Brief
 /// Experimental cut to select two muons
 /// \EndBrief
Bool_t ExpParameters::twoMuonCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2]){
	
	if(expNum == 0){ //charm: 1 GeV mimimum: Mip should be seen throughout detector material (this is still generous)
		if(iOnMuonDetector>1 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 1 || expNum == 19){ // BEBC
		if( pCA[0].P() > 1 && pCA[1].P() > 1 && iOnECal == 2 && iOnMuonDetector == 2 ) {
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			if((p_vis.Pt() + p_vis.Mt()) < (MDCh - MMu)) return kTRUE;
		}
		return kFALSE;
	}
	else if(expNum == 2){ // according to bluemlein
		if(iOnECal==2 && (pCA[0].E()+pCA[1].E())>10) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 3){ // NuTeV
		// FV cuts (assert decay vertex is 1m away from material)
		if(posDecay.Perp()>2.3)return kFALSE; // R=2.3m helium bags
		for(Double_t deltaZMaterial:{15.4, 25.2}){
			Double_t zMaterial = ZFVIn + deltaZMaterial;
			if( posDecay.Z() > zMaterial - 1. && posDecay.Z() < zMaterial + 1. ) return kFALSE;
		}
		// kinematic cuts
		if( pCA[0].E() > 2 && pCA[1].E() > 2  && iOnTracker1 == 2 && iOnTracker4 == 2 && iOnECal == 2 && iOnMuonDetector == 2) {
			TLorentzVector p_vis = (pCA[0]+pCA[1]), pHarder = pCA[0].E()>pCA[1].E() ? pCA[0]:pCA[1];
			Double_t pNu  =  pCA[1].Z() + pCA[0].Z();
			Double_t Q2_vis = (TLorentzVector(0,0,pNu,pNu) - pHarder).M2(), nu_vis = TMath::Sqrt(MP*MP+p_vis.Perp2());
			if(  (Q2_vis < 0.2 * MP*nu_vis ) && ( MP*MP + 2*MP*nu_vis - Q2_vis < 4.) && (p_vis.Pt() + p_vis.Mt() < 3. )) return kTRUE;
		}
		return kFALSE;
	}
	else if(expNum == 4|| expNum == 17 || expNum == 16 ){ // NA62, shipecn4 and SHADOWS 
		if (iOnECal ||  iOnTracker4) return kTRUE;
		// if( pCA[0].E() > 5 && pCA[1].E() > 5  && iOnTracker1 == 2 && iOnTracker4 == 2) return kTRUE; //similar spectrometers  && iOnMuonDetector == 2
		return kFALSE;
	}
	else if(expNum == 5 || expNum == 6){ // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"  https://www.sciencedirect.com/science/article/abs/pii/S016890021930347X
	// also there is a 1m iron absorber which takes 1GeV off the muons
		if(iOnTracker1 == 2 && iOnTracker4 == 2 && (pCA[0].E()+pCA[1].E())>4.2 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 7){ //seeing in https://arxiv.org/pdf/2103.13910.pdf page 2-71, muon energies below 1GeV are fully OK. Above the angle wr.t. to the beam axis has to be smaller than 40degrees (700mrad)
		if(iOnMuonDetector == 2){
			if((pCA[0].E()-pCA[0].M()) < 1 && (pCA[1].E()-pCA[1].M()) < 1) return kTRUE;
			if(pCA[0].Theta()<0.7 && pCA[1].Theta()<0.7) return kTRUE;
		}
		return kFALSE;
	} // exp 5
	else if(expNum == 8){ //SHiP
		if(  iOnTracker1 == 2 && iOnTracker4 == 2 && iOnECal == 2&& iOnMuonDetector==2 && pCA[0].P() > 1 && pCA[1].P() > 1 ) {
			if (!(this->GetDecayModeOpen() )) return kTRUE;
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			Double_t delta_z_FV = (posDecay.Z() - ZFVIn) / (ZFVEnd - ZFVIn);
			return ((posDecay.Z() - ZFVIn) > 1 )  &&  ( fabs(posDecay.X()) < (0.645 + delta_z_FV *  1.5) ) && ( fabs(posDecay.Y()) < (1.495 + delta_z_FV *  1.65 ) ) &&  fabs((posDecay.Perp() - p_vis.Z() / p_vis.Pt() * (posDecay.Z() - ZTarget )) < 2.5); //
		}
		return kFALSE;
	}
	return kFALSE;
}

 /// \fn twoHadronCondition
 /// \Brief
 /// Experimental cut to select two hadrons
 /// \EndBrief
Bool_t ExpParameters::twoHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2]){
	
	if(expNum == 0){
		// if(iChOnECal==0) return kTRUE; //Charm vetos hadronic showers in (https://inspirehep.net/literature/214233)
		return kFALSE;
	}
	else if(expNum == 2){
		// if(iChOnECal==0) return kTRUE;  //Nucal explicitely comments on vetoing hadronic showers above 1.5GeV
		return kFALSE;
	}
	else if(expNum == 3){ //NuTeV
		return 0; // Only final states with a muon were considered
	}
	else if(expNum == 4|| expNum == 17 || expNum == 16){ // NA62, shipecn4 and SHADOWS
		if(pCA[0].E() > 5 && pCA[1].E() > 5  && iOnTracker1 == 2 && iOnTracker4 == 2 && iChOnECal == 2) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 5 || expNum == 6){ // DarkQuest. in phase 1 0 sensi to hadrons because of background. will chance for phase 2, see paper 2008.08108
		//return kTRUE;
		//if(iOnTracker1 == 2 && iOnTracker4 == 2 && (pCA[0].E()+pCA[1].E())>4.2) return kTRUE; // mimic the muon conditions from twoMuonCondition
		//if(iOnMuonDetector == 2 && iChOnECal>=2 && iOnTracker1 == 2 && iOnTracker4 == 2 && totalEnergyAcceptance > energyCut[experiment]) kTRUE; // DarkQuest <-??
		return kFALSE;
	}
	else if(expNum == 7){
		Double_t zDecay = posDecay.Z();
		if((zDecay > ZFVIn) && (zDecay < ZECal )){ //LAr
			if(iChOnECal==2) return kTRUE;
			return kFALSE;
		} else if((zDecay > ZECal) && (zDecay < ZFVEnd)){ //GAr
			if(iChOnECal==2 && (pCA[0].E() - finState.at(0)->GetMass()) > 0.005 && (pCA[1].E() - finState.at(1)->GetMass()) > 0.005) return kTRUE; //kinetic energy above 5MeV
			return kFALSE;
		}
		return kFALSE;
	}
	else if(expNum == 8){
		if(  iOnTracker1 == 2 && iOnTracker4 == 2 && iChOnECal==2 && pCA[0].P() > 1 && pCA[1].P() > 1 ) {
			if (!(this->GetDecayModeOpen() )) return kTRUE;
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			Double_t delta_z_FV = (posDecay.Z() - ZFVIn) / (ZFVEnd - ZFVIn);
			return ((posDecay.Z() - ZFVIn) > 1 )  &&  ( fabs(posDecay.X()) < (0.645 + delta_z_FV *  1.5) ) && ( fabs(posDecay.Y()) < (1.495 + delta_z_FV *  1.65 ) ) &&  fabs((posDecay.Perp() - p_vis.Z() / p_vis.Pt() * (posDecay.Z() - ZTarget )) < 2.5); //
		}
		return kFALSE;
	}
	return kFALSE;
}

 /// \fn muonHadronCondition
 /// \Brief
 /// Experimental cut to select muon+hadron
 /// \EndBrief
Bool_t ExpParameters::muonHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2]){
	
	if(expNum == 0){ // CHARM
		return 0; //Charm vetos hadronic showers
	}
	else if(expNum == 1||expNum == 19){ // BEBC
		if( pCA[0].P() > 1 && pCA[1].P() > 1 && iChOnECal == 2 && iOnMuonDetector == 1) return kTRUE;
		return kFALSE;
	}
	if(expNum == 3){//NuTeV
		// FV cuts (assert decay vertex is 1m away from material)
		if(posDecay.Perp()>2.3)return kFALSE; // R=2.3m helium bags
		for(Double_t deltaZMaterial:{15.4, 25.2}){
			Double_t zMaterial = ZFVIn + deltaZMaterial;
			if( posDecay.Z() > zMaterial - 1. && posDecay.Z() < zMaterial + 1. ) return kFALSE;
		}
		// kinematic cuts
		if( pCA[0].E() > 10. && pCA[1].E() > 2. && iOnTracker1 == 2 && iOnTracker4 == 2 &&  iChOnECal == 2  && iOnMuonDetector == 1) return kTRUE;
		return kFALSE;
	}
	if(expNum == 4 || expNum == 16){ // NA62 and Shadows
		if( pCA[0].E() > 5. && pCA[1].E() > 5.  && iOnTracker1 == 2 && iOnTracker4 == 2 &&  iChOnECal == 2  && iOnMuonDetector == 1) return kTRUE;//similar spectrometers
		return kFALSE;
	}
	if(expNum == 5 || expNum == 6){ // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"  https://www.sciencedirect.com/science/article/abs/pii/S016890021930347X
	// also there is a 1m iron absorber which takes 1GeV off the muons
		if(iOnTracker1 == 2 && iOnTracker4 == 2 && (pCA[0].E()+pCA[1].E())>4.2 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	if(expNum == 7){ // DUNE 
		if(iOnMuonDetector == 1  && iChOnECal == 2 && ((pCA[1].E()-pCA[1].M()) < 1 || pCA[1].Theta()<0.7)){ //if muon is detectable and one of the Calos detected both particles
			Double_t zDecay = posDecay.Z();
			if((zDecay > ZFVIn) && (zDecay < ZECal )) return kTRUE; //LAr
			if((zDecay > ZECal) && (zDecay < ZFVEnd) && (pCA[0].E() - finState.at(0)->GetMass()) > 0.005) return kTRUE; //GAr //kinetic energy above 5MeV
		}
		return kFALSE;
	}
	if(expNum == 8){ // SHiP
		if( iOnTracker1 == 2 && iOnTracker4 == 2 &&  iChOnECal == 2  && iOnMuonDetector == 1 && pCA[0].P() > 1. && pCA[1].P() > 1.) return kTRUE;
		return kFALSE;
	}
	std::cout<<"\n[Info] : muonHadron condition not implemented for selected experiment."<<std::endl;//..
	return 0;
	}

 /// \fn electronMuonCondition
 /// \Brief
 /// Experimental cut to select muon+electron
 /// \EndBrief
Bool_t ExpParameters::electronMuonCondition(Double_t totalEnergyAcceptance, Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TVector3 posECal[2], TLorentzVector pCA[2]){
	Double_t radiuspos0 = posECal[0].XYvector().Mod();
	TVector3 diffECal = posECal[1]-posECal[0];
	Double_t distance = diffECal.Mag();
	if(expNum == 0){ //charm: 1 GeV mimimum: Mip should be seen throughout detector material (this is still generous)
		if(iChOnECal == 2 && iOnMuonDetector==1 && 5 < pCA[0].E() && pCA[0].E() < 50 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 1 || expNum == 19){ // BEBC
		if( pCA[0].P()>1 && pCA[1].P() > 1 && iChOnECal == 2 ) {
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			if((p_vis.Pt() + p_vis.Mt()) < (MDCh - MMu)) return kTRUE; // This is too strict for BC6 as it should only require MDCh - MEl -> conservative
		}
		return kFALSE;
	}
	else if(expNum == 3){ //NuTeV
		// FV cuts (assert decay vertex is 1m away from material)
		if(posDecay.Perp()>2.3)return kFALSE; // R=2.3m helium bags
		for(Double_t deltaZMaterial:{15.4, 25.2}){
			Double_t zMaterial = ZFVIn + deltaZMaterial;
			if(posDecay.Z() > zMaterial - 1. && posDecay.Z() < zMaterial + 1.) return kFALSE;
		}
		// // kinematic cuts
		if( pCA[0].E() > 2 && pCA[1].E() > 2  && iOnTracker1 == 2 && iOnTracker4 == 2 && iChOnECal == 2 && iOnMuonDetector == 2) {
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			Double_t pNu  = p_vis.Z();
			Double_t Q2_vis = (TLorentzVector(0,0,pNu,pNu) - pCA[1]).M2(), nu_vis = TMath::Sqrt(MP*MP+p_vis.Perp2());
			if(  (Q2_vis < 0.2 * MP*nu_vis ) && ( MP*MP + 2*MP*nu_vis - Q2_vis < 4.) && (p_vis.Pt() + p_vis.Mt() < 3. )) return kTRUE;
		}
		return kFALSE;
	}
	else if(expNum == 4){ //NA62
		if(pCA[0].E() > 5. && pCA[1].E() > 5. && iOnTracker1 == 2 && iOnTracker4 == 2 && iChOnECal==2 && radiuspos0 > 0.15 && iOnMuonDetector == 1 ) return kTRUE; 
		return kFALSE;
	}
	else if(expNum == 7){ // DUNE 
		if(iOnMuonDetector == 1  && iChOnECal == 2 && ((pCA[1].E()-pCA[1].M()) < 1 || pCA[1].Theta()<0.7)){ //if muon is detectable and one of the Calos detected both particles
			Double_t zDecay = posDecay.Z();
			if((zDecay > ZFVIn) && (zDecay < ZECal )) return kTRUE; //LAr
			if((zDecay > ZECal) && (zDecay < ZFVEnd) && (pCA[0].E() - finState.at(0)->GetMass()) > 0.005 && (pCA[1].E() - finState.at(1)->GetMass()) > 0.005) return kTRUE; //GAr //kinetic energy above 5MeV
		}
		return kFALSE;
	}
	else if(expNum == 8){ //SHiP@ECN3
		if(  iOnTracker1 == 2 && iChOnECal==2 && iOnMuonDetector==1 && pCA[0].P() > 1 && pCA[1].P() > 1 ){
			if (!(this->GetDecayModeOpen() )) return kTRUE;
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			Double_t delta_z_FV = (posDecay.Z() - ZFVIn) / (ZFVEnd - ZFVIn);
			return ((posDecay.Z() - ZFVIn) > 1 )  &&  ( fabs(posDecay.X()) < (0.645 + delta_z_FV *  1.5) ) && ( fabs(posDecay.Y()) < (1.495 + delta_z_FV *  1.65 ) ) &&  fabs((posDecay.Perp() - p_vis.Z() / p_vis.Pt() * (posDecay.Z() - ZTarget )) < 2.5); //
		}
		return kFALSE;
	}
	else if(expNum == 16){ //SHADOWS
		if(iChOnECal==2 && distance > 0.1 && iOnMuonDetector==1 && pCA[0].E() > 1 && pCA[1].E() > 1 && totalEnergyAcceptance >3) return kTRUE; // SHADOWS
		return kFALSE;
	}
	std::cout<<"\n[Info] : electronMuon condition not implemented for selected experiment."<<std::endl;//..
	return kFALSE;
	}

 /// \fn electronHadronCondition
 /// \Brief
 /// Experimental cut to select hadron+electron
 /// \EndBrief
Bool_t ExpParameters::electronHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, TVector3 posDecay, TVector3 posECal[2], TLorentzVector pCA[2]){
	
	Double_t radiuspos0 = posECal[0].XYvector().Mod();
	
	if(expNum == 1||expNum == 19){ // BEBC
		if( pCA[0].P() > 1 && pCA[1].P() > 1 && iChOnECal == 2 ) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 3){ //NuTeV
		return kFALSE; // Only final states with a muon were considered
	}
	else if(expNum == 4||expNum == 16){ //NA62
		if(pCA[0].E() > 5. && pCA[1].E() > 5. && iOnTracker1 == 2 && iOnTracker4 == 2 && iChOnECal==2 && radiuspos0 > 0.15 ) return kTRUE; 
		return kFALSE;
	}
	else if(expNum == 5 || expNum == 6){ // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"  https://www.sciencedirect.com/science/article/abs/pii/S016890021930347X
		if(iOnTracker1 == 2 && iOnTracker4 == 2 && (pCA[0].E()+pCA[1].E())>4.2 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 7){ // DUNE
		if(iChOnECal == 2){ //
			Double_t zDecay = posDecay.Z();
			if((zDecay > ZFVIn) && (zDecay < ZECal )) return kTRUE; //LAr
			if((zDecay > ZECal) && (zDecay < ZFVEnd) && (pCA[0].E() - finState.at(0)->GetMass()) > 0.005 && (pCA[1].E() - finState.at(1)->GetMass()) > 0.005) return kTRUE; //GAr //kinetic energy above 5MeV
		}
		return kFALSE;
	}
	else if(expNum == 8){ //SHiP@ECN3
		if(iOnTracker1==2  && iOnTracker4==2 && iChOnECal==2 && pCA[0].P() > 1 && pCA[1].P() > 1 ) return kTRUE;
		return kFALSE;
	}
	std::cout<<"\n[Info] : electronHadron condition not implemented for selected experiment."<<std::endl;//..
	return 0;
	}

 /// \fn twoElectronCondition
 /// \Brief
 /// Experimental cut to select two electrons
 /// \EndBrief
Bool_t ExpParameters::twoElectronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iOnECal, Double_t totalEnergyAcceptance,TVector3 posDecay, TVector3 posFVEnd[2], TVector3 posECal[2], TLorentzVector pCA[2]){

	if(expNum == 1||expNum == 19){ // BEBC
		if( pCA[0].P() > 1 && pCA[1].P() > 1 && iOnECal == 2 ) {
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			if((p_vis.Pt() + p_vis.Mt()) < (MDCh - MEl)) return kTRUE;
		}
		return kFALSE;
	}
	else if(expNum == 3){ //NuTeV
		return kFALSE; // Only final states with a muon were considered
	}
	else if(expNum == 5 || expNum == 6){ // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"  https://www.sciencedirect.com/science/article/abs/pii/S016890021930347X
	// also there is a 1m iron absorber which takes 1GeV off the muons
		if(iOnTracker1 == 2 && iOnTracker4 == 2 && (pCA[0].E()+pCA[1].E())>4.2 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 8){ //SHiP@ECN3
		if(iOnTracker1==2  && iOnTracker4==2 && iOnECal==2 && pCA[0].P() > 1 && pCA[1].P() > 1 ) {
			if (!(this->GetDecayModeOpen() )) return kTRUE;
			TLorentzVector p_vis = (pCA[0]+pCA[1]);
			Double_t delta_z_FV = (posDecay.Z() - ZFVIn) / (ZFVEnd - ZFVIn);
			return ((posDecay.Z() - ZFVIn) > 1 )  &&  ( fabs(posDecay.X()) < (0.645 + delta_z_FV *  1.5) ) && ( fabs(posDecay.Y()) < (1.495 + delta_z_FV *  1.65 ) ) &&  fabs((posDecay.Perp() - p_vis.Z() / p_vis.Pt() * (posDecay.Z() - ZTarget )) < 2.5); //
		}
		return kFALSE;
	}
	else if(expNum == 18){ //ORCA

		Double_t energy0 = pCA[0].E();
		Double_t energy1 = pCA[1].E();
		if(iOnECal == 2 &&( 1.257*MEl < energy0) && (1.257*MEl < energy1)  ) return 1 ; // Minimum distance between 2 stations should be considered and energies above cherenkov threshold in seawater
	}
	return this->twoPhotonCondition(iOnECal,totalEnergyAcceptance, posDecay, posFVEnd, posECal, pCA);
}