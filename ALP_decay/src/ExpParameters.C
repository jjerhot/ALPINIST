#include "ExpParameters.h"

//ClassImp(ExpParameters)

ExpParameters::ExpParameters(){
	//Intitializing common values:
	energyMinInFile = 5.;
	energyMaxInFile = 295.;
	energyValuesInFile = 30;// check: one more than steps in ALP Mathematica output
	energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

	thetaMinInFile = 0.00018;
	thetaMaxInFile = 0.01098;
	thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
	thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); // used to be log scale

	nWidths = 101;
	widthMin = log10(1.*pow(10,-25));
	widthMax = log10(1.*pow(10,-10));
	wdWidth = (widthMax-widthMin)/(nWidths-1);

	massXMin = {log10(0.0001),log10(0.01)}; //min mass in massfile[i] //log scale output
	massXMax = {log10(0.01),log10(3.01)}; //max mass in massfile[i] //log scale output
	nMX = {101,101}; //number of bins for output
	nMassFiles = 2;
	mAMinInFile = {0.0001,0.01}; //min mass in massfile[i] //lin scale input
	mAMaxInFile = {0.01,3.01}; //max mass in massfile[i] //lin scale input
	mAValuesInFile = {101,601};//number of bins for input //check: one more than steps in ALP Mathematica output
	minMassFile = {"01","10"};
	maxMassFile = {"10","3010"};

	for (Int_t ifile=0; ifile < nMassFiles; ifile++){
		massXWid[ifile] = (massXMax[ifile]-massXMin[ifile])/(nMX[ifile]-1); //log scale output
		mAWid[ifile] = (mAMaxInFile[ifile]-mAMinInFile[ifile])/(mAValuesInFile[ifile]-1); //lin scale input
	}

	rndm = new TRandom3(0);
	rndm->SetSeed(0);

	PhiBeamEuler = 0;
	ThetaBeamEuler = 0;
	PsiBeamEuler = 0;

	nRegions = 1;
}

ExpParameters::ExpParameters(Int_t Experiment_, Int_t ProductionMode_, Int_t DecayMode_) : ExpParameters::ExpParameters(){
	//Intitializing experiment specific values
	expNum = Experiment_;
	switch(DecayMode_){ // so far "2Gamma", "2El","2Mu", "3Pi0", "3Pi", "2PiGamma", "2Pi0Eta", "2PiEta", "2Pi0EtaPrim", "2PiEtaPrim",
		case 0: mFinState={0., 0., 0.}; chargeFinState={0,0,0}; break;
		case 1: mFinState={MEl, MEl, 0.}; chargeFinState={1,-1,0}; break;
		case 2: mFinState={MMu, MMu, 0.}; chargeFinState={1,-1,0}; break;
		case 3:	mFinState={MPi0, MPi0, MPi0}; chargeFinState={0,0,0};	break;
		case 4: mFinState={MPiCh, MPiCh, MPi0}; chargeFinState={1,-1,0}; break;
		case 5: mFinState={MPiCh, MPiCh, 0.}; chargeFinState={1,-1,0}; break;
		case 6: mFinState={MPi0, MPi0, MEta}; chargeFinState={0,0,0}; break;
		case 7: mFinState={MPiCh, MPiCh, MEta}; chargeFinState={1,-1,0}; break;
		case 8: mFinState={MPi0, MPi0, MEtaPrim}; chargeFinState={0,0,0}; break;
		case 9: mFinState={MPiCh, MPiCh, MEtaPrim}; chargeFinState={1,-1,0}; break;
		default: mFinState={0., 0., 0.}; chargeFinState={0,0,0}; break;
	}
	switch(Experiment_){
		case 0 :
			expLabel = expLabels[Experiment_];
			expName = "NA62";
			beamEnergy = 400;
			ATarget = 63.5; //ACu
			ZTarget = 29; //ZCu
			ZTAX = 23.070; //XTAX101024
			X0 = 0.;  // m
			Y0 = -0.022;  // m
			Z0 = ZTAX; //was 20 m
			ZFVIn = 102; //102.425; // 105, comparing with NA62MC
			zStraw1 = 183. ; // m, position of 1st straw chamber
			zStraw4 = 218.7; // m, position of 4th straw chamber
			ZMagnet     = 0.5*(196.350+197.650);
			MagnetZLength       = 1.3; // m
			MagnetFieldStrength = 0.6928; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 241.495;
			ZFVEnd = zStraw1;
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; // neutral decay (we allow decay downstream of spectrometer), comparing with NA62MC
			ZMUV3 = 246.850;

			Sigacceptance = 1; // better call that "efficiency?"
			Sigacceptancemumu = 1;

			acceptanceHole = 0.1 ; // m, half the hole length
			acceptanceSide = 1.; // m
			//POT = 1.3*1.E16; // 1 day
			//POT = 1.*1.E18; // 2021-2023
			POT = 1.*1.E19; //2026-2030 period (possible x4 update beyond)
			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(6.*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 3.601*10^6 [pb]: sigma_cc @ 400 GeV (from Pythia)
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "[Error] Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 1 :
			expLabel = expLabels[Experiment_];
			expName = "CHARM";
			beamEnergy = 400;
			ATarget = 63.5;//ACu
			ZTarget = 29.; //ZCu
			ZTAX = 0;
			X0 = 0.;  // m
			Y0 = -4.8;// m since the angle is ~10 mrad, but the detector is parallel to the beam line
			Z0 = 0.; //was 20 m
			ZFVIn = 480.; //m
			ZLKR = 515.; //m
			ZFVEnd = ZLKR;
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; // ?????
			ZMUV3 = 515.; //  dummy, same as "LKr"

			Sigacceptance = 0.51; // for photons
			Sigacceptancemumu  = 0.85; // fof muons TO CHECK !!!

			POT = 2.4*1.E18;

			thetaMinInFile = 0.0025; 
			thetaMaxInFile = 0.0175; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 5.5;
			energyMaxInFile = 324.5;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(6.*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 3.601*10^6 [pb]: sigma_cc @ 400 GeV (from Pythia)
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 2:
			expLabel = expLabels[Experiment_];
			expName = "NuCal";
			beamEnergy = 70;
			ATarget = 55.8; //AFe
			ZTarget = 26.; //ZFe
			ZTAX = 0;
			X0 = 0.;  // m
			Y0 = 0.;  // m
			Z0 = ZTAX; //was 20 m
			ZFVIn = 64.; //m
			ZLKR = 23.+64.;
			ZFVEnd = ZLKR;
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR;
			ZMUV3 = 87.;

			Sigacceptance =0.7; // see Bluemlein for photons https://arxiv.org/pdf/1311.3870.pdf
			Sigacceptancemumu =0.8; // see Bluemlein for muons https://arxiv.org/pdf/1311.3870.pdf

			POT = 1.7*1.E18;

			thetaMinInFile = 0.00025; 
			thetaMaxInFile = 0.01525; // factor 5 w.r.t. shadows
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); //same binning as SHADOWS
			energyMinInFile = 1.1;
			energyMaxInFile = 64.9;
			energyWid = (energyMaxInFile-energyMinInFile)/(energyValuesInFile-1); // used to be log scale

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(5.1*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 1.15*1E-7/(38.38*1E9)*TMath::Power(ATarget,1./3); break; //38.38*10^9 [pb]: sigma_pp @ 70GeV (from sx file), 1.15*10^-7 [pb]: sigma_bb @ 70 GeV (from Pythia)
				case 6: normCrossSec = 1.15*1E-7/(38.38*1E9)*TMath::Power(ATarget,1./3); break; //38.38*10^9 [pb]: sigma_pp @ 70GeV (from sx file), 1.15*10^-7 [pb]: sigma_bb @ 70 GeV (from Pythia)
				case 7: normCrossSec = 0.1551*1E6/(38.38*1E9)*TMath::Power(ATarget,1./3); break; //38.38*10^9 [pb]: sigma_pp @ 70GeV (from sx file), 1.551*10^5 [pb]: sigma_cc @ 70 GeV (from Pythia)
	//			case 8: normCrossSec = 1/(38.38*1e9); break; //38.38*10^9 [pb]: sigma_pp @ 70GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 3:
			expLabel = expLabels[Experiment_];
			expName = "SHiP";
			beamEnergy = 400;
			ATarget = 95;//molybdenum?
			ZTarget = 42.; //molybdenum?
			ZTAX = 0;
			X0 = 0.;  // m
			Y0 = 0.;  // m straight for mathematica comparison
			Z0 = 0.; //was 20 m
			ZFVIn = 45.; //m
			zStraw1 = ZFVIn+50.76 ; // m, position of 1st straw chamber
			zStraw4 = ZFVIn+60.16; // m, position of 4th straw chamber
			ZMagnet     = ZFVIn + 53.45;//
			MagnetZLength       = 5.; // m // in figure 1b https://cds.cern.ch/record/2005715/files/main.pdf they also state 0.75 Tm
			MagnetFieldStrength = 0.15; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; //
			ZLKR = zStraw4 + 0.5; //m
			ZFVEnd = zStraw1;
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR;
			ZMUV3 = 1E10; // not condition

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

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(12.66*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 3.601*10^6 [pb]: sigma_cc @ 400 GeV (from Pythia)
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 4:
			expLabel = expLabels[Experiment_]; // DarkQuest as in 1804.00661v1
			expName = "DarkQuest";
			beamEnergy = 120;
			ATarget = 56;//IRON
			ZTarget = 26.; //IRON
			ZTAX = 0;
			X0 = 0.;  // m
			Y0 = 0.;  // m straight for mathematica comparison
			Z0 = 0.; //was 20 m
			ZFVIn = 5.; //m it used to be 7 according to text of 1804.00661v1, but now it is 5 accoring to 2008.08108
			zStraw1 = 6.16 ; // m, position of 1st tracker of DarkQuest (x-view), https://www.sciencedirect.com/science/article/pii/S016890021930347X
			zStraw4 =  18.79; // position of 3nd tracker x-view,
			ZMagnet     = 9;// // follow fig 1 in 2008.08108, check
			MagnetZLength       = 3.; // m //
			MagnetFieldStrength = 0.4; // Tesla, see https://www.sciencedirect.com/science/article/pii/S016890021930347X
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; //
			// assume that the second tracker is also hit when the third is hit, so 2nd is not included explcitye
			// 13.47; // m, position of 2nd tracker x-view, https://www.sciencedirect.com/science/article/pii/S016890021930347X
			ZLKR = 18.5;// 25.; //m (total length is 25 m but Calo needs to go before wall)
			// THIS IS PHASE 1, in phase 2 the fiducial region would be larger
			ZFVEnd = zStraw1; // it used to be 8, see fig 1 in 1804.00661v1 [hep-ph], but now it is 6 accoring to 2008.08108
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; /// gamma gamma
			ZMUV3 = 21.; // according to 2008.08108

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

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(5.1*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 3.8/(38.54*1E9)*TMath::Power(ATarget,1./3); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 3.8 [pb]: sigma_bb @ 120 GeV (from Pythia)
				case 6: normCrossSec = 3.8/(38.54*1E9)*TMath::Power(ATarget,1./3); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 3.8 [pb]: sigma_bb @ 120 GeV (from Pythia)
				case 7: normCrossSec = 0.518*1E6/(39.54*1E9)*TMath::Power(ATarget,1./3); break; //39.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 0.518*10^6 [pb]: sigma_cc @ 120 GeV (from Pythia)
//				case 8: normCrossSec = 1/(38.54*1e9); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 5:
			// based on https://arxiv.org/pdf/2103.13910.pdf and https://arxiv.org/pdf/2011.05995.pdf
			expLabel = expLabels[Experiment_];
			expName = "DUNE";
			beamEnergy = 120;
			ATarget = 12.;//Graphite
			ZTarget = 6.; //Graphite
			ZTAX = 0;
			X0 = 0.;
			Y0 = 0.; //by default but detector can move off-axis up to 30.5m (see DUNE-PRISM)
			Z0 = 0.;
			ZFVIn = 574.5; // front of LAr
			ZLKR = 580;// front of GAr
			ZFVEnd = 584.; // end of GAr
			ZMUV3 = ZFVEnd;

			Sigacceptance = 1;
			Sigacceptancemumu =1; //??

			POT = 1.1*1.E21; //1 year of data taking
//				POT = 1.1*1.E22; //10 years of data taking?
			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(0.35*1e12); break; //from https://www.nist.gov/pml/xcom-photon-cross-sections-database
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 3.8/(38.54*1E9)*TMath::Power(ATarget,1./3); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 3.8 [pb]: sigma_bb @ 120 GeV (from Pythia)
				case 6: normCrossSec = 3.8/(38.54*1E9)*TMath::Power(ATarget,1./3); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 3.8 [pb]: sigma_bb @ 120 GeV (from Pythia)
				case 7: normCrossSec = 0.518*1E6/(39.54*1E9)*TMath::Power(ATarget,1./3); break; //39.54*10^9 [pb]: sigma_pp @ 120GeV (from sx file), 0.518*10^6 [pb]: sigma_cc @ 120 GeV (from Pythia)
//				case 8: normCrossSec = 1/(38.54*1e9); break; //38.54*10^9 [pb]: sigma_pp @ 120GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 6:
			//based on PBC proposal https://indico.cern.ch/event/1002356/contributions/4229628/attachments/2201710/3781990/shadows_PBC.pdf
			expLabel = expLabels[Experiment_];
			expName = "SHADOWS";
			beamEnergy = 400;
			PhiBeamEuler = TMath::Pi()/2; //SHADOWS is offaxis in X axis
			ThetaBeamEuler = 0.0624188; //(~62mrad) for calorimeter located 36m downstream the target and 2.25m shifted center w.r.t. beam axis
			PsiBeamEuler = 0;
			ATarget = 63.5; //ACu
			ZTarget = 29; //ZCu
			ZTAX = 23.070; //NA62
			X0 = 0.;//-(1+1.25);  //shifted 1m wrt NA62, 2.5m calo size (also option with 3m size in consideration)
			Y0 = -0.022; //as for NA62
			Z0 = ZTAX; //NA62
			ZFVIn = ZTAX + 10.; //10m from target
			zStraw1 = ZFVIn + 20. + 6. - 3. ; // m, position of 1st straw chamber (decay volume end + 6m Calo location - 3m spectrometer length)
			zStraw4 = zStraw1 + 2.5; // m, position of 4th straw chamber (0.5m between 1,2 and 3,4 station and 1.5m magnet)
			ZMagnet     = 0.5*(zStraw1+zStraw4);
			MagnetZLength       = 1.5; // m
			MagnetFieldStrength = 1.; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = zStraw4 + 0.5; //~50cm behind 4th tracker station
			ZFVEnd = zStraw1; //20(+3?) m decay volume for charged particles
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR + 1.; //1m distance between calo and muon detector as measured from the image in the SHADOWS paper

			Sigacceptance = 1;
			Sigacceptancemumu = 1;

			POT = 1.*1.E19; //2026-2030 period (possible x4 update beyond)

			thetaMinInFile = 0.012;
			thetaMaxInFile = 0.12;
			thetaValuesInFile = 31;// check: one more than steps in ALP Mathematica output
			thetaWid = (thetaMaxInFile-thetaMinInFile)/(thetaValuesInFile-1); // used to be log scale
			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(6.*1e12); break;
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
				case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 1.866*10^3 [pb]: sigma_bb @ 400 GeV (from Pythia)
				case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from sx file), 3.601*10^6 [pb]: sigma_cc @ 400 GeV (from Pythia)
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 7:
			// based on presentation https://indico.cern.ch/event/1055867/contributions/4437601/attachments/2280165/3874009/darkSectorKOTOv1.pdf
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			ATarget = 197; // Gold (AU)
			ZTarget = 79; //
			ZTAX = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56;  //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTAX;
			ZFVIn = ZTAX + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			zStraw1 = 0; // no straw
			zStraw4 = 0; // no straw
			ZMagnet   = 0.5*(zStraw1+zStraw4); // no straw
			MagnetZLength       = 0; // m
			MagnetFieldStrength = 0; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 27;
			ZFVEnd = ZLKR; //20(+3?) m decay volume for charged particles
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR + 1.; //summy

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
			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(38.1*1e12); break; // for Au, 4GeV photons
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
	//			case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
	//			case 8: normCrossSec = 1; break;
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 8:
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			ATarget = 197; // Gold (AU)
			ZTarget = 79; //
			ZTAX = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56;  //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTAX;
			ZFVIn = ZTAX + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			zStraw1 = 0; // no straw
			zStraw4 = 0; // no straw
			ZMagnet   = 0.5*(zStraw1+zStraw4); // no straw
			MagnetZLength       = 0; // m
			MagnetFieldStrength = 0; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 27;
			ZFVEnd = ZLKR; //20(+3?) m decay volume for charged particles
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR + 1.; //summy

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

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(38.1*1e12); break; // for Au, 4GeV photons
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
	//			case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 8: normCrossSec = 1; break;
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 9:
			expLabel = expLabels[Experiment_];
			expName = "KOTO2";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.087266; //5∘ = 87.266 mrad
			PsiBeamEuler = 0;
			ATarget = 197; // Gold (AU)
			ZTarget = 79; //
			ZTAX = 0; //presume that they dumped on T1 target
			X0 = 0.;//
			Y0 = 0.;
			Z0 = ZTAX;
			ZFVIn = ZTAX + 43 + 1 + 1.75; // 3-m long beam line and 1-m long space; The Front Barrel Counter is 1.75-m
			zStraw1 = 0; // no straw
			zStraw4 = 0; // no straw
			ZMagnet   = 0.5*(zStraw1+zStraw4); // no straw
			MagnetZLength       = 0; // m
			MagnetFieldStrength = 0; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 44 + 20;
			ZFVEnd = ZLKR; 
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR; //summy

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
			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(38.1*1e12); break; // for Au, 4GeV photons
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
	//			case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
	//			case 8: normCrossSec = 1; break;
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 10:
			expLabel = expLabels[Experiment_];
			expName = "KOTO2";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.087266; //5∘ = 87.266 mrad
			PsiBeamEuler = 0;
			ATarget = 197; // Gold (AU)
			ZTarget = 79; //
			ZTAX = 0; //presume that they dumped on T1 target
			X0 = 0.;//
			Y0 = 0.;
			Z0 = ZTAX;
			ZFVIn = ZTAX + 43 + 1 + 1.75; // 3-m long beam line and 1-m long space; The Front Barrel Counter is 1.75-m
			zStraw1 = 0; // no straw
			zStraw4 = 0; // no straw
			ZMagnet   = 0.5*(zStraw1+zStraw4); // no straw
			MagnetZLength       = 0; // m
			MagnetFieldStrength = 0; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 44 + 20;
			ZFVEnd = ZLKR; 
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR; //summy

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

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(38.1*1e12); break; // for Au, 4GeV photons
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
	//			case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
	//			case 8: normCrossSec = 1; break;
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;
		case 11:
			expLabel = expLabels[Experiment_];
			expName = "KOTO";
			beamEnergy = 30;
			PhiBeamEuler = 0;
			ThetaBeamEuler = 0.279253; //16∘ = 279.253 mrad
			PsiBeamEuler = 0;
			ATarget = 197; // Gold (AU)
			ZTarget = 79; //
			ZTAX = 0; //presume that they dumped on T1 target
			X0 = 0.;//7.56;  //The KL beam line is off-axis by an angle of 16∘ = 280mrad with respect to the primary proton beam line. And the whole thing is l=27m long, we compute x=l\times theta. should we put a -sign here?
			Y0 = 0.;
			Z0 = ZTAX;
			ZFVIn = ZTAX + 21. + 2.9; // The front of the KOTO detector is located at the end of the beam line, about 21 m from the T1 target. + 2.9 front barrel
			zStraw1 = 0; // no straw
			zStraw4 = 0; // no straw
			ZMagnet   = 0.5*(zStraw1+zStraw4); // no straw
			MagnetZLength       = 0; // m
			MagnetFieldStrength = 0; // Tesla
			MagKick = MagnetFieldStrength*MagnetZLength*0.3; // GeV
			ZLKR = 27;
			ZFVEnd = ZLKR; //20(+3?) m decay volume for charged particles
			if(chargeFinState[0] == 0)	ZFVEnd = ZLKR; //for photons
			ZMUV3 = ZLKR + 1.; //summy

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

			switch(ProductionMode_){
				case 0: normCrossSec = 1/(53.*TMath::Power(ATarget,0.77)*1E9); break;
				case 1: normCrossSec = 1/(38.1*1e12); break; // for Au, 4GeV photons
				case 2: normCrossSec = 1; break;
				case 3: normCrossSec = 1; break;
				case 4: normCrossSec = 1; break;
	//			case 5: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 6: normCrossSec = 1.866*1E3/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 7: normCrossSec = 3.601*1E6/(39.85*1E9)*TMath::Power(ATarget,1./3); break; //TO BE ADAPTED
	//			case 8: normCrossSec = 1; break;
	//			case 8: normCrossSec = 1/(39.85*1e9); break; //39.85*10^9 [pb]: sigma_pp @ 400GeV (from meson sx file)
				default: std::cout << "Invalid production mode: " <<  ProductionMode_ << " for experiment " << expName << std::endl; std::exit(1);
			}
			break;

		default:
			std::cout << "Invalid experiment label: " <<  Experiment_ << ". Use values ";
			for(Int_t iExp=0;iExp<expLabels.size();iExp++) std::cout << iExp << " (" << expLabels[iExp] << "), ";
			std::cout << std::endl;
			exit(1);

	}
	ZDist = ZFVIn-ZTAX;
	BeamDecayLength = ZFVEnd - ZFVIn; // m
}

ExpParameters::~ExpParameters(){
	if(rndm) delete rndm;
}

Bool_t ExpParameters::inCaloAcceptance(TVector3 posFVEnd, TVector3 posLKr, Double_t zDecay){
	// experiment = 0 --> NA62
	// experiment = 1 --> CHARM
	// experiment = 2 --> nucal U70
	// experiment = 3 --> SHiP
	// experiment = 4 --> DarkQuest 2m x 2m
    // experiment = 5 --> DUNE
    // experiment = 6 --> SHADOWS
    // experiment = 7 --> KOTO
	Double_t xa = fabs(posLKr.X()*1000.); //in mm
	Double_t ya = fabs(posLKr.Y()*1000.); //in mm

	if (expNum == 0) {
		if (xa > 1130)return kFALSE;
		if (ya > 1130)return kFALSE;
		if (xa+ya > 1598) return kFALSE;
		if (xa > 632 && ya > 847) return kFALSE; // !83.7?
		if (ya > 947 && xa > 522) return kFALSE;
		if(pow(xa*xa+ya*ya,0.5)<80) return kFALSE;
	}
	else if (expNum == 1) { // this works since you offset the proton beam
		if (xa > 1500.) return kFALSE; //for CHARM 3 m by 3 m
		if (ya > 1500.) return kFALSE; //for CHARM 3 m by 3 m
	}
	else if (expNum ==2) { //nucal
		if(pow(xa*xa+ya*ya,0.5)>1300) return kFALSE;
	}
	else if (expNum == 3) { // ship // new setup has no ellipse anymore
		// Double_t xellip = 2500;
		// Double_t yellip = 5000;
		// Double_t rad= pow(xa/xellip,2) +pow(ya/yellip,2);
		// if (rad > 1.) return kFALSE; //for ship
		if (xa > 5000.) return kFALSE; // 10 m by
		if (ya > 2500.) return kFALSE; // 5 m

	}
	else if (expNum == 4) {
		if (xa > 1000.) return kFALSE; // 2 m by
		if (ya > 1000.) return kFALSE; // 2 m
	}
	else if (expNum == 5) { //DUNE
		if((zDecay > ZFVIn) && (zDecay < ZLKR )){ //LAr
			if (abs(posLKr.X()) > 3.) return kFALSE;
			if (abs(posLKr.Y()) > 1.) return kFALSE;
		} else if((zDecay > ZLKR) && (zDecay < ZFVEnd)) { //DUNE GAr: 2x4.8 m
			if (abs(posFVEnd.X()) > 2.4) return kFALSE;
			if (abs(posFVEnd.Y()) > 1.) return kFALSE;
		} else return kFALSE;
	}
	else if (expNum == 6) { //SHADOWS
		if (xa > 1250.) return kFALSE; // 2.5 m by
		if (ya > 1250.) return kFALSE; // 2.5 m
	}
	else if (expNum == 7 || expNum == 8 || expNum == 11) { //KOTO
		if(pow(xa*xa+ya*ya,0.5)>1000) return kFALSE; // The CsI calorimeter consists of 2716 undoped CsI crystals stacked in a cylindrical shape of 2-m diameter
		if (abs(xa) < 75. && abs(ya) < 75.) return kFALSE; //A 15 x 15 cm2 center region of the stacked CSI permitted the beam particles to pass through the detector.
	}
	else if (expNum == 9 || expNum == 10) { //KOTO2
		if(pow(xa*xa+ya*ya,0.5)>1500) return kFALSE; //3m diameter
		if (abs(xa) < 100. && abs(ya) < 100.) return kFALSE; //20x20cm central hole
	}
	// else if (expNum == 7) {
	// 	if (xa > 100000.) return kFALSE; // 200 m by
	// 	if (ya > 100000.) return kFALSE; // 200 m
	// }
	return kTRUE;
}

Bool_t ExpParameters::inCaloOuterEdgeAcceptance(TVector3 posLKr){
	// experiment = 0 --> NA62
	// experiment = 1 --> CHARM
	// experiment = 2 --> nucal U70
	// experiment = 3 --> SHiP
	// experiment = 4 --> DarkQuest 2m x 2m
    // experiment = 5 --> DUNE
    // experiment = 6 --> SHADOWS
    // experiment = 7 --> KOTO
	Double_t xa = fabs(posLKr.X()*1000.); //in mm
	Double_t ya = fabs(posLKr.Y()*1000.); //in mm

	if (expNum == 0) {
		if (xa > 1130)return kFALSE;
		if (ya > 1130)return kFALSE;
		if (xa+ya > 1598) return kFALSE;
		if (xa > 632 && ya > 847) return kFALSE; // !83.7?
		if (ya > 947 && xa > 522) return kFALSE;
	}
	else if (expNum == 1) { // this works since you offset the proton beam
		if (xa > 1500.) return kFALSE; //for CHARM 3 m by 3 m
		if (ya > 1500.) return kFALSE; //for CHARM 3 m by 3 m
	}
	else if (expNum ==2) { //nucal
		if(pow(xa*xa+ya*ya,0.5)>1300) return kFALSE;
	}
	else if (expNum == 3) { // ship // new setup has no ellipse anymore
		// Double_t xellip = 2500;
		// Double_t yellip = 5000;
		// Double_t rad= pow(xa/xellip,2) +pow(ya/yellip,2);
		// if (rad > 1.) return kFALSE; //for ship
		if (xa > 5000.) return kFALSE; // 10 m by
		if (ya > 2500.) return kFALSE; // 5 m

	}
	else if (expNum == 4) {
		if (xa > 1000.) return kFALSE; // 2 m by
		if (ya > 1000.) return kFALSE; // 2 m
	}
	else if (expNum == 5) { //DUNE
		if (abs(posLKr.X()) > 3.) return kFALSE;
		if (abs(posLKr.Y()) > 1.) return kFALSE;
	}
	else if (expNum == 6) { //SHADOWS
		if (xa > 1250.) return kFALSE; // 2.5 m by
		if (ya > 1250.) return kFALSE; // 2.5 m
	}
	else if (expNum == 7 || expNum == 8 || expNum == 11) { //KOTO
		if(pow(xa*xa+ya*ya,0.5)>1000) return kFALSE; // The CsI calorimeter consists of 2716 undoped CsI crystals stacked in a cylindrical shape of 2-m diameter
	}
	else if (expNum == 9 || expNum == 10) { //KOTO2
		if(pow(xa*xa+ya*ya,0.5)>1500) return kFALSE; //3m diameter
	}
	// else if (expNum == 7) {
	// 	if (xa > 100000.) return kFALSE; // 200 m by
	// 	if (ya > 100000.) return kFALSE; // 200 m
	// }
	return kTRUE;
}

Bool_t ExpParameters::inSpectrometerAcceptance(TVector3 posStraw){

	Double_t strawRadius =pow(pow(posStraw.X(),2) + pow(posStraw.Y(),2),0.5);
	if(expNum == 0){ // NA62 is round
		if((acceptanceHole < strawRadius) && (strawRadius<acceptanceSide)) return kTRUE;
	}
	if(expNum == 1) return kFALSE; //no straws
	if(expNum == 2) return kFALSE; //no straws
	if(expNum == 3){ // ship is a rectangle 5 by 10
		if((abs(posStraw.X())<2.5) && (abs(posStraw.Y())<5))	return kTRUE;
	}
	if(expNum == 4){ // DarkQuest is a rectangle, 2 by 2
		if((abs(posStraw.X())<1) && (abs(posStraw.Y())<1))	return kTRUE;
	}
	if(expNum == 5) return kFALSE; //no straws
	if(expNum == 6){ // SHADOWS is a rectangle, 2.5 by 2.5
		if((abs(posStraw.X())<1.25) && (abs(posStraw.Y())<1.25))	return kTRUE;
	}

	return kFALSE;
}

Bool_t ExpParameters::inMuonVetoAcceptance(TVector3 posFVEnd, TVector3 posMUV, Double_t zDecay){
	Double_t muonRadius =pow(pow(posMUV.X(),2) + pow(posMUV.Y(),2),0.5);
	if(expNum == 0){ // NA62 MUV3: 2.64 x 2.64 m with 0.22 diameter (central tile) hole
		if((abs(posMUV.X())<1.32) && (abs(posMUV.Y())<1.32) && (muonRadius>0.11)) return kTRUE;
	}
	else if(expNum == 1) { // CHARM has no straws. Sinmply ask both muons within 3x3 m at end of decay volume?
		if((abs(posMUV.X())<1.5) && (abs(posMUV.Y())<1.5)) return kTRUE; // both Muons within 3x3m. Distance requirement?
	}
	else if(expNum == 2) { // nucal, see https://arxiv.org/abs/1311.3870
		if(pow(pow(posMUV.X(),2) + pow(posMUV.Y(),2),0.5)<1.3) return kTRUE; // follow bluemlein
	}
	if(expNum == 3){ // ship, assume the same size as calo: 5 by 10
		if((abs(posMUV.X())<2.5) && (abs(posMUV.Y())<5))	return kTRUE;
	}
	if(expNum == 4){ // DarkQuest, assume the same size as calo: 2 by 2
		if((abs(posMUV.X())<1) && (abs(posMUV.Y())<1))	return kTRUE;
	}
	else if(expNum == 5) { // DUNE
		if((zDecay > ZFVIn) && (zDecay < ZLKR )){ //decays in DUNE LAr
			if((abs(posMUV.X())<3.) && (abs(posMUV.Y())<1.))	return kTRUE;
		} else if((zDecay > ZLKR) && (zDecay < ZFVEnd)){ //decays in DUNE GAr
			if((abs(posFVEnd.X())<2.4) && (abs(posFVEnd.Y())<1.))	return kTRUE;
		}
	}
	if(expNum == 6){ // SHADOWS, assume the same size as calo: 2.5 by 2.5
		if((abs(posMUV.X())<1.25) && (abs(posMUV.Y())<1.25))	return kTRUE;
	}

	return kFALSE;
}

Int_t ExpParameters::twoPhotonCondition(Int_t iOnLKr, Double_t totalEnergyAcceptance, Double_t zDecay, TVector3 posFVEnd[2], TVector3 posLKr[2], TLorentzVector pgA[2]){
	TVector3 diffLKr = posLKr[1]-posLKr[0];
	Double_t distance = diffLKr.Mag();
	Double_t radiuspos0 = posLKr[0].XYvector().Mod();
	Double_t radiuspos1 = posLKr[1].XYvector().Mod();
	Double_t energy0 = pgA[0].E();
	Double_t energy1 = pgA[1].E();

	if(expNum == 0){ //NA62
		if(iOnLKr==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3 && radiuspos1> 0.15 && radiuspos0> 0.15) return 1; //update
		return 0;
	}
	else if(expNum == 1){ //CHARM
		if(iOnLKr==1 && energy0>5. && energy0<50.) return 1;
		else if(iOnLKr==2 && energy0>5. && energy1 > 5. && energy0<50. && energy1 < 50.) return 1; //iOnLKr>=1 && ; improved CHARM
		return 0;
	}
	else if(expNum == 2){ //NuCal
		if((iOnLKr==1 && (totalEnergyAcceptance>10)) || (iOnLKr ==2 && ((energy0>10.) || (energy1>10.)))) return 1; // NuCal
		return 0;
	}
	else if(expNum == 3){ //SHiP
		if(iOnLKr==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; //ship
		return 0;
	}
	else if(expNum == 4){ //DarkQuest
		if(iOnLKr==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; // DarkQuest
		return 0;
	}
	else if(expNum == 5){ //DUNE
		if((zDecay > ZFVIn) && (zDecay < ZLKR )){ //LAr
			if(iOnLKr==2) return 1;
			return 0;
		} else if((zDecay > ZLKR) && (zDecay < ZFVEnd)){ //GAr
			if(iOnLKr==2 && energy0>0.02 && energy1 > 0.02 && pgA[0].Angle(pgA[1].Vect()) > 0.14) return 1; //at least 20MeV and 0.8deg separation angle
			return 0;
		}
		return 0;
	}
	else if(expNum == 6){ //SHADOWS
		if(iOnLKr==2 && distance > 0.1 && energy0>1. && energy1 > 1. && totalEnergyAcceptance >3) return 1; // SHADOWS
		return 0;
	}
	else if(expNum == 7){ //KOTOdump must resolve lower gamma energies
		// `dream cuts` similar to first attempt, see #34 on GITHUB.
		// for the photon cluster distance cut we follow the thesis by Melissa Hutcheson (30cm)
		if(iOnLKr==2 && distance > 0.3 && energy0>0.05 && energy1 > 0.05) return kTRUE; // 50MeV

		return 0;
	}
	else if(expNum == 8){ //KOTO
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnLKr==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posLKr[0].X()),abs(posLKr[0].Y())) > 0.15 && TMath::Min(abs(posLKr[1].X()),abs(posLKr[1].Y())) > 0.15)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posLKr[0].X()*energy0+posLKr[1].X()*energy1)/(energy0+energy1),
			(posLKr[0].Y()*energy0+posLKr[1].Y()*energy1)/(energy0+energy1),
			posLKr[0].Z()
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
		Double_t zDecayPnnInDV = ZLKR-zVtxDistPnn-21; //21m - beginning of DV in pnn analysis

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
	else if(expNum == 9){ //KOTO2dump
		if(iOnLKr==2 && distance > 0.3 && energy0>0.1 && energy1 > 0.1 && (energy0 + energy1) > 0.5) return 1; // cuts 1-5 (p234) from Pnn ensuring photon cluster quality
		return 0;
	}
	else if(expNum == 10){ //KOTO2
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnLKr==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0; //stays the same as for KOTO1
		if(!(radiuspos0 < 1.35 && radiuspos1 < 1.35 && distance > 0.3 && TMath::Min(abs(posLKr[0].X()),abs(posLKr[0].Y())) > 0.175 && TMath::Min(abs(posLKr[1].X()),abs(posLKr[1].Y())) > 0.175)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.5 && energy0 > 0.1 && energy1 > 0.1)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posLKr[0].X()*energy0+posLKr[1].X()*energy1)/(energy0+energy1),
			(posLKr[0].Y()*energy0+posLKr[1].Y()*energy1)/(energy0+energy1),
			posLKr[0].Z()
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

		Double_t zDecayPnnInDV = ZLKR-zVtxDistPnn-44; //the Pnn decay volume starts after 44m beam line

		TVector3 vectCOE(posCOE.X(),posCOE.Y(),zVtxDistPnn);
		TLorentzVector pCOE;
		Double_t momPi0 = TMath::Sqrt(TMath::Power(energy0 + energy1, 2) - MPi0*MPi0);
		pCOE.SetPxPyPzE(vectCOE.X()/vectCOE.Mag()*momPi0,vectCOE.Y()/vectCOE.Mag()*momPi0,vectCOE.Z()/vectCOE.Mag()*momPi0,energy0+energy1);

		if( zDecayPnnInDV < 15 && zDecayPnnInDV > 1.75 && pCOE.Pt() > 0.4) //exclude the 1.75m long front barrel counter
			return 1;
		else
			return 0;

	}
	else if(expNum == 11){ //KOTO Pnn region excluded
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(!(iOnLKr==2)) return 0; //2 clusters passing inCaloAcceptance()
		if(angleSepXY*180/TMath::Pi() > 150) return 0;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posLKr[0].X()),abs(posLKr[0].Y())) > 0.15 && TMath::Min(abs(posLKr[1].X()),abs(posLKr[1].Y())) > 0.15)) return 0; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return 0; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posLKr[0].X()*energy0+posLKr[1].X()*energy1)/(energy0+energy1),
			(posLKr[0].Y()*energy0+posLKr[1].Y()*energy1)/(energy0+energy1),
			posLKr[0].Z()
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
		Double_t zDecayPnnInDV = ZLKR-zVtxDistPnn-21; //21m - beginning of DV in pnn analysis

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
	return 0;
}

Bool_t ExpParameters::multiplePhotonCondition(Int_t nGammas, Double_t totalEnergyAcceptance, Double_t zDecay, TVector3 posFVEnd[6], TVector3 posLKr[6], TLorentzVector pgA[6]){

// evaluate number of distinct gammas seen
	double clusterMerginDistance = 0.1;
	if(expNum == 7 || expNum == 9)
		clusterMerginDistance = 0.3; //different threshold for KOTO
	int flagMerged[6] = {-1,-1,-1,-1,-1,-1};
	for (Int_t k1=0; k1<nGammas; k1++){
		if (!inCaloAcceptance(posFVEnd[k1],posLKr[k1],zDecay)) continue;
		if (flagMerged[k1] >= 0) continue;

		for (Int_t k2=k1+1; k2<nGammas; k2++){
			if (!inCaloAcceptance(posFVEnd[k2],posLKr[k2],zDecay)) continue;
			if (flagMerged[k2] >= 0) continue;

			if ((posLKr[k1]-posLKr[k2]).Mag() < clusterMerginDistance) {
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

	int iOnLKr = 0;
	int i3OnLKr = 0;
	int i5to50OnLKr = 0;
	int indexAbove3[6] = {0,0,0,0,0,0};
	int index5to50[6] = {0,0,0,0,0,0};
	for (Int_t k1=0; k1<nGammas; k1++){
		if (!inCaloAcceptance(posFVEnd[k1],posLKr[k1],zDecay)) continue;
		if (pgA[k1].E() > 3){
			indexAbove3[i3OnLKr] = k1;
			i3OnLKr++; // for nuCal
		}
		if (pgA[k1].E() > 5 && pgA[k1].E()<50){
			index5to50[i5to50OnLKr] = k1;
			i5to50OnLKr++; // for charm
		}
		if (flagMerged[k1] >= 0) continue;
  		if (expNum == 5 && zDecay > ZLKR && zDecay < ZFVEnd && pgA[k1].E() < 0.02) continue; //DUNE GAr energy threshold
		if ((expNum == 7 && pgA[k1].E() < 0.05) || (expNum == 9 && pgA[k1].E() < 0.1)) continue; //KOTO(2)dump single cluster energy threshold
		iOnLKr++;
	}

	//exp. conditions
	if(expNum == 0 || expNum == 3 || expNum == 4 || expNum == 6){ //NA62, SHiP, SHADOWS, DarkQuest
		if(nGammas==iOnLKr && totalEnergyAcceptance>3) return kTRUE; // a minimum requirement on the individual gamma energy should be probably added here ?
		return kFALSE;
	}
	else if(expNum == 1){
		if(i5to50OnLKr==1 || i5to50OnLKr==2) return kTRUE;//CHARM, as in digamma, but not more than two photons allowed
		return kFALSE;
	}
	else if(expNum == 2){ //NuCal
		Int_t i3OnLKrdetected=0;
		Double_t totalEnergyDetected=0;
		for(Int_t n3onLKr=0; n3onLKr< i3OnLKr; n3onLKr++){
			// take into account efficiency of photon detection
			if(rndm->Rndm()>0.7){
				i3OnLKrdetected++;
				totalEnergyDetected+=pgA[indexAbove3[n3onLKr]].E();
			}
		}
		// mimic the following statement:
		//Requiring events with an electromagnetic shower of Eel m > 3 GeV, 106 events remain. A large fraction of them has an energetic hadron shower. This part can be removed by excluding all events with hadron energies above 1.5 GeV. For the remaining 21 events the electromagnetic shower energy is shown in Fig. 3a.
		if(i3OnLKrdetected>=1 &&totalEnergyDetected >10) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 5){ //DUNE
		if(nGammas==iOnLKr) return kTRUE; //we drop the angle conditions
		return kFALSE;
	}
	else if(expNum == 7){ //KOTOdump
		if(nGammas==iOnLKr && nMerged == 0 && totalEnergyAcceptance>0.5) return kTRUE; // to be checked
		return kFALSE;
	}
	else if(expNum == 8){ //KOTO - same conditions as for two photons but require 2 clusters (or more photons in acceptance if clusters overlap)
		if(nClusters != 2) return kFALSE;
		std::vector<Int_t> realIndex;
		for (Int_t k1=0; k1<nGammas; k1++) if(flagMerged[k1]==-1) realIndex.push_back(k1);
		if(realIndex.size()!=nClusters){
			std::cout << "[KOTO multiplePhotonCondition] Error: unexpected number of calorimeter clusters" << std::endl;
			exit(1);
		}
		TVector3 diffLKr = posLKr[realIndex[1]]-posLKr[realIndex[0]];
		Double_t distance = diffLKr.Mag();
		Double_t radiuspos0 = posLKr[realIndex[0]].XYvector().Mod();
		Double_t radiuspos1 = posLKr[realIndex[1]].XYvector().Mod();
		Double_t energy0 = pgA[realIndex[0]].E();
		Double_t energy1 = pgA[realIndex[1]].E();
		Double_t energyH = TMath::Max(energy0,energy1);
		Double_t energyL = TMath::Min(energy0,energy1);
		Double_t angleSepXY = TMath::ACos((radiuspos0*radiuspos0+radiuspos1*radiuspos1-distance*distance)/(2*radiuspos0*radiuspos1));
		///// KOTO cuts aplicable directly for ALPs
		if(angleSepXY*180/TMath::Pi() > 150) return kFALSE;
		if(!(radiuspos0 < 0.85 && radiuspos1 < 0.85 && distance > 0.3 && TMath::Min(abs(posLKr[realIndex[0]].X()),abs(posLKr[realIndex[0]].Y())) > 0.15 && TMath::Min(abs(posLKr[realIndex[1]].X()),abs(posLKr[realIndex[1]].Y())) > 0.15)) return kFALSE; //additional geometry cuts
		if(!(totalEnergyAcceptance > 0.65 && energy0 > 0.1 && energy0 < 2. && energy1 > 0.1 && energy1 < 2. && energyL/energyH > 0.2)) return kFALSE; //energy cuts

		///// KOTO cuts requiring reconstruction of Pnn-like event
		TVector3 posCOE(
			(posLKr[realIndex[0]].X()*energy0+posLKr[realIndex[1]].X()*energy1)/(energy0+energy1),
			(posLKr[realIndex[0]].Y()*energy0+posLKr[realIndex[1]].Y()*energy1)/(energy0+energy1),
			posLKr[realIndex[0]].Z()
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
		Double_t zDecayPnnInDV = ZLKR-zVtxDistPnn-21;

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

	else if(expNum == 9){ //KOTO2dump
		if(nGammas==iOnLKr && nMerged == 0 && totalEnergyAcceptance>0.5) return kTRUE; // to be checked
		return kFALSE;
	}
	return kFALSE;
}

Bool_t ExpParameters::twoMuonCondition(Int_t iOnStraw1, Int_t iOnStraw4, Int_t iOnLKr, Int_t iOnMUV3, Double_t zDecay, TLorentzVector pCA[2]){
	if(expNum == 0|| expNum == 3 || expNum == 6){ // NA62, ship and SHADOWS
		if( pCA[0].E() > 5 && pCA[1].E() > 5  && iOnStraw1 == 2 && iOnStraw4 == 2) return kTRUE; //similar spectrometers
		return kFALSE;
	}
	else if(expNum == 1){ //charm: 1 GeV mimimum: Mip should be seen throughout detector material (this is still generous)
		if(iOnMUV3>1 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 2){ // according to bluemlein
		if(iOnLKr==2 && (pCA[0].E()+pCA[1].E())>10) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 4){ // DarkQuest, follow Fig 14 in "The DarkQuest spectrometer at Fermilab"  https://www.sciencedirect.com/science/article/abs/pii/S016890021930347X
	// also there is a 1m iron absorber which takes 1GeV off the muons
		if(iOnStraw1 == 2 && iOnStraw4 == 2 && (pCA[0].E()+pCA[1].E())>4.2 && pCA[0].E() > 1 && pCA[1].E() > 1) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 5){ //seeing in https://arxiv.org/pdf/2103.13910.pdf page 2-71, muon energies below 1GeV are fully OK. Above the angle wr.t. to the beam axis has to be smaller than 40degrees (700mrad)
		if(iOnMUV3 == 2){
			if((pCA[0].E()-pCA[0].M()) < 1 && (pCA[1].E()-pCA[1].M()) < 1) return kTRUE;
			if(pCA[0].Theta()<0.7 && pCA[1].Theta()<0.7) return kTRUE;
		}
		else return kFALSE;
	} // exp 5
	return kFALSE;
}

Bool_t ExpParameters::twoHadronCondition(Int_t iOnStraw1, Int_t iOnStraw4, Int_t iChOnLKr, Int_t iOnMUV3, Double_t zDecay, TLorentzVector pCA[6]){
	if(expNum == 0|| expNum == 3 || expNum == 6){ // NA62, ship and SHADOWS
		if(pCA[0].E() > 5 && pCA[1].E() > 5  && iOnStraw1 == 2 && iOnStraw4 == 2 && iChOnLKr == 2) return kTRUE;
		return kFALSE;
	}
	else if(expNum == 1){
		if(iChOnLKr==0) return kTRUE; //Charm vetos hadronic showers in (https://inspirehep.net/literature/214233)
		return kFALSE;
	}
	else if(expNum == 2){
		if(iChOnLKr==0) return kTRUE;  //Nucal explicitely comments on vetoing hadronic showers above 1.5GeV
		return kFALSE;
	}
	else if(expNum == 4){ // DarkQuest. in phase 1 0 sensi to hadrons because of background. will chance for phase 2, see paper 2008.08108
		//return kTRUE;
		//if(iOnStraw1 == 2 && iOnStraw4 == 2 && (pCA[0].E()+pCA[1].E())>4.2) return kTRUE; // mimic the muon conditions from twoMuonCondition
		//if(iOnMUV3 == 2 && iChOnLKr>=2 && iOnStraw1 == 2 && iOnStraw4 == 2 && totalEnergyAcceptance > energyCut[experiment]) kTRUE; // DarkQuest <-??
		return kFALSE;
	}
	else if(expNum == 5){
		if((zDecay > ZFVIn) && (zDecay < ZLKR )){ //LAr
			if(iChOnLKr==2) return kTRUE;
			return kFALSE;
		} else if((zDecay > ZLKR) && (zDecay < ZFVEnd)){ //GAr
			if(iChOnLKr==2 && (pCA[0].E() - mFinState[0]) > 0.005 && (pCA[1].E() - mFinState[1]) > 0.005) return kTRUE; //kinetic energy above 5MeV
			return kFALSE;
		}
		return kFALSE;
	}
	return kFALSE;
}
