#include <iostream>
#include <array>
#include <fstream>
#include "TError.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include <TROOT.h>
#include "TGenPhaseSpace.h"

#include "ExpParameters.h"
//#include "DecayMCGlobal.h"

using namespace std;

Bool_t allInAcceptance;
Bool_t verbose;

//functions
void decayTo2Body(ExpParameters &, AxionParameters &, TH1D &, std::vector<TH1D*> &, TRandom3 *);
void decayTo3Body(ExpParameters &, AxionParameters &, TH1D &, std::vector<TH1D*> &, TH2D*, TRandom3 *, TH2D**);

void DecayMCProcess(Int_t experiment, Int_t productionmode, Int_t decaymode, Int_t numberOfMCevents, Bool_t acc = true, Bool_t verb = false){
	// put default values for optional arguments for running with interpreter
	allInAcceptance = acc;
	verbose = verb;
	// check argument validity when running with interpreter
	if(productionmode<0||productionmode>=prodmodeNames.size()){
		std::cout << "[Error] Invalid production mode: " <<  productionmode << ". Use values ";
		for(Int_t imode=0;imode<prodmodeNames.size();imode++) std::cout << imode << " (" << prodmodeNames[imode] << "), ";
		std::cout << std::endl;
		exit(1);
	}
	TString prodmodeName = prodmodeNames[productionmode];

	if(decaymode<0||decaymode>=decaymodeNames.size()){
		std::cout << "[Error] Invalid decay mode: " <<  decaymode << ". Use values ";
		for(Int_t imode=0;imode<decaymodeNames.size();imode++) std::cout << imode << " (" << decaymodeNames[imode] << "), ";
		std::cout << std::endl;
		exit(1);
	}
	TString decaymodeName = decaymodeNames[decaymode];
	ExpParameters genParam(experiment,productionmode,decaymode);

	std::cout << "[Info] Simulating " << genParam.GetExpName() << " as " << genParam.GetExpLabel() << " in " << prodmodeName << " production mode and " << decaymodeName << " decay mode with POT " << genParam.GetPOT() << " and normCrossSec " << genParam.GetNormCrossSec() << std::endl;

	//variables for 3-body decays
	TH1D** rescaleDalitz1d;
	TH2D** rescaleDalitz;
	TH2D*** originalDalitz;
	std::vector<Double_t> mALPList; //list of allowed mass bins
	std::vector<Double_t> valGamma,m12List,m23List;
	std::vector<std::array<std::vector<Double_t>,3>> ALPlist; // each mALP has unique m12 and m23 axis + m12*m23 values of diff. width  (contains m12List and m23List and valGamma)

	if (decaymode>2){ //width tables only for 3-body decays
		std::cout << "[Info] Simulating 3-body decay, will read differential width table" << std::endl;
		// reading diff. width tables for 3-body decays
		ifstream widthFile;
		widthFile.open(Form("%s/%s.dat",widthPath.Data(),decaymodeName.Data()));
		if (!widthFile.is_open()) {
			std::cout << "[Error] Differential width table not found" << std::endl;
			exit(1);
		}
		Double_t mALPCon, m12Con, m23Con, valCon; //containers
		widthFile >> mALPCon >> m12Con >> m23Con >> valCon;
		mALPList.push_back(mALPCon);
		m12List.push_back(m12Con);
		m23List.push_back(m23Con);
		valGamma.push_back(valCon);
		while (widthFile >> mALPCon >> m12Con >> m23Con >> valCon) {
			if(mALPCon == mALPList.back()){
				valGamma.push_back(valCon);
				if(m12Con == m12List.back()) m23List.push_back(m23Con);
				else{
					m23List.clear();
					m23List.push_back(m23Con);
					m12List.push_back(m12Con);
				}
			} else{
				//updating ALP mass
				mALPList.push_back(mALPCon);

				//// FOR DEBUGGING /////
				if(verbose)
					std::cout << "Printing number of m12,m23,value points for mass " << mALPCon << " : " << m12List.size() << " " << m23List.size() << " " << valGamma.size() << std::endl;
				///////////
				
				//updating list
				std::array<std::vector<Double_t>,3> mALPaxis = {m12List,m23List,valGamma};
				ALPlist.push_back(mALPaxis);
				valGamma.clear();
				valGamma.push_back(valCon);
				m12List.clear();
				m12List.push_back(m12Con);
				m23List.clear();
				m23List.push_back(m23Con);
			}
		}
		if (widthFile.eof()){
			std::array<std::vector<Double_t>,3> mALPaxis = {m12List,m23List,valGamma};
			ALPlist.push_back(mALPaxis);
			valGamma.clear();
			m12List.clear();
			m23List.clear();
			std::cout << "[Info] Finished reading differential width table" << std::endl << std::endl;
		} else{
			std::cout << "[Error] Did not reach EoF when reading differential width table" << std::endl;
			exit(1);
		}

		widthFile.close();

		std::cout << "[Info] Plots for axion widths declaration and filling:" << std::endl;

		rescaleDalitz = new TH2D*[ALPlist.size()];
		rescaleDalitz1d = new TH1D*[ALPlist.size()];
		originalDalitz = new TH2D**[ALPlist.size()]; // dalitz plot before and after acceptance criteria
		//// FOR DEBUGGING /////
		if(verbose){
			std::cout << "m12 begin: " << ALPlist[0][0].front() << " end: " << ALPlist[0][0].back() << std::endl;
			std::cout << "m23 begin: " << ALPlist[0][1].front() << " end: " << ALPlist[0][1].back() << std::endl;
		}
		////////////////
		for (Int_t mb=0; mb < ALPlist.size(); mb++){
			Double_t wdM12 = (ALPlist[mb][0].back()-ALPlist[mb][0].front())/ALPlist[mb][0].size();
			Double_t wdM23 = (ALPlist[mb][1].back()-ALPlist[mb][1].front())/ALPlist[mb][1].size();
			rescaleDalitz[mb] =
				new TH2D(Form("DalitzDens%f",mALPList[mb]),Form("rescaleDalitz for ALP mass %f GeV;m12^2;m23^2;weight",TMath::Power(10.,mALPList[mb])),
				ALPlist[mb][0].size(),ALPlist[mb][0].front()-0.5*wdM12,ALPlist[mb][0].back()+0.5*wdM12,
				ALPlist[mb][1].size(),ALPlist[mb][1].front()-0.5*wdM23,ALPlist[mb][1].back()+0.5*wdM23);
			originalDalitz[mb] = new TH2D*[4];
			for (int ib = 0; ib<4; ib++){
			  originalDalitz[mb][ib] =
			    new TH2D(Form("OriginalDalitz%f_befAft%d",mALPList[mb],ib),Form("originalDalitz before/after %d for ALP mass %f GeV;m12^2;m23^2",ib,TMath::Power(10.,mALPList[mb])),
				     ALPlist[mb][0].size(),ALPlist[mb][0].front()-0.5*wdM12,ALPlist[mb][0].back()+0.5*wdM12,
				     ALPlist[mb][1].size(),ALPlist[mb][1].front()-0.5*wdM23,ALPlist[mb][1].back()+0.5*wdM23);
			}
		}
		std::cout << "[Info] Plots declared." << std::endl;

		for (Int_t mb=0; mb < ALPlist.size(); mb++){
			Int_t linecounter=0;
			for (Int_t m12b=0; m12b < ALPlist[0][0].size(); m12b++){
				for (Int_t m23b=0; m23b < ALPlist[0][1].size(); m23b++){
					rescaleDalitz[mb]->SetBinContent(m12b+1, m23b+1, ALPlist[mb][2][linecounter]);
					linecounter+=1;
				}
			}
			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << "linecounter reached: " << linecounter << std::endl;
			/////////
			rescaleDalitz1d[mb] = rescaleDalitz[mb]->ProjectionX(Form("%s_m12",rescaleDalitz[mb]->GetName()));
		}

		std::cout << "[Info] Plots filled." << std::endl;

		genParam.SetMinMassA(TMath::Power(10.,mALPList.front()));
		genParam.SetMaxMassA(TMath::Power(10.,mALPList.back()));

	} else if (decaymode == 1){
		genParam.SetMinMassA(2*MEl);
		genParam.SetMaxMassA(genParam.GetMaxMassAInFile()[1]); //3 GeV
	} else if (decaymode == 2){
		genParam.SetMinMassA(2*MMu);
		genParam.SetMaxMassA(genParam.GetMaxMassAInFile()[1]); //3 GeV
	} else{
		genParam.SetMinMassA(0.);
		genParam.SetMaxMassA(genParam.GetMaxMassAInFile()[1]);
	}
	std::cout << "[Info] For " << decaymodeName << " decay mode min ALP mass = " << genParam.GetMinMassA() << " [GeV] and max ALP mass = " << genParam.GetMaxMassA() << " [GeV] NmassFiles = " << genParam.GetNMassFiles() << std::endl << std::endl;

	//plots declaration
	std::cout << "[Info] Plots for axion distributions declaration:" << std::endl;

	std::vector<std::vector<TH1D*>> beforeLog;
	std::vector<std::vector<std::vector<TH1D*>>> afterLogMC;
	std::vector<std::vector<TH2D*>> ExpectedNum;
	for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){
		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << " declaration for mass file: #" << ifile << " MALP in range " << genParam.GetMinMassFile()[ifile] << "-" << genParam.GetMaxMassFile()[ifile] << " Npoints = " << genParam.GetNMassX()[ifile] << std::endl;
		////////////////
		std::vector<TH1D*> beforeLogFile;
		std::vector<std::vector<TH1D*>> afterLogMCFile;
		for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
			beforeLogFile.push_back(new TH1D(Form("beforeLog[%d]%d",ifile,mb),"beforeLog;weight in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
			std::vector<TH1D*> afterLogMCbin;
			if(genParam.GetNRegions() < 1)
				afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d",ifile,mb),"afterLogMC;weight in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
			else
				for (Int_t iReg=0; iReg < genParam.GetNRegions(); iReg++)
					afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d-reg%d",ifile,mb,iReg+1),"afterLogMC;weight in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
			afterLogMCFile.push_back(afterLogMCbin);
		} //massBins
		beforeLog.push_back(beforeLogFile);
		afterLogMC.push_back(afterLogMCFile);

		std::vector<TH2D*> ExpectedNumFile;
		if(genParam.GetNRegions() < 1)
			ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]",ifile),"Expected Num",genParam.GetNMassX()[ifile],genParam.GetMinMassX()[ifile]-0.5*genParam.GetWdMassX()[ifile],genParam.GetMaxMassX()[ifile]+0.5*genParam.GetWdMassX()[ifile],genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
		else
			for (Int_t iReg=0; iReg < genParam.GetNRegions(); iReg++)
				ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]-reg%d",ifile,iReg+1),"Expected Num",genParam.GetNMassX()[ifile],genParam.GetMinMassX()[ifile]-0.5*genParam.GetWdMassX()[ifile],genParam.GetMaxMassX()[ifile]+0.5*genParam.GetWdMassX()[ifile],genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
		ExpectedNum.push_back(ExpectedNumFile);
	} //massFiles

	std::cout << "[Info] Plots declared." << std::endl << std::endl << "[Info] Reading axion distribution now:" << std::endl;
	// read in axion spectrum
	TH2D*** fHAxionETheta = new TH2D**[genParam.GetNMassFiles()];
	for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){
		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << std::endl << " reading mass file: #" << ifile << " MALP in range " << genParam.GetMinMassFile()[ifile] << "-" << genParam.GetMaxMassFile()[ifile] << " Npoints = " << genParam.GetNMassX()[ifile] << std::endl;
		///////////////

		fHAxionETheta[ifile] = new TH2D*[genParam.GetNMassX()[ifile]];

		ifstream infile;
		string infilePath = Form("%s/%s/alp_%s_beam%dGeV_%sto%sMeV_%s.dat",sourcePath.Data(),genParam.GetExpName().Data(),prodmodeName.Data(),genParam.GetBeamEnergy(),genParam.GetMinMassFile()[ifile].Data(),genParam.GetMaxMassFile()[ifile].Data(),genParam.GetExpName().Data());
		infile.open(infilePath);
		if (!infile.is_open()) {
			std::cout << "[Error] Axion momentum distribution table not found" << std::endl;
			std::cout << "[Error] with path: " << infilePath << std::endl;
			exit(1);
		}

		Double_t Theta, Energy, Mass, Val;
		std::vector<Double_t> valV;
		while (infile >> Theta >> Energy >> Mass >> Val) {
			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << "\r" << "  ReadAxionETheta >> reading from file to memory. Line: " << valV.size() << std::flush;
			///////////////////////
			if (Energy<0 || Theta <0 || Mass <0 ){
				std::cout << "[Error] Negative value encountered in dat file " << std::endl;
				exit(1);
			}
			valV.push_back(Val);
		}
		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << std::endl << "  ReadAxionETheta >> reading from file to memory completed" << std::endl;
		/////////////////////
		infile.close();

		for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){ //loop over output mass bins //mapping linear input on log output
			Double_t logMassX = genParam.GetMinMassX()[ifile] + genParam.GetWdMassX()[ifile]*mb;
			Double_t massX = TMath::Power(10.,logMassX);

			//reading axion files
			Int_t mbIn = TMath::Floor((massX-genParam.GetMinMassAInFile()[ifile])/genParam.GetWdMassA()[ifile] + 0.5);//corresponding input bin for output bin mb //avoid issues with C++ not correctly rounding
			if (mbIn < 0) {
				std::cout << "[Error] Mass specified for Axion is out of range. Input mass bin is " << mbIn << std::endl;
				std::cout << "[Error] problem occured with masswid= " << genParam.GetWdMassA()[ifile] << " mAMinInFile " << genParam.GetMinMassAInFile()[ifile] << " and mass: " << massX << std::endl;
				exit(1);
			}
			if (mbIn >= genParam.GetValuesMassAInFile()[ifile]){
				std::cout << "[Error] Mass specified for Axion is out of range. Input mass bin is " << mbIn << std::endl;
				std::cout << "[Error] problem occured with masswid= " << genParam.GetWdMassA()[ifile] << " mAMinInFile " << genParam.GetMinMassAInFile()[ifile] << " and mass: " << massX << std::endl;
				exit(1);
			}

			fHAxionETheta[ifile][mb] = new TH2D(Form("AxionETheta[%d]%d",ifile,mb), Form("AxionETheta for ALP mass %f GeV (bin:%d)",massX,mb),
					genParam.GetValuesThetaInFile(), genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta(),genParam.GetMaxThetaInFile()+0.5*genParam.GetWdTheta(),
					genParam.GetValuesEnergyInFile(), genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy(),genParam.GetMaxEnergyInFile()+0.5*genParam.GetWdEnergy()); // Assure that Bin spacing does not overshoot to negative bin edge!
			if((genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy())<0){std::cout << "WARNING " << (genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy()) << std::endl;}
			if((genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta())<0){std::cout << "WARNING " << (genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta()) << std::endl;}

			Int_t offset = genParam.GetValuesThetaInFile()*genParam.GetValuesEnergyInFile();

			Double_t massVal = genParam.GetMinMassAInFile()[ifile] + genParam.GetWdMassA()[ifile]*mbIn;
			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << "\r" << "  ReadAxionETheta >> reading from line " << offset*mbIn + 1 << " to " << offset*(mbIn+1) << " for mb: " << mb << " (mass " << massX << " GeV) selected massBin: " << mbIn << " (original mass " << massVal << " GeV)" << "          " << std::flush;
			////////////
			for (Int_t lineCounter = offset*mbIn+1; lineCounter < offset*(mbIn+1); lineCounter++){
				if(lineCounter>valV.size()){
					std::cout << "[Error] Line " << lineCounter << " is out of buffer size." << std::endl;
					exit(1);
				}
				Int_t iEnergy = (lineCounter - offset*mbIn)/genParam.GetValuesThetaInFile();
				Int_t iTheta = lineCounter - offset*mbIn - iEnergy*genParam.GetValuesThetaInFile();
				//// FOR DEBUGGING /////
				if(verbose)
					cout << "Line: " << lineCounter << " corresponds to theta bin: " << iTheta << " | energy bin: " << iEnergy << " | mass bin: " << mb << " | and value: " << valV[lineCounter] << std::endl;
				//////
				fHAxionETheta[ifile][mb]->SetBinContent(iTheta+1, iEnergy+1, valV[lineCounter-1]);
			}
		}
		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << std::endl << "  ReadAxionETheta >> reading completed " << std::endl;
		//////////////
		valV.clear();
	}

	std::cout << std::endl << "[Info] Finished reading axion momentum distribution" << std::endl;

	Double_t massXprevisous = 0;
	for (Int_t ifile=1; ifile < genParam.GetNMassFiles(); ifile++){

		std::cout << std::endl << "[Info] Entering MC for mass range " << genParam.GetMinMassAInFile()[ifile] << " - " << genParam.GetMaxMassAInFile()[ifile] << std::endl;

		TRandom3 *rndmGen = new TRandom3(0);
		rndmGen->SetSeed(0);

		AxionParameters alpParam;

		//mass loop
		Int_t mbTrue = -1;
		Int_t mbTot = genParam.GetNMassX()[ifile];
		for (Int_t mb=0; mb < mbTot; mb++){
			//progress bar:
			std::cout << "\r" << " Simulating evts for mass bin: " << mb << "    " << std::flush;

			Double_t logMassX = genParam.GetMinMassX()[ifile] + genParam.GetWdMassX()[ifile]*mb;
			alpParam.massA = TMath::Power(10.,logMassX);
			if (alpParam.massA < genParam.GetMinMassA() || alpParam.massA > genParam.GetMaxMassA()) continue; //ALP mass constraint
			mbTrue+=1;

			for(Int_t iEv=0; iEv< numberOfMCevents; iEv++){ //MC loop

				//Pick random from input distributions
				fHAxionETheta[ifile][mb]->GetRandom2(alpParam.thetaA,alpParam.energyA);
				alpParam.phiA = rndmGen->Rndm()*2*TMath::Pi(); //random uniform

				Double_t pA2 = alpParam.energyA*alpParam.energyA - alpParam.massA*alpParam.massA;
				if(pA2<0) pA2 = 0.; //can happen since we pick randomly from the energy table
				alpParam.pA = TMath::Sqrt(pA2);
				alpParam.betaA = alpParam.pA/alpParam.energyA;

				int iEnergy = (alpParam.energyA-genParam.GetMinEnergyInFile())/genParam.GetWdEnergy(); //bin
				int iTheta  = (alpParam.thetaA-genParam.GetMinThetaInFile())/genParam.GetWdTheta(); //bin
				alpParam.crossSecA = fHAxionETheta[ifile][mb]->GetBinContent(iTheta+1,iEnergy+1);

				if (genParam.GetPhiBeamEuler() != 0 || genParam.GetThetaBeamEuler() != 0 || genParam.GetPsiBeamEuler() != 0){
					TVector3 oldAlp(TMath::Sin(alpParam.thetaA)*TMath::Cos(alpParam.phiA),
							TMath::Sin(alpParam.thetaA)*TMath::Sin(alpParam.phiA),
							TMath::Cos(alpParam.thetaA));
					TRotation a;
					//// FOR DEBUGGING /////
					if(verbose)
					  	cout << "Original X:" << oldAlp.X() << " | Y:" << oldAlp.Y() << " | Z:" << oldAlp.Z() << std::endl;
					///////////////////////
					a.SetXEulerAngles(genParam.GetPhiBeamEuler(),genParam.GetThetaBeamEuler(),genParam.GetPsiBeamEuler());
					TVector3 newAlp = a*oldAlp; // active rotation
					//// FOR DEBUGGING /////
					if(verbose)
					  	cout << "New X:" << newAlp.X() << " | Y:" << newAlp.Y() << " | Z:" << newAlp.Z() << std::endl;
					//////////////////////
					
					alpParam.thetaA = newAlp.Theta();
					alpParam.phiA = newAlp.Phi();
				}

				//increase efficiency by skipping ALPs out of calorimeter acceptance
				Bool_t ALPInAcceptance = true;
				if(allInAcceptance){
					Double_t dZ = genParam.GetZLKR() -genParam.GetZ0();
					TVector3 ALPPosAtLKr(genParam.GetX0() + dZ*TMath::Sin(alpParam.thetaA)*TMath::Cos(alpParam.phiA),
										 genParam.GetY0() + dZ*TMath::Sin(alpParam.thetaA)*TMath::Sin(alpParam.phiA),
										 genParam.GetZLKR());
					ALPInAcceptance = genParam.inCaloOuterEdgeAcceptance(ALPPosAtLKr);
				}

				//decay width loop
				for (Int_t nW=0; nW < genParam.GetNWidths(); nW++){
					alpParam.widthExpA = genParam.GetMinWidth() + genParam.GetWdWidth()*nW;  //y-axis
					if(allInAcceptance && !ALPInAcceptance){
						beforeLog[ifile][mb]->Fill(alpParam.widthExpA,1);
						continue;
					}
					alpParam.decayLengthA = alpParam.pA/alpParam.massA*hc* 1./TMath::Power(10.,alpParam.widthExpA); // [m] beta*gamma/(Gamma/hc)
					if (decaymode <= 2) decayTo2Body(genParam,alpParam,*beforeLog[ifile][mb],afterLogMC[ifile][mb],rndmGen);
					else decayTo3Body(genParam,alpParam,*beforeLog[ifile][mb],afterLogMC[ifile][mb],rescaleDalitz[mbTrue],rndmGen,originalDalitz[mbTrue]);
				} // width loop

			} // ievt

		} // mass loop

		std::cout << std::endl << "[Info] Event loop finished! Will now plot efficiencies" << std::endl;

		std::vector<std::vector<TH1D*>> eff1Log;

		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << "new event " << std::endl;
		//////////////
		for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
			std::vector<TH1D*> eff1LogBin;
			TH1D* deno = (TH1D*) beforeLog[ifile][mb]->Clone();
			for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
				TH1D* nume = (TH1D*) afterLogMC[ifile][mb].at(iReg)->Clone();
				eff1LogBin.push_back((TH1D*) nume->Clone(Form("eff1Log_bin_%d_region_%d",mb,iReg+1)));
				for(Int_t jb=0; jb<eff1LogBin[iReg]->GetNbinsX(); jb++)	if(deno->GetBinContent(jb)>0) eff1LogBin[iReg]->SetBinContent(jb,nume->GetBinContent(jb)/deno->GetBinContent(jb));
				delete nume;
			}
			delete deno;
			eff1Log.push_back(eff1LogBin);
		} // mass


		for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
			Double_t logMassX = genParam.GetMinMassX()[ifile] + genParam.GetWdMassX()[ifile]*mb;
			//// FOR DEBUGGING /////
				if(verbose)
				cout << "Output mass:" << TMath::Power(10,logMassX) << " bin:" << mb << std::endl;
			////////////////////////
			for (Int_t i=0; i< beforeLog[ifile][mb]->GetNbinsX(); i++) {
				Double_t widthExp = beforeLog[ifile][mb]->GetBinCenter(i+1);
				for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++) ExpectedNum[ifile].at(iReg)->Fill(logMassX,widthExp, eff1Log[mb][iReg]->GetBinContent(i+1));
			} // massbins

		} // mass

		std::cout << "[Info] Writing root output " << std::endl;

		TString filename;
		filename = Form("%s/%s/ALPsfromPi0_experiment%dprodmode%d_MCevents%d_%sto%sMeV_decaymode%s.root",outPath.Data(),genParam.GetExpName().Data(),genParam.GetExpNum(),productionmode,numberOfMCevents,genParam.GetMinMassFile()[ifile].Data(),genParam.GetMaxMassFile()[ifile].Data(),decaymodeName.Data());
		TFile *fileOut = new TFile(filename.Data(),"RECREATE");
		fileOut->cd();
		for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++)
			ExpectedNum[ifile].at(iReg)->Write();

		for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
			beforeLog[ifile][mb]->Write();
			delete beforeLog[ifile][mb];
			
			fHAxionETheta[ifile][mb]->Write();
			delete fHAxionETheta[ifile][mb];
			for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
				afterLogMC[ifile][mb].at(iReg)->Write();
				delete afterLogMC[ifile][mb].at(iReg);

				eff1Log[mb][iReg]->Write();
				delete eff1Log[mb][iReg]; 
			}
		}

		fileOut->Close();
		delete fileOut;
		if(rndmGen) delete rndmGen;


		//// FOR DEBUGGING - Dalitz plot rescale check/////
		if(verbose){
			TString testFilename = Form("%s/%s/DalitzStudy_experiment%dprodmode%d_MCevents%d_%sto%sMeV_decaymode%s.root",outPath.Data(),genParam.GetExpName().Data(),genParam.GetExpNum(),productionmode,numberOfMCevents,genParam.GetMinMassFile()[ifile].Data(),genParam.GetMaxMassFile()[ifile].Data(),decaymodeName.Data());//_%s.root",outPath.Data(),genParam.expName.Data(),decaymodeName.Data());
			TFile *testFileOut = new TFile(testFilename.Data(),"RECREATE");
			testFileOut->cd();
			TGraph* acceptance3BodyFlat = new TGraph();
			acceptance3BodyFlat->SetName("acceptance3BodyFlat");
			TGraph* acceptance3BodyDalitzDensity = new TGraph();
			acceptance3BodyDalitzDensity->SetName("acceptance3BodyDalitzDensity");
			for (Int_t mb=0; mb < ALPlist.size(); mb++) {
			rescaleDalitz[mb]->Write();

			for (int ib=0; ib<4; ib++){
				originalDalitz[mb][ib]->Write();
			}
			if (originalDalitz[mb][0]->GetEntries()) acceptance3BodyFlat->SetPoint(acceptance3BodyFlat->GetN(),TMath::Power(10.,mALPList[mb]),originalDalitz[mb][1]->GetSum()/originalDalitz[mb][0]->GetSum());
			if (originalDalitz[mb][2]->GetEntries()) acceptance3BodyDalitzDensity->SetPoint(acceptance3BodyDalitzDensity->GetN(),TMath::Power(10.,mALPList[mb]),originalDalitz[mb][3]->GetSum()/originalDalitz[mb][2]->GetSum());
			}
			acceptance3BodyFlat->Write();
			acceptance3BodyDalitzDensity->SetMarkerStyle(24);
			acceptance3BodyDalitzDensity->SetMarkerColor(2);
			acceptance3BodyDalitzDensity->SetLineColor(2);
			acceptance3BodyDalitzDensity->Write();

			for (Int_t mb=0; mb < ALPlist.size(); mb++) {
			rescaleDalitz1d[mb]->Write();
			}

			testFileOut->Close();
			delete testFileOut;
		}
		////////////////



		std::cout << "[Info] Data written" << std::endl;

	} //massFiles loop

	std::cout << "[Info] Writing dat files " << std::endl;

	Double_t x,y,z;
	ofstream RunOutput; // Export data points
	TString outFileName = Form("%s/%s/%s_%s_%s.dat",outPath.Data(),genParam.GetExpName().Data(),genParam.GetExpLabel().Data(),prodmodeName.Data(),decaymodeName.Data());

	for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
		if (genParam.GetNRegions() > 1) outFileName = Form("%s/%s/%s_%s_%s_reg%d.dat",outPath.Data(),genParam.GetExpName().Data(),genParam.GetExpLabel().Data(),prodmodeName.Data(),decaymodeName.Data(),iReg+1);

		/// open and close to clear the file
		RunOutput.open(outFileName.Data());
		RunOutput.close();

		// open again
		RunOutput.open(outFileName.Data(), ios::app);
		// extract correct limit plot

		// extract exactly "paperstyle"", to have compareable ranges easily
		Double_t xPrevious = 0;
		for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){
			for( int i=1; i<=ExpectedNum[ifile].at(iReg)->GetNbinsX(); i++){ //0 is underflow bin
				x=TMath::Power(10.,ExpectedNum[ifile].at(iReg)->GetXaxis()->GetBinCenter(i)); //GeV!!
				if(TMath::Abs(x-xPrevious)<1E-6) continue; //skipping repeating mass bins

				for( int j=1; j<=ExpectedNum[ifile].at(iReg)->GetNbinsY(); j++){
					y=TMath::Power(10.,ExpectedNum[ifile].at(iReg)->GetYaxis()->GetBinCenter(j));
					z=ExpectedNum[ifile].at(iReg)->GetBinContent(i,j);
					RunOutput << x << " " << y << " " << z << "\n"; //
				} // loop over y
				xPrevious = x;

			} // loop over x
			delete ExpectedNum[ifile].at(iReg);
		} //massFiles loop
		RunOutput.close();
	}


	exit(0);

}  // main

/// fct defs
void decayTo2Body(ExpParameters &gen, AxionParameters &alp, TH1D &before, std::vector<TH1D*> &after, TRandom3 *rndm){ // a -> mumu decay
	Double_t randoms[5];
	for (Int_t i=0; i<5; i++) randoms[i] = rndm->Rndm();

	/////// from gamma ->
	// sampling of the a decay position
	Double_t w1 = (1-TMath::Exp(-gen.GetBeamDecayLength()/alp.decayLengthA)); // decay

	//// FOR DEBUGGING /////
	if(verbose)
		std::cout << " BeamDecayLength " << gen.GetBeamDecayLength() << " lamb " << alp.decayLengthA << " ZDist " << gen.GetZDist() << "w1" << w1 << std::endl;
	//////

	TVector3 FVentrance( gen.GetX0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA),gen.GetY0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA),gen.GetZ0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Cos(alp.thetaA));
	TVector3 entrypoint(gen.GetX0(),gen.GetY0(),gen.GetZ0());
	Double_t tau = (FVentrance-entrypoint).Mag()/alp.decayLengthA; // basically always (105-23)/decaylength
	Double_t decayweight = TMath::Exp(-tau); // reach

	// here comes a test to check with NA62MC //
	Double_t rProduction_a = -alp.decayLengthA*TMath::Log(1.-w1*randoms[0]); // now OK, see documentation for NA62 analysis

	// rProduction_a is within BeamDecayLength, i.e. in between ZFVIn and the ZFVEnd, so we need to add the distance travelled //
	Double_t xProduction_a = gen.GetX0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA); // exponential distribution
	Double_t yProduction_a = gen.GetY0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA); // exponential distribution
	Double_t zProduction_a = gen.GetZ0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Cos(alp.thetaA); // exponential distribution?


	// for the flat decay we have to remember which random was chosen and reqeight in the end
	TVector3 decaypoint(xProduction_a,yProduction_a,zProduction_a);
	Double_t deltar= (decaypoint-FVentrance).Mag();
	///// <- from gamma

	Double_t pStar_a = TMath::Sqrt(alp.massA*alp.massA-4.*gen.GetMassFinState()[0]*gen.GetMassFinState()[0])/2;
	TVector3 bb2(alp.betaA*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA), alp.betaA*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA), alp.betaA*TMath::Cos(alp.thetaA));
	Double_t cTheta = -1. + 2.*randoms[3]; // flat distribution on costheta
	Double_t thetaProduction_a = TMath::ACos(cTheta);  // theta value in (0,pi);
	Double_t phiProduction_a = 2.*TMath::Pi()*randoms[4]; // flat distribution

	TLorentzVector pfA[2];
	pfA[0].SetPxPyPzE(pStar_a*TMath::Sin(thetaProduction_a)*TMath::Cos(phiProduction_a),
			pStar_a*TMath::Sin(thetaProduction_a)*TMath::Sin(phiProduction_a),
			pStar_a*TMath::Cos(thetaProduction_a),
			TMath::Sqrt(pStar_a*pStar_a + gen.GetMassFinState()[0]*gen.GetMassFinState()[0]));
	pfA[1].SetPxPyPzE(-pStar_a*TMath::Sin(thetaProduction_a)*TMath::Cos(phiProduction_a),
			-pStar_a*TMath::Sin(thetaProduction_a)*TMath::Sin(phiProduction_a),
			-pStar_a*TMath::Cos(thetaProduction_a),
			TMath::Sqrt(pStar_a*pStar_a + gen.GetMassFinState()[1]*gen.GetMassFinState()[1]));
	TVector3 posFVEnd[2];
	TVector3 posStraw1[2];
	TVector3 posMagnet[2];
	TVector3 posStraw4[2];
	TVector3 posLKr[2];
	TVector3 posMUV3[2];
	Int_t iOnMUV3 = 0;
	Int_t iChOnLKr = 0;
	Int_t iOnLKr = 0;
	Int_t iOnStraw1 = 0;
	Int_t iOnStraw4 = 0;
	Double_t totalEnergyAcceptance = 0.;

	Double_t thetaC[2];
	Double_t strawRadius1[2];
	Double_t strawRadius4[2];

	for (Int_t k=0; k<2; k++) {
		pfA[k].Boost(bb2);
		Double_t pTf = pow(pow(pfA[k].X(),2)+pow(pfA[k].Y(),2),0.5);
		thetaC[k] = pTf/pfA[k].Z();
		//// FOR DEBUGGING /////
		if(verbose)
			std::cout << "thetaC " << thetaC[k] << "pfA[k].P() " << pfA[k].P() << std::endl;
		///////////////

		posFVEnd[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZFVEnd()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZFVEnd()-zProduction_a),gen.GetZFVEnd());
		if(gen.GetChargeFinState()[k] !=0 ){
			posStraw1[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZStraw1()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZStraw1()-zProduction_a),gen.GetZStraw1());
			if (gen.GetMagKick() != 0){
				posMagnet[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZMagnet()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZMagnet()-zProduction_a),gen.GetZMagnet());
				posStraw4[k].SetXYZ(posMagnet[k].X() + (pfA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pfA[k].Z()*(gen.GetZStraw4()-gen.GetZMagnet()),posMagnet[k].Y() + pfA[k].Y()/pfA[k].Z()*(gen.GetZStraw4()-gen.GetZMagnet()),gen.GetZStraw4());
				posLKr[k].SetXYZ(posMagnet[k].X() + (pfA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pfA[k].Z()*(gen.GetZLKR()-gen.GetZMagnet()),posMagnet[k].Y() + pfA[k].Y()/pfA[k].Z()*(gen.GetZLKR()-gen.GetZMagnet()),gen.GetZLKR());
				posMUV3[k].SetXYZ(posMagnet[k].X() + (pfA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pfA[k].Z()*(gen.GetZMUV3()-gen.GetZMagnet()),posMagnet[k].Y() + pfA[k].Y()/pfA[k].Z()*(gen.GetZMUV3()-gen.GetZMagnet()),gen.GetZMUV3());
			} else{
				posStraw4[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZStraw4()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZStraw4()-zProduction_a),gen.GetZStraw4());
				posLKr[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZLKR()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZLKR()-zProduction_a),gen.GetZLKR());
				posMUV3[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZMUV3()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZMUV3()-zProduction_a),gen.GetZMUV3());
			}

			if (gen.inSpectrometerAcceptance(posStraw1[k])) iOnStraw1++;
			if (gen.inSpectrometerAcceptance(posStraw4[k])) iOnStraw4++;
			if (gen.inCaloAcceptance(posFVEnd[k], posLKr[k], zProduction_a)){
				iChOnLKr++;
				totalEnergyAcceptance += pfA[k].E();
			}
			if (gen.inMuonVetoAcceptance(posFVEnd[k], posMUV3[k],zProduction_a)) iOnMUV3++;
		} else{
			posLKr[k].SetXYZ(xProduction_a + pfA[k].X()/pfA[k].Z()*(gen.GetZLKR()-zProduction_a),yProduction_a + pfA[k].Y()/pfA[k].Z()*(gen.GetZLKR()-zProduction_a),gen.GetZLKR());
			if (gen.inCaloAcceptance(posFVEnd[k], posLKr[k], zProduction_a)){
				iOnLKr++;
				totalEnergyAcceptance += pfA[k].E();
			}
		}
	} // both muons

	Int_t condition = 0;
	Double_t sigacc = gen.GetSigacceptance();
	if (gen.GetMassFinState()[0] + gen.GetMassFinState()[1] < MEl){ //2gamma
		condition = gen.twoPhotonCondition(iOnLKr, totalEnergyAcceptance, zProduction_a, posFVEnd, posLKr, pfA);
	} else if (gen.GetMassFinState()[0] + gen.GetMassFinState()[1] == 2*MEl){ //2el, same conditions as for photons?
		condition = gen.twoPhotonCondition(iChOnLKr, totalEnergyAcceptance, zProduction_a, posFVEnd, posLKr, pfA);
	} else if (gen.GetMassFinState()[0] + gen.GetMassFinState()[1] == 2*MMu){ //2mu
		sigacc = gen.GetSigacceptancemumu();
		condition = gen.twoMuonCondition(iOnStraw1, iOnStraw4, iChOnLKr, iOnMUV3, zProduction_a, pfA);
	} else if ((gen.GetMassFinState()[0] + gen.GetMassFinState()[1] == 2*MPiCh ) || (gen.GetMassFinState()[0] + gen.GetMassFinState()[1] == 2*MKCh )){ //2hadrons
		sigacc = gen.GetSigacceptancemumu();
		condition = gen.twoHadronCondition(iOnStraw1, iOnStraw4, iChOnLKr, iOnMUV3, zProduction_a, pfA);
	}
	if(condition > 0) after.at(condition-1)->Fill(alp.widthExpA, w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
	before.Fill(alp.widthExpA,1); //  decayweight in or out? This is firstweightlogsteps in NA62MC
} // decayTo2Body

void decayTo3Body(ExpParameters &gen, AxionParameters &alp, TH1D &before, std::vector<TH1D*> &after, TH2D* rescaleDal, TRandom3 *rndm, TH2D** originalDal){

  TGenPhaseSpace event;

	Double_t randoms[7];
	for (Int_t i=0; i<6; i++) randoms[i] = rndm->Rndm();

	// sampling of the a decay position
	Double_t w1 = (1-TMath::Exp(-gen.GetBeamDecayLength()/alp.decayLengthA)); // decay

	//// FOR DEBUGGING /////
	if(verbose)
		std::cout << " BeamDecayLength " << gen.GetBeamDecayLength() << " lamb " << alp.decayLengthA << " ZDist " << gen.GetZDist() << "w1" << w1 << std::endl;
	/////////////////

	TVector3 FVentrance( gen.GetX0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA),gen.GetY0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA),gen.GetZ0() + (gen.GetZFVIn() -gen.GetZ0())*TMath::Cos(alp.thetaA));
	TVector3 entrypoint(gen.GetX0(),gen.GetY0(),gen.GetZ0());
	Double_t tau = (FVentrance-entrypoint).Mag()/alp.decayLengthA;
	Double_t decayweight = TMath::Exp(-tau); // reach

	// here comes a test to check with NA62MC //
	Double_t rProduction_a = -alp.decayLengthA*TMath::Log(1.-w1*randoms[0]); // now OK, see documentation for NA62 analysis

	// rProduction_a is within BeamDecayLength, i.e. in between ZFVIn and the ZFVEnd, so we need to add the distance travelled //
	Double_t xProduction_a = gen.GetX0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA); // exponential distribution
	Double_t yProduction_a = gen.GetY0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA); // exponential distribution
	Double_t zProduction_a = gen.GetZ0() + (gen.GetZFVIn() -gen.GetZ0() + rProduction_a)*TMath::Cos(alp.thetaA); // exponential distribution?


	// for the flat decay we have to remember which random was chosen and reqeight in the end
	TVector3 decaypoint(xProduction_a,yProduction_a,zProduction_a);
	Double_t deltar= (decaypoint-FVentrance).Mag();

	// a -> 3body decay
	Double_t mds[3]={gen.GetMassFinState()[0],gen.GetMassFinState()[1],gen.GetMassFinState()[2]};
	TLorentzVector pgA[3]; // to include 3pi0->6gamma decays
	// a-> 1 2 3
	TVector3 bbA(alp.betaA*TMath::Sin(alp.thetaA)*TMath::Cos(alp.phiA), alp.betaA*TMath::Sin(alp.thetaA)*TMath::Sin(alp.phiA), alp.betaA*TMath::Cos(alp.thetaA));
	TLorentzVector palp(alp.energyA*bbA.X(),alp.energyA*bbA.Y(),alp.energyA*bbA.Z(),alp.energyA);
	event.SetDecay(palp, 3, mds);
	Double_t weight = event.Generate();
	for (int i=0; i<3; i++) pgA[i] = (*event.GetDecay(i));

	double m12 = mds[0]*mds[0]+mds[1]*mds[1]+2.*pgA[0].Dot(pgA[1]);
	double m23 = mds[1]*mds[1]+mds[2]*mds[2]+2.*pgA[1].Dot(pgA[2]);
	int m12bin = rescaleDal->GetXaxis()->FindBin(m12);
	int m23bin = rescaleDal->GetYaxis()->FindBin(m23);
	double densityWeight = rescaleDal->GetBinContent(m12bin,m23bin)*rescaleDal->GetNbinsX()*rescaleDal->GetNbinsY()/rescaleDal->GetSum();

	//increase efficiency by checking if all decay product projections are in calorimeter acceptance
	if(allInAcceptance){
		for (int i=0; i<3; i++){
			Double_t dZ = gen.GetZLKR() - decaypoint.Z();
			TVector3 PosAtLKr(decaypoint.X() + dZ*pgA[i].Px()/pgA[i].P(),
							  decaypoint.Y() + dZ*pgA[i].Py()/pgA[i].P(),
							  gen.GetZLKR());
			if(!gen.inCaloOuterEdgeAcceptance(PosAtLKr)){
				before.Fill(alp.widthExpA,1);
				//// FOR DEBUGGING /////
				if(verbose){
					originalDal[0]->Fill(m12,m23,w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*gen.GetSigacceptance());
					originalDal[2]->Fill(m12,m23,densityWeight*w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*gen.GetSigacceptance());
				}
				///////////////
				return;
			}
			
		}
	}


// assign charge to the particles

	int charges[3] = {gen.GetChargeFinState()[0],gen.GetChargeFinState()[1],gen.GetChargeFinState()[2]};
	int ncharged = 0;
	int chargedIndex[2] ={-999, -999};
	int charge =999;
	for (Int_t k=0; k<3; k++) {
	  if ((fabs(mds[k]-MKCh) < 0.001) || (fabs(mds[k]-MPiCh) < 0.001) || (fabs(mds[k]-MMu) < 0.001) || (fabs(mds[k]-MEl) < 0.0001)) {
	    if (ncharged == 0) charge = 1;
	    else charge = -1;
	    chargedIndex[ncharged] = k;
	    ncharged++;
	  }
	  else charge = 0;
	  charges[k] = charge;
	}

	if (ncharged%2 != 0) {
	  std::cerr << "[Error] Wrong mass setting " << mds[0] << " " << mds[1] << " " << mds[2] << std::endl;
	  return;
	}

// extrapolate the charged particles

//	TVector3 rStart[3];
	TVector3 posFVEnd[3];
	TVector3 posStraw1[3];
	TVector3 posMagnet[3];
	TVector3 posStraw4[3];
	TVector3 posLKr[3];
	TVector3 posMUV3[3];
	Int_t iOnStraw1 = 0;
	Int_t iOnStraw4 = 0;
	Int_t iOnMUV3 = 0;
	Int_t iChOnLKr = 0;
	double dz;
	Double_t totalEnergyAcceptance =0;


	for (int k=0;k<3;k++) {
		if (fabs(mds[k]-MPi0) < 0.001 || fabs(mds[k]-MEta) < 0.001 || fabs(mds[k]-MEtaPrim) < 0.001 || mds[k] < 0.0001) continue;// use gammas, instead
		if (charges[k]==0) continue;

		posFVEnd[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZFVEnd()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZFVEnd()-zProduction_a),gen.GetZFVEnd());
		posStraw1[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZStraw1()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZStraw1()-zProduction_a),gen.GetZStraw1());
		if (gen.GetMagKick() != 0){
			posMagnet[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZMagnet()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZMagnet()-zProduction_a),gen.GetZMagnet());
			posStraw4[k].SetXYZ(posMagnet[k].X() + (pgA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pgA[k].Z()*(gen.GetZStraw4()-gen.GetZMagnet()),posMagnet[k].Y() + pgA[k].Y()/pgA[k].Z()*(gen.GetZStraw4()-gen.GetZMagnet()),gen.GetZStraw4());
			posLKr[k].SetXYZ(posMagnet[k].X() + (pgA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pgA[k].Z()*(gen.GetZLKR()-gen.GetZMagnet()),posMagnet[k].Y() + pgA[k].Y()/pgA[k].Z()*(gen.GetZLKR()-gen.GetZMagnet()),gen.GetZLKR());
			posMUV3[k].SetXYZ(posMagnet[k].X() + (pgA[k].X() + gen.GetChargeFinState()[k]*gen.GetMagKick())/pgA[k].Z()*(gen.GetZMUV3()-gen.GetZMagnet()),posMagnet[k].Y() + pgA[k].Y()/pgA[k].Z()*(gen.GetZMUV3()-gen.GetZMagnet()),gen.GetZMUV3());
		} else{
			posStraw4[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZStraw4()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZStraw4()-zProduction_a),gen.GetZStraw4());
			posLKr[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZLKR()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZLKR()-zProduction_a),gen.GetZLKR());
			posMUV3[k].SetXYZ(xProduction_a + pgA[k].X()/pgA[k].Z()*(gen.GetZMUV3()-zProduction_a),yProduction_a + pgA[k].Y()/pgA[k].Z()*(gen.GetZMUV3()-zProduction_a),gen.GetZMUV3());
		}

		if (gen.inSpectrometerAcceptance(posStraw1[k])) iOnStraw1++;
		if (gen.inSpectrometerAcceptance(posStraw4[k])) iOnStraw4++;
		if (gen.inCaloAcceptance(posFVEnd[k], posLKr[k], zProduction_a)){
			iChOnLKr++;
			if (fabs(pgA[k].M()-MMu)>0.0001) totalEnergyAcceptance += pgA[k].E(); // consider pi, e, K+- to release all of their energy in the LKr <-- not actually true
		}
		if (gen.inMuonVetoAcceptance(posFVEnd[k], posMUV3[k],zProduction_a)) iOnMUV3++;
	}

// evaluate decays to gg and consider the original gammas, if any

	TLorentzVector pgamma[6];
	int ngamma = 0;
	for (Int_t k=0; k<3; k++) {
		if (fabs(mds[k]-MPi0) < 0.001 || fabs(mds[k]-MEta) < 0.001 || fabs(mds[k]-MEtaPrim) < 0.001){
			// sample angle of 1,2 photons wrt their total momentum in the lab
			Double_t cThetaProduction_gg = -1. + 2.*rndm->Rndm(); // flat distribution on costheta
			Double_t thetaProduction_gg = TMath::ACos(cThetaProduction_gg);  // theta value in (0,pi);
			Double_t phiProduction_gg = 2.*TMath::Pi()*rndm->Rndm(); // flat distribution
			Double_t pStar = mds[k]*0.5;
			Double_t beta = pgA[k].P()/pgA[k].E();
			TVector3 bbmeson = (beta/pgA[k].P())*pgA[k].Vect();
			pgamma[ngamma].SetXYZT(
					pStar*TMath::Sin(cThetaProduction_gg)*TMath::Cos(phiProduction_gg),
					pStar*TMath::Sin(cThetaProduction_gg)*TMath::Sin(phiProduction_gg),
					pStar*TMath::Cos(cThetaProduction_gg),
					pStar);
			pgamma[ngamma].Boost(bbmeson);
			ngamma++;

			pgamma[ngamma].SetXYZT(
					-pStar*TMath::Sin(cThetaProduction_gg)*TMath::Cos(phiProduction_gg),
					-pStar*TMath::Sin(cThetaProduction_gg)*TMath::Sin(phiProduction_gg),
					-pStar*TMath::Cos(cThetaProduction_gg),
					pStar);
			pgamma[ngamma].Boost(bbmeson);
			ngamma++;
		}
		else if (fabs(mds[k]) < 0.0001) { // photon
			pgamma[ngamma].SetXYZT(pgA[k].X(),pgA[k].Y(),pgA[k].Z(),pgA[k].T());
			ngamma++;
		}
	}

//  extrapolate gammas to LKr

	TVector3 posgamma[6];
	TVector3 posGammaFVEnd[6];
	for (Int_t k=0; k<ngamma; k++){
		dz = (gen.GetZLKR()-zProduction_a);
		posgamma[k].SetXYZ(xProduction_a + pgamma[k].X()/pgamma[k].Z()*dz,yProduction_a + pgamma[k].Y()/pgamma[k].Z()*dz,gen.GetZLKR());
		dz = (gen.GetZFVEnd()-zProduction_a);
		posGammaFVEnd[k].SetXYZ(xProduction_a + pgamma[k].X()/pgamma[k].Z()*dz,yProduction_a + pgamma[k].Y()/pgamma[k].Z()*dz,gen.GetZFVEnd());
		if (gen.inCaloAcceptance(posGammaFVEnd[k],posgamma[k],zProduction_a)) totalEnergyAcceptance += pgamma[k].E();
	}

	//// FOR DEBUGGING /////
	// printout to test conditions
	if(verbose){
  		if (ncharged > 0) std::cout << "iOnStraw1 " << iOnStraw1 << " iOnStraw4 " << iOnStraw4 << std::endl;
		if (ncharged == 2) std::cout << "pgA[chargedIndex[0]].E() " << pgA[chargedIndex[0]].E() << " pgA[chargedIndex[1]].E() " << pgA[chargedIndex[1]].E() << std::endl;
		std::cout << "ngamma " << ngamma << " totalEnergyAcceptance  " << totalEnergyAcceptance << std::endl;
	}
	//// DEBUGGING /////


// experimental conditions
	Bool_t condition = 1;
	double sigacc = gen.GetSigacceptance();

	if(ncharged == 2){
		sigacc = gen.GetSigacceptancemumu();
		condition = gen.twoHadronCondition(iOnStraw1, iOnStraw4, iChOnLKr, iOnMUV3, zProduction_a, pgA);
	}

	if(ngamma) condition = condition && gen.multiplePhotonCondition(ngamma, totalEnergyAcceptance, zProduction_a, posGammaFVEnd, posgamma, pgamma);

	if(condition) {
		after.at(condition-1)->Fill(alp.widthExpA, densityWeight*w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
		originalDal[1]->Fill(m12,m23,w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
		originalDal[3]->Fill(m12,m23,densityWeight*w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
	}
	before.Fill(alp.widthExpA,1); //  decayweight in or out? This is firstweightlogsteps in NA62MC
	originalDal[0]->Fill(m12,m23,w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
	originalDal[2]->Fill(m12,m23,densityWeight*w1*decayweight*gen.GetPOT()*alp.crossSecA*gen.GetNormCrossSec()*sigacc);
} // close decayTo3Body
