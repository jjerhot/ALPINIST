 /// \file DecayMCProcess.C
 /// \Brief
 /// main part of the decay module running simulation
 /// \EndBrief
 /// \Detailed
 /// Loads input tables, simulates decay, checks if event passes selection, writes output.
 /// Can be run separately as a ROOT macro. Example:
 /// \code
 /// .include ../include
 /// .L ExpParameters.C
 /// .x DecayMCProcess.C(0,1,0,0,1000000)
 /// \endcode
 /// \EndDetailed

#include <iostream>
#include <array>
#include <fstream>
#include "TError.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h" //for TTree
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include <TROOT.h>
#include "TGenPhaseSpace.h"

#include "ExpParameters.h"

/// Parameter declaration
Bool_t allInAcceptance;
Bool_t verbose;
Bool_t flatDecay;
Bool_t flatDalitz;

/// Global class dec and init
ExoticParticleProperty exoProp("alp");
Particle *exotic;
TRandom3 rndmGen(0);

 /// \fn inDalitzBoundary
 /// \Brief
 /// Checks if event is within boundary of the Dalitz plot
 /// \EndBrief
Bool_t inDalitzBoundary(Double_t m122, Double_t m232, Double_t m122wd, Double_t m232wd,
                        Double_t m1, Double_t m2, Double_t m3, Double_t m){
   if((m122+0.5*m122wd) > (m-m3)*(m-m3) || (m122-0.5*m122wd) < (m1+m2)*(m1+m2))
      return false;
   if((m232+0.5*m232wd) > (m-m1)*(m-m1) || (m232-0.5*m232wd) < (m2+m3)*(m2+m3))
      return false;
   Double_t m23max2 = (pow(m*m - m1*m1 + m2*m2 - m3*m3,2)
                        -pow(TMath::Sqrt(pow(m122-m1*m1+m2*m2,2)-4*m2*m2*m122)-TMath::Sqrt(pow(m122-m*m+m3*m3,2)-4*m3*m3*m122),2))/(4*m122);
   if((m232+0.5*m232wd) > m23max2)
      return false;
   Double_t m23min2 = (pow(m*m - m1*m1 + m2*m2 - m3*m3,2)
                        -pow(TMath::Sqrt(pow(m122-m1*m1+m2*m2,2)-4*m2*m2*m122)+TMath::Sqrt(pow(m122-m*m+m3*m3,2)-4*m3*m3*m122),2))/(4*m122);
   if((m232-0.5*m232wd) < m23min2)
      return false;
   return true;
}


/// Function declaration
Int_t decayTo2Body(ExpParameters &, AxionParameters &, Double_t &, Particle *);
Int_t decayTo3Body(ExpParameters &, AxionParameters &, Double_t &, Particle *, TH2D*, TH2D**);
Int_t decayTo4Body(ExpParameters &, AxionParameters &, Double_t &, Particle *);
void GetAllFinalStates(Particle*, std::vector<Particle*>&);
void FillCDAzTAX(std::vector<Particle*>, TH2D*, const TVector3, const Double_t);

 /// \fn DecayMCProcess
 /// \Brief
 /// Main loop over masses and widths, simulating events
 /// \EndBrief
void DecayMCProcess(Int_t exoticIndex, Int_t experiment, Int_t productionmode, Int_t decaymode, Int_t numberOfMCevents,
					Bool_t acc = true, Int_t nAttempts = 1000, Bool_t verb = false, Int_t seed = 0,
					Bool_t flat = false, Bool_t dalitz = false, Int_t nBinsX = 101, Int_t nBinsY = 101,
					Bool_t xIsLin = false, Bool_t yIsLin = false, Int_t yVar = 0){
	std::cout << "xIsLin:" << xIsLin << " yIsLin:" << yIsLin << " yVar:" << yVar <<  " yName:" << yVars[yVar].Data() << std::endl;
	// put default values for optional arguments for running with interpreter
	TString exoName = exoLabels[exoticIndex];

	rndmGen.SetSeed(seed);

	allInAcceptance = acc;
	verbose = verb;
	flatDecay = flat;
	flatDalitz = dalitz;
	// check argument validity when running with interpreter
	if(productionmode<0||productionmode>=prodmodeNames[exoticIndex].size()){
		std::cout << "[Error] Invalid production mode: " <<  productionmode << ". Use values ";
		for(Int_t imode=0;imode<prodmodeNames[exoticIndex].size();imode++) std::cout << imode << " (" << prodmodeNames[exoticIndex][imode] << "), ";
		std::cout << std::endl;
		exit(1);
	}
	TString prodmodeName = prodmodeNames[exoticIndex][productionmode];

	if(decaymode<0||decaymode>=decaymodeNames[exoticIndex].size()){
		std::cout << "[Error] Invalid decay mode: " <<  decaymode << ". Use values ";
		for(Int_t imode=0;imode<decaymodeNames[exoticIndex].size();imode++) std::cout << imode << " (" << decaymodeNames[exoticIndex][imode] << "), ";
		std::cout << std::endl;
		exit(1);
	}
	TString decaymodeName = decaymodeNames[exoticIndex][decaymode];
	ExpParameters genParam(exoticIndex,experiment,prodmodeName,decaymode,nBinsX,nBinsY,xIsLin,yIsLin,yVar);
	// ExpParameters genParam(exoticIndex,experiment,prodmodeName,decaymode);

	std::cout << "[Info] Simulating "<< exoName << "s at " << genParam.GetExpLabel() << " in " << prodmodeName << " production mode and " << decaymodeName << " decay mode with POT " << genParam.GetPOT() << " and normCrossSec " << genParam.GetNormCrossSec() << std::endl;


	//Ensuring out_directory exists
	std::filesystem::path out_dir = Form("%s/%s/%s/",outPath.Data(),genParam.GetExpName().Data(),exoName.Data());
	if (verbose) std::cout<< "[Info] Trying to create "<< out_dir <<std::endl;
	std::filesystem::create_directories(out_dir);


	TString activeCoupling; // loop to cover the different possible hnl couplings
	for (Int_t activeCouplingMode = 0; activeCouplingMode<activeCouplings[exoticIndex].size(); ++activeCouplingMode){
		activeCoupling = activeCouplings[exoticIndex][activeCouplingMode];
		if(activeCoupling) std::cout << "\n[Info] Starting simulation for "<< activeCoupling << ":" <<std::endl;

		//variables for 3-body decays
		TH1D** rescaleDalitz1d;
		TH2D** rescaleDalitz;
		TH2D*** originalDalitz;
		std::vector<Double_t> mALPList; //list of allowed mass bins
		std::vector<Double_t> valGamma,m12List,m23List;
		std::vector<std::array<std::vector<Double_t>,3>> ALPlist; // each mALP has unique m12 and m23 axis + m12*m23 values of diff. width  (contains m12List and m23List and valGamma)

		if (genParam.finState.size() == 3 && !flatDalitz && (strcmp(exoName.Data(), "alp")||strcmp(exoName.Data(), "hnl"))){ //width tables only available for 3-body decays of alps and hnls
			std::cout << "[Info] Simulating 3-body decay, will read differential width table" << std::endl;
			// reading diff. width tables for 3-body decays
			std::ifstream widthFile;
			widthFile.open(Form("%s/%s%s.dat",widthPath.Data(),decaymodeName.Data(),activeCoupling.Data()));
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

			std::cout << "[Info] Plots for "<< exoName <<" widths declaration and filling:" << std::endl;

			rescaleDalitz = new TH2D*[ALPlist.size()];
			rescaleDalitz1d = new TH1D*[ALPlist.size()];
			originalDalitz = new TH2D**[ALPlist.size()]; // dalitz plot before and after acceptance criteria
			//// FOR DEBUGGING /////
			if(verbose){
				std::cout << "m12 begin: " << ALPlist[0][0].front() << " end: " << ALPlist[0][0].back() << std::endl;
				std::cout << "m23 begin: " << ALPlist[0][1].front() << " end: " << ALPlist[0][1].back() << std::endl;
			}
			////////////////
			for (Int_t mb=0; mb < ALPlist.size(); mb++){ // m12 and m23 steps are linear 
				Double_t wdM12 = (ALPlist[mb][0].back()-ALPlist[mb][0].front())/ALPlist[mb][0].size();
				Double_t wdM23 = (ALPlist[mb][1].back()-ALPlist[mb][1].front())/ALPlist[mb][1].size();
				rescaleDalitz[mb] =
					new TH2D(Form("DalitzDens%d",mb),Form("rescaleDalitz for %s mass %f GeV (bin %d);m12^2;m23^2;weight",exoName.Data(),TMath::Power(10.,mALPList[mb]),mb),
					ALPlist[mb][0].size(),ALPlist[mb][0].front()-0.5*wdM12,ALPlist[mb][0].back()+0.5*wdM12,
					ALPlist[mb][1].size(),ALPlist[mb][1].front()-0.5*wdM23,ALPlist[mb][1].back()+0.5*wdM23);
				originalDalitz[mb] = new TH2D*[4];
				for (int ib = 0; ib<4; ib++){
				originalDalitz[mb][ib] =
					new TH2D(Form("OriginalDalitz%d_befAft%d",mb,ib),Form("originalDalitz before/after %d for %s mass %f GeV;m12^2;m23^2",ib,exoName.Data(),TMath::Power(10.,mALPList[mb])),
						ALPlist[mb][0].size(),ALPlist[mb][0].front()-0.5*wdM12,ALPlist[mb][0].back()+0.5*wdM12,
						ALPlist[mb][1].size(),ALPlist[mb][1].front()-0.5*wdM23,ALPlist[mb][1].back()+0.5*wdM23);
				}
			}
			std::cout << "[Info] Plots declared." << std::endl;

		for (Int_t mb=0; mb < ALPlist.size(); mb++){
			Int_t linecounter=0;
			Int_t nBinsInsideBoundary=0;
			for (Int_t m12b=0; m12b < ALPlist[0][0].size(); m12b++){
				Double_t m122 = rescaleDalitz[mb]->GetXaxis()->GetBinCenter(m12b+1);
				Double_t m122wd = 2*(m122-rescaleDalitz[mb]->GetXaxis()->GetBinLowEdge(m12b+1));
				for (Int_t m23b=0; m23b < ALPlist[0][1].size(); m23b++){
					Double_t m232 = rescaleDalitz[mb]->GetYaxis()->GetBinCenter(m23b+1);
					Double_t m232wd = 2*(m232-rescaleDalitz[mb]->GetYaxis()->GetBinLowEdge(m23b+1));
					// set 0 if outside of kinematically allowed region
					if(inDalitzBoundary(m122,m232,m122wd,m232wd,genParam.finState.at(0)->GetMass(),genParam.finState.at(1)->GetMass(),genParam.finState.at(2)->GetMass(),TMath::Power(10,mALPList[mb]))){
						rescaleDalitz[mb]->SetBinContent(m12b+1, m23b+1, ALPlist[mb][2][linecounter]);
						nBinsInsideBoundary++;
					}
					else
						rescaleDalitz[mb]->SetBinContent(m12b+1, m23b+1, 0);
					linecounter+=1;
				}
			}
			rescaleDalitz[mb]->Scale((Double_t)nBinsInsideBoundary/rescaleDalitz[mb]->GetSum());
			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << "linecounter reached: " << linecounter << std::endl;
			/////////
			rescaleDalitz1d[mb] = rescaleDalitz[mb]->ProjectionX(Form("%s_m12",rescaleDalitz[mb]->GetName()));
		}

			std::cout << "[Info] Plots filled." << std::endl;

			genParam.SetMinMassA(TMath::Power(10.,mALPList.front()));
			genParam.SetMaxMassA(TMath::Power(10.,mALPList.back()));

		} else{
			Double_t minMass = 0.;
			for(auto finState : genParam.finState)
				minMass += finState->GetMass();
			genParam.SetMinMassA(minMass);
			if(prodmodeName == "Dmeson") genParam.SetMaxMassA(MDs);
			else genParam.SetMaxMassA(genParam.GetMaxMassAInFile()[1]); //5.310 GeV
			// genParam.SetMaxMassA()
		}
		std::cout << "[Info] For " << decaymodeName << " decay mode min "<<exoName<<" mass = " << genParam.GetMinMassA() << " [GeV] and max "<<exoName<<" mass = " << genParam.GetMaxMassA() << " [GeV] NmassFiles = " << genParam.GetNMassFiles() << '\n' << std::endl;

		//plots declaration
		std::cout << "[Info] Plots for "<< exoName << " distributions declaration:" << std::endl;

		std::vector<std::vector<TH1D*>> beforeLog;
		std::vector<std::vector<std::vector<TH1D*>>> afterLogMC;
		std::vector<std::vector<TH2D*>> ExpectedNum;
		//storing z and momentum
		std::vector<std::vector<TH2D*>> beforeLogFlat;
		std::vector<std::vector<std::vector<TH2D*>>> afterLogMCFlat;
		//input distribution
		std::vector<std::vector<TH2D*>> fHExoETheta;
		//TGenPhaseSpaceNormalisation
		std::vector<std::vector<Double_t>> tGenCorrection;
		//look-up-table of output mass bins to the Dalitz plot associated with next larger mass
		std::vector<std::vector<Int_t>> mb_to_Dalitz_mbs;
		Int_t curr_Dalitz_mb {};
		//store header information of the input file 
		TString input_header, input_coupling;

		for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){
			// //read in exotic spectrum
			std::cout << "[Info] Reading "<<exoName<<" distribution for file " << ifile << std::endl;
			std::ifstream infile;
			std::string infilePath = Form("%s/%s/%s_%s%s_beam%dGeV_%sto%sMeV_%s.dat", sourcePath.Data(), genParam.GetExpName().Data(),exoName.Data(),prodmodeName.Data(), activeCoupling.Data(), 
																					  genParam.GetBeamEnergy(), genParam.GetMinMassFile()[ifile].Data(), genParam.GetMaxMassFile()[ifile].Data(), genParam.GetExpName().Data());
			infile.open(infilePath);
			if (!infile.is_open()) {
				std::cout << "[Error] "<<exoName<<" momentum distribution table not found" << std::endl;
				std::cout << "[Error] with path: " << infilePath << std::endl;
				exit(1);
			}
			std::vector<std::array<std::vector<Double_t>,3>> ParticleList;
			std::vector<Double_t> MassList,EnergyList,ThetaList,ValList;
			Double_t Mass, Energy, Theta, Val;
			std::string line;
			while (std::getline(infile, line)) { // initial scan of input file
				if (!line.empty() && line[0] == '#') { //reading header information
					std::istringstream iss(line);
					if (!input_header.Length()) {
						std::string token;
						while(std::getline(iss, token, '|')) continue;
						input_header = token.substr(1);
					} else {
						while(iss >> input_coupling) continue;
						input_coupling = input_coupling(0,input_coupling.Length()-2);
					}
				} else break;
			}
			infile >> Theta >> Energy >> Mass >> Val; //first data line
			MassList.push_back(Mass);
			EnergyList.push_back(Energy);
			ThetaList.push_back(Theta);
			ValList.push_back(Val);
			while (infile >> Theta >> Energy >> Mass >> Val) {
				if(Mass == MassList.back()){ //same mass bin
					ValList.push_back(Val);
					if(Energy == EnergyList.back()) //repeating theta sequence
						ThetaList.push_back(Theta);
					else{ //
						EnergyList.push_back(Energy);
						ThetaList.clear();
						ThetaList.push_back(Theta);
					}
				} else{
					//updating ALP mass
					MassList.push_back(Mass);

					//// FOR DEBUGGING /////
					if(verbose)
						std::cout << "Printing number of theta,energy,value points for mass " << Mass << " : " << ThetaList.size() << " " << EnergyList.size() << " " << ValList.size() << std::endl;
					///////////
					
					//updating list
					ParticleList.push_back({ThetaList,EnergyList,ValList});
					ValList.clear();
					ValList.push_back(Val);
					ThetaList.clear();
					ThetaList.push_back(Theta);
					EnergyList.clear();
					EnergyList.push_back(Energy);
				}
			}
			if (infile.eof()){
				ParticleList.push_back({ThetaList,EnergyList,ValList});
				genParam.SetValuesEnergyInFile(EnergyList.size());
				genParam.SetMinEnergyInFile(EnergyList.front());
				genParam.SetMaxEnergyInFile(EnergyList.back());
				genParam.SetWdEnergy((EnergyList.back()-EnergyList.front())/(EnergyList.size()-1));
				genParam.SetValuesThetaInFile(ThetaList.size());
				genParam.SetMinThetaInFile(ThetaList.front());
				genParam.SetMaxThetaInFile(ThetaList.back());
				genParam.SetWdTheta((ThetaList.back()-ThetaList.front())/(ThetaList.size()-1));

				ValList.clear();
				EnergyList.clear();
				ThetaList.clear();
			} else{
				std::cout << "[Error] Did not reach EoF when reading input table" << std::endl;
				exit(1);
			}

			infile.close();

			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << " declaration for mass file: #" << ifile << " M_"<<exoName<<" in range " << genParam.GetMinMassFile().at(ifile) << "-" << genParam.GetMaxMassFile().at(ifile) << " Npoints = " << genParam.GetNMassX().at(ifile) << std::endl;
			////////////////
			if(flatDecay){
				std::vector<TH2D*> beforeLogFile;
				std::vector<std::vector<TH2D*>> afterLogMCFile;
				for (Int_t mb=0; mb < genParam.GetNMassX().at(ifile); mb++){
					beforeLogFile.push_back(new TH2D(Form("beforeLogFlat[%d]%d",ifile,mb),"beforeLog (flat);momentum;zvtx",100,0,genParam.GetBeamEnergy(),100,genParam.GetZFVIn(),genParam.GetZFVEnd()));
					std::vector<TH2D*> afterLogMCbin;
					if(genParam.GetNRegions() < 1)
						afterLogMCbin.push_back(new TH2D(Form("afterLogMCFlat[%d]%d",ifile,mb),"afterLogMC (flat);momentum;zvtx",100,0,genParam.GetBeamEnergy(),100,genParam.GetZFVIn(),genParam.GetZFVEnd()));
					else
						for (Int_t iReg=0; iReg < genParam.GetNRegions(); iReg++)
							afterLogMCbin.push_back(new TH2D(Form("afterLogMCFlat[%d]%d-reg%d",ifile,mb,iReg+1),"afterLogMC (flat);momentum;zvtx",100,0,genParam.GetBeamEnergy(),100,genParam.GetZFVIn(),genParam.GetZFVEnd()));
					afterLogMCFile.push_back(afterLogMCbin);
				} //massBins
				beforeLogFlat.push_back(beforeLogFile);
				afterLogMCFlat.push_back(afterLogMCFile);
				beforeLogFile.clear();
				afterLogMCFile.clear();
			}
			std::vector<TH1D*> beforeLogFile;
			std::vector<std::vector<TH1D*>> afterLogMCFile;

			mb_to_Dalitz_mbs.push_back(std::vector<Int_t> (genParam.GetNMassX().at(ifile),-1));

			for (Int_t mb=0; mb < genParam.GetNMassX().at(ifile); mb++){
				// beforeLogFile.push_back(new TH1D(Form("beforeLog[%d]%d",ifile,mb),"beforeLog;width in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
				beforeLogFile.push_back(new TH1D(Form("beforeLog[%d]%d",ifile,mb),Form("beforeLog;%s;weight",yVars[yVar].Data()),  genParam.nY,genParam.minY-0.5*genParam.wdY,genParam.maxY+0.5*genParam.wdY));
				std::vector<TH1D*> afterLogMCbin;
				if(genParam.GetNRegions() < 1)
					// afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d",ifile,mb),"afterLogMC;width in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
					afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d",ifile,mb),Form("afterLogMC;%s;weight",yVars[yVar].Data()), genParam.nY,genParam.minY-0.5*genParam.wdY,genParam.maxY+0.5*genParam.wdY));
				else
					for (Int_t iReg=0; iReg < genParam.GetNRegions(); iReg++)
						// afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d-reg%d",ifile,mb,iReg+1),"afterLogMC;width in log;weight", genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
						afterLogMCbin.push_back(new TH1D(Form("afterLogMC[%d]%d-reg%d",ifile,mb,iReg+1),Form("afterLogMC;%s;weight",yVars[yVar].Data()), genParam.nY,genParam.minY-0.5*genParam.wdY,genParam.maxY+0.5*genParam.wdY));
				afterLogMCFile.push_back(afterLogMCbin);
			} //massBins
			beforeLog.push_back(beforeLogFile);
			afterLogMC.push_back(afterLogMCFile);
			beforeLogFile.clear();
			afterLogMCFile.clear();

			std::vector<TH2D*> ExpectedNumFile;
			if(genParam.GetNRegions() < 1)
				// ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]",ifile),"Expected Num",genParam.GetNMassX().at(ifile),genParam.GetMinMassX().at(ifile)-0.5*genParam.GetWdMassX().at(ifile),genParam.GetMaxMassX().at(ifile)+0.5*genParam.GetWdMassX().at(ifile),genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
				ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]",ifile),Form("Expected Num;m;%s",yVars[yVar].Data()),genParam.GetNMassX().at(ifile),genParam.GetMinMassX().at(ifile)-0.5*genParam.GetWdMassX().at(ifile),genParam.GetMaxMassX().at(ifile)+0.5*genParam.GetWdMassX().at(ifile),genParam.nY,genParam.minY-0.5*genParam.wdY,genParam.maxY+0.5*genParam.wdY));
			else
				for (Int_t iReg=0; iReg < genParam.GetNRegions(); iReg++)
					// ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]-reg%d",ifile,iReg+1),"Expected Num",genParam.GetNMassX().at(ifile),genParam.GetMinMassX().at(ifile)-0.5*genParam.GetWdMassX().at(ifile),genParam.GetMaxMassX().at(ifile)+0.5*genParam.GetWdMassX().at(ifile),genParam.GetNWidths(),genParam.GetMinWidth()-0.5*genParam.GetWdWidth(),genParam.GetMaxWidth()+0.5*genParam.GetWdWidth()));
					ExpectedNumFile.push_back(new TH2D(Form("ExpectedNum[%d]-reg%d",ifile,iReg+1),Form("Expected Num;m;%s",yVars[yVar].Data()),genParam.GetNMassX().at(ifile),genParam.GetMinMassX().at(ifile)-0.5*genParam.GetWdMassX().at(ifile),genParam.GetMaxMassX().at(ifile)+0.5*genParam.GetWdMassX().at(ifile), genParam.nY,genParam.minY-0.5*genParam.wdY,genParam.maxY+0.5*genParam.wdY));
			ExpectedNum.push_back(ExpectedNumFile);
			ExpectedNumFile.clear();

			genParam.SetValuesMassAInFile(ifile,MassList.size());
			genParam.SetMinMassAInFile(ifile,MassList.front());
			genParam.SetMaxMassAInFile(ifile,MassList.back());
			genParam.SetWdMassA(ifile,(MassList.back()-MassList.front())/(MassList.size()-1));
			genParam.SetLogWdMassA(ifile,(log10(MassList.back())-log10(MassList.front()))/(MassList.size()-1));



			//Input scale test
			Double_t linHypothesis = log10(genParam.GetMinMassAInFile().at(ifile) + ((Double_t) ((MassList.size()-1)/2)) * genParam.GetWdMassA().at(ifile));
			Double_t logHypothesis = log10(genParam.GetMinMassAInFile().at(ifile)) + ((Double_t) ((MassList.size()-1)/2)) * genParam.GetLogWdMassA().at(ifile);
			Double_t actual = log10(MassList.at((MassList.size()-1)/2));
			if( 1e-4 < TMath::Abs(actual-logHypothesis)/actual   && 1e-4 < TMath::Abs(actual-linHypothesis)/actual){
				std::cout << "[Error] could not interpret mass input scale with linear hypothesis " << linHypothesis << " and log hypothesis "<< logHypothesis <<" for target "<< actual <<"."<< std::endl; 
				exit(1);
			}
			Bool_t inputMassesAreLog = (TMath::Abs(actual-logHypothesis) < TMath::Abs(actual-linHypothesis)) ? kTRUE: kFALSE;
			if(inputMassesAreLog) std::cout << "[Info] Interpreted input masses as log10 scale"<< std::endl; 
			else std::cout << "[Info] Interpreted input masses as lin scale"<< std::endl; 

			std::cout << "[Info] Writing "<<exoName<<" distribution for file " << ifile << std::endl;
			std::vector<TH2D*> fHExoEThetaFile;

			//for correction
			TGenPhaseSpace normEvent;
			std::vector<Double_t> tGenCorrectionMF;
			Double_t minMass = 0.0;
			std::vector<Double_t> finMasses;
			for(auto finState : genParam.finState){
				minMass += finState->GetMass();
				finMasses.push_back(finState->GetMass());
			}
			for (Int_t mb=0; mb < genParam.GetNMassX().at(ifile); mb++){ //loop over output mass bins mapping identified input scale on output
				Double_t massX = genParam.GetMinMassX().at(ifile) + genParam.GetWdMassX().at(ifile)*mb;
				Double_t logMassX;
				if(!xIsLin){
					logMassX = massX;
					massX = TMath::Power(10.,massX);
				} else logMassX = TMath::Log10(massX);

				//reading exotic files
				Int_t mbIn = inputMassesAreLog? TMath::Floor((logMassX-log10(genParam.GetMinMassAInFile().at(ifile)))/genParam.GetLogWdMassA().at(ifile) + 0.5) : 
												TMath::Floor((massX-genParam.GetMinMassAInFile().at(ifile))/genParam.GetWdMassA().at(ifile) + 0.5);
				//corresponding input bin for output bin mb 
			
				// sanity checks for bin allocation 
				if(inputMassesAreLog && ((MassList[mbIn]<TMath::Power(10.,logMassX - 0.5* genParam.GetWdMassX().at(ifile) )) || (TMath::Power(10.,logMassX + 0.5* genParam.GetWdMassX().at(ifile) ))<MassList[mbIn] ) ){
					std::cout << "[Error] Mass specified for "<<exoName<<" is out of range. Input mass bin is " << mbIn <<"." << std::endl;
					std::cout << "[Error] " << MassList[mbIn]<< " outside of allocated masswindow ("<< (TMath::Power(10.,logMassX - 0.5* genParam.GetWdMassX().at(ifile) ))<<", "<< (TMath::Power(10.,logMassX + 0.5* genParam.GetWdMassX().at(ifile) )) <<")"<< std::endl;
					exit(1);
				}else if(!inputMassesAreLog && ( (MassList[mbIn]< massX - 0.5 * genParam.GetWdMassA().at(ifile) )|| (massX + 0.5* genParam.GetWdMassA().at(ifile) <MassList[mbIn] ) )){
					std::cout << "[Error] Mass specified for "<<exoName<<" is out of range. Input mass bin is " << mbIn <<"." << std::endl;
					std::cout << "[Error] " << MassList[mbIn]<< " outside of allocated masswindow ("<< (massX - 0.5 * genParam.GetWdMassA().at(ifile) )<<", "<< (massX + 0.5 * genParam.GetWdMassA().at(ifile) ) <<")"<< std::endl;
					exit(1);
				}
				if (mbIn < 0) {
					std::cout << "[Error] Mass specified for "<<exoName<<" is out of range. Input mass bin is " << mbIn << std::endl;
					std::cout << "[Error] problem occured with log masswid= " << genParam.GetLogWdMassA().at(ifile) << " mAMinInFile " << genParam.GetMinMassAInFile().at(ifile) << " and mass: " << massX << std::endl;
					exit(1);
				}
				if (mbIn >= genParam.GetValuesMassAInFile().at(ifile)){
					std::cout << "[Error] Mass specified for "<<exoName<<" is out of range. Input mass bin is " << mbIn << std::endl;
					std::cout << "[Error] problem occured with log masswid= " << genParam.GetLogWdMassA().at(ifile) << " mAMinInFile " << genParam.GetMinMassAInFile().at(ifile) << " and mass: " << massX << std::endl;
					exit(1);
				}

				fHExoEThetaFile.push_back(new TH2D(Form("%s_ETheta[%d]%d",exoName.Data(),ifile,mb), Form("%s ETheta for %s mass %f GeV (bin:%d)",exoName.Data(),exoName.Data(),massX,mb),
						genParam.GetValuesThetaInFile(), genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta(),genParam.GetMaxThetaInFile()+0.5*genParam.GetWdTheta(),
						genParam.GetValuesEnergyInFile(), genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy(),genParam.GetMaxEnergyInFile()+0.5*genParam.GetWdEnergy())); // Assure that Bin spacing does not overshoot to negative bin edge!
				if((genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy())<0){std::cout << "WARNING Min Energy out of bounds : " << (genParam.GetMinEnergyInFile()-0.5*genParam.GetWdEnergy()) << std::endl;}
				if((genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta())<0){std::cout << "WARNING Min Theta out of bounds : " << (genParam.GetMinThetaInFile()-0.5*genParam.GetWdTheta()) << std::endl;}
				Int_t linecounter=0;
				for (Int_t iEnergy=1; iEnergy <= ParticleList[mbIn][1].size(); iEnergy++){
					for (Int_t iTheta=1; iTheta <= ParticleList[mbIn][0].size(); iTheta++){
						fHExoEThetaFile.back()->SetBinContent(iTheta, iEnergy, ParticleList[mbIn][2][linecounter]);
						linecounter+=1;
						//// FOR DEBUGGING /////
						if(verbose)
							std::cout << "Line: " << linecounter << " corresponds to theta bin: " << iTheta << " | energy bin: " << iEnergy << " | mass bin: " << mb << " | and value: " << ParticleList[mb][2][linecounter] << std::endl;
						//////
					}
				}
				//TGenPhaseSpace correction
				if(genParam.finState.size() > 2){
					TLorentzVector exo(0,0,0,massX);
					if(massX >= minMass){
						normEvent.SetDecay(exo,finMasses.size(), finMasses.data());
						Double_t normalization = 0;
						for (Int_t n=0;n<1000000;n++)
							normalization += normEvent.Generate()/1000000;
						if(normalization)
							tGenCorrectionMF.push_back(1./normalization);
						else
							tGenCorrectionMF.push_back(0.);
					}
					else
						tGenCorrectionMF.push_back(0.);
					if(!flatDalitz){//DalitzMatching
			 		while((logMassX - mALPList[curr_Dalitz_mb]) > ( (fabs(mALPList[curr_Dalitz_mb]) < fabs(logMassX) ? fabs(logMassX) : fabs(mALPList[curr_Dalitz_mb])) * 1.E-5)){ // evaluates for mALPList[curr_Dalitz_mb] definitelyLessThan logMassX 
						++curr_Dalitz_mb;
						if (curr_Dalitz_mb == mALPList.size()){
							std::cout << "[Error] Mass bin with mass "<< massX <<"GeV above pregenerated Dalitz table mass range(max="<<pow(10,mALPList.back())<<"GeV). Please extend Dalitz table at least up to "<< pow(10,genParam.GetMaxMassX().back())<<"GeV or use the \"--flat-Dalitz\" flag " << std::endl;
							exit(1);
						}
					}
					mb_to_Dalitz_mbs.back()[mb] = curr_Dalitz_mb;
				}	
				} else
					tGenCorrectionMF.push_back(1.); //no correction
			}

			fHExoETheta.push_back(fHExoEThetaFile);
			fHExoEThetaFile.clear();
			tGenCorrection.push_back(tGenCorrectionMF);
			tGenCorrectionMF.clear();
			ParticleList.clear();

		} //massFiles

		for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){

			std::cout << std::endl << "[Info] Entering MC for mass range " << genParam.GetMinMassAInFile()[ifile] << " - " << genParam.GetMaxMassAInFile()[ifile] << std::endl;

			TRandom3 *rndmGen = new TRandom3(0);
			rndmGen->SetSeed(0);

			AxionParameters alpParam;
			
			Double_t rProduction_a = 0;
			Double_t eventWeight = 0;
			Double_t sigAcc = genParam.GetSigacceptance();
			Int_t signalRegion = 0; // 0 if not in acceptance, otherwise number of signal region
			if (genParam.GetDecayModeName() != "2Gamma" || genParam.GetDecayModeName() != "2El") //2 charged hadron or mu efficiency
				sigAcc = genParam.GetSigacceptancemumu();

			//mass loop
			for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
				//progress bar:
				std::cout << "\r" << " Simulating evts for mass bin: " << mb << "    "<< std::flush;

				alpParam.massA = genParam.GetMinMassX()[ifile] + genParam.GetWdMassX()[ifile]*mb;
				if(!xIsLin)	alpParam.massA = TMath::Power(10.,alpParam.massA);
				if (alpParam.massA < genParam.GetMinMassA()) continue; //ALP mass constraint
				if (alpParam.massA < genParam.GetMinMassAInFile()[ifile] || alpParam.massA > genParam.GetMaxMassAInFile()[ifile]) continue;
				if(genParam.GetMaxMassA() != 0 && alpParam.massA > genParam.GetMaxMassA()) continue;
				alpParam.crossSecA = fHExoETheta[ifile][mb]->Integral("width");
				if(alpParam.crossSecA <= 0) continue;

				exoProp.Mass = alpParam.massA;
				exoProp.Position.SetXYZ(genParam.GetX0(),genParam.GetY0(),genParam.GetZ0());
				exotic = new Particle(&exoProp);

				Int_t iEv=0;
				Int_t iAtt=0;
				while(iEv < numberOfMCevents){

					//Pick random from input distributions
					fHExoETheta[ifile][mb]->GetRandom2(alpParam.thetaA,alpParam.energyA);
					alpParam.phiA = rndmGen->Rndm()*2*TMath::Pi(); //random uniform

					Double_t pA2 = alpParam.energyA*alpParam.energyA - alpParam.massA*alpParam.massA;
					if(pA2<0) pA2 = 0.; //can happen since we pick randomly from the energy table
					alpParam.pA = TMath::Sqrt(pA2);
					alpParam.betaA = alpParam.pA/alpParam.energyA;
					// if(alpParam.pA == 0){ //skip, nothing to simulate
					// 	// tMC->Fill();
					// 	continue; //?? do something?
					// }

					if (genParam.GetPhiBeamEuler() != 0 || genParam.GetThetaBeamEuler() != 0 || genParam.GetPsiBeamEuler() != 0){
						TVector3 oldAlp(TMath::Sin(alpParam.thetaA)*TMath::Cos(alpParam.phiA),
								TMath::Sin(alpParam.thetaA)*TMath::Sin(alpParam.phiA),
								TMath::Cos(alpParam.thetaA));
						TRotation a;
						//// FOR DEBUGGING /////
						if(verbose)
							std::cout << "Original X:" << oldAlp.X() << " | Y:" << oldAlp.Y() << " | Z:" << oldAlp.Z() << std::endl;
						///////////////////////
						a.SetXEulerAngles(genParam.GetPhiBeamEuler(),genParam.GetThetaBeamEuler(),genParam.GetPsiBeamEuler());
						TVector3 newAlp = a*oldAlp; // active rotation
						//// FOR DEBUGGING /////
						if(verbose)
							std::cout << "New X:" << newAlp.X() << " | Y:" << newAlp.Y() << " | Z:" << newAlp.Z() << std::endl;
						//////////////////////
						
						alpParam.thetaA = newAlp.Theta();
						alpParam.phiA = newAlp.Phi();
					}

					exoProp.Momentum.SetXYZ(alpParam.pA*TMath::Sin(alpParam.thetaA)*TMath::Cos(alpParam.phiA),
											alpParam.pA*TMath::Sin(alpParam.thetaA)*TMath::Sin(alpParam.phiA),
											alpParam.pA*TMath::Cos(alpParam.thetaA));
					exotic->Reset(exoProp.Momentum);

					//increase efficiency by skipping ALPs out of calorimeter acceptance
					Bool_t ALPInAcceptance = true;
					if(allInAcceptance){
						exotic->PropagateToZ(genParam.GetZECal());
						ALPInAcceptance = genParam.inCaloOuterEdgeAcceptance(exotic->GetPosition());
					}
					eventWeight = 1.;
					if(flatDecay){
						// for flat reweight we will need to store distance travelled in FV and exotic momentum
						rProduction_a = genParam.GetBeamDecayLength()*rndmGen->Rndm(); // flat
						// rProduction_a = genParam.GetBeamDecayLength()/2; // flat
						exoProp.ProductionCrossSection = alpParam.crossSecA;
						exotic->PropagateToZ(genParam.GetZFVIn());
						exotic->Propagate(rProduction_a);
						// std::cout << "rProduction_a:" << rProduction_a << std::endl;
						// std::cout << "genParam.GetBeamDecayLength():" << genParam.GetBeamDecayLength() << std::endl;
						// std::cout << "genParam.GetZFVIn():" << genParam.GetZFVIn() << std::endl;
						// std::cout << "genParam.GetZFVEnd():" << genParam.GetZFVEnd() << std::endl;
						// std::cout << "exotic->GetPosition().Z():" << exotic->GetPosition().Z() << std::endl;
						beforeLogFlat[ifile][mb]->Fill(alpParam.pA,exotic->GetPosition().Z(),eventWeight);
						// tMC->Fill(); // fill TTree
						if(allInAcceptance && !ALPInAcceptance || alpParam.pA == 0){
							// tMC->Fill();
							iAtt++; // count attempts to simulate event
							if(iAtt > nAttempts)
								iEv++;
							continue;
						}
						Double_t decayWeight = 1.;
						// eventWeight = tGenCorrection[ifile][mb]*decayWeight*alpParam.crossSecA*genParam.GetNormCrossSec()*sigAcc; //fill TTrees
						eventWeight = tGenCorrection[ifile][mb]*decayWeight*genParam.GetPOT()*alpParam.crossSecA*genParam.GetNormCrossSec()*sigAcc;
						// std::cout << "eventWeight:" << eventWeight/genParam.GetPOT() << std::endl;
						// std::cout << "alpParam.crossSecA:" << alpParam.crossSecA << std::endl;
						// std::cout << "corrected alpParam.crossSecA:" << alpParam.crossSecA*tGenCorrection[ifile][mb] << std::endl;
						// std::cout << "alpParam.massA:" << alpParam.massA << std::endl;


						if (genParam.finState.size() == 2) signalRegion = decayTo2Body(genParam,alpParam,eventWeight,exotic);
						else if(genParam.finState.size() == 3){
							if(!flatDalitz)
								signalRegion = decayTo3Body(genParam,alpParam,eventWeight,exotic,rescaleDalitz[mb_to_Dalitz_mbs[ifile][mb]],originalDalitz[mb_to_Dalitz_mbs[ifile][mb]]);
							else signalRegion = decayTo3Body(genParam,alpParam,eventWeight,exotic,nullptr,nullptr);
						} else if (genParam.finState.size() == 4) signalRegion = decayTo4Body(genParam,alpParam,eventWeight,exotic);

						// std::cout << "set weight:" << eventWeight << std::endl;
						// std::cout << "set weight no POT:" << eventWeight/genParam.GetPOT() << std::endl;
						if(signalRegion){
								afterLogMCFlat[ifile][mb].at(signalRegion-1)->Fill(alpParam.pA,exotic->GetPosition().Z(),eventWeight);
								// tMCAcc->Fill();//fill TTree
							iEv++;
							iAtt=0; // reset when successful
						} else{
							iAtt++; // count attempts to simulate event
							if(iAtt > nAttempts)
								iEv++;
						}
						// tMC->Fill();
					} else{
						//decay width loop
						// for (Int_t nW=0; nW < genParam.GetNWidths(); nW++){
						for (Int_t iY=0; iY < genParam.nY; iY++){ //ctau
							// std::cout << "iY:" << iY << " genParam.nY:" << genParam.nY << " genParam.minY:" << genParam.minY << " genParam.wdY:" << genParam.wdY << std::endl;
							Double_t valY = genParam.minY + genParam.wdY*iY;
							Double_t valYTrue = valY;
							if(!yIsLin) valYTrue = TMath::Power(10.,valYTrue);
							exoProp.DecayWidth = valYTrue; //default
							if(yVar == 1) exoProp.DecayWidth = hc/ (valYTrue*c*1e15);//tau [fs]
							else if(yVar == 2) exoProp.DecayWidth = hc/ (valYTrue);//ctau [m]
							alpParam.widthExpA = TMath::Log10(exoProp.DecayWidth); //obsolete
							// std::cout << "valY:" << valY << "valYTrue:" << valYTrue << " exoProp.DecayWidth" << exoProp.DecayWidth << std::endl;
							// alpParam.widthExpA = genParam.GetMinWidth() + genParam.GetWdWidth()*nW;  //y-axis
								eventWeight = 1.;
							// beforeLog[ifile][mb]->Fill(alpParam.widthExpA,eventWeight);
								beforeLog[ifile][mb]->Fill(valY,eventWeight);
							// tMC->Fill(); // fill TTree
							if(allInAcceptance && !ALPInAcceptance || alpParam.pA == 0){
								// tMC->Fill();
								iAtt++; // count attempts to simulate event
								if(iAtt > nAttempts)
									iEv++;
								continue;
							}
							alpParam.decayLengthA = alpParam.pA/alpParam.massA*hc* 1./exoProp.DecayWidth; // [m] beta*gamma/(Gamma/hc)
							
							Double_t w1 = (1-TMath::Exp(-genParam.GetBeamDecayLength()/alpParam.decayLengthA)); // probability to decay in FV
							rProduction_a = -alpParam.decayLengthA*TMath::Log(1.-w1*rndmGen->Rndm()); // now OK, see documentation for NA62 analysis

							exoProp.ProductionCrossSection = alpParam.crossSecA;

							exotic->Reset(exoProp.Momentum);
							exotic->PropagateToZ(genParam.GetZFVIn());
							Double_t decayWeight = TMath::Exp(-(exotic->GetPosition()-exotic->GetInitPosition()).Mag()/alpParam.decayLengthA); // probability to reach the FV and decay therein
							exotic->Propagate(rProduction_a);

							// eventWeight = tGenCorrection[ifile][mb]*w1*decayWeight*alpParam.crossSecA*genParam.GetNormCrossSec()*sigAcc; //for ttree, no POT
							eventWeight = tGenCorrection[ifile][mb]*w1*decayWeight*genParam.GetPOT()*alpParam.crossSecA*genParam.GetNormCrossSec()*sigAcc;

							if (genParam.finState.size() == 2) signalRegion = decayTo2Body(genParam,alpParam,eventWeight,exotic);
							else if(genParam.finState.size() == 3){
								if(!flatDalitz) 
									signalRegion = decayTo3Body(genParam,alpParam,eventWeight,exotic, rescaleDalitz[mb_to_Dalitz_mbs[ifile][mb]], originalDalitz[mb_to_Dalitz_mbs[ifile][mb]]);
								else signalRegion = decayTo3Body(genParam,alpParam,eventWeight,exotic,nullptr, nullptr);
							} else if (genParam.finState.size() == 4) signalRegion = decayTo4Body(genParam,alpParam,eventWeight,exotic);

							if(signalRegion){
								afterLogMC[ifile][mb].at(signalRegion-1)->Fill(valY, eventWeight); //ctau
								// afterLogMC[ifile][mb].at(signalRegion-1)->Fill(alpParam.widthExpA, eventWeight);
									// tMCAcc->Fill(); //fill TTree
								iEv++;
								iAtt=0; // reset when successful
							} else{
								iAtt++; // count attempts to simulate event
								if(iAtt > nAttempts)
									iEv++;
							}
						// tMC->Fill();
						} // width loop
					}
				} // ievt
				delete exotic;

			} // mass loop

			std::cout << std::endl << "[Info] Event loop finished! Will now plot efficiencies" << std::endl;

			std::vector<std::vector<TH1D*>> eff1Log;
			// std::vector<std::vector<TH2D*>> eff1LogFlat;
			//// FOR DEBUGGING /////
			if(verbose)
				std::cout << "new event " << std::endl;
			//////////////

			for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
				Double_t massX = genParam.GetMinMassX()[ifile] + genParam.GetWdMassX()[ifile]*mb;
				Double_t massXTrue = massX;
				if(!xIsLin)	massXTrue = TMath::Power(10.,massX);
				//// FOR DEBUGGING /////
					if(verbose)
					std::cout << "Output mass:" << massXTrue << " bin:" << mb << std::endl;
				////////////////////////

				if(flatDecay){ //fill nume and deno with reweight
					// for (Int_t nW=0; nW < genParam.GetNWidths(); nW++){
					for (Int_t iY=0; iY < genParam.nY; iY++){ //ctau
					Double_t valY = genParam.minY + genParam.wdY*iY;
					Double_t valYTrue = valY;
					if(!yIsLin) valYTrue = TMath::Power(10.,valYTrue);
					Double_t width = valYTrue; //default
					if(yVar == 1) width = hc/ (valYTrue*c*1e15);//tau [fs]
					else if(yVar == 2) width = hc/ (valYTrue);//ctau [m]
						for (Int_t iMom=1; iMom<=beforeLogFlat[ifile][mb]->GetNbinsX(); iMom++) {
							Double_t decayLength = beforeLogFlat[ifile][mb]->GetXaxis()->GetBinCenter(iMom)/massXTrue*hc/width;
							
							for (Int_t iZ=1; iZ<=beforeLogFlat[ifile][mb]->GetNbinsY(); iZ++) {
								Double_t decayWeight = TMath::Exp(-(genParam.GetZFVIn()-genParam.GetZ0())/decayLength)*(1-TMath::Exp(-genParam.GetBeamDecayLength()/decayLength));
								Double_t extraReweight = TMath::Exp(-(beforeLogFlat[ifile][mb]->GetYaxis()->GetBinCenter(iZ)-genParam.GetZFVIn())/decayLength)/((1-TMath::Exp(-genParam.GetBeamDecayLength()/decayLength))*decayLength);
								
								beforeLog[ifile][mb]->Fill(valY,extraReweight*beforeLogFlat[ifile][mb]->GetBinContent(iMom,iZ));
								// beforeLog[ifile][mb]->Fill(widthExp,beforeLogFlat[ifile][mb]->GetBinContent(iMom,iZ)); //for TTree
							
								for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++)
									// afterLogMC[ifile][mb].at(iReg)->Fill(widthExp,afterLogMCFlat[ifile][mb].at(iReg)->GetBinContent(iMom,iZ)); //for TTree
								afterLogMC[ifile][mb].at(iReg)->Fill(valY,decayWeight*extraReweight*afterLogMCFlat[ifile][mb].at(iReg)->GetBinContent(iMom,iZ));

							}
						}
					}
				}
				std::vector<TH1D*> eff1LogBin;
				TH1D* deno = (TH1D*) beforeLog[ifile][mb]->Clone();
				for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
					TH1D* nume = (TH1D*) afterLogMC[ifile][mb].at(iReg)->Clone();
					eff1LogBin.push_back((TH1D*) nume->Clone(Form("eff1Log_bin_%d_region_%d",mb,iReg+1)));
					for(Int_t jb=1; jb<=eff1LogBin.back()->GetNbinsX(); jb++)
						if(deno->GetBinContent(jb)>0)
							eff1LogBin.back()->SetBinContent(jb,nume->GetBinContent(jb)/deno->GetBinContent(jb));
					delete nume;
				}
				delete deno;
				eff1Log.push_back(eff1LogBin);
				eff1LogBin.clear();
				for (Int_t i=1; i<=beforeLog[ifile][mb]->GetNbinsX(); i++) {
					Double_t valY = beforeLog[ifile][mb]->GetBinCenter(i);
					for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
						ExpectedNum[ifile].at(iReg)->Fill(massX,valY, eff1Log[mb][iReg]->GetBinContent(i));
					}
					
				} // massbins
			} // mass

			std::cout << "[Info] Writing root output " << std::endl;

			TString filename;
			filename = Form("%s/%s/%s/experiment%dprodmode%d_MCevents%d_%sto%sMeV_decaymode%s.root",outPath.Data(),genParam.GetExpName().Data(),exoName.Data(),genParam.GetExpNum(),productionmode,numberOfMCevents,genParam.GetMinMassFile()[ifile].Data(),genParam.GetMaxMassFile()[ifile].Data(),decaymodeName.Data());
			TFile *fileOut = new TFile(filename.Data(),"RECREATE");
			fileOut->cd();
			// tMC->Write(); //for TTree
			// tMCAcc->Write(); //for TTree
			for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++)
				ExpectedNum[ifile].at(iReg)->Write();

			for (Int_t mb=0; mb < genParam.GetNMassX()[ifile]; mb++){
				if(flatDecay){
					beforeLogFlat[ifile][mb]->Write();
					delete beforeLogFlat[ifile][mb];
				}
				beforeLog[ifile][mb]->Write();
				delete beforeLog[ifile][mb];
				
				fHExoETheta[ifile][mb]->Write();
				delete fHExoETheta[ifile][mb];

				for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
					if(flatDecay){
						afterLogMCFlat[ifile][mb].at(iReg)->Write();
						delete afterLogMCFlat[ifile][mb].at(iReg);
					}
					afterLogMC[ifile][mb].at(iReg)->Write();
					delete afterLogMC[ifile][mb].at(iReg);

					eff1Log[mb][iReg]->Write();
					delete eff1Log[mb][iReg]; 
				}
			}

			fileOut->Close();
			delete fileOut;

			//// FOR DEBUGGING - Dalitz plot rescale check/////
			if(verbose){
				TString testFilename = Form("%s/%s/%s/DalitzStudy_experiment%dprodmode%d_MCevents%d_%sto%sMeV_decaymode%s.root",outPath.Data(),genParam.GetExpName().Data(),exoName.Data(),genParam.GetExpNum(),productionmode,numberOfMCevents,genParam.GetMinMassFile()[ifile].Data(),genParam.GetMaxMassFile()[ifile].Data(),decaymodeName.Data());//_%s.root",outPath.Data(),genParam.expName.Data(),decaymodeName.Data());
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
		std::ofstream RunOutput; // Export data points
		TString outFileNameBase = Form("%s/%s/%s/%s_%s%s_%s",outPath.Data(),genParam.GetExpName().Data(),exoName.Data(),genParam.GetExpLabel().Data(),prodmodeName.Data(),activeCoupling.Data(),decaymodeName.Data());
		if(xIsLin) outFileNameBase += "_xLin";
		if(yIsLin) outFileNameBase += "_yLin";
		if(yVar) outFileNameBase += "_" + yVars[yVar];
			
		for(Int_t iReg = 0; iReg < genParam.GetNRegions(); iReg++){
			TString outFileName = outFileNameBase;
			if (genParam.GetNRegions() > 1) outFileName += Form("_reg%d.dat",iReg+1);
			else outFileName += ".dat";
			// if (genParam.GetNRegions() > 1) outFileName = Form("%s/%s/%s/%s_%s_%s_reg%d.dat",outPath.Data(),genParam.GetExpName().Data(),exoName.Data(),genParam.GetExpLabel().Data(),prodmodeName.Data(),decaymodeName.Data(),iReg+1);

			/// open and close to clear the file
			RunOutput.open(outFileName.Data());
			RunOutput.close();

			// open again
			RunOutput.open(outFileName.Data(), std::ios::app);
			RunOutput << Form("# Expected yield for %s%s at %s experiment produced in %s and decaying to %s | %s ",exoName.Data(),activeCoupling.Data(),genParam.GetExpLabel().Data(),prodmodeName.Data(),decaymodeName.Data(),exoName.Data()) <<"decays were generated with statistics of "<< numberOfMCevents<< " | production was "<< input_header <<"\n";
			RunOutput << Form("# m_X[GeV] Gamma_X[GeV] Y[per(%s Br(%s->%s))]\n", input_coupling.Data(),exoName.Data(),decaymodeName.Data());
			// extract correct limit plot
			// extract exactly "paperstyle"", to have compareable ranges easily
			Double_t xPrevious = 0;
			for (Int_t ifile=0; ifile < genParam.GetNMassFiles(); ifile++){
				for( int i=1; i<=ExpectedNum[ifile].at(iReg)->GetNbinsX(); i++){ //0 is underflow bin
					x=ExpectedNum[ifile].at(iReg)->GetXaxis()->GetBinCenter(i); //GeV!! // ctau
					if(!xIsLin) x=TMath::Power(10.,x);
					if(TMath::Abs(x-xPrevious)<1E-6) continue; //skipping repeating mass bins

					for( int j=1; j<=ExpectedNum[ifile].at(iReg)->GetNbinsY(); j++){
						y=ExpectedNum[ifile].at(iReg)->GetYaxis()->GetBinCenter(j); // ctau
						if(!yIsLin) y=TMath::Power(10.,y);
						z=ExpectedNum[ifile].at(iReg)->GetBinContent(i,j);
						RunOutput << x << " " << y << " " << z << "\n"; //
					} // loop over y
					xPrevious = x;

				} // loop over x
				delete ExpectedNum[ifile].at(iReg);
			} //massFiles loop
			RunOutput.close();
		}

		// .. clearing Dalitz Plots to avoid memory leaks
		if(!flatDalitz){
			for(Int_t mb=0; mb<ALPlist.size(); ++mb){
				delete rescaleDalitz[mb];
				for(Int_t ib=0; ib<4; ++ib)
					delete originalDalitz[mb][ib];
		}}
	}
	exit(0);
	
}  // main

 /// \fn decayTo2Body
 /// \Brief
 /// Performes a 2-body decay and checks if event passes cuts
 /// \EndBrief
Int_t decayTo2Body(ExpParameters &gen, AxionParameters &alp, Double_t &eventWeight, Particle* exotic){ 
	exotic->Decay2Body(gen.finState.at(0),gen.finState.at(1));

	std::vector<TVector3> posFVEndAll; // store final state projections
	std::vector<TVector3> posECalAll; // store final state projections
	std::vector<TLorentzVector> lorMomAll; // store final staxrte 4-mom
	TLorentzVector missingMom(0.,0.,0.,0.);
	Int_t iOnMuonDetector = 0;
	Int_t iChOnECal = 0;
	// Int_t iChOnTrigger = 0;
	Int_t iOnECal = 0;
	Int_t iOnTracker1 = 0;
	Int_t iOnTracker4 = 0;
	Double_t totalEnergyInAcceptance = 0.;
	/// loop over all objects in the tree
	std::vector<Particle*> finalStates;
	finalStates.clear();
	GetAllFinalStates(exotic, finalStates);
	if(exotic->GetName()=="hnl" && rndmGen.Integer(2))
		for(Particle* finalState:finalStates) finalState->ChargeConjugate();
	for(auto finalState : finalStates){
		if (finalState->GetName() == "nu"){
			missingMom += TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy());
			continue;
		} //don't care about neutrinos
		lorMomAll.push_back(TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy()));
		finalState->PropagateToZ(gen.GetZFVEnd());
		posFVEndAll.push_back(finalState->GetPosition());
		if(finalState->GetCharge()!=0){//check acceptances for charged particles
			finalState->PropagateToZ(gen.GetZTracker1());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker1++;
			for(Magnet magnet:gen.Get_Magnets()){
				if ((magnet.GetZ()+magnet.GetZLength())<finalState->GetInitPosition().Z()) continue;
				if (magnet.GetZ()<finalState->GetInitPosition().Z()){
					finalState->PropagateToZ(finalState->GetInitPosition().Z());
					finalState->Kick((finalState->GetInitPosition().Z()-magnet.GetZ())/magnet.GetZLength()*magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				} else {
					finalState->PropagateToZ(magnet.GetZ());
					finalState->Kick(magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				}
			}
			finalState->PropagateToZ(gen.GetZTracker4());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker4++;
			// finalState->PropagateToZ(gen.GetZTrigger());
			// if (gen.inTriggerAcceptance(finalState->GetPosition()) && finalState->GetCharge()) iChOnTrigger++;
			finalState->PropagateToZ(gen.GetZECal());
			posECalAll.push_back(finalState->GetPosition());
			if(!finalState->GetMomentum().Z()) {
				posECalAll.push_back(TVector3(TMath::Infinity(),TMath::Infinity(),gen.GetZECal()));
				continue;} // ensuring soft particle is not considered but does not affect evaluation of other final states
			if (gen.inCaloAcceptance(posFVEndAll.back(),posECalAll.back(),finalState->GetInitPosition().Z())){
				iChOnECal++;
				totalEnergyInAcceptance += finalState->GetEnergy();
			}
			if(!(finalState->GetName() == "mu+" || finalState->GetName() == "mu-")) 
				continue;
			finalState->PropagateToZ(gen.GetZMuonDetector());
			if (gen.inMuonDetectorAcceptance(posFVEndAll.back(), finalState->GetPosition(),finalState->GetInitPosition()))
				iOnMuonDetector++;
		} else{ //check only calorimeter for neutral particle
			finalState->PropagateToZ(gen.GetZECal());
			posECalAll.push_back(finalState->GetPosition());
			if (gen.inCaloAcceptance(posFVEndAll.back(),posECalAll.back(),finalState->GetInitPosition().Z())){
				iOnECal++;
				totalEnergyInAcceptance += finalState->GetEnergy();
			}
				// std::cout << " particle:" << finalState->GetName() << " inAccCalo:1" << std::endl;
			// } else std::cout << " particle:" << finalState->GetName() << " inAccCalo:0" << std::endl;
			// std::cout << " X:" << finalState->GetPosition().X() << " Y:" << finalState->GetPosition().Y() << " Z:" << finalState->GetPosition().Z() << std::endl;
		}
		finalState->EndParticle(); // assume it has been absorbed
	}
	
	Int_t condition = 0;
	Double_t sigacc = gen.GetSigacceptance();
	if(finalStates.size()==2){
		if(gen.GetDecayModeName() == "2Gamma") //for 2gamma decay check twoPhotonCondition
			condition = gen.twoPhotonCondition(iOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndAll.front(), &posECalAll.front(), &lorMomAll.front());
		else if(gen.GetDecayModeName() == "2El") //for 2el decay check twoPhotonCondition for charged
			condition = gen.twoElectronCondition(iOnTracker1, iOnTracker4, iChOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndAll.front(), &posECalAll.front(), &lorMomAll.front(), missingMom);
		else if (gen.GetDecayModeName() == "2Mu"){ //2mu
			sigacc = gen.GetSigacceptancemumu();
			condition = gen.twoMuonCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomAll.front(), missingMom);
		} else if (gen.GetDecayModeName() == "2Pi" || gen.GetDecayModeName() == "2K"){ //2hadrons
			sigacc = gen.GetSigacceptancemumu();
			condition = gen.twoHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomAll.front());
		}else if ((gen.GetDecayModeName() == "PiMu" ) || (gen.GetDecayModeName() == "KMu")){ //1hadron 1mu
			sigacc = gen.GetSigacceptancemumu();
			condition = gen.muonHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomAll.front());
		}
		else if ((gen.GetDecayModeName() == "PiEl" ) || (gen.GetDecayModeName() == "KEl")){ //1hadron 1el
			condition = gen.electronHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, exotic->GetEndPosition(), &posECalAll.front(), &lorMomAll.front());
		}
		else{
			std::cout<<"\r"<<"[Info] unrecognised final state treated as invisible decay"<<std::flush;
		}
	}else if (finalStates.size()==3) {
		if(gen.GetDecayModeName() == "RhoNu") {//nupipi
			condition = gen.twoHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomAll.front());
		} else { //nugammagamma from PiNu or EtaNu
			condition = gen.twoPhotonCondition(iOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndAll.front(), &posECalAll.front(), &lorMomAll.front());
		}

	}else if( finalStates.size()==4){
		std::vector<TVector3> posFVEndGammaSystem { {posFVEndAll.at(2), posFVEndAll.at(3)} }; // store final state projections
		std::vector<TVector3> posECalGammaSystem   { {posECalAll.at(2), posECalAll.at(3)} }; // store final state projections
		std::vector<TLorentzVector> lorMomGammaSystem { {lorMomAll.at(2), lorMomAll.at(3)} }; // store final state 4-mom
		if (gen.GetDecayModeName() == "RhoMu"  ){ // 1mu 1hadron 2 gammas 
			sigacc = gen.GetSigacceptancemumu();
			condition = gen.muonHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomAll.front())
						&& gen.twoPhotonCondition(iOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndGammaSystem.front(), &posECalGammaSystem.front(), &lorMomGammaSystem.front());
		} else if (gen.GetDecayModeName() == "RhoEl" ){
			condition = gen.electronHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, exotic->GetEndPosition(), &posECalAll.front(), &lorMomAll.front())
						&& gen.twoPhotonCondition(iOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndGammaSystem.front(), &posECalGammaSystem.front(), &lorMomGammaSystem.front());
		}
	}
	finalStates.clear();
	return condition;
} // decayTo2Body

 /// \fn decayTo3Body
 /// \Brief
 /// Performes a 3-body decay and checks if event passes cuts
 /// \EndBrief
Int_t decayTo3Body(ExpParameters &gen, AxionParameters &alp, Double_t &eventWeight, Particle* exotic, TH2D* rescaleDal, TH2D** originalDal){
	Double_t m12 = -1.;
	Double_t m23 = -1.;
	eventWeight *= exotic->Decay3Body(gen.finState.at(0),gen.finState.at(1),gen.finState.at(2),m12,m23,rescaleDal);

	std::vector<TVector3> posFVEndPhotons, posFVEndRest; // store photon projections
	std::vector<TVector3> posECalPhotons, posECalRest; // store photon projections
	std::vector<TLorentzVector> lorMomPhotons, lorMomRest; // store photon 4-mom
	TLorentzVector missingMom(0.,0.,0.,0.);
	Int_t iOnMuonDetector = 0;
	Int_t iChOnECal = 0;
	Int_t iOnECal = 0;
	Int_t iOnTracker1 = 0;
	Int_t iOnTracker4 = 0;
	// Int_t iChOnTrigger = 0;
	Double_t totalEnergyInAcceptance = 0.;
	Int_t nCharged = 0;

	/// loop over all objects in the tree
	std::vector<Particle*> finalStates;
	finalStates.clear();
	GetAllFinalStates(exotic, finalStates);
	if(exotic->GetName()=="hnl" && rndmGen.Integer(2))
		for(Particle* finalState:finalStates) finalState->ChargeConjugate();
	for(auto& finalState : finalStates){
		if (finalState->GetName() == "nu"){
			missingMom += TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy());
			continue;
		} 
		finalState->PropagateToZ(gen.GetZFVEnd());
		if(finalState->GetName() == "gamma"){
			lorMomPhotons.push_back(TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy()));
			posFVEndPhotons.push_back(finalState->GetPosition());
		}
		else{
			lorMomRest.push_back(TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy()));
			posFVEndRest.push_back(finalState->GetPosition());
		}
		if(finalState->GetCharge()!=0){//check acceptances for charged particles
			nCharged++;
			finalState->PropagateToZ(gen.GetZTracker1());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker1++;
			for(Magnet magnet:gen.Get_Magnets()){
				if ((magnet.GetZ()+magnet.GetZLength())<finalState->GetInitPosition().Z()) continue;
				if (magnet.GetZ()<finalState->GetInitPosition().Z()){
					finalState->PropagateToZ(finalState->GetInitPosition().Z());
					finalState->Kick((finalState->GetInitPosition().Z()-magnet.GetZ())/magnet.GetZLength()*magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				} else {
					finalState->PropagateToZ(magnet.GetZ());
					finalState->Kick(magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				}
			}
			// finalState->PropagateToZ(gen.GetZTrigger());
			// if (gen.inTriggerAcceptance(finalState->GetPosition()) && finalState->GetCharge()) iChOnTrigger++;
			if(!finalState->GetMomentum().Z()) {
				posECalRest.push_back(TVector3(TMath::Infinity(),TMath::Infinity(),gen.GetZECal()));
				continue;} // ensuring soft particle is not considered but does not affect evaluation of other final states
			finalState->PropagateToZ(gen.GetZTracker4());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker4++;
			finalState->PropagateToZ(gen.GetZECal());
			posECalRest.push_back(finalState->GetPosition());
			if (gen.inCaloAcceptance(posFVEndRest.back(),posECalRest.back(),finalState->GetInitPosition().Z())){
				iChOnECal++;
				totalEnergyInAcceptance += finalState->GetEnergy();
			}
			if(!(finalState->GetName() == "mu+" || finalState->GetName() == "mu-")) 
				continue;
			finalState->PropagateToZ(gen.GetZMuonDetector());
			if (gen.inMuonDetectorAcceptance(posFVEndRest.back(), finalState->GetPosition(),finalState->GetInitPosition()))
				iOnMuonDetector++;
		} else{ //check only calorimeter for neutral particle
			finalState->PropagateToZ(gen.GetZECal());
			if(finalState->GetName() == "gamma"){
				posECalPhotons.push_back(finalState->GetPosition());
				if (gen.inCaloAcceptance(posFVEndPhotons.back(),posECalPhotons.back(),finalState->GetInitPosition().Z())){
					iOnECal++;
					totalEnergyInAcceptance += finalState->GetEnergy();
				}
			}
			else{
				posECalRest.push_back(finalState->GetPosition());
				if (gen.inCaloAcceptance(posFVEndRest.back(),posECalRest.back(),finalState->GetInitPosition().Z())){
					iOnECal++;
					totalEnergyInAcceptance += finalState->GetEnergy();
				}
			}
		}
		// std::cout << "Particle:" << finalState->GetName() << std::endl;
		// std::cout << "PosX:" << finalState->GetPosition().X() << std::endl;
		// std::cout << "PosY:" << finalState->GetPosition().Y() << std::endl;
		// std::cout << "PosZ:" << finalState->GetPosition().Z() << std::endl;
		// std::cout << "InCaloAcc:" << gen.inCaloAcceptance(posFVEndRest.back(),posECalRest.back(),finalState->GetInitPosition().Z()) << std::endl;
		
		finalState->EndParticle(); // assume it has been absorbed if not a neutrino
	}

	Bool_t condition = true;
	Double_t sigacc = gen.GetSigacceptance();;

	if ((gen.GetDecayModeName() == "NuMuMu" ) ){ 
		sigacc = gen.GetSigacceptancemumu();
		condition = iOnMuonDetector == 2 && gen.twoMuonCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomRest.front(), missingMom);
	}else if ((gen.GetDecayModeName() == "NuElEl" )  ){ 
		condition = gen.twoElectronCondition(iOnTracker1, iOnTracker4, iChOnECal, totalEnergyInAcceptance, exotic->GetEndPosition(), &posFVEndRest.front(), &posECalRest.front(), &lorMomRest.front(), missingMom);
	}else if ((gen.GetDecayModeName() == "NuElMu" )  ){ 
		sigacc = gen.GetSigacceptancemumu();
		condition = gen.electronMuonCondition(totalEnergyInAcceptance, iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &posECalRest.front(), &lorMomRest.front(), missingMom);
	} else {
		if(nCharged == 2){
			sigacc = gen.GetSigacceptancemumu();
			condition = gen.twoHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomRest.front());
		}

		if(lorMomPhotons.size()) 
			condition = condition && gen.multiplePhotonCondition(lorMomPhotons.size(), totalEnergyInAcceptance, exotic->GetEndPosition(),
															 &posFVEndPhotons.front(), &posECalPhotons.front(), &lorMomPhotons.front());
	}
	if(verbose){
		if(condition) {
			originalDal[1]->Fill(m12,m23,eventWeight); //should be adapted
			originalDal[3]->Fill(m12,m23,eventWeight); //should be adapted
		}
		originalDal[0]->Fill(m12,m23,eventWeight); //should be adapted
		originalDal[2]->Fill(m12,m23,eventWeight); //should be adapted
	}
	finalStates.clear();
	return condition;

} //decayTo3Body

 /// \fn decayTo4Body
 /// \Brief
 /// Performes a 4-body decay and checks if event passes cuts
 /// \EndBrief
Int_t decayTo4Body(ExpParameters &gen, AxionParameters &alp, Double_t &eventWeight, Particle* exotic){
	eventWeight *= exotic->Decay4Body(gen.finState.at(0),gen.finState.at(1),gen.finState.at(2),gen.finState.at(3));

	std::vector<TVector3> posFVEndPhotons, posFVEndRest; // store photon projections
	std::vector<TVector3> posECalPhotons, posECalRest; // store photon projections
	std::vector<TLorentzVector> lorMomPhotons, lorMomRest; // store photon 4-mom
	TLorentzVector missingMom(0.,0.,0.,0.);
	Int_t iOnMuonDetector = 0;
	Int_t iChOnECal = 0;
	Int_t iOnECal = 0;
	Int_t iOnTracker1 = 0;
	Int_t iOnTracker4 = 0;
	Double_t totalEnergyInAcceptance = 0.;
	Int_t nCharged = 0;

	/// loop over all objects in the tree
	std::vector<Particle*> finalStates;
	finalStates.clear();
	GetAllFinalStates(exotic, finalStates);
	if(exotic->GetName()=="hnl" && rndmGen.Integer(2))
		for(Particle* finalState:finalStates) finalState->ChargeConjugate();
	for(auto finalState : finalStates){
		if (finalState->GetName() == "nu"){
			missingMom += TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy());
			continue;
		} 
		finalState->PropagateToZ(gen.GetZFVEnd());
		if(finalState->GetName() == "gamma"){
			lorMomPhotons.push_back(TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy()));
			posFVEndPhotons.push_back(finalState->GetPosition());
		}
		else{
			lorMomRest.push_back(TLorentzVector(finalState->GetMomentum().x(),finalState->GetMomentum().y(),finalState->GetMomentum().z(),finalState->GetEnergy()));
			posFVEndRest.push_back(finalState->GetPosition());
		}
		if(finalState->GetCharge()!=0){//check acceptances for charged particles
			nCharged++;
			finalState->PropagateToZ(gen.GetZTracker1());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker1++;
			for(Magnet magnet:gen.Get_Magnets()){
				if ((magnet.GetZ()+magnet.GetZLength())<finalState->GetInitPosition().Z()) continue;
				if (magnet.GetZ()<finalState->GetInitPosition().Z()){
					finalState->PropagateToZ(finalState->GetInitPosition().Z());
					finalState->Kick((finalState->GetInitPosition().Z()-magnet.GetZ())/magnet.GetZLength()*magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				} else {
					finalState->PropagateToZ(magnet.GetZ());
					finalState->Kick(magnet.GetMagKick(), magnet.GetMagFieldInfo(), magnet.GetIsToroidal());
				}
			}
			if(!finalState->GetMomentum().Z()) {
				posECalRest.push_back(TVector3(TMath::Infinity(),TMath::Infinity(),gen.GetZECal()));
				continue;} // ensuring soft particle is not considered but does not affect evaluation of other final states
			finalState->PropagateToZ(gen.GetZTracker4());
			if (gen.inSpectrometerAcceptance(finalState->GetPosition())) iOnTracker4++;
			finalState->PropagateToZ(gen.GetZECal());
			posECalRest.push_back(finalState->GetPosition());
			if (gen.inCaloAcceptance(posFVEndRest.back(),posECalRest.back(),finalState->GetInitPosition().Z())){
				iChOnECal++;
				totalEnergyInAcceptance += finalState->GetEnergy();
			}
			finalState->PropagateToZ(gen.GetZMuonDetector());
			if (gen.inMuonDetectorAcceptance(posFVEndRest.back(), finalState->GetPosition(),finalState->GetInitPosition()))
				iOnMuonDetector++;
		} else{ //check only calorimeter for neutral particle
			finalState->PropagateToZ(gen.GetZECal());
			if(finalState->GetName() == "gamma"){
				posECalPhotons.push_back(finalState->GetPosition());
				if (gen.inCaloAcceptance(posFVEndPhotons.back(),posECalPhotons.back(),finalState->GetInitPosition().Z())){
					iOnECal++;
					totalEnergyInAcceptance += finalState->GetEnergy();
				}
			}
			else{
				posECalRest.push_back(finalState->GetPosition());
				if (gen.inCaloAcceptance(posFVEndRest.back(),posECalRest.back(),finalState->GetInitPosition().Z())){
					iOnECal++;
					totalEnergyInAcceptance += finalState->GetEnergy();
				}
			}
		}
		finalState->EndParticle(); // assume it has been absorbed if not a neutrino
	}

	Bool_t condition = true;
	Double_t sigacc = gen.GetSigacceptance();;

	if(nCharged == 2){
		sigacc = gen.GetSigacceptancemumu();
		condition = gen.twoHadronCondition(iOnTracker1, iOnTracker4, iChOnECal, iOnMuonDetector, exotic->GetEndPosition(), &lorMomRest.front());
	}

	if(lorMomPhotons.size()) 
		condition = condition && gen.multiplePhotonCondition(lorMomPhotons.size(), totalEnergyInAcceptance, exotic->GetEndPosition(),
															 &posFVEndPhotons.front(), &posECalPhotons.front(), &lorMomPhotons.front());

	finalStates.clear();
	return condition;

} //decayTo4Body

 /// \fn GetAllFinalStates
 /// \Brief
 /// Iterates through the decay tree and returns a pointer to all stable final states
 /// \EndBrief
void GetAllFinalStates(Particle* particle, std::vector<Particle*>& finalStates){
	//recursive function to dig out all the stable final states from the tree of decays
	// auto nDaughters = particle->GetChildren().size();
	auto nDaughters = particle->GetNChildren();
	if(nDaughters!=0){
		for(UInt_t daughter = 0; daughter < nDaughters; daughter++)
			// GetAllFinalStates(particle->GetChildren().at(daughter), finalStates);
			GetAllFinalStates(particle->GetDaughter(daughter), finalStates);
	} else {
		if(!particle->GetUnstable())
			finalStates.push_back(particle);
	}
}

 /// \fn FillCDAzTAX
 /// \Brief
 /// Fills a CDA vs Z histo for given event
 /// \EndBrief
void FillCDAzTAX(std::vector<Particle*> finalStates, TH2D* CDAzTAXs, const TVector3 TAX, const Double_t eventWeight=1.){
	TVector3 decayVertex{finalStates[0]->GetInitPosition()};
	TVector3 motherMom(0.,0.,0.), protonPath(0,0,1);
	for(auto finalState : finalStates) 
		if(finalState->GetName() != "nu")
			motherMom += finalState->GetMomentum();
	TVector3 connectingPath {motherMom.Cross(protonPath)};

	if(connectingPath==TVector3(0,0,0)) return;
	Double_t CDA 	{	TMath::Abs(connectingPath.Dot(decayVertex-TAX) / connectingPath.Mag())	};
	Double_t ZTarget	{	decayVertex.Z()+(protonPath.Cross(connectingPath)).Dot(TAX-decayVertex) / connectingPath.Mag2() * motherMom.Z() };

	CDAzTAXs->Fill(ZTarget,CDA*1000,eventWeight);
}