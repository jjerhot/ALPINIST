#ifndef EXPPARAMETERS_H
#define EXPPARAMETERS_H

#include "Particle.h"

#include <iostream>
#include <array>
#include <fstream>
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TROOT.h>
#include "Rtypes.h"

class Magnet;
class ExpParameters{

public:

	ExpParameters(Int_t nBinsX, Int_t nBinsY, Bool_t xIsLin, Bool_t yIsLin, Int_t yVar);
	ExpParameters(Int_t ExoticType_, Int_t Experiment_, TString ProductionMode_, Int_t DecayMode_,
					Int_t nBinsX, Int_t nBinsY, Bool_t xIsLin, Bool_t yIsLin, Int_t yVar);
	virtual ~ExpParameters();

//	ClassDef(ExpParameters,1);

	//parameters input
	Int_t GetNMassFiles() const             		{return nMassFiles; }
	std::array<TString,2> GetMinMassFile() const    {return minMassFile; }
	std::array<TString,2> GetMaxMassFile() const    {return maxMassFile; }
	
	Int_t GetValuesEnergyInFile() const   	{return energyValuesInFile; }
	Double_t GetMinEnergyInFile() const   	{return energyMinInFile; }
	Double_t GetMaxEnergyInFile() const   	{return energyMaxInFile; }
	Double_t GetWdEnergy() const			{return energyWid; }
	void SetValuesEnergyInFile (Int_t n)    { energyValuesInFile = n; }
	void SetMinEnergyInFile (Double_t min)  { energyMinInFile = min; }
	void SetMaxEnergyInFile (Double_t max)  { energyMaxInFile = max; }
	void SetWdEnergy (Double_t wd)     		{ energyWid = wd; }

	Int_t GetValuesThetaInFile() const    	{return thetaValuesInFile; }
	Double_t GetMinThetaInFile() const    	{return thetaMinInFile; }
	Double_t GetMaxThetaInFile() const    	{return thetaMaxInFile; }
	Double_t GetWdTheta() const             {return thetaWid; }
	void SetValuesThetaInFile (Int_t n)    	{ thetaValuesInFile = n; }
	void SetMinThetaInFile (Double_t min)  	{ thetaMinInFile = min; }
	void SetMaxThetaInFile (Double_t max)  	{ thetaMaxInFile = max; }
	void SetWdTheta (Double_t wd)     		{ thetaWid = wd; }

	std::array<Int_t,2> GetValuesMassAInFile() const	{ return mAValuesInFile; }
	std::array<Double_t,2> GetMinMassAInFile() const    { return mAMinInFile; }
	std::array<Double_t,2> GetMaxMassAInFile() const    { return mAMaxInFile; }
	std::array<Double_t,2> GetWdMassA() const      		{ return mAWid; }
	std::array<Double_t,2> GetLogWdMassA() const 		{ return logMAWid; }
	void SetValuesMassAInFile (Int_t iFile,Int_t n)     { mAValuesInFile.at(iFile) = n; }
	void SetMinMassAInFile (Int_t iFile,Double_t min)   { mAMinInFile.at(iFile) = min; }
	void SetMaxMassAInFile (Int_t iFile,Double_t max)   { mAMaxInFile.at(iFile) = max; }
	void SetWdMassA (Int_t iFile,Double_t wd)     		{ mAWid.at(iFile) = wd; }
	void SetLogWdMassA (Int_t iFile, Double_t logwd)	{ logMAWid.at(iFile) = logwd; }

	//parameters output
	Int_t GetNWidths() const              	{return nWidths; }
	Double_t GetMinWidth() const          	{return widthMin; }
	Double_t GetMaxWidth() const          	{return widthMax; }
	Double_t GetWdWidth() const           	{return wdWidth; }

	std::array<Int_t,2> GetNMassX() const               {return nMX; }
	std::array<Double_t,2> GetMinMassX() const        	{return massXMin; }
	std::array<Double_t,2> GetMaxMassX() const        	{return massXMax; }
	std::array<Double_t,2> GetWdMassX() const  			{return massXWid; }	


	Int_t GetExoNum() const					{return exoNum; }
	// TString GetExoName() const 				{return exoName; }
	Int_t GetExpNum() const                	{return expNum; }
	TString GetExpName() const              {return expName; }
	TString GetExpLabel() const             {return expLabel; }
	Int_t GetBeamEnergy() const             {return beamEnergy; }
	Double_t GetZECal() const                {return ZECal; }
	Double_t GetZFVEnd() const              {return ZFVEnd; }
	Double_t GetZFVIn() const               {return ZFVIn; }
	Double_t GetZTarget() const                {return ZTarget; }
	Double_t GetX0() const                  {return X0; }
	Double_t GetY0() const                  {return Y0; }
	Double_t GetZ0() const                  {return Z0; }
	Double_t GetPhiBeamEuler() const        {return PhiBeamEuler;}
	Double_t GetThetaBeamEuler() const      {return ThetaBeamEuler;}
	Double_t GetPsiBeamEuler() const        {return PsiBeamEuler;}
	Double_t GetSigacceptance() const       {return Sigacceptance; } // used for gamma gamma
	Double_t GetSigacceptancemumu() const   {return Sigacceptancemumu; } // can be different from the above
	Double_t GetPOT() const                 {return POT; }
	Double_t GetNormCrossSec() const        {return normCrossSec; }
	Double_t GetZTracker1() const             {return ZTracker1; }
	Double_t GetZTracker4() const             {return ZTracker4; }
	Double_t GetAcceptanceHole() const      {return acceptanceHole; }
	Double_t GetAcceptanceSide() const      {return acceptanceSide; }
	Double_t GetMuonDetectorSize() const            {return MuonDetectorSize; }
	Double_t GetZTrigger() const               {return ZTrigger; }
	Double_t GetZMuonDetector() const               {return ZMuonDetector; }
	std::vector<Magnet> Get_Magnets() const	{return Magnets;}

	Double_t GetZDist() const               {return ZDist; }
	Double_t GetBeamDecayLength() const     {return BeamDecayLength; }

	Double_t GetMinMassA() const            {return mAMin; }
	Double_t GetMaxMassA() const            {return mAMax; } //when there is a kinematic restriction on ALP mass

	void SetMinMassA     (Double_t	ma)     { mAMin = ma; }
	void SetMaxMassA     (Double_t	ma)     { mAMax = ma; }

	Int_t GetDecayMode() const {return decayMode; }
	TString GetDecayModeName() const {return decayModeName;}
	Bool_t GetDecayModeOpen() const {return decayModeOpen; }

	std::array<Double_t,3> GetMassFinState() const   	{return mFinState; } //final state particles masses
	std::array<Int_t,3> GetChargeFinState() const     	{return chargeFinState; } //final state particles charges 

	Int_t GetNRegions() const 				{return nRegions;}

	//detector geometry methods
	Bool_t inCaloAcceptance(TVector3 posFVEnd, TVector3 posECal, Double_t zDecay); //detector calorimeter geometries
	Bool_t inTriggerAcceptance(TVector3 posTrigger); //detector Trigger geometries
	Bool_t inCaloOuterEdgeAcceptance(TVector3); //detector calorimeter geometries
	Bool_t inSpectrometerAcceptance(TVector3 posTracker); //detector spectrometer geometries
	Bool_t inMuonDetectorAcceptance(TVector3 posFVEnd, TVector3 posMuDet, TVector3 posDecay); //detector muon veto geometries

    //experimental conditions methods
    Int_t   twoPhotonCondition(Int_t iOnECal, Double_t totalEnergyAcceptance, TVector3 posDecay, TVector3 posFVEnd[2], TVector3 posECal[2], TLorentzVector pgA[2]);
    Bool_t  multiplePhotonCondition(Int_t nGammas, Double_t totalEnergyAcceptance, TVector3 posDecay, TVector3 posFVEnd[6], TVector3 posECal[6], TLorentzVector pgA[6]);	
    Bool_t  twoMuonCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2], TLorentzVector pMiss);    
	Bool_t  twoHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2]);
	Bool_t  muonHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector, TVector3 posDecay, TLorentzVector pCA[2]);
	Bool_t  electronMuonCondition(Double_t totalEnergyAcceptance, Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal, Int_t iOnMuonDetector,  TVector3 posDecay, TVector3 posECal[2], TLorentzVector pCA[2], TLorentzVector pMiss);
	Bool_t	electronHadronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iChOnECal,  TVector3 posDecay, TVector3 posECal[2], TLorentzVector pCA[2]);
	Bool_t	twoElectronCondition(Int_t iOnTracker1, Int_t iOnTracker4, Int_t iOnECal, Double_t totalEnergyAcceptance, TVector3 posDecay, TVector3 posFVEnd[2], TVector3 posECal[2], TLorentzVector pCA[2], TLorentzVector pMiss);
	std::vector<SMParticleProperty*> finState; //final state particles masses

	// Int_t ncTau = 0;
	// Double_t cTauMin = 0;
	// Double_t cTauMax = 0;
	// Double_t cTauWidth = 0;
	
	Int_t nY = 0;
	Double_t minY = 0;
	Double_t maxY = 0;
	Double_t wdY = 0;
	
private:

	//Declarations
	Int_t decayMode = -1;
	TString decayModeName = "";
	Bool_t decayModeOpen = 0;

	Double_t energyMinInFile = 0;
	Double_t energyMaxInFile = 0;
	Int_t energyValuesInFile = 0;

	Double_t thetaMinInFile = 0;
	Double_t thetaMaxInFile = 0;
	Int_t thetaValuesInFile = 0;

	Int_t nWidths = 0;
	Double_t widthMin = 0;
	Double_t widthMax = 0;
	Double_t wdWidth = 0;

	std::array<Int_t,2> mAValuesInFile = {0,0};

	Int_t nMassFiles = 0;
	std::array<Int_t,2> nMX = {0,0};
	std::array<Double_t,2> massXMin, mAMinInFile = {0,0};
	std::array<Double_t,2> massXMax, mAMaxInFile = {0,0};
	std::array<TString,2> minMassFile = {"",""};
	std::array<TString,2> maxMassFile = {"",""};

	Double_t energyWid = 0;
	Double_t thetaWid = 0;

	std::array<Double_t,2> massXWid = {0,0};
	std::array<Double_t,2> logMassXWid = {0,0};
	std::array<Double_t,2> mAWid = {0,0};
	std::array<Double_t,2> logMAWid = {0,0};

	Int_t exoNum = 0;
	Int_t expNum = 0;
	// TString exoName = "";
	TString expName = "";
	TString expLabel = "";
	Int_t beamEnergy = 0;
	std::array<Double_t,3> target;
	Double_t ZECal = 0;
	Double_t ZFVEnd = 0;
	Double_t ZFVIn = 0;
	Double_t ZTarget = 0;
	Double_t X0 = 0;
	Double_t Y0 = 0;
	Double_t Z0 = 0;
	Double_t PhiBeamEuler = 0;   // euler phi angle of beam: rotation of beam axis about the z-axis (i.e., tranforms {x,y,z} in {x',y',z})
	Double_t ThetaBeamEuler = 0; // euler theta angle of beam: rotation about the x'-axis (i.e., tranforms {x',y',z} in {x',y'',z'})
	Double_t PsiBeamEuler = 0; // euler psi angle of beam: rotation about the z'-axis (i.e., tranforms {x',y'',z'} in {x'',y''',z'})
	Double_t Sigacceptance = 0; // used for gamma gamma
	Double_t Sigacceptancemumu = 0; // can be different from the above
	Double_t POT = 0;
	Double_t normCrossSec = 0;
	Double_t ZTracker1 = 0;
	Double_t ZTracker4 = 0;
	Double_t acceptanceHole = 0;
	Double_t acceptanceSide = 0;
	Double_t MuonDetectorSize;
	Double_t ZTrigger = 0;
	Double_t ZMuonDetector = 0;
	std::vector<Magnet> Magnets = {};

	Double_t ZDist;
	Double_t BeamDecayLength;

	Double_t mAMin, mAMax; //when there is a kinematic restriction on ALP mass

	std::array<Double_t,3> mFinState; //final state particles masses
	std::array<Int_t,3> chargeFinState; //final state particles charges

	Int_t nRegions;

};

 /// \class Magnet
 /// \Brief
 /// Magnet object with field strengh and length
 /// \EndBrief
class Magnet{
	
private:
	Bool_t   _IsTorroidal; // refers to a torroidal magnet with a coil in the center
	Double_t _FieldInfo; //rad
	Double_t _ZMagnet; //m
	Double_t _ZLengthMagnet; //m
	Double_t _MagKick; //GeV
public:
	Magnet() : _IsTorroidal(0), _FieldInfo(TMath::PiOver2()), _ZMagnet(0.), _MagKick(0.){}
	Magnet(Double_t ZMagnet, Double_t MagnetZLength, Double_t MagnetFieldStrength): _IsTorroidal(0), _FieldInfo(TMath::PiOver2()){
		_ZMagnet		= ZMagnet;
		_ZLengthMagnet  = MagnetZLength;
		_MagKick		= c*1E-9*MagnetZLength*MagnetFieldStrength;
	}
	Magnet(Double_t ZMagnet, Double_t MagnetZLength, Double_t MagnetFieldStrength, Double_t MagnetFieldPhi) : Magnet(ZMagnet, MagnetZLength, MagnetFieldStrength) {
		_FieldInfo	= MagnetFieldPhi;
		_IsTorroidal= 0;
	}
	Magnet(Double_t ZMagnet, Double_t MagnetZLength, Double_t MagnetFieldStrength, Bool_t IsToroidal, Double_t DetectorRadius) : Magnet(ZMagnet, MagnetZLength, MagnetFieldStrength) {
		_IsTorroidal	= IsToroidal;
		_FieldInfo		= DetectorRadius;
	}
	virtual ~Magnet(){}
	Double_t GetZ()				const {return _ZMagnet;}
	Double_t GetZLength()		const {return _ZLengthMagnet;}
	Double_t GetMagKick()		const {return _MagKick;}
	Double_t GetMagFieldInfo()	const {return _FieldInfo;}
	Double_t GetMagRadius()		const {return _FieldInfo;}
	Double_t GetIsToroidal()	const {return _IsTorroidal;}
};

#endif
