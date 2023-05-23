#ifndef EXPPARAMETERS_H
#define EXPPARAMETERS_H

#include "DecayMCGlobal.h"

#include <iostream>
#include <array>
#include <fstream>
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TROOT.h>
#include "Rtypes.h"

class ExpParameters{

public:

	ExpParameters();
	ExpParameters(Int_t Experiment_, Int_t ProductionMode_, Int_t DecayMode_);
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

	std::array<Int_t,2> GetValuesMassAInFile() const	{return mAValuesInFile; }
    std::array<Double_t,2> GetMinMassAInFile() const    {return mAMinInFile; }
    std::array<Double_t,2> GetMaxMassAInFile() const    {return mAMaxInFile; }
	std::array<Double_t,2> GetWdMassA() const      		{return mAWid; }
	void SetValuesMassAInFile (Int_t iFile,Int_t n)     { mAValuesInFile.at(iFile) = n; }
	void SetMinMassAInFile (Int_t iFile,Double_t min)   { mAMinInFile.at(iFile) = min; }
	void SetMaxMassAInFile (Int_t iFile,Double_t max)   { mAMaxInFile.at(iFile) = max; }
	void SetWdMassA (Int_t iFile,Double_t wd)     		{ mAWid.at(iFile) = wd; }

    //parameters output
	Int_t GetNWidths() const              	{return nWidths; }
	Double_t GetMinWidth() const          	{return widthMin; }
	Double_t GetMaxWidth() const          	{return widthMax; }
	Double_t GetWdWidth() const           	{return wdWidth; }

	std::array<Int_t,2> GetNMassX() const               {return nMX; }
	std::array<Double_t,2> GetMinMassX() const        	{return massXMin; }
	std::array<Double_t,2> GetMaxMassX() const        	{return massXMax; }
	std::array<Double_t,2> GetWdMassX() const  			{return massXWid; }

	Int_t GetExpNum() const                	{return expNum; }
	TString GetExpName() const              {return expName; }
	TString GetExpLabel() const              {return expLabel; }
	Int_t GetBeamEnergy() const             {return beamEnergy; }
	Double_t GetATarget() const             {return ATarget; }
	Double_t GetZTarget() const             {return ZTarget; }
	Double_t GetZLKR() const                {return ZLKR; }
	Double_t GetZFVEnd() const              {return ZFVEnd; }
	Double_t GetZFVIn() const               {return ZFVIn; }
	Double_t GetZTAX() const                {return ZTAX; }
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
	Double_t GetZMagnet() const             {return ZMagnet; }
	Double_t GetMagnetZLength() const       {return MagnetZLength; }
	Double_t GetMagnetFieldStrength() const {return MagnetFieldStrength; }
	Double_t GetMagKick() const             {return MagKick; }
	Double_t GetZStraw1() const             {return zStraw1; }
	Double_t GetZStraw4() const             {return zStraw4; }
	Double_t GetAcceptanceHole() const      {return acceptanceHole; }
	Double_t GetAcceptanceSide() const      {return acceptanceSide; }
	Double_t GetMUV3Size() const            {return MUV3Size; }
	Double_t GetZMUV3() const               {return ZMUV3; }

	Double_t GetZDist() const               {return ZDist; }
	Double_t GetBeamDecayLength() const     {return BeamDecayLength; }

	Double_t GetMinMassA() const            {return mAMin; }
    Double_t GetMaxMassA() const            {return mAMax; } //when there is a kinematic restriction on ALP mass

	void SetMinMassA     (Double_t	ma)     { mAMin = ma; }
	void SetMaxMassA     (Double_t	ma)     { mAMax = ma; }

	std::array<Double_t,3> GetMassFinState() const   	{return mFinState; } //final state particles masses
	std::array<Int_t,3> GetChargeFinState() const     	{return chargeFinState; } //final state particles charges    

	Int_t GetNRegions() const 				{return nRegions;}

    //detector geometry methods
    Bool_t  inCaloAcceptance(TVector3 posFVEnd, TVector3 posLKr, Double_t zDecay); //detector calorimeter geometries
	Bool_t  inCaloOuterEdgeAcceptance(TVector3); //detector calorimeter geometries
    Bool_t  inSpectrometerAcceptance(TVector3 posStraw); //detector spectrometer geometries
    Bool_t  inMuonVetoAcceptance(TVector3 posFVEnd, TVector3 posMUV, Double_t zDecay); //detector muon veto geometries

    //experimental conditions methods
    Int_t  twoPhotonCondition(Int_t iOnLKr, Double_t totalEnergyAcceptance, Double_t zDecay, TVector3 posFVEnd[2], TVector3 posLKr[2], TLorentzVector pgA[2]);
    Bool_t  multiplePhotonCondition(Int_t nGammas, Double_t totalEnergyAcceptance, Double_t zDecay, TVector3 posFVEnd[6], TVector3 posLKr[6], TLorentzVector pgA[6]);	
    Bool_t  twoMuonCondition(Int_t iOnStraw1, Int_t iOnStraw4, Int_t iOnLKr, Int_t iOnMUV3, Double_t zDecay, TLorentzVector pCA[2]);    
	Bool_t  twoHadronCondition(Int_t iOnStraw1, Int_t iOnStraw4, Int_t iOnLKr, Int_t iOnMUV3, Double_t zDecay, TLorentzVector pCA[2]);

private:

	//Declarations
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
	std::array<Double_t,2> mAWid = {0,0};

	Int_t expNum = 0;
	TString expName = "";
	TString expLabel = "";
	Int_t beamEnergy = 0;
	Double_t ATarget = 0;
	Double_t ZTarget = 0;
	Double_t ZLKR = 0;
	Double_t ZFVEnd = 0;
	Double_t ZFVIn = 0;
	Double_t ZTAX = 0;
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
	Double_t ZMagnet = 0;
	Double_t MagnetZLength = 0;
	Double_t MagnetFieldStrength = 0;
	Double_t MagKick = 0;
	Double_t zStraw1 = 0;
	Double_t zStraw4 = 0;
	Double_t acceptanceHole = 0;
	Double_t acceptanceSide = 0;
	Double_t MUV3Size;
	Double_t ZMUV3 = 0;


	Double_t ZDist;
	Double_t BeamDecayLength;

	Double_t mAMin, mAMax; //when there is a kinematic restriction on ALP mass

	std::array<Double_t,3> mFinState; //final state particles masses
	std::array<Int_t,3> chargeFinState; //final state particles charges

	TRandom3 *rndm;

	Int_t nRegions;

};

#endif
