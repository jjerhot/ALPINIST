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

    //parameters
	Double_t GetMinEnergyInFile() const   	{return energyMinInFile; }
	Double_t GetMaxEnergyInFile() const   	{return energyMaxInFile; }
	Int_t GetValuesEnergyInFile() const   	{return energyValuesInFile; }

	Double_t GetMinThetaInFile() const    	{return thetaMinInFile; }
	Double_t GetMaxThetaInFile() const    	{return thetaMaxInFile; }
	Int_t GetValuesThetaInFile() const    	{return thetaValuesInFile; }

	Int_t GetNWidths() const              	{return nWidths; }
	Double_t GetMinWidth() const          	{return widthMin; }
	Double_t GetMaxWidth() const          	{return widthMax; }
	Double_t GetWdWidth() const           	{return wdWidth; }

	std::array<Int_t,2> GetValuesMassAInFile() const	{return mAValuesInFile; }

	Int_t GetNMassFiles() const                       	{return nMassFiles; }
	std::array<Int_t,2> GetNMassX() const               {return nMX; }
	std::array<Double_t,2> GetMinMassX() const        	{return massXMin; }
    std::array<Double_t,2> GetMinMassAInFile() const    {return mAMinInFile; }
	std::array<Double_t,2> GetMaxMassX() const        	{return massXMax; }
    std::array<Double_t,2> GetMaxMassAInFile() const    {return mAMaxInFile; }
	std::array<TString,2> GetMinMassFile() const      	{return minMassFile; }
	std::array<TString,2> GetMaxMassFile() const      	{return maxMassFile; }

	Double_t GetWdEnergy() const         	{return energyWid; }
	Double_t GetWdTheta() const             {return thetaWid; }

	std::array<Double_t,2> GetWdMassX() const  			{return massXWid; }
	std::array<Double_t,2> GetWdMassA() const      		{return mAWid; }

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
	Double_t energyMinInFile;
	Double_t energyMaxInFile;
	Int_t energyValuesInFile;

	Double_t thetaMinInFile;
	Double_t thetaMaxInFile;
	Int_t thetaValuesInFile;

	Int_t nWidths;
	Double_t widthMin;
	Double_t widthMax;
	Double_t wdWidth;

	std::array<Int_t,2> mAValuesInFile;

	Int_t nMassFiles;
	std::array<Int_t,2> nMX;
	std::array<Double_t,2> massXMin, mAMinInFile;
	std::array<Double_t,2> massXMax, mAMaxInFile;
	std::array<TString,2> minMassFile;
	std::array<TString,2> maxMassFile;

	Double_t energyWid;
	Double_t thetaWid;

	std::array<Double_t,2> massXWid;
	std::array<Double_t,2> mAWid;

	Int_t expNum;
	TString expName;
	TString expLabel;
	Int_t beamEnergy;
	Double_t ATarget;
	Double_t ZTarget;
	Double_t ZLKR;
	Double_t ZFVEnd;
	Double_t ZFVIn;
	Double_t ZTAX;
	Double_t X0;
	Double_t Y0;
	Double_t Z0;
	Double_t PhiBeamEuler;   // euler phi angle of beam: rotation of beam axis about the z-axis (i.e., tranforms {x,y,z} in {x',y',z})
	Double_t ThetaBeamEuler; // euler theta angle of beam: rotation about the x'-axis (i.e., tranforms {x',y',z} in {x',y'',z'})
	Double_t PsiBeamEuler; // euler psi angle of beam: rotation about the z'-axis (i.e., tranforms {x',y'',z'} in {x'',y''',z'})
	Double_t Sigacceptance; // used for gamma gamma
	Double_t Sigacceptancemumu; // can be different from the above
	Double_t POT;
	Double_t normCrossSec;
	Double_t ZMagnet = 0;
	Double_t MagnetZLength = 0;
	Double_t MagnetFieldStrength = 0;
	Double_t MagKick = 0;
	Double_t zStraw1 = 0;
	Double_t zStraw4 = 0;
	Double_t acceptanceHole;
	Double_t acceptanceSide;
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
