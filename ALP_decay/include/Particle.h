//////////////////////////////////////////////////////////////////
// ALPINIST/DecayMC
// Author: J.Jerhot 20/6/2023
// Classes for manipulating particles and storing their properties
//////////////////////////////////////////////////////////////////

#ifndef PARTICLE_H
#define PARTICLE_H

#include "DecayMCGlobal.h"

#include <iostream>
#include <TROOT.h>
#include "TGenPhaseSpace.h"
#include "TH2D.h"

class ExoticParticleProperty{
private:
	TString name;

public:
	ExoticParticleProperty() :  name(""), Mass(0.), Position(0.,0.,0.), Momentum(0.,0.,0.),
								DecayWidth(0.), ProductionCrossSection(0.) {}

	ExoticParticleProperty(TString name_) : ExoticParticleProperty(){
		name = name_;
	}

	TString GetName() const { return name; }
	Double_t Mass;	//GeV
	TVector3 Position;
	TVector3 Momentum;	//GeV/c
	Double_t DecayWidth; //GeV
	Double_t ProductionCrossSection;
};

class SMParticleProperty{
private:
	void Init(TString name_, Bool_t unstable_, Double_t mass_, Int_t charge_){
		name = name_;
        unstable = unstable_;
        mass = mass_;
        charge = charge_;
	}
		TString name;
		Bool_t unstable; //if not stable decay immediately
		Double_t mass;
		Int_t charge;
public:
	SMParticleProperty() :  name(""), unstable(false), mass(0.), charge(0) {}
	SMParticleProperty(TString name_) : SMParticleProperty(){
		if(name_=="gamma")
			Init(name_,false,0.,0);
		else if(name_=="el+")
			Init(name_,false,MEl,1);
		else if(name_=="el-")
			Init(name_,false,MEl,-1);
		else if(name_=="mu+")
			Init(name_,false,MMu,1);
		else if(name_=="mu-")
			Init(name_,false,MMu,-1);
		else if(name_=="pi+")
			Init(name_,false,MPiCh,1);
		else if(name_=="pi-")
			Init(name_,false,MPiCh,-1);
		else if(name_=="pi0")
			Init(name_,true,MPi0,0);
		else if(name_=="eta")
			Init(name_,true,MEta,0);
		else if(name_=="eta_prime")
			Init(name_,true,MEtaPrim,0);
		else if(name_=="kaon+")
			Init(name_,false,MKCh,1);
		else if(name_=="kaon-")
			Init(name_,false,MKCh,-1);
		else if(name_=="rho+")
			Init(name_,true,MRho,1);
		else if(name_=="rho-")
			Init(name_,true,MRho,-1);
		else if(name_=="rho0")
			Init(name_,true,MRho,0);
		else if(name_=="nu")
			Init(name_,false,0.,0); //assume 0 mass for nu
		else
			Init("notSM",false,0,0); //default (not a SM particle)
	}
	TString GetName() const { return name; }
	Bool_t GetUnstable() const { return unstable; }
	Double_t GetMass() const { return mass; }
	Double_t GetCharge() const { return charge; }
};

//simplified particles .. maybe add momentum, width, etc. and apply also for exotics?
class Particle{

private:
	void Init(TString, Bool_t, Double_t, Int_t);

	void CheckStability();

	//particle properties
	TString name;
	Bool_t unstable; //if not stable decay immediately
	Double_t mass;
	Int_t charge;
	TVector3 initpos;
	TVector3 endpos;
	TVector3 initmom;
	TVector3 endmom;
	TVector3 pos;
	TVector3 mom;
	std::vector<Particle*> children;
	Particle *parent;
	Bool_t alive;

public:
	Particle();
	Particle(Particle*);
	Particle(TString, Bool_t, Double_t, Int_t, TVector3, TVector3, Particle* parent_ = nullptr);
	Particle(SMParticleProperty*, TVector3, TVector3, Particle* parent_ = nullptr);
	Particle(ExoticParticleProperty*, Particle* parent_ = nullptr);
	~Particle();

	Bool_t operator==(const Particle& par) const{
		if(name == par.name) return true;
		return false;
	}

	Particle& operator=(const Particle& orig){ //to be fixed
		name = orig.name;
		unstable = orig.unstable;
		mass = orig.mass;
		charge = orig.charge;
		pos = orig.pos;
		mom = orig.mom;
		initpos = orig.initpos;
		endpos = orig.endpos;
		initmom = orig.initmom;
		endmom = orig.endmom;
		children = orig.children;
		parent = orig.parent;
		alive = orig.alive;
		return *this;
	}

	//get particle properties
	TString GetName() const { return name; }
	Bool_t GetUnstable() const { return unstable; }
	Double_t GetMass() const { return mass; }
	Double_t GetCharge() const { return charge; }
	TVector3 GetPosition() const { return pos;}
	TVector3 GetMomentum() const { return mom;}
	TVector3 GetInitPosition() const { return initpos;}
	TVector3 GetInitMomentum() const { return initmom;}
	TVector3 GetEndPosition() const { return endpos;}
	TVector3 GetEndMomentum() const { return endmom;}
	
	UInt_t GetNChildren() const { return children.size();}
	Particle* GetDaughter(Int_t i) const { return children.at(i); }

	Particle* GetParent() const { return parent; }
	Bool_t GetAlive() const { return alive; }
	Double_t GetEnergy() const { return TMath::Sqrt(mom.Mag2() + mass*mass); }

	//modify particle properties
	void Propagate(Double_t);
	void PropagateToZ(Double_t);
	void Kick(Double_t, Double_t, Bool_t);
	void KickX(Double_t);
	void KickY(Double_t);
	void Decay2Body(TString, TString);
	void Decay2Body(SMParticleProperty*, SMParticleProperty*);
	Double_t Decay3Body(TString, TString, TString, Double_t &, Double_t &, TH2D* = nullptr);
	Double_t Decay3Body(SMParticleProperty*, SMParticleProperty*, SMParticleProperty*, Double_t &, Double_t &, TH2D* = nullptr);
	Double_t Decay4Body(TString, TString,TString, TString);
	Double_t Decay4Body(SMParticleProperty*, SMParticleProperty*,SMParticleProperty*, SMParticleProperty*);
	void EndParticle();
	void Reset(TVector3);
	void ChargeConjugate(){ charge *= -1;}
};

#endif