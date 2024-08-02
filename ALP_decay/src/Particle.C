 /// \class Particle
 /// \Brief
 /// Contains particle properties and pointers to decay products
 /// \EndBrief
 /// \Detailed
 /// Base class for both exotic and SM particles
 /// Automatically decays rhos and light neutral pseudoscalars.
 /// \EndDetailed

#include "Particle.h"

Particle::Particle() : parent(nullptr), name(""), unstable(false), mass(0.), charge(0), pos(0.,0.,0.), mom(0.,0.,0.),
				 initpos(0.,0.,0.), endpos(0.,0.,0.), initmom(0.,0.,0.),endmom(0.,0.,0.), alive(true){}

Particle::Particle(Particle* orig) : Particle(){ //copy all properties
	name = orig->name;
	unstable = orig->unstable;
	mass = orig->mass;
	charge = orig->charge;
	pos = orig->pos;
	mom = orig->mom;
	initpos = orig->initpos;
	endpos = orig->endpos;
	initmom = orig->initmom;
	endmom = orig->endmom;
	children = orig->children;
	parent = orig->parent;
	alive = orig->alive;
}

Particle::Particle(TString name_, Bool_t unstable_, Double_t mass_, Int_t charge_, TVector3 position_, TVector3 momentum_, Particle* parent_) : Particle() {
	initpos = position_;
	pos = position_;
	initmom = momentum_;
	mom = momentum_;
	parent = parent_; //check
	Init(name_,unstable_,mass_,charge_);
	CheckStability();
}

Particle::Particle(SMParticleProperty *prop, TVector3 position_, TVector3 momentum_, Particle* parent_) : Particle() {
	initpos = position_;
	pos = position_;
	initmom = momentum_;
	mom = momentum_;
	parent = parent_; //check
	if(prop->GetName()=="notSM")
		std::cerr << "[Particle] Warning: not a SM particle; maybe try to initialize manually?" << std::endl;
	Init(prop->GetName(),prop->GetUnstable(),prop->GetMass(),prop->GetCharge());
	CheckStability();
}

Particle::Particle(ExoticParticleProperty *prop, Particle* parent_) : Particle() {
	Init(prop->GetName(),false,prop->Mass,0);
	initpos = prop->Position;
	pos = prop->Position;
	initmom = prop->Momentum;
	mom = prop->Momentum;
	parent = parent_; //check
}

Particle::~Particle(){
    while(!children.empty()) {
        delete children.back();
        children.pop_back();
    }
	children.clear();
}

void Particle::Init(TString name_, Bool_t unstable_, Double_t mass_, Int_t charge_){
	name = name_;
	unstable = unstable_;
	mass = mass_;
	charge = charge_;
}

void Particle::CheckStability(){
	if(unstable){
		if(name =="pi0" || name =="eta" || name =="eta_prime")
			Decay2Body("gamma","gamma");
		else if(name =="rho+")
			Decay2Body("pi+","pi0");
		else if(name =="rho-")
			Decay2Body("pi-","pi0");
		else if(name =="rho0")
			Decay2Body("pi0","pi0");
	}
}

 /// \fn Propagate
 /// \Brief
 //// Propagate particle by distance
 /// \EndBrief
void Particle::Propagate(Double_t dist){
	if(alive)
		pos.SetXYZ(pos.X() + mom.X()/mom.Mag()*dist,pos.Y() + mom.Y()/mom.Mag()*dist,pos.Z() + mom.Z()/mom.Mag()*dist);
	else
		std::cerr << "[Particle] Warning: cannot propagate, particle doesn't exist" << std::endl;
}

 /// \fn Propagate
 /// \Brief
 //// Propagate particle to Z position
 /// \EndBrief
void Particle::PropagateToZ(Double_t zFin){
	if(alive)
		pos.SetXYZ(pos.X() + mom.X()/mom.Z()*(zFin-pos.Z()),pos.Y() + mom.Y()/mom.Z()*(zFin-pos.Z()),zFin);
	else
		std::cerr << "[Particle] Warning: cannot propagate, particle doesn't exist" << std::endl;
}

 /// \fn Kick
 /// \Brief
 //// Apply momentum kick
 /// \EndBrief
void Particle::Kick(Double_t kick, Double_t fieldInfo, Bool_t isToroidal){
	if(alive){
		if(charge){
			if(isToroidal){ 
				Double_t scale  = (fieldInfo - pos.Perp()) / fieldInfo ; //using effective field radius to scale maximal kick appropriately 
				kick *= (0 < scale && scale < 0.975) ? scale : 0.; 		 //coil at the centre of the magnet
				fieldInfo = pos.Phi()+TMath::PiOver2();	  //making FieldInfo appropriate MagnetFiledPhi for consequent calculation
			}
			Double_t mom2Before = mom.Mag2();
			mom.SetX(mom.X() - TMath::Sin(fieldInfo)*charge*kick); //adjust in x direction, assume momentum in GeV
			mom.SetY(mom.Y() + TMath::Cos(fieldInfo)*charge*kick);
			mom.SetZ(TMath::Sqrt(TMath::Max(mom2Before - mom.Perp2(), 0.))); //conserve momentum
		}
	} else
		std::cerr << "[Particle] Warning: cannot kick, particle doesn't exist" << std::endl;
}

// void Particle::KickX(Double_t kick){ //magnetic kick in x-direction (kick [T.m])
// 	if(alive){
// 		if(charge){
// 			Double_t momBefore = mom.Mag();
// 			mom.SetX(mom.X() + charge*kick*c*1E-9); //adjust in x direction, assume momentum in GeV
// 			(momBefore/mom.Mag())*mom; //conserve momentum
// 		}
// 	} else
// 		std::cerr << "[Particle] Warning: cannot kick, particle doesn't exist" << std::endl;
// }

// void Particle::KickY(Double_t kick){ //magnetic kick in y-direction (kick [T.m])
// 	if(alive){
// 		if(charge){
// 			// Double_t momBefore = mom.r();
// 			Double_t momBefore = mom.Mag();
// 			mom.SetY(mom.Y() + charge*kick*c*1E-9); //adjust in y-direction, assume momentum in GeV
// 			(momBefore/mom.Mag())*mom; //conserve momentum
// 		}
// 	} else
// 		std::cerr << "[Particle] Warning: cannot kick, particle doesn't exist" << std::endl;
// }

 /// \fn Decay2Body
 /// \Brief
 //// Decay particle via 2-body decay to given final state objects
 /// \EndBrief
void Particle::Decay2Body(SMParticleProperty* prop1, SMParticleProperty* prop2){
	if(prop1->GetMass()+prop2->GetMass() > mass){
		std::cerr << "[Particle] Warning: decay not kinematically allowed, skipping" << std::endl;
		return;
	}
	Double_t pStar = TMath::Sqrt(mass*mass-4.*prop1->GetMass()*prop2->GetMass())/2;
	Double_t betaMag = mom.Mag()/GetEnergy();
	TVector3 beta(betaMag*TMath::Sin(mom.Theta())*TMath::Cos(mom.Phi()),
					betaMag*TMath::Sin(mom.Theta())*TMath::Sin(mom.Phi()), 
					betaMag*TMath::Cos(mom.Theta()));
	Double_t thetaDecay = TMath::ACos(-1. + 2.*rndmGen.Rndm());  // theta value in (0,pi); (flat in costheta)
	Double_t phiDecay = 2.*TMath::Pi()*rndmGen.Rndm(); // flat distribution
	// Double_t thetaDecay = 0.;
	// Double_t phiDecay = 0.;
	TLorentzVector mom1(pStar*TMath::Sin(thetaDecay)*TMath::Cos(phiDecay),
			pStar*TMath::Sin(thetaDecay)*TMath::Sin(phiDecay),
			pStar*TMath::Cos(thetaDecay),
			TMath::Sqrt(pStar*pStar + prop1->GetMass()*prop1->GetMass()));
	TLorentzVector mom2(-pStar*TMath::Sin(thetaDecay)*TMath::Cos(phiDecay),
			-pStar*TMath::Sin(thetaDecay)*TMath::Sin(phiDecay),
			-pStar*TMath::Cos(thetaDecay),
			TMath::Sqrt(pStar*pStar + prop2->GetMass()*prop2->GetMass()));

	mom1.Boost(beta);
	mom2.Boost(beta);

	children.push_back(new Particle(prop1,pos,mom1.Vect(),this));
	children.push_back(new Particle(prop2,pos,mom2.Vect(),this));

	EndParticle();
}

 /// \fn Decay2Body
 /// \Brief
 //// 2-body decay given only final state names
 /// \EndBrief
void Particle::Decay2Body(TString particle1, TString particle2){
	if(alive){
		SMParticleProperty prop1(particle1);
		SMParticleProperty prop2(particle2);

		if(prop1.GetName()=="notSM" || prop2.GetName()=="notSM") //maybe add in the future decays into exotics
			std::cerr << "[Particle] Warning: decaying to undefined type of particle" << std::endl;
		
		Decay2Body(&prop1, &prop2);
	} else
		std::cerr << "[Particle] Warning: cannot decay, particle doesn't exist" << std::endl;
}

 /// \fn Decay3Body
 /// \Brief
 //// Decay particle via 3-body decay to given final state objects assuming some Dalitz distribution - returns corresponding weight
 /// \EndBrief
Double_t Particle::Decay3Body(SMParticleProperty* prop1, SMParticleProperty* prop2, SMParticleProperty* prop3, Double_t &m12, Double_t &m23, TH2D* dalitzPlot){
	Double_t weight = 0.;
	if(alive){
		if(prop1->GetMass()+prop2->GetMass()+prop3->GetMass() > mass){
			std::cerr << "[Particle] Warning: decay not kinematically allowed, skipping" << std::endl;
			return weight;
		}
		Double_t daughterMasses[3]={prop1->GetMass(),prop2->GetMass(),prop3->GetMass()};
		TGenPhaseSpace event;
		TLorentzVector lorentzMom(mom,GetEnergy());
		event.SetDecay(lorentzMom, 3, daughterMasses);
		weight = event.Generate();

		children.push_back(new Particle(prop1,pos,(*event.GetDecay(0)).Vect(),this));
		children.push_back(new Particle(prop2,pos,(*event.GetDecay(1)).Vect(),this));
		children.push_back(new Particle(prop3,pos,(*event.GetDecay(2)).Vect(),this));

		Double_t m12sq = daughterMasses[0]*daughterMasses[0]+daughterMasses[1]*daughterMasses[1]+2.*(*event.GetDecay(0)).Dot(*event.GetDecay(1));
		Double_t m23sq = daughterMasses[1]*daughterMasses[1]+daughterMasses[2]*daughterMasses[2]+2.*(*event.GetDecay(1)).Dot(*event.GetDecay(2));
		
		if(dalitzPlot != nullptr) {//weight from DalitzPlot
			if (m12sq < dalitzPlot->GetXaxis()->GetXmin() || dalitzPlot->GetXaxis()->GetXmax() < m12sq )
				std::cout << "M_12^2 value " << m12sq << "out of bounds " << dalitzPlot->GetXaxis()->GetXmin()<< ", "<<dalitzPlot->GetXaxis()->GetXmax()<<std::endl;
			else if(m23sq < dalitzPlot->GetYaxis()->GetXmin() || dalitzPlot->GetYaxis()->GetXmax() < m23sq)
				std::cout << "M_23^2 value " << m23sq << "out of bounds " << dalitzPlot->GetYaxis()->GetXmin()<< ", "<<dalitzPlot->GetYaxis()->GetXmax()<<std::endl;
			else {
				weight *= dalitzPlot->GetBinContent(dalitzPlot->GetXaxis()->FindFixBin(m12sq),dalitzPlot->GetYaxis()->FindFixBin(m23sq));
				// if (dalitzPlot->GetSum()) weight *= (dalitzPlot->GetNbinsX())*(dalitzPlot->GetNbinsY())/(dalitzPlot->GetSum()); 
				// else weight = 0.;
			}
		}
		//add sampling for nonflat
		EndParticle();
		return weight;
	} else
		std::cerr << "[Particle] Warning: cannot decay, particle doesn't exist" << std::endl;

	return weight;
}

 /// \fn Decay3Body
 /// \Brief
 //// Decay particle via 3-body decay to given final state names assuming some Dalitz distribution - returns corresponding weight
 /// \EndBrief
Double_t Particle::Decay3Body(TString particle1, TString particle2, TString particle3, Double_t &m12, Double_t &m23, TH2D* dalitzPlot){
	Double_t weight = 0.;
	if(alive){
		SMParticleProperty prop1(particle1);
		SMParticleProperty prop2(particle2);
		SMParticleProperty prop3(particle3);

		if(prop1.GetName()=="notSM" || prop2.GetName()=="notSM" || prop3.GetName()=="notSM") //maybe add in the future decays into exotics
			std::cerr << "[Particle] Warning: decaying to undefined type of particle" << std::endl;
		else
			weight = Decay3Body(&prop1, &prop2, &prop3, m12, m23, dalitzPlot);
	} else
		std::cerr << "[Particle] Warning: cannot decay, particle doesn't exist" << std::endl;
	return weight;
}

 /// \fn Decay4Body
 /// \Brief
 //// Decay particle via 4-body decay to given final state objects
 /// \EndBrief
Double_t Particle::Decay4Body(SMParticleProperty* prop1, SMParticleProperty* prop2, SMParticleProperty* prop3, SMParticleProperty* prop4){
	Double_t weight = 0.;
	if(alive){
		if(prop1->GetMass()+prop2->GetMass()+prop3->GetMass()+prop4->GetMass() > mass){
			std::cerr << "[Particle] Warning: decay not kinematically allowed, skipping" << std::endl;
			return weight;
		}
		Double_t daughterMasses[4]={prop1->GetMass(),prop2->GetMass(),prop3->GetMass(),prop4->GetMass()};
		TGenPhaseSpace event;
		TLorentzVector lorentzMom(mom,GetEnergy());
		event.SetDecay(lorentzMom, 4, daughterMasses);
		weight = event.Generate();

		children.push_back(new Particle(prop1,pos,(*event.GetDecay(0)).Vect(),this));
		children.push_back(new Particle(prop2,pos,(*event.GetDecay(1)).Vect(),this));
		children.push_back(new Particle(prop3,pos,(*event.GetDecay(2)).Vect(),this));
		children.push_back(new Particle(prop4,pos,(*event.GetDecay(3)).Vect(),this));

		//add sampling for nonflat
		EndParticle();
		return weight;
	} else
		std::cerr << "[Particle] Warning: cannot decay, particle doesn't exist" << std::endl;

	return weight;
}

 /// \fn Decay4Body
 /// \Brief
 //// 4-body decay given only final state names
 /// \EndBrief
Double_t Particle::Decay4Body(TString particle1, TString particle2, TString particle3, TString particle4){
	Double_t weight = 0.;
	if(alive){
		SMParticleProperty prop1(particle1);
		SMParticleProperty prop2(particle2);
		SMParticleProperty prop3(particle3);
		SMParticleProperty prop4(particle3);

		if(prop1.GetName()=="notSM" || prop2.GetName()=="notSM" || prop3.GetName()=="notSM" || prop4.GetName()=="notSM") //maybe add in the future decays into exotics
			std::cerr << "[Particle] Warning: decaying to undefined type of particle" << std::endl;
		else
			weight = Decay4Body(&prop1, &prop2, &prop3, &prop4);
	} else
		std::cerr << "[Particle] Warning: cannot decay, particle doesn't exist" << std::endl;
	return weight;
}

 /// \fn EndParticle
 /// \Brief
 //// absorb or end by some other means than Decay2Body or Decay3Body
 /// \EndBrief
void Particle::EndParticle(){ 
	endpos = pos;
	endmom = mom;
	alive = false;
}

 /// \fn Reset
 /// \Brief
 //// Put back in initial state
 /// \EndBrief
void Particle::Reset(TVector3 momentum_) {
	alive = true;
    while(!children.empty()) {
        delete children.back();
        children.pop_back();
    }
	pos = initpos;
	mom = momentum_;
	initmom = momentum_;
	endpos.SetXYZ(0., 0., 0.);
	endmom.SetXYZ(0., 0., 0.);
	CheckStability();
}