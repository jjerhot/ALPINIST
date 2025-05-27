 /// \file DecayMC.C
 /// \Brief
 /// ALP_decay main with argument parser
 /// \EndBrief
 /// \Detailed
 /// After compilation to get available mandatory and optional arguments:
 /// \code
 /// ./DecayMC -h
 /// \endcode
 /// \EndDetailed

#include <iostream>
#include "DecayMCProcess.C"
//#include "DecayMCGlobal.h"

 /// \struct InputParser
 /// \Brief
 /// Parser for input command line arguments
 /// \EndBrief
class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    virtual ~InputParser(){
        tokens.clear();
    }
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
                != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

 /// \fn main
 /// \Brief
 /// Calls parser and runs the DecayMCProcess() with the obtained arguments
 /// \EndBrief
int main(int argc, char** argv) {
    Int_t exotic, experiment,  productionmode, decaymode, nEvents;
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        std::cout<< "usage: ALPs_BD [-h] [-x EXOTIC] [-e EXP] [-p PROD] [-d DECAY] [-n NEVTS]" << std::endl << std::endl;
        std::cout<< "ALP MC decay module. Select the experiment, production and decay modes and number of MC events" << std::endl << std::endl;
        std::cout<< "required arguments:" << std::endl;
        std::cout<< "\t -h\t\tshow this help message and exit" << std::endl;
        std::cout<< "\t -x EXOTIC\t\toptions available: " << std::endl << "\t\t\t";
        for(auto iexp : exoLabels) std::cout << iexp.Data() << ", ";

        std::cout <<"\b\b "<<std::endl<<"\t -e EXP\t\toptions available: " << std::endl << "\t\t\t";
        for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";   

        std::cout <<"\b\b  "<<std::endl<<"\t -p PROD\toptions available for  ";
        for(Int_t iexo{}; iexo<exoLabels.size();++iexo){
            std::cout << "\b\b  \n\t\t" <<   exoLabels[iexo] << ":\t";
            for(auto imode : prodmodeNames[iexo]) std::cout << imode.Data() << ", ";
        }

        std::cout << "\b\b  "<<std::endl<<"\t -d DECAY\toptions available for  ";
        for(Int_t iexo{}; iexo<exoLabels.size();++iexo){
            std::cout << "\b\b  \n\t\t" <<   exoLabels[iexo] << ":\t";
            for(auto imode : decaymodeNames[iexo]) std::cout << imode.Data() << ", ";
        }std::cout << "\b\b  ";
        std::cout << std::endl<<"\t -n NEVTS\tnumber of MC events to generate per bin (not POT! only changing precision)" << std::endl;
        std::cout << "optional arguments:" << std::endl;
        std::cout << std::endl << "\t --natt NEVTS\tnumber of attempts to simulate an MC event in acceptance (default = 1000)" << std::endl;
        std::cout << std::endl << "\t -v\t\tverbosity flag" << std::endl;
        std::cout << std::endl << "\t -a 0/1\t\tall final states in calorimeter acceptance (true by default)" << std::endl;
        std::cout << std::endl << "\t -s SEED\t(system time by default)" << std::endl;
        std::cout << std::endl << "\t --flat-decay\t\tflat decay and reweight flag" << std::endl;
        std::cout << std::endl << "\t --nbins-x\t\tnumber of bins for x-axis (default = 101)" << std::endl;
        std::cout << std::endl << "\t --nbins-y\t\tnumber of bins for y-axis (default = 101)" << std::endl;
        std::cout << std::endl << "\t --linear-x\t\tuse linear scale for x-axis instead of logarithmic" << std::endl;
        std::cout << std::endl << "\t --linear-y\t\tuse linear scale for y-axis instead of logarithmic" << std::endl;
        std::cout << std::endl << "\t --variable-y\t\tchoose variable for y-axis, default is Width [GeV]; options available: " << std::endl << "\t\t\t";
            for(auto ivar : yVars) std::cout << ivar.Data() << ", ";
        std::cout << std::endl << "\t --active-coupling\t\tchoose active coupling for which to simulate decays (only applicable to hnl, default all); options available: El, Mu, Tau, all";
        std::cout << std::endl;
        return 0;
    }

    if(!input.cmdOptionExists("-x")) {
		std::cout << "[Error] -x is required parameter, options available: ";
        for(auto iexo : exoLabels) std::cout << iexo.Data() << ", ";
		std::cout << "\b\b " << std::endl;
        return 1;
	} else{
        const std::string &exo = input.getCmdOption("-x");
        auto itr = std::find(exoLabels.begin(), exoLabels.end(), exo);
        if(itr == exoLabels.end()){
			std::cout << "[Error] Invalid exotic label: " << exo << ". Use values ";
			for(auto iexo : exoLabels) std::cout << iexo.Data() << ", ";
			std::cout << "\b\b " << std::endl;
            return 1;
        }
        exotic = std::distance(exoLabels.begin(), itr);
    }

	if(!input.cmdOptionExists("-e")) {
		std::cout << "[Error] -e is required parameter, options available: ";
        for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";
		std::cout << "\b\b " << std::endl;
        return 1;
	} else{
        const std::string &exp = input.getCmdOption("-e");
        auto itr = std::find(expLabels.begin(), expLabels.end(), exp);
        if(itr == expLabels.end()){
			std::cout << "[Error] Invalid experiment label: " <<  exp << ". Use values ";
			for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";
			std::cout << "\b\b " << std::endl;
            return 1;
        }
        experiment = std::distance(expLabels.begin(), itr);
    }

	if(!input.cmdOptionExists("-p")) {
		std::cout << "[Error] -p is required parameter, options available: " << std::endl;
        for(auto imode : prodmodeNames[exotic]) std::cout << imode.Data() << ", ";
		std::cout << "\b\b " << std::endl;
        return 1;
	} else{
        const std::string &prod = input.getCmdOption("-p");
        auto itr = std::find(prodmodeNames[exotic].begin(), prodmodeNames[exotic].end(), prod);
        if(itr == prodmodeNames[exotic].end()){
			std::cout << "[Error] Invalid production mode: " <<  prod << ". Use values ";
			for(auto imode : prodmodeNames[exotic]) std::cout << imode.Data() << ", ";
			std::cout << "\b\b " << std::endl;
            return 1;
        }
        productionmode = std::distance(prodmodeNames[exotic].begin(), itr);
    }

	if(!input.cmdOptionExists("-d")) {
		std::cout << "[Error] -d is required parameter, options available: " << std::endl;
        for(auto imode : decaymodeNames[exotic]) std::cout << imode.Data() << ", ";
		std::cout << "\b\b " << std::endl;
        return 1;
	} else{
        const std::string &dec = input.getCmdOption("-d");
        auto itr = std::find(decaymodeNames[exotic].begin(), decaymodeNames[exotic].end(), dec);
        if(itr == decaymodeNames[exotic].end()){
			std::cout << "[Error] Invalid decay mode: " <<  dec << ". Use values ";
			for(auto imode : decaymodeNames[exotic]) std::cout << imode.Data() << ", ";
			std::cout << std::endl;
            return 1;
        }
        decaymode = std::distance(decaymodeNames[exotic].begin(), itr);
    }

	if(!input.cmdOptionExists("-n")) {
		std::cout << "[Error] -n is required parameter, exiting" << std::endl;
		return 1;
	} else{
        const std::string &nevt = input.getCmdOption("-n");
        if(std::all_of(nevt.cbegin(), nevt.cend(), ::isdigit))
            nEvents = std::stoi(nevt);
        else{
            std::cout << "[Error] -n has to be integer, exiting" << std::endl;
            return 1;
        }
    }

    //optional arguments
    Bool_t allInAcceptance = true;
    Int_t nAttempts = 1000;
    Bool_t verbose = false;
    Int_t seed = 0;
    Bool_t flatDecay = false;
    Bool_t flatDalitz = false;
    Int_t nBinsX = 101;
    Int_t nBinsY = 101;
    Bool_t xIsLin = false;
    Bool_t yIsLin = false;
    Int_t yVar = 0; //0-Decay width (GeV), 1-tau (fs), 2-ctau (m)
    Int_t actCoup;

	if(input.cmdOptionExists("--natt")){
        const std::string &nbins = input.getCmdOption("--natt");
        if(std::all_of(nbins.cbegin(), nbins.cend(), ::isdigit)){
            nAttempts = std::stoi(nbins);
            if(nAttempts < 1){
                std::cout << "[Error] --natt has to be positive" << std::endl;
                return 1;
            }
        }
        else{
            std::cout << "[Error] --natt has to be positive integer, exiting" << std::endl;
            return 1;
        }
    }

	if(input.cmdOptionExists("-v"))
        verbose = true;

	if(input.cmdOptionExists("-a")){
        if(std::stoi(input.getCmdOption("-a")) == 0)
            allInAcceptance = false;
        else if(std::stoi(input.getCmdOption("-a")) == 1)
            allInAcceptance = true;
        else{
            std::cout << "[Error] -a invalid value use 0 or 1" << std::endl;
            return 1;
        }
    }

	if(input.cmdOptionExists("-s")) {
        const std::string &seedstr = input.getCmdOption("-s");
        if(std::all_of(seedstr.cbegin(), seedstr.cend(), ::isdigit))
            seed = std::stoi(seedstr);
        else{
            std::cout << "[Error] -s has to be integer, exiting" << std::endl;
            return 1;
        }
    }

	if(input.cmdOptionExists("--flat-decay"))
        flatDecay = true;

	if(input.cmdOptionExists("--flat-dalitz"))
        flatDalitz = true;

	if(input.cmdOptionExists("--nbins-x")){
        const std::string &nbins = input.getCmdOption("--nbins-x");
        if(std::all_of(nbins.cbegin(), nbins.cend(), ::isdigit))
            nBinsX = std::stoi(nbins);
        else{
            std::cout << "[Error] --nbins-x has to be integer, exiting" << std::endl;
            return 1;
        }
    }

	if(input.cmdOptionExists("--nbins-y")){
        const std::string &nbins = input.getCmdOption("--nbins-y");
        if(std::all_of(nbins.cbegin(), nbins.cend(), ::isdigit))
            nBinsY = std::stoi(nbins);
        else{
            std::cout << "[Error] --nbins-y has to be integer, exiting" << std::endl;
            return 1;
        }
    }

	if(input.cmdOptionExists("--linear-x"))
        xIsLin = true;

	if(input.cmdOptionExists("--linear-y"))
        yIsLin = true;

	if(!input.cmdOptionExists("--variable-y")) {
        yVar = 0; //0-Decay width (GeV)
	} else{
        const std::string &yVarChosen = input.getCmdOption("--variable-y");
        auto itr = std::find(yVars.begin(), yVars.end(), yVarChosen);
        if(itr == yVars.end()){
			std::cout << "[Error] Invalid y-axis variable: " << yVarChosen << ". Use values ";
			for(auto ivar : yVars) std::cout << ivar.Data() << ", ";
			std::cout << "\b\b " << std::endl;
            return 1;
        }
        yVar = std::distance(yVars.begin(), itr);
    }

    if(!input.cmdOptionExists("--active-coupling")) {
        actCoup = -1; //0-Decay width (GeV)
	} else if (exotic != 1){
        std::cout << "[Error] Flag --active-coupling only applicable to exotic hnl, exotic chosen was '";
        std::cout << exoLabels.at(exotic) << "'."<< std::endl;
        return 1;
    }else{
        std::string actCoupChosen = input.getCmdOption("--active-coupling");
        if (actCoupChosen != "all"){
            actCoupChosen.insert(0, 1 , '-'); 
            actCoupChosen.append("Mixing");
            auto itr = std::find(activeCouplings[exotic].begin(), activeCouplings[exotic].end(), actCoupChosen);
            if(itr == activeCouplings[exotic].end()){
                std::cout << "[Error] Invalid active coupling: '" << actCoupChosen.substr(1,actCoupChosen.length()-7) << "'. Use values ";
                for(auto coupling : activeCouplings[exotic]) std::cout << coupling.Remove(0,1).Remove(coupling.Length()-6,6) << ", ";
                std::cout << "all " << std::endl;
                return 1;
            }
            actCoup = std::distance(activeCouplings[exotic].begin(), itr);
        } else actCoup = -1;
    }

    DecayMCProcess(exotic, experiment, productionmode, decaymode, nEvents, allInAcceptance, nAttempts, verbose, seed, flatDecay, flatDalitz,nBinsX,nBinsY,xIsLin,yIsLin,yVar,actCoup);

    return 0;
}