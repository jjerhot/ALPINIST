#include <iostream>
#include "DecayMCProcess.C"
//#include "DecayMCGlobal.h"

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
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
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

int main(int argc, char** argv) {
    Int_t experiment,  productionmode,  decaymode, nEvents;
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        std::cout<< "usage: ALPs_BD [-h] [-e EXP] [-p PROD] [-d DECAY] [-n NEVTS]" << std::endl << std::endl;
        std::cout<< "ALP MC decay module. Select the experiment, production and decay modes and number of MC events" << std::endl << std::endl;
        std::cout<< "required arguments:" << std::endl;
        std::cout<< "\t -h\t\tshow this help message and exit" << std::endl;
        std::cout<< "\t -e EXP\t\toptions available: " << std::endl << "\t\t\t";
            for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";
        std::cout << std::endl << "\t -p PROD\toptions available: " << std::endl << "\t\t\t";
            for(auto imode : prodmodeNames) std::cout << imode.Data() << ", ";
        std::cout << std::endl << "\t -d DECAY\toptions available: " << std::endl << "\t\t\t";
            for(auto imode : decaymodeNames) std::cout << imode.Data() << ", ";
        std::cout << std::endl << "\t -n NEVTS\tnumber of MC events to generate per bin (not POT! only changing precision)" << std::endl;
        std::cout<< "optional arguments:" << std::endl;
        std::cout << std::endl << "\t -v\t\tverbosity flag" << std::endl;
        std::cout << std::endl << "\t -a 0/1\t\tall final states in calorimeter acceptance (true by default)" << std::endl;

        return 0;
    }

	if(!input.cmdOptionExists("-e")) {
		std::cout << "[Error] -e is required parameter, options available: " << std::endl;
        for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";
		std::cout << std::endl;
        return 1;
	} else{
        const std::string &exp = input.getCmdOption("-e");
        auto itr = std::find(expLabels.begin(), expLabels.end(), exp);
        if(itr == expLabels.end()){
			std::cout << "[Error] Invalid experiment label: " <<  exp << ". Use values ";
			for(auto iexp : expLabels) std::cout << iexp.Data() << ", ";
			std::cout << std::endl;
            return 1;
        }
        experiment = std::distance(expLabels.begin(), itr);
    }

	if(!input.cmdOptionExists("-p")) {
		std::cout << "[Error] -p is required parameter, options available: " << std::endl;
        for(auto imode : prodmodeNames) std::cout << imode.Data() << ", ";
		std::cout << std::endl;
        return 1;
	} else{
        const std::string &prod = input.getCmdOption("-p");
        auto itr = std::find(prodmodeNames.begin(), prodmodeNames.end(), prod);
        if(itr == prodmodeNames.end()){
			std::cout << "[Error] Invalid production mode: " <<  prod << ". Use values ";
			for(auto imode : prodmodeNames) std::cout << imode.Data() << ", ";
			std::cout << std::endl;
            return 1;
        }
        productionmode = std::distance(prodmodeNames.begin(), itr);
    }

	if(!input.cmdOptionExists("-d")) {
		std::cout << "[Error] -d is required parameter, options available: " << std::endl;
        for(auto imode : decaymodeNames) std::cout << imode.Data() << ", ";
		std::cout << std::endl;
        return 1;
	} else{
        const std::string &dec = input.getCmdOption("-d");
        auto itr = std::find(decaymodeNames.begin(), decaymodeNames.end(), dec);
        if(itr == decaymodeNames.end()){
			std::cout << "[Error] Invalid decay mode: " <<  dec << ". Use values ";
			for(auto imode : decaymodeNames) std::cout << imode.Data() << ", ";
			std::cout << std::endl;
            return 1;
        }
        decaymode = std::distance(decaymodeNames.begin(), itr);
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
    Bool_t verbose = false;

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

    DecayMCProcess(experiment,  productionmode,  decaymode, nEvents, allInAcceptance, verbose);

    return 0;
}