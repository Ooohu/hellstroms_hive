#include "bdt_var.h"

std::string gadget_MakeSafeName(std::string safe_name){
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '('), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), ')'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '\\'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '/'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '['), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), ']'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '+'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '-'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '*'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '.'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), ' '), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), ','), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), '|'), safe_name.end());
	safe_name.erase(std::remove(safe_name.begin(), safe_name.end(), ':'), safe_name.end());

	return safe_name;
};


void gadget_loadV( std::vector<int>& out_v, std::string in){ out_v.push_back( std::stod(in));}
void gadget_loadV( std::vector<double>& out_v, std::string in){ out_v.push_back( std::stod(in));}
void gadget_loadV( std::vector<std::string>& out_v, std::string in){ out_v.push_back( in);}
void gadget_loadV( std::vector<TString>& out_v, std::string in){ out_v.push_back( in.c_str());}
