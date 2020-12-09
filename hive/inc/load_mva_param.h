#ifndef LOADMVAPARAM_H
#define LOADMVAPARAM_H

#include "tinyxml.h"

//#include "method_struct.h"
#include "bdt_var.h"
#include "bdt_file.h"
#include "bdt_systematics.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "TColor.h"
#include "TMVA/Types.h"


class MVALoader {

	protected:

//	std::vector<method_struct> vec_methods;
    std::map<std::string,std::string> aliasMap; 
	std::vector< std::string > stage_cuts;
	std::vector< std::string > pre_cuts;

	std::vector< bdt_file* > all_files;

	std::vector< bdt_sys > all_systematics;

	std::vector< bdt_info > all_bdtinfos;

	std::vector< bdt_variable > all_variables;

	public:
	
	//environmental parameters;
	TiXmlDocument* doc;
	int verbosity;//0 - minimal, 1 - normal message, 2 - debug message;
	bool isSys;
	std::string whichxml;//xml file name;
    std::string analysis_tag;
    std::string topo_name;
    std::string inputdir;

	//essential funcitons
	MVALoader(std::string xmlname)
		: 
		isSys(false),
		whichxml(xmlname),
		verbosity(1){
			LoadEnv();
			LoadBDTfiles();//this turns on isSys to be ture, if any systematic file is defined;
			LoadBDTSettings();
			LoadVariables();
		}


	void LoadEnv();
	void LoadBDTfiles();
	void LoadBDTSettings();
	void LoadVariables();

	std::string AliasParse(std::string in);

//	std::vector<method_struct> GetMethods();
    size_t GetNFiles(){return 0;}
    
	std::vector< bdt_file* > GetFiles(){ return all_files;};
	std::vector< bdt_sys > GetSys(){ return all_systematics;};
	std::vector< bdt_variable > GetVar(){return all_variables;};
	std::vector< bdt_info > GetBDTInfo(){return all_bdtinfos;};
	std::vector<std::string> GetPreCuts( ){return pre_cuts;};
	std::string GetCuts( int nth_cut){return stage_cuts[nth_cut];};
	int GetStagesCount( ){return stage_cuts.size();};
	
};

#endif
