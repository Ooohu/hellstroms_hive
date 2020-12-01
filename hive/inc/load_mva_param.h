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

	std::vector< bdt_file > all_files;

	std::vector< bdt_variable > all_variables;

	std::vector< bdt_sys > all_systematics;

	std::vector< bdt_info > all_bdtinfos;


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
    
	std::vector< bdt_file > GetFiles(){ return all_files;};
	std::vector< bdt_sys > GetSys(){ return all_systematics;};
	std::vector< bdt_variable > GetVar(){return all_variables;};
	std::vector< bdt_info > GetInfo(){return all_bdtinfos;};
	std::string GetCuts( int nth_cut){return stage_cuts[nth_cut];};
	
    //BDT file info

	//bdt_files
    size_t n_bdt_files;
    std::vector<std::string> bdt_filenames;
    std::vector<std::string> bdt_tags;
    std::vector<std::string> bdt_hist_styles;
    std::vector<std::string> bdt_dirs;
    std::vector<std::string> bdt_plotnames;
    std::vector<std::string> bdt_weight;
    std::vector<bool> bdt_is_validate_file;

    std::vector<TColor*> bdt_cols;
    std::vector<int> bdt_linecols;
    std::vector<int> bdt_linestyles;
    std::vector<int> bdt_group;//Keng

    std::vector<int> bdt_fillstyles;
    std::vector<double> bdt_scales;
    std::vector<std::vector<std::string>> bdt_definitions;
    std::vector<std::vector<std::string>> bdt_training_cuts;
    std::vector<bool> bdt_is_onbeam_data;
    std::vector<double> bdt_fixpot;
    std::vector<double> bdt_onbeam_pot;
    std::vector<double> bdt_onbeam_spills;
    std::vector<bool> bdt_is_offbeam_data;
    std::vector<double> bdt_offbeam_spills;

    std::vector<bool> bdt_is_signal;
    std::vector<bool> bdt_on_top;
    std::vector<bool> bdt_is_training_signal;



    std::vector<int> v_eff_denom_stage;
    std::vector<int> v_eff_numer_stage;
    std::vector<std::string> v_eff_denom_cut;
    std::vector<std::string> v_eff_numer_cut;
    


    std::vector<double> bdt_cuts;

    std::vector<std::string> recomc_names;
	std::vector<std::string> recomc_defs;
	std::vector<TColor*> recomc_cols;

	//these are for systematics;
	int n_sys;//keep track with the number of systematics
	std::vector< TString>  sys_tag;
	std::vector< TString>  sys_dir;
	std::vector< TString>  sys_filename;
	std::vector< TString>  sys_treename;

	std::vector< std::vector< TString >> sys_for_files;
	std::vector< std::vector< TString >> sys_vars;
	std::vector< std::vector< TString >> sys_vars_name;

	std::vector< double> sys_pot;
	std::vector< double> sys_throws;

	std::vector< bool> sys_its_CV;
	std::vector< bool> sys_its_multithrows;
	std::vector< bool> sys_its_OpticalModel;


};

#endif
