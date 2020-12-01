#ifndef BDT_INFO_H
#define BDT_INFO_H

#include <vector>
#include <string>
#include <iostream>
#include "bdt_var.h"
#include "method_struct.h"
/******** Our includes *****/

/******** Root includes *****/

struct bdt_info{

	public:
	
	//some convientant labels
	TMVA::Types::EMVA  type;
	std::string identifier;
	std::string name;	
	std::string binning;
	std::string str;//read from type; =  XGBoost;
    std::vector<std::pair<std::string,std::string>> xg_config;

	std::string base_cuts ;
 	std::string pre_cuts;
	std::string bdt_cosmic_cuts;
	std::string mid_cuts;
	std::string bdt_bnb_cuts;
	std::string post_cuts;

	std::string signal_definition;
	std::string background_definition;

    method_struct TMVAmethod;
    std::vector<bdt_variable> train_vars;
//    std::vector<bdt_variable> spec_vars;
	std::string topo_name;
	
	//CHECK who need method_struct??
	std::string option;


    std::string bdt_name;
    std::string bdt_binning;
    std::string bdt_tag;

    std::vector<bdt_variable> bdt_all_vars;
    std::vector<bdt_variable> bdt_train_vars;
//    std::vector<bdt_variable> bdt_spec_vars;

    std::string filename;
    std::string foldername;
    std::string training_cut;
    double training_fraction;

    //tags for training v2.2
    std::string bkg_test_tag;
    std::string bkg_train_tag;
    std::string bkg_test_cut;
    std::string sig_train_tag;
    std::string sig_test_tag;
    std::string sig_test_cut;

    double scan_max;
    double scan_min;
    double scan_steps;

    std::string topological_name;
    std::string topological_definition;
    std::vector<std::string> v_topological_definition;

    std::vector<std::string> precuts;



	bdt_info(){identifier = "null"; name = "null"; binning = "null";};
	
//    bdt_info(std::string in_identifier, std::string in_name, std::string in_bin) : identifier(in_identifier), name(in_name), binning(in_bin), base_cuts("1"), pre_cuts("1"), mid_cuts("1"),post_cuts("1"), bdt_cosmic_cuts("1"), bdt_bnb_cuts("1"),signal_definition("1"),background_definition("1"){
//	topo_name = "2#gamma1p";	
//	};
	
		//to be removed, CHECK
		bdt_info(std::string analysis_tag, method_struct inmethod):  TMVAmethod(inmethod), binning(inmethod.bdt_binning){ 
			topo_name = TMVAmethod.topological_name;	
			train_vars = TMVAmethod.bdt_train_vars;
			//    spec_vars = TMVAmethod.bdt_spec_vars;
		}

	void SetDefaultAttributes(){
		base_cuts = "1";
		pre_cuts = "1";
		mid_cuts = "1";
		post_cuts = "1";
		bdt_cosmic_cuts = "1";
		bdt_bnb_cuts = "1";
		signal_definition = "1";
		background_definition = "1";
	}

	int setName(std::string in){ identifier = in;return 0;};
	int setBaseCuts(std::string in){ base_cuts = in; return 0;}; 
	int setPreCuts(std::string in){ pre_cuts = in; return 0;}; 
	int setMidCuts(std::string in){ mid_cuts = in; return 0;}; 
	int setPostCuts(std::string in){ post_cuts = in; return 0;}; 
	int setBdtCosmicCuts(std::string in){ bdt_cosmic_cuts = in; return 0;}; 
	int setBdtBnbCuts(std::string in){ bdt_bnb_cuts = in; return 0;}; 

	int setTopoName(std::string in){ topo_name = in; return 0;};
	
	int setSignalDefinition(std::string in){ signal_definition = in; return 0;}; 
	int setBackgroundDefinition(std::string in){ background_definition = in; return 0;}; 

};

#endif
