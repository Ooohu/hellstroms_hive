#ifndef BDT_VAR_H
#define BDT_VAR_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
/******** Our includes *****/

/******** Root includes *****/

#include "TTreeFormula.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMVA/Types.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TFriendElement.h"


struct bdt_variable{

	public:
		std::string name;//variable name
		TString mininame;
        int id;
		int cat;
        std::string safe_name;
		std::string binning;
		std::string unit;
		std::string safe_unit;
		bool is_track;
		std::string type;
        bool is_logplot;
        bool has_covar;
//        std::string covar_name;
        std::string covar_file;
        std::string covar_legend_name;

        double plot_min;
        double plot_max;
		
        int n_bins;
		int int_n_bins;//use this as initial bin# when using variable binning;

        std::vector<double> edges;

		bool is_custombin;
		

//        int addCovar(std::string name, std::string file){}
        int addCovar(std::string file){
            has_covar=true;
//            covar_name = name;
            covar_file = file;
            return 0;
        }

		bdt_variable(std::string inname, std::string inbin, std::string inunit,bool intrack, std::string intype) : bdt_variable(inname,inbin,inunit,intrack,intype,-1){}; 
		
        bdt_variable(std::string inname, std::string inbin, std::string inunit,bool intrack, std::string intype,int in_id) : 
			name(inname), 
			binning(inbin),
			unit(inunit),
			is_track(intrack),
			type(intype),
            id(in_id),
            is_logplot(false)
	{//WARNING!! This is the orignal version, use the overflow one below!
		plot_min =-999;
		plot_max =-999;
		safe_name = name;
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

		safe_unit = unit;
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), ' '), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '('), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), ')'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '\\'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '/'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '['), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), ']'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '+'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '-'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '*'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '|'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), ':'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '#'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '^'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '{'), safe_unit.end());
		safe_unit.erase(std::remove(safe_unit.begin(), safe_unit.end(), '}'), safe_unit.end());

		has_covar = false;

		std::string bins = binning;
		edges.clear();


		bins.erase(std::remove(bins.begin(), bins.end(), '('), bins.end());
		bins.erase(std::remove(bins.begin(), bins.end(), ')'), bins.end());

		size_t pos = 0;
		std::string delim = ",";
		bins = bins + delim+"0";//Keng yea, need this to read the plot_max below
		std::string token;
		n_bins = -1;

		int ith_number = 1;
		while ((pos = bins.find(delim)) != std::string::npos) {
			token = bins.substr(0, pos);
//			if(n_bins<0) n_bins = (int)std::stod(token);//Keng commented this
			switch(ith_number++){
				case 1:
					n_bins = (int)std::stod(token);
					break;
				case 2:
					plot_min = (int)std::stod(token);
					break;
				case 3:
					plot_max = (int)std::stod(token);
					break;
			}

			edges.push_back(std::stod(token));
			bins.erase(0, pos + delim.length());
		}

		//edges.push_back(std::stod(bins));

		//            cat = 0;
	};



		bdt_variable(){};

		bdt_variable(std::string inname, std::string inbin, std::string inunit,bool intrack) : 
			name(inname), 
			binning(inbin),
			is_track(intrack),
			unit(inunit)
	{ 
		plot_min =-999;
		plot_max =-999;
		//            cat = 0;
	};

};





#endif
