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

std::string gadget_MakeSafeName(std::string safe_name);

struct bdt_variable{

	public:
		//raw from xml files;
		std::string name;//variable name
		TString mininame;//CHECK, has disable this;
		std::string binning;//default binning string
		std::string unit;//axis unit
		std::string type;//CHECK, useless
        std::string covar_file;//systematic cov. location

		//optional
        std::string covar_legend_name;//for adjusting legend name
        bool is_logplot;//make a log plot if true

		//derivated members
        int id;
		int cat;//group id
        std::string safe_name;//tag for output
		std::string safe_unit;//tag for output
		bool is_track;//CHECK, useless
        bool has_covar;//has (not) systematics
		bool is_custombin;//true - use variable binnings
//      std::string covar_name;

        double plot_min;//left edge;
        double plot_max;//right edge;
		double plot_height;//NEW, constrainied max bin height of plots;
		
        int n_bins;//bins for equal binnings;
		int int_n_bins;//use this as initial bin# when using variable binning;

        std::vector<double> edges;//update this;
        std::vector<std::string> edges_str;//NEW, edges in strings;

		

//        int addCovar(std::string name, std::string file){}
        int addCovar(std::string file){
            has_covar=true;
//            covar_name = name;
            covar_file = file;
            return 0;
        }

//		bdt_variable(std::string inname, std::string inbin, std::string inunit,bool intrack, std::string intype) : bdt_variable(inname,inbin,inunit,intrack,intype,-1){}; 
		
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
		safe_name = gadget_MakeSafeName(name);
		safe_unit = gadget_MakeSafeName(unit);

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
					plot_min = (double)std::stod(token);
					break;
				case 3:
					plot_max = (double)std::stod(token);
					break;
			}

			edges.push_back(std::stod(token));
			edges_str.push_back(token);
			bins.erase(0, pos + delim.length());
		}

		//edges.push_back(std::stod(bins));

		//            cat = 0;
	};



//		bdt_variable(){};
		//short-handed constructor?
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
