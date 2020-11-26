#ifndef BDT_VAR_H
#define BDT_VAR_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <ctype.h>
//#include <typeinfo>
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

void gadget_loadV( std::vector<int>& out_v, std::string in);
void gadget_loadV( std::vector<double>& out_v, std::string in);
void gadget_loadV( std::vector<std::string>& out_v, std::string in);
void gadget_loadV( std::vector<TString>& out_v, std::string in);

//Execption: template has to be in the header, but other funciton cant be defined here..
//see: https://stackoverflow.com/questions/1111440/undefined-reference-error-for-template-method
//--
//template usage reference: https://stackoverflow.com/questions/17649136/function-which-is-able-to-return-different-types
//T (for now) doule be std::string, TString, int, double
//
//This function return string to a vector of T. usage (for int) : gadget_Tokenizer< int> ("str");
template<typename T > 
std::vector< T> gadget_Tokenizer( std::string input){
	//exp:
	//(0,1, spit out {"0","1"}
	//var1[0],var2[0] spit out {"var1[0]", "var2[0]"}
	
	input.erase(std::remove(input.begin(), input.end(), '('), input.end());
	input.erase(std::remove(input.begin(), input.end(), ')'), input.end());
	std::replace(input.begin(), input.end(), ',', ' ');

	std::istringstream iss(input);

	std::vector<T> output;

	std::string one_str;
	//use overloaded function to fill the vector;
	while( iss>> one_str) gadget_loadV( output, one_str);
	
	return output;
};

struct bdt_variable{

	public:
		//raw from xml files;
		std::string name;//variable name
		std::string binning;//default binning string
		std::string unit;//axis unit

		//optional
		int id;
		std::string covar_file;//systematic cov. location
		std::string covar_legend_name;//for adjusting legend name
		std::vector<int> cats;//group id
		bool is_train;
		bool is_logplot;//make a log plot if true
		bool is_custombin;//true - use variable binnings
		double plot_height;//NEW, constrainied max bin height of plots;

		//derivated members
		std::string safe_name;//tag for output
		std::string safe_unit;//tag for output
		bool has_covar;//has (not) systematics

		double plot_min;//left edge;
		double plot_max;//right edge;

		int int_n_bins;//use this as initial bin# when using variable binning;
		double bin_gap;//final bin numbers
		int n_bins;//final bin numbers

		std::vector<double> edges;//update this;
		std::vector<std::string> edges_str;//NEW, edges in strings;



		bdt_variable(std::string inname, std::string inbinning, std::string inunit) 
		{
			reset_var();
			name = inname;
			binning = inbinning;
			unit = inunit;

			safe_name = gadget_MakeSafeName(name);
			safe_unit = gadget_MakeSafeName(unit);
		}



		void load_bininfo()
		{

			std::vector<std::string> bin_info = gadget_Tokenizer<std::string> ( binning);
			std::vector<double> bin_info_num = gadget_Tokenizer<double> ( binning);

			if(bin_info.size() < 3){ 
				std::cout<<"Binning "<< binning<<" is incorrect, please check "<<std::endl;
				exit(EXIT_FAILURE);
			}

			plot_min = bin_info_num[1];//(double) std::stod( bin_info[1] );
			plot_max = bin_info_num.back();

			if(is_custombin){//bin_info is ( bingap, edge1, edge2,...)
				int_n_bins =  (plot_max-plot_min) / bin_info_num.front() ;
				n_bins = bin_info.size() - 2;
//				std::cout<<"Variable binning "<<std::endl;

			} else{//bin_info is ( nbins, left_edge,right_edge)
				int_n_bins = (int) bin_info_num.front();
				n_bins = int_n_bins;
			}
			bin_gap = (double) ( (plot_max-plot_min)/n_bins);

//			std::cout<<" Bins "<<n_bins <<" max "<<plot_max<<" min "<<plot_min<<" gap "<< bin_gap<<std::endl;

			edges.clear();
			edges_str.clear();

			for(int index = 0; index < n_bins+1; index++){
				if(is_custombin){//( bingap, edge1, edge2,...)
					edges_str.push_back( bin_info[index+1] );

				} else{//( bingap, edge1, edge2,...)
					double edge = plot_min + bin_gap * index;
					std::ostringstream edge_str;
					edge_str << edge;
					edges_str.push_back( edge_str.str());
				}
				edges.push_back( (double) std::stod( edges_str.back() ) );
//				std::cout<<index<<" add edge "<<edges.back()<<std::endl;
			}
//			std::cout<<" finish "<<std::endl;
		};



		void reset_var(){
			name = "EMPTY";//variable name
			binning = "EMPTY";//default binning string
			unit = "EMPTY";//axis unit
			covar_file = "EMPTY";//systematic cov. location
			id = -1;

			covar_legend_name = "EMPTY";//for adjusting legend name
			cats = {-1};//group id
			is_train = false;
			is_logplot = false;//make a log plot if true
			is_custombin = false;//true - use variable binnings
			plot_height = 0;//NEW, constrainied max bin height of plots;

			safe_name = "EMPTY";//tag for output
			safe_unit = "EMPTY";//tag for output
			has_covar = false;//has (not) systematics

			plot_min = -999;//left edge;
			plot_max = -999;//right edge;

			int_n_bins = 0;//use this as initial bin# when using variable binning;
			bin_gap = 0;//final bin numbers
			n_bins = 0;//final bin numbers

			edges = {0};//update this;
			edges_str = {"0"};//NEW, edges in strings;
		}

};





#endif
