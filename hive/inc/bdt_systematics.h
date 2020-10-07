#ifndef BDT_SYSTEMATICS_H
#define BDT_SYSTEMATICS_H

#include <iostream>

#include "bdt_var.h"
#include "bdt_file.h"
#include "load_mva_param.h"

#include "TStyle.h"
#include "TSystem.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TF1.h"

//NOTE, train/app part of the hive.cxx  need to be updated, such that the BDT cuts can be used on systematics;

class bdt_sys;//more detail follows later

class sys_env{

	public:

		/*
		 * Initialize Systematics of given variable 1dhistogram of weights;
		 * This will add contents into syss;
		 */
		TString top_dir;
		TString drawn_dir;
		TString root_dir;
		double out_POT;

		int fstage;
		std::vector<double> fbdt_cuts;
		TString fcut_hash;

		TString hist_prefix;//raw histograms for making covaraicne matrices;
		TString cov_prefix;//all (fractional) covariane maatrices;
		TString final_prefix;//fractional matrices for plotting
		
		bool verbose;
		bool debug_verbose;

		//filled on the flight
		std::vector<TString> tag_collection;//tags = systags[?]+"_"+bdtftags[?];
		std::multimap< TString, bdt_sys* > tag2SWmap;
		std::map< TString, bdt_sys* > tag2CVmap;

		//set directories and POT;
		void setEnv(TString input_top_dir, TString subdir_root, TString subdir_drawn){
				top_dir = (input_top_dir+"systematics/");
				root_dir =  top_dir+"/"+subdir_root+"/";
				drawn_dir = top_dir+"/"+subdir_drawn+"/";
		};

		void catchupEnv(TString in_top_dir, TString indir_root, TString indir_drawn, double cur_pot, int cur_stage, TString cur_hash, std::vector<double> cur_bdt_cuts){//for son class, i.e. bdt_sys, to catch up info.
			top_dir = in_top_dir;
			root_dir = indir_root;
			drawn_dir = indir_drawn;
			out_POT = cur_pot;
			fstage = cur_stage;
			fcut_hash = cur_hash;
			fbdt_cuts = cur_bdt_cuts;
		};

		//get fstage & cut_hash filled.
		void setStageHash(int cur_stage, bdt_file* file, std::vector<double> input_bdt_cuts){
			fstage = cur_stage;
			if(fstage < 0 ){
				std::cout<<"Please select a stage using -s. Abort."<<std::endl;
				exit(EXIT_FAILURE);
			}
			std::string cuts =  file->getStageCuts(fstage, input_bdt_cuts);
			//see the function at line 951 of bdt_file.cxx
			unsigned long calHash = file->jenkins_hash(cuts);
			fcut_hash = "h"+ to_string_prec(calHash,0);

			fbdt_cuts = input_bdt_cuts;
		};


		int checkEnv();


		TString getFileName(TString plot_type, std::vector<bdt_variable> var);

		//make 1dhist or load 1dhist; then make cov matrices;
		void InitSys(std::vector<bdt_variable> var, std::vector<bdt_sys*> syss);
		/*
		 * make covariance matrix according to the histograms;
		 */
		void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss);


		void setVerbose(int lv){	

			verbose = false;
			debug_verbose = false;

			switch(lv){
				case 2:
					debug_verbose = true;
				case 1:
					verbose = true;
			}
		};
};


/*
 * set-ups for each systematic files (1 bdt_file x1 particulara systematics)
 *
 */
class bdt_sys : public sys_env, public bdt_file{
	//so bdt_sys is now have the public functions(), elements from bdtfile;
	//we need the following functions from bdt_sys:
	//-- TString bdt_file::getStageCutsIndex(int fstage, std::vector<double> bdt_cuts, int vec_index);
	//
	//and elements:
	// tag
	// vec_entry_lists
	
	public:
		//From the XML;
		TString tag;// BkgMC, dirt | MCUnism, Multism,OpticalModel _ dirt, pi0misd, delta..
		TString systag;//BkgMC, dirt | MCUnism, Multism,OpticalModel

//		TString dir;//for multi files, need to know the parent directory.
		TString dir;//directory of the input systematic file
		TString filename;
		TString treename;

		std::vector< TString > vars;
		std::vector< TString > vars_name;

		double pot;
		double throws;
		
		bool its_CV;
		bool its_multithrows;
		bool its_OM;
		bool fullyloaded;
		std::vector< bool > loads;

		//deduced elements
		std::vector< TString> TdirNames;
		std::vector< std::vector<TH1F*> > hists;
		std::vector< std::vector<TString> > histNames;
		std::vector< std::vector<TH2D*> > twodhists;
		std::vector< std::vector<TString> > twodhistNames;
		std::vector< TMatrixD > CVhists;
		std::vector< TMatrixD > FracMatrices;

		int num_vars;

//		int start_with_weight = 0;//a label for which weight to start, this is for the case that the code was interrupted

		//constructor, now each bdt_sys contains one bdt_file; i.e. 2 systematics will have 2xfiles# bdt_sys classes;
		bdt_sys(int index, MVALoader XMLconfig, int bdtfile_index, bdt_flow inflow);
		~bdt_sys();


		/*
		 * Make 1d histograms;
		 * this fills hist;
		 */
//		void Make1dhist(TString histfilename, bdt_variable* var);//,  double plot_pot,std::vector<double> bdt_cuts);
		void Make1dhist(bdt_variable *var);//,  double plot_pot,std::vector<double> bdt_cuts);

//		void Load1dhist(TString histfilename, bdt_variable* var);
		bool Load1dhist(TFile* cur_file);


};

/*
 * make covaraince matrix according to input histograms, hist - weights and cv - cv;
 *
 */
TH2D* MakeCov(TString name,TH1F* hist, TH1F* cv);

/*
 * covariane matrix propagation, create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TString name,TH2D* cov_temp, TH1F* oldcv, bool, TH2D*);

/*
 * Propagate matrix
 */
TH2D* PropagateCov(TString name,TH2D* frac_cov, TH1F* newcv);



/*
 * Smooth the matrix
 */
TH1F SmoothSW(TH1F* sw, TH1F* cv, bool special);

/*
 * Make Covaraince matrix from 2d histograms
 *
 */
TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv);

#endif
