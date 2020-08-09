#ifndef BDT_SYSTEMATICS_H
#define BDT_SYSTEMATICS_H

#include <iostream>

#include "bdt_var.h"
#include "bdt_file.h"
#include "load_mva_param.h"

#include "TStyle.h"
#include "TSystem.h"

struct bdt_sys{
	public:
		TString tag;// BkgMC, dirt | MCUnism, Multism,OpticalModel
		TString dir;//for multi files, need to know the parent directory.
		TString filename;
		TString treename;


		std::vector< TString > vars;
		std::vector< TString > vars_name;


		double pot;
		double throws;
		
		bool its_CV;
		bool its_multithrows;
		bool its_OM;

		//deduced elements
		std::vector< std::vector<TH1F*> > hist;
		std::vector< std::vector<TH2D*> > twodhist;
		int num_vars;
		int start_with_weight = 0;//a label for which weight to start, this is for the case that the code was interrupted

		//constructor
		bdt_sys(int index, MVALoader XMLconfig);
		~bdt_sys();
};

/*
 * Initialize Systematics of given variable 1dhistogram of weights;
 * This will add contents into syss;
 */

void InitSys(std::vector<bdt_variable> var, std::vector<std::string> precuts, std::vector<bdt_sys*> syss, double plot_pot, TString dir_root, TString dir_drawn);

/*
 * Make 1d histograms;
 */
void Make1dhist(TFile* hist_root, bdt_variable* var, TString event_cuts, bdt_sys* temps, double plot_pot);


/*
 * make covariance matrix according to the histograms;
 *
 */

void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss, TString dir_root, TString dir_drawn, double plot_pot);

/*
 * make covaraince matrix according to input histograms, hist - weights and cv - cv;
 *
 */
TH2D* MakeCov(TString name,TH1F* hist, TH1F* cv);

/*
 * covariane matrix propagation, create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TString name,TH2D* cov_temp, TH1F* oldcv, TH1F* newcv);

/*
 * Make Covaraince matrix from 2d histograms
 *
 */


TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv);

#endif
