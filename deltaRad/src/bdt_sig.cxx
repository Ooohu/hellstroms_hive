#include "bdt_sig.h"



/*
How do I normalize..

Two cuts( on this single cut)

INTIME
RECO2 has 1013198 events in 31497 files
DETSIM has 1030055 events in 31980 files
I ended with 989875 events


BNBCOSMIC
DETSIM has 2412300 in 48246 files
RECO2 should have 2360950 events
I ended with 2360950 events //ooohh all of them 

Weight each intime cosmic event with a factor 10.279*N_gen_BNB/(N_gen_cosmic*my_rate)
= 10.279*2412300/(1030055*989875/1013198) = 24.639718178714663

times whatever POT scaling we need to put on the BNB events to get to 6.6e20
which for v3.0_with calo is 2.38091e+21
= 24.6397*6.6e20/2.38091e21  = 6.830246418386248

	precut--
	for loop over BDT1
	for loop over BDT2
		calc S/sqrt(S+S+BKG)

*/


std::vector<double> scan_significance(TFile * fout, std::vector<bdt_file*> sig_files, std::vector<bdt_file*> bkg_files, bdt_cuts cosmic_cut, bdt_cuts bnb_cut){
	std::cout<<"Starting to Scan Significance"<<std::endl;
	double best_significance = 0;
	double best_mva_cut = DBL_MAX;
	double best_mva_cut2 = DBL_MAX;

	double plot_pot = 6.6e20;
	

	//for nice plots make the 50, 25 is quicker tho
	int nsteps_cosmic = 20;//50
	double cut_min_cosmic = 999;
	double cut_max_cosmic = -999;

	int nsteps_bnb = 20;//50
	double cut_min_bnb = 999;//0.52;
	double cut_max_bnb = -999;



	for(size_t i = 0; i < sig_files.size(); ++i) {
	//	double tmin_cos = sig_files.at(i)->tvertex->GetMinimum( (sig_files.at(i)->getBDTVariable(cosmic_cut).name + ">0").c_str()    );
		double tmax_cos = sig_files.at(i)->tvertex->GetMaximum( sig_files.at(i)->getBDTVariable(cosmic_cut).name.c_str()    );
		double tmax_bnb = sig_files.at(i)->tvertex->GetMaximum( sig_files.at(i)->getBDTVariable(bnb_cut).name.c_str()    );

		//if( tmin_cos <= cut_min_cosmic) cut_min_cosmic=tmin_cos;
		if( tmax_cos >= cut_max_cosmic) cut_max_cosmic=tmax_cos;
		if( tmax_bnb >= cut_max_bnb) cut_max_bnb=tmax_bnb;

	}
	cut_min_cosmic = cut_max_cosmic*0.82;
	cut_min_bnb = cut_max_bnb*0.85;

	std::cout<<"BNB sig scan from: "<<cut_min_bnb<<" to "<<cut_max_bnb<<std::endl;
	std::cout<<"COSMIC sig scan from: "<<cut_min_cosmic<<" to "<<cut_max_cosmic<<std::endl;

	double step_cosmic = (cut_max_cosmic-cut_min_cosmic)/((double)nsteps_cosmic);
	double step_bnb = (cut_max_bnb-cut_min_bnb)/((double)nsteps_bnb);


	TH2D * h2_sig_cut = new TH2D( "significance_2D",  "significance_2D",nsteps_cosmic, cut_min_cosmic, cut_max_cosmic, nsteps_bnb, cut_min_bnb, cut_max_bnb);
	std::vector<double> vec_sig;//some vectors to store TGraph info;
	std::vector<double> vec_cut;	

	for(int di=1; di<=nsteps_cosmic; di++) {
		double d  = (double)(di-1.0)*step_cosmic + cut_min_cosmic; ;	

		for(int di2=1; di2<=nsteps_bnb; di2++) {
			double d2  = (double)(di2-1.0)*step_bnb + cut_min_bnb ;	

			double signal = 0;
			double background = 0;
			std::vector<double> bkg;	

			for(size_t i = 0; i < sig_files.size(); ++i) {
				double pot_scale = (plot_pot/sig_files.at(i)->pot )*sig_files.at(i)->scale_data;
			
				std::string bnbcut = sig_files.at(i)->getStageCuts(3,d,d2); 
				signal += sig_files.at(i)->tvertex->GetEntries(bnbcut.c_str())*pot_scale;

			}

			for(size_t i = 0; i < bkg_files.size(); ++i) {
				double pot_scale = (plot_pot/bkg_files.at(i)->pot)*bkg_files.at(i)->scale_data;
		
	
				std::string bnbcut = bkg_files.at(i)->getStageCuts(3,d,d2); 
				bkg.push_back(	bkg_files.at(i)->tvertex->GetEntries(bnbcut.c_str())*pot_scale);			

				background += bkg.back();
			}
			double significance =0;
			if(signal==0){
				 significance =0;
			}else if(background !=0){
				significance = signal/sqrt(background);
			}else{
				std::cout<<"method_best_significane_seperate || signal2+background2 == 0, so significance  = nan @ cut1: "<<d<<", cut2: "<<d2<<std::endl;
				break;
			}


			if(significance > best_significance) {
				best_significance = significance;
				best_mva_cut = d;
				best_mva_cut2 = d2;
			}


			std::cout<<d<<" "<<d2<<" "<<significance<<" #signal: "<<signal<<" #bkg: "<<background<<" || "<<" bnb: "<<bkg.at(0)<<" cos: "<<bkg.at(1)<<std::endl;
			vec_sig.push_back(significance);
			vec_cut.push_back(d2);
			h2_sig_cut->SetBinContent(di,di2, significance);
		}

	}
		

	h2_sig_cut->SetStats(false);
	TCanvas * c_sig_cuts =  new TCanvas( "significance_cuts_colz", "significance_cuts_colz", 2000,1600 );
	c_sig_cuts->Divide(2,1);
	TPad *p1 = (TPad*)c_sig_cuts->cd(1);
	p1->SetRightMargin(0.13);
	h2_sig_cut->Draw("colz");
	h2_sig_cut->GetXaxis()->SetTitle("Cosmic Cut");
	h2_sig_cut->GetYaxis()->SetTitle("BNB Cut");
	
   	std::vector<double> vec_bf_cut1 = {best_mva_cut};
   	std::vector<double> vec_bf_cut2 = {best_mva_cut2};
	TGraph *graph_bf = new TGraph(vec_bf_cut1.size(), &vec_bf_cut1[0], &vec_bf_cut2[0]);
	graph_bf->SetMarkerStyle(29);
	graph_bf->SetMarkerSize(2);
	graph_bf->SetMarkerColor(kBlack);
	graph_bf->Draw("same p");
 

	TGraph * graph_cut = new TGraph(vec_sig.size(), &vec_cut[0], &vec_sig[0]);
	graph_cut->SetTitle("1D slices");
	c_sig_cuts->cd(2);
	graph_cut->Draw("alp");

	h2_sig_cut->Write();
	c_sig_cuts->Write();

	return std::vector<double>{best_mva_cut, best_mva_cut2, best_significance};

}
