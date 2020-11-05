#include "bdt_varplot.h"

double gadget_getchi2(TH1* obser, TH1* hypo, TMatrixD cov){

	bool debug_message = true;
	bool include_sys = false;
	double mychi = 0;
	
	int bins = cov.GetNrows();


		std::vector<std::stringstream> chiMap(bins);

		TMatrixD covInverse = cov.Invert();
		for(int ib=1; ib<bins+1; ib++){
			for(int jb=1; jb<bins+1; jb++){//use Matrices_set to evaluate the chi^2;
				double dav_i = obser->GetBinContent(ib);//the obser
				double mcv_i = hypo->GetBinContent(ib);

				double dav_j = obser->GetBinContent(jb);//the obser
				double mcv_j = hypo->GetBinContent(jb);

				double numerator = (dav_i - mcv_i)*(dav_j - mcv_j);
				double inv_M = covInverse(ib-1,jb-1);
				double curchi = numerator*inv_M;
				if(ib==1 && debug_message){
//					std::cout<<"Chi2 Bini "<<ib<<", excess "<<dav_i<<", err "<<obser->GetBinError(ib)<<",MC "<<mcv_i<<",err "<<hypo->GetBinError(ib)<<",inv_Matrix "<<covInverse(ib-1,jb-1)<<std::endl;
					std::cout<<"Chi2 Binj "<<jb<<", excess "<<dav_j<<", err "<<obser->GetBinError(jb)<<",MC "<<mcv_j<<",err "<<hypo->GetBinError(jb)<<",inv_Matrix "<<inv_M<<std::endl;
					std::cout<<"\t\t ---> Resulting chi2 "<<curchi<<std::endl;
				}
				chiMap[ib-1]<<std::setw(12)<<curchi<<" ";
				mychi+=curchi;
			} 
		}

		if(debug_message){
			std::cout<<"ChiSq values are "<<std::endl;
			for(int printdex = 0; printdex < bins; ++printdex){
				std::cout<<chiMap[printdex].rdbuf()<<std::endl;
			}
		}
	return mychi;
}


double gadget_getlikelihood(TH1* obser, TH1* hypo, TMatrixD cov){

	bool debug_message = true;
	bool include_sys = false;
	double ML = 1;
	double log_norm = 0;
	double temp_chi2 = 0;
	
	int bins = cov.GetNrows();

//	if(include_sys){

//		TMatrixD covInverse = cov.Invert();
	for(int ib=1; ib<bins+1; ib++){
		double ErrSq = cov(ib-1,ib-1);
		double curml = 1.0/sqrt(2*3.14*ErrSq);
		double chi2i = pow(obser->GetBinContent(ib)-hypo->GetBinContent(ib),2)/(ErrSq);
		temp_chi2+=chi2i;
		curml *= exp( - chi2i/2);
		ML*=curml;	
		log_norm += log(2*3.14*ErrSq);
	}

	if(debug_message){
		std::cout<<"ChiSq from Log(L): "<<2*(-0.5*log_norm-log(ML))<<std::endl;
		std::cout<<"Should be equal to :"<<temp_chi2<<std::endl;
	}
	return ML;
}

//void gadget_addDoubleAxis(v, "pow(x,3)" , maxy*yaxis_factor);

void gadget_addDoubleAxis(bdt_variable var, TString express, double y_coord){
//	std::cout<<"Draw extra x-axis"<<std::endl;
	double minv = var.plot_min;
	double maxv = var.plot_max;

	//						double mapped_maxv = pow(maxv,(1/3.0))*500;
	TF1 *temp_axis = new TF1("temp_axis",experss,0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;
	//						TF1 *temp_axis = new TF1("temp_axis","pow(x,3.0)",0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;

	TGaxis *new_axis = new TGaxis(minv, y_coord, maxv, y_coord, "temp_axis", 403, "-");//This says where to draw axis (xmin,ymin,xmax,ymax,TF1,num_division =N1 + 100*N2 + 10000*N3 [by default, ROOT will take numbers not larger than these to optimize the axis] ,option"-");https://root.cern.ch/doc/master/classTGaxis.html
	new_axis->SetTitle("Radius [cm]");
	new_axis->SetLabelSize(0.03);
	new_axis->SetTitleSize(0.03);
	new_axis->SetTitleOffset(1.3);
	new_axis->Draw();
	//std::cout<<"min "<<minv<<" max "<<maxv<<" maxy "<<maxy<<std::endl;
}

std::vector<double> gadget_getYranges(std::vector<TH1*> hists, TH1* hist){

	double maxy = 0;
	double miny = 0;

	hists.push_back(hist);

	for( auto hs : hists){

		double maxy_of_hs = hs->GetBinContent(hs->GetMaximumBin())+hs->GetBinError(hs->GetMaximumBin());
		if( maxy < maxy_of_hs ) maxy = maxy_of_hs;

		double miny_of_hs = hs->GetBinContent(hs->GetMinimumBin())-hs->GetBinError(hs->GetMaximumBin());
		if( miny > miny_of_hs ) miny = miny_of_hs;

	}
	std::vector<double> output = {miny, maxy};
	return output;

}


//int plot_bdt_variables(bdt_file * signal_pure, bdt_file * background_pure, std::vector<bdt_variable> vars, bdt_info input_bdt_info, bool isSpectator,int stage,std::vector<double> cuts);
int  plot_var_allF(std::vector< bdt_file *> MCfiles, bdt_file* datafile, std::vector<bdt_variable> vars, std::vector<bdt_sys*> syss, int stage, std::vector<double> bdtcuts, bool isNorm){
//    std::vector<std::string> title = {"Topological Selection","Pre-Selection Cuts","Cosmic","BNB","Pi0","Other"};
		
		//isNorm true -  scale MC to the excess;
		
		int num_MC = MCfiles.size();
		double is_bestfit_in = false;//false - not include signal in the output;
//		bool do_pair_plot = true;//true - excess vs each MC;
		bool draw_estimator = true;
		bool debug_message = true;

		double yaxis_factor = 1.2;
		TString legend_title = "MiniBooNE Preliminary";

		std::cout<<"on stage : "<<stage<<"\n"<<"Setting sig stage entry lists."<<std::endl;
        if(stage>1){
            datafile->calcBDTEntryList(stage,bdtcuts);
        }

		for(size_t index = 0; index <num_MC; ++index){
			if(stage>1) MCfiles[index]->calcBDTEntryList(stage,bdtcuts);
			MCfiles[index]->setStageEntryList(stage);
		}
        std::cout<<" Setting back stage entry lists."<<std::endl;
        datafile->setStageEntryList(stage);
        std::cout<<" Set stage entry lists."<<std::endl;

        int nv = 0;

		for(auto &v: vars){

			std::cout<<"On variable: "<<v.name<<std::endl;

			//Prepare data first, modified it according to MC and draw 
			TH1* data_original = datafile->getTH1(v,"1",v.safe_name+"_data_var" ,datafile->pot);
			TH1* excess = (TH1*) data_original->Clone();

			//STEP 1: Prepare MC
			std::vector<TH1*> MC(num_MC);
			TH1* AllNue_MC;
			bool first_hist = true;
			bool using_sys = v.has_covar;
			int nbins = v.n_bins;

			for(size_t index = 0; index <num_MC; ++index){//loop through each NueMC; do the substraction;

				MC[index] = MCfiles[index]->getTH1(v,"1",v.safe_name+"_MC_var" , datafile->pot);
				
				MC[index]->SetLineColor((MCfiles[index]->is_signal)? kBlack : MCfiles[index]->col);
				MC[index]->SetLineWidth(2);

//				MC[index]->SetFillColor(MCfiles[index]->col);
				MC[index]->SetFillStyle(3445);


				if(is_bestfit_in || !MCfiles[index]->is_signal ){//NueMC must go through, Best-fit may not;

					if(first_hist){//copy the first MC hist
						AllNue_MC = (TH1*) MC[index]->Clone("AllNuemc");
						first_hist = false;
					} else{//add to the current MC hist;
						AllNue_MC->Add(MC[index],1.0);

					}
				}
			}//next MC[index]; 

			excess->Add(AllNue_MC,-1.0);//Calculate Data - MC


			//STEP 2: Adjust errors & MC[index], prepare estimator;
			//STEP 2.1: Introduce error matrix, include systematic error (or not) in the excess;
			std::vector< TMatrixD > Matrices_set(3, TMatrixD(nbins,nbins));
			if(using_sys){
				TFile *covar_f = new TFile(v.covar_file.c_str(),"read");
				TMatrixD * covar_collapsed = new TMatrixD(v.n_bins,v.n_bins);
				*covar_collapsed = gadget_PrepareMatrix(syss, covar_f, AllNue_MC, datafile->pot, datafile->tag);
				Matrices_set = gadget_SeparateMatrix( covar_collapsed, AllNue_MC, "vars/");//0: shape-only, 1: mixed, 2: normalization-only
				covar_f->Close();
			}

			TMatrixD M_MixedShape_orig = Matrices_set[0]+Matrices_set[1];
			TMatrixD M_orig = Matrices_set[0]+Matrices_set[1] + Matrices_set[2];

			//add AllNue statError to the diagonal;
			for(int ib=1; ib<nbins+1; ib++){
//				if(debug_message)std::cout<<"\tAll NueMC bin: "<<ib<<" has "<<AllNue_MC->GetBinContent(ib)<< " with statE "<< AllNue_MC->GetBinError(ib) << " sys "<< sqrt(mii)<<" data statE "<<excess->GetBinError(ib)<<std::endl;
				double cur_error = AllNue_MC->GetBinError(ib);
				M_MixedShape_orig(ib-1, ib-1) +=  pow(cur_error,2);
				M_orig(ib-1, ib-1) +=  pow(cur_error,2);
				if(debug_message) std::cout<<"\t Add AllNue Stat Error to the diagonal (Bin,Error) = ("<<ib<<","<<cur_error<<")"<<std::endl;
			}


			//STEP 2.2: Scale MC[index] and Prepare estimators;
			std::vector< double > mc_scale_factor(num_MC);
			std::vector<std::string > estimators(num_MC);

			for(size_t index = 0; index <num_MC; ++index){

				mc_scale_factor[index] = excess->Integral()/MC[index]->Integral();
				if(isNorm) MC[index]->Scale( mc_scale_factor[index] );

				if(debug_message) std::cout<<"\n Evaluate statistical estimators of "<<MCfiles[index]->tag<<std::endl;
				TH1* cur_th1 = MC[index];
				//get ks
				TString ks = "KS: "+to_string_prec(cur_th1->KolmogorovTest(excess), 3);
				TMatrixD Mms = M_MixedShape_orig;
				TMatrixD Mall = M_orig;

				for(int ib=1; ib<nbins+1; ib++){
					double addonvalues = pow(cur_th1->GetBinError(ib),2);
//					double addonvalues = (using_sys)? pow(cur_th1->GetBinError(ib),2) : pow(data_original->GetBinError(ib),2) + pow(cur_th1->GetBinError(ib), 2);
					Mms(ib-1,ib-1) += addonvalues;
					Mall(ib-1,ib-1) += addonvalues;
				}
				double mychi= gadget_getchi2(excess, cur_th1, Mms);
				double ML =  gadget_getlikelihood(excess, cur_th1, Mall);

				estimators[index] = ks +"#chi^{2}/n#it{DOF}: "+to_string_prec(mychi, 2) +"/"+to_string_prec(nbins-1)
					+", #chi^{2} P^{val}: "+to_string_prec(TMath::Prob(mychi, nbins - 1 ),3) //-1 is due to the normalization factor
					+ " Log(L): " + to_string_prec(log(ML),2);
				std::cout<<MCfiles[index]->tag<<" summary: "<<estimators[index]<<" Scale factor "<<mc_scale_factor[index]<<std::endl;
//				std::cout<<"ori data bin1 "<<data_original->GetBinError(1)<<std::endl;
//				std::cout<<"ori data bin1 "<<AllNue_MC->GetBinError(1)<<std::endl;
			}


			//STEP 3: Draw, only 1 - 1 comparison plots;
			//set max&min, and scale MC evnets if isNorm is true;

			//Canvas to be Drawn
			TCanvas *c_var = new TCanvas(("cvar_"+v.name).c_str(), ("cvar_"+v.name).c_str(),1200,1200);
			std::vector<double>  Yrange = gadget_getYranges(MC, excess);
			
			for(int index = 0; index < num_MC; ++index){

//				if(debug_message) std::cout<<"Plot "<<MCfiles[index]->tag<<std::endl;
				c_var->Clear();
				//adjust excess error bar
				TH1* cur_th1 = MC[index];
				for(int ib=1; ib<nbins+1; ib++){
					//AllMC Total Error +  MCi Error
					double datapErrSq = M_orig(ib-1,ib-1) + pow(cur_th1->GetBinError(ib),2);
					if(!using_sys) datapErrSq += pow(data_original->GetBinError(ib),2);
					
					excess->SetBinError(ib, sqrt(datapErrSq));//now data has AllNue sys + AllNue stat (+ data Stat)
//					std::cout<<"Bin "<<ib<<" statErr "<<cur_th1->GetBinError(ib)<<" data point Error "<<excess->GetBinError(ib)<<std::endl;
				}
				excess->SetTitle(" ");
				excess->SetMarkerSize(3);
				excess->SetMarkerStyle(20);
				excess->SetLineColor(kBlack);

				excess->GetYaxis()->SetTitle((isNorm)? "Events [MC Rescaled]" : "Events");
				excess->GetYaxis()->SetTitleOffset(1.5);
				excess->SetMaximum(Yrange[1]*yaxis_factor);
				excess->SetMinimum(Yrange[0]);

				excess->Draw("E1");

				cur_th1->Draw("hist same");
				

				if((v.unit).find("R/500")!= std::string::npos) gadget_addDoubleAxis(v, "pow(x,3)" , Yrange[1]*yaxis_factor);

				//Legend
				TLegend *legend = new TLegend(0.52,0.56,0.88,0.80);//x1,y1,x2,y2

				TString legend_label = (is_bestfit_in)? "Data - (NueMC+Best Fit)": "Data - NueMC";
				legend_label += using_sys? " w. Stat & Sys. Error" : "";//" w. Stat & Sys. Error";

				legend->AddEntry(excess, legend_label, "lp" );
				legend->AddEntry(MC[index], (MCfiles[index]->plot_name+" scaled:x"+to_string_prec(mc_scale_factor[index],2)).c_str(),"lf");

				legend->SetHeader(legend_title,"C");
				legend->SetBorderSize(0);
				legend->SetLineColor(kWhite);
				legend->SetFillStyle(0);
				legend->SetNColumns(1);

				legend->Draw();

				//Estimator
				TLatex *estimator = new TLatex(0.12, 0.85,estimators[index].c_str());
				estimator->SetNDC();
				estimator->SetTextColor(kRed-7); 
				estimator->SetTextSize(0.03);
				if(draw_estimator) estimator->Draw("same");


				//Conclude
				c_var->Print(("vars/"+std::to_string(nv)+"_"+v.safe_unit+"_"+MCfiles[index]->tag+"_stage_"+std::to_string(stage)+".pdf").c_str(),"pdf");

			}


			if(debug_message){
				std::vector< TH1*> pHists(MC);
				pHists.push_back(excess);
				std::vector<bdt_file*> pFiles(MCfiles);
				pFiles.push_back(datafile);
				std::cout<<"------------------- Summary -------------------"<<std::endl;
				gadget_summary( v, pFiles, pHists);
			}

			delete c_var;
			nv++;

		}

    return 0;
}


int  plot_bdt_variables(bdt_file * signal_pure, bdt_file * background_pure, std::vector<bdt_variable> vars, bdt_info input_bdt_info, bool isSpectator, int stage, std::vector<double> bdtcuts){

    std::vector<std::string> title = {"Topological Selection","Pre-Selection Cuts","Cosmic","BNB","Pi0","Other"};

    {
        int j = stage;

                std::cout<<"on stage : "<<j<<std::endl;
        std::cout<<" Setting sig stage entry lists."<<std::endl;
        if(j>1){
            signal_pure->calcBDTEntryList(j,bdtcuts);
            background_pure->calcBDTEntryList(j,bdtcuts);
        }


        signal_pure->setStageEntryList(j);
        std::cout<<" Setting back stage entry lists."<<std::endl;
        background_pure->setStageEntryList(j);
        std::cout<<" Set stage entry lists."<<std::endl;

        int nv = 0;
        for(auto &v: vars){

            //std::string cut_signal = signal_pure->getStageCuts(j,-9,-9); 
            //std::string cut_background_pure = background_pure->getStageCuts(j,-9,-9); 

            //		TH1* sig = signal_pure->getTH1(v,cut_signal.c_str(),v.safe_name+"_sig_var" ,1.0);
            //		TH1* bkg = background_pure->getTH1(v,cut_background_pure.c_str(),v.safe_name+"_bkg_var" ,1.0);
            TCanvas *c_var = new TCanvas(("cvar_"+v.name+"_"+input_bdt_info.identifier).c_str(), ("cvar_"+v.name+"_"+input_bdt_info.identifier).c_str(),1200,1200);
            c_var->cd();

            std::cout<<"On variable: "<<v.name<<std::endl;

            TH1* sig = signal_pure->getTH1(v,"1",v.safe_name+"_sig_var" ,1.0);
            TH1* bkg = background_pure->getTH1(v,"1",v.safe_name+"_bkg_var" ,1.0);

            std::cout<<"Integrals: "<<sig->Integral()<<" "<<bkg->Integral()<<std::endl;
    
            if(!TMath::Finite(sig->Integral())){
                std::cout<<"scanning"<<std::endl;
                signal_pure->tvertex->Scan(v.name.c_str());
            }else{
                std::cout<<"So its fnite"<<std::endl;
            }



            sig->Scale(1.0/sig->Integral());			
            bkg->Scale(1.0/bkg->Integral());			
            sig->SetLineColor(signal_pure->col);
            bkg->SetLineColor(background_pure->col);
            sig->SetLineWidth(2);
            bkg->SetLineWidth(2);
            c_var->cd();

            sig->SetFillColor(signal_pure->col);
            bkg->SetFillColor(background_pure->col);
            sig->SetFillStyle(3445);
            bkg->SetFillStyle(3454);

            if(j!=1){
                sig->SetTitle(title.at(j).c_str());
            }else{
                sig->SetTitle(" ");
            }

            c_var->cd();			

            sig->Draw("hist");
            sig->SetMinimum(0);
            bkg->Draw("hist same");
            //sig->GetXaxis()->SetTitle(v.unit.c_str());
            sig->GetYaxis()->SetTitle("Events [Area Normalized]");
            sig->GetYaxis()->SetTitleOffset(1.5);

            TLegend *l = new TLegend(0.11,0.75,0.89,0.89);
            l->SetLineColor(kWhite);
            l->SetFillStyle(0);
            l->SetNColumns(2);


            l->AddEntry(sig, signal_pure->plot_name.c_str(),"lf");	
            l->AddEntry(bkg, background_pure->plot_name.c_str(),"lf");	
            l->Draw();

            TText *pre;
            if (isSpectator) {
                pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulaton In-Progress - Spectator Variable");
            }else {
                pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulaton In-Progress - Training Variable");

            }


            //TText *pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulation - In Progress");
            pre->Draw();

            TLatex latex;
            latex.SetTextSize(0.06);
            latex.SetTextAlign(13);  //align at top
            latex.SetNDC();
            latex.DrawLatex(.7,.71, input_bdt_info.topo_name.c_str());


            double max_height = std::max( sig->GetMaximum(), bkg->GetMaximum());
            sig->SetMaximum(max_height*1.3);


            c_var->Print(("vars/"+std::to_string(nv)+"_"+input_bdt_info.identifier+"_"+v.safe_unit+"_stage_"+std::to_string(j)+".pdf").c_str(),"pdf");


            delete sig;
            delete bkg;
            delete c_var;
            nv++;
        }

    }


    return 0;
}


int plot_bdt_variable(bdt_file * signal_pure, bdt_file * background_pure, bdt_variable v,bdt_info input_bdt_info, bool isSpectator,int stage, std::vector<double> bdtcuts){

    std::vector<bdt_variable> vars = {v};
    return plot_bdt_variables(signal_pure, background_pure, vars, input_bdt_info, isSpectator,stage, bdtcuts);
}


int  plot_bdt_variables(bdt_file * signal_pure, bdt_file *background_pure, std::vector<bdt_variable> vars, std::vector<double> bdt_cuts, int stage){
//Old
//
std::cout<<"Outdated function "<<__func__<<", see "<<__FILE__<<" Line: "<<__LINE__<<std::endl;

/*
    std::vector<std::string> title = {"Topological Selection","Pre-Selection Cuts","cos","bnb","ncpi0","nue","tm"};

        if(stage>1) {
        signal_pure->calcBDTEntryList(stage,bdt_cuts);
        background_pure->calcBDTEntryList(stage,bdt_cuts);
        }
        signal_pure->setStageEntryList(stage);
        background_pure->setStageEntryList(stage);

        for(auto &v: vars){

            //std::string cut_signal = signal_pure->getStageCuts(j,-9,-9); 
            //std::string cut_background_pure = background_pure->getStageCuts(j,-9,-9); 

            //		TH1* sig = signal_pure->getTH1(v,cut_signal.c_str(),v.safe_name+"_sig_var" ,1.0);
            //		TH1* bkg = background_pure->getTH1(v,cut_background_pure.c_str(),v.safe_name+"_bkg_var" ,1.0);
            TCanvas *c_var = new TCanvas(("cvar_"+v.name).c_str(), ("cvar_"+v.name).c_str(),1200,1200);
            c_var->cd();

            std::cout<<"On variable: "<<v.name<<std::endl;

            TH1* sig = signal_pure->getTH1(v,"1",v.safe_name+"_sig_var" ,1.0);
            TH1* bkg = background_pure->getTH1(v,"1",v.safe_name+"_bkg_var" ,1.0);

            std::cout<<"Integrals: "<<sig->Integral()<<" "<<bkg->Integral()<<std::endl;
            sig->Scale(1.0/sig->Integral());			
            bkg->Scale(1.0/bkg->Integral());			
            sig->SetLineColor(signal_pure->col);
            bkg->SetLineColor(background_pure->col);
            sig->SetLineWidth(2);
            bkg->SetLineWidth(2);
            c_var->cd();

            sig->SetFillColor(signal_pure->col);
            bkg->SetFillColor(background_pure->col);
            sig->SetFillStyle(3445);
            bkg->SetFillStyle(3454);

            if(stage!=1){
                sig->SetTitle(title.at(stage).c_str());
            }else{
                sig->SetTitle(" ");
            }

            c_var->cd();			

            sig->Draw("hist");
            sig->SetMinimum(0);
            bkg->Draw("hist same");
            //sig->GetXaxis()->SetTitle(v.unit.c_str());
            sig->GetYaxis()->SetTitle("Events [Area Normalized]");
            sig->GetYaxis()->SetTitleOffset(1.5);

            TLegend *l = new TLegend(0.11,0.75,0.89,0.89);
            l->SetLineColor(kWhite);
            l->SetFillStyle(0);
            l->SetNColumns(2);


            l->AddEntry(sig, signal_pure->plot_name.c_str(),"lf");	
            l->AddEntry(bkg, background_pure->plot_name.c_str(),"lf");	
            l->Draw();

            TText *pre;
            if (false) {
                pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulaton In-Progress - Spectator Variable");
            }else {
                pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulaton In-Progress - Training Variable");

            }


            //TText *pre = drawPrelim(0.1,0.915,0.03,"MicroBooNE Simulation - In Progress");
            pre->Draw();

            TLatex latex;
            latex.SetTextSize(0.06);
            latex.SetTextAlign(13);  //align at top
            latex.SetNDC();
            latex.DrawLatex(.7,.71, signal_pure->flow.bdt_vector[0].topo_name.c_str());


            double max_height = std::max( sig->GetMaximum(), bkg->GetMaximum());
            sig->SetMaximum(max_height*1.3);


            c_var->Print(("var/"+v.safe_unit+"_stage_"+std::to_string(stage)+".pdf").c_str(),"pdf");


            delete sig;
            delete bkg;
            delete c_var;
        }


*/
    return 0;

}
