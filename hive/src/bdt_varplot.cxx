#include "bdt_varplot.h"

//int plot_bdt_variables(bdt_file * signal_pure, bdt_file * background_pure, std::vector<bdt_variable> vars, bdt_info input_bdt_info, bool isSpectator,int stage,std::vector<double> cuts);
int  plot_var_allF(std::vector< bdt_file *> MCfiles, bdt_file* datafile, std::vector<bdt_variable> vars, std::vector<bdt_sys*> syss, int stage, std::vector<double> bdtcuts, bool isNorm){
//    std::vector<std::string> title = {"Topological Selection","Pre-Selection Cuts","Cosmic","BNB","Pi0","Other"};
    {
		int num_MC = MCfiles.size();
		double is_bestfit_in = false;//false - not include signal in the output;
		bool do_pair_plot = true;//true - excess vs each MC;
		bool draw_estimator = true;
		bool debug_message = true;

		double yaxis_factor = 1.6;

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

			double maxy = 0;
			double miny = 0;

			std::cout<<"On variable: "<<v.name<<std::endl;

			//Prepare data first, modified it according to MC and draw 
			TH1* data = datafile->getTH1(v,"1",v.safe_name+"_data_var" ,datafile->pot);
			TLegend *legend = new TLegend(0.13,0.73,0.89,0.88);

			data->SetMarkerSize(3);
			data->SetMarkerStyle(20);
			data->SetLineColor(kBlack);
//			data->SetBinErrorOption(TH1::kPoisson);

			//Prepare MC
			std::vector<TH1*> MC(num_MC);
			TH1* AllNue_MC;
			bool first_hist = true;

			for(size_t index = 0; index <num_MC; ++index){//loop through each NueMC; do the substraction;

				MC[index] = MCfiles[index]->getTH1(v,"1",v.safe_name+"_MC_var" , datafile->pot);
				
				//adjust MC error using the Matrix;
				TFile *covar_f = new TFile(v.covar_file.c_str(),"read");
				TMatrixD * covar_collapsed = new TMatrixD(v.n_bins,v.n_bins);
				*covar_collapsed = gadget_PrepareMatrix(syss, covar_f, MC[index] , datafile->pot, MCfiles[index]->tag);
				//			std::vector< TH2D* > gadget_SeparateMatrix( covar_collapsed, MC[index], "vars/");
				for(int jndex = 0; jndex < MC[index]->GetNbinsX() ; ++jndex){

					double temp_binErrSq = pow(MC[index]->GetBinError(jndex+1),2) + pow(MC[index]->GetBinContent(jndex+1),2)*(*covar_collapsed)(jndex,jndex);
					std::cout<<"Update "<<MCfiles[index]->tag<<" Bin "<<1+jndex<<"/"<<MC[index]->GetNbinsX()<<" Error "<<MC[index]->GetBinError(jndex+1);
					std::cout<<" to ErrSq: "<<temp_binErrSq<<" MC: "<<MC[index]->GetBinContent(jndex+1)<<" FM_ii:"<<(*covar_collapsed)(jndex,jndex)<< std::endl;
					MC[index]->SetBinError(jndex+1, sqrt(temp_binErrSq));
				}
				MC[index]->SetLineColor((MCfiles[index]->is_signal)? kBlack:MCfiles[index]->col);
				MC[index]->SetLineWidth(2);

//				MC[index]->SetFillColor(MCfiles[index]->col);
				MC[index]->SetFillStyle(3445);

				//addd legend
				legend->AddEntry(MC[index], (MCfiles[index]->plot_name).c_str(),"lf");	
//				legend->AddEntry(MC[index], (MCfiles[index]->plot_name+" w. Stat & Sys. Error").c_str(),"lf");	

				if(is_bestfit_in || !MCfiles[index]->is_signal ){//NueMC must go through, Best-fit may not;

					if(first_hist){//copy the first MC hist
						AllNue_MC = (TH1*) MC[index]->Clone("AllNuemc");
						first_hist = false;
					} else{//add to the current MC hist;
						AllNue_MC->Add(MC[index],1.0);

					}
					std::cout<<" data - "<< MCfiles[index]->tag <<std::endl;
					data->Add(MC[index],-1.0);//Data - MC
				}
			}

			double maxy_of_data = data->GetBinContent(data->GetMaximumBin());
			//Canvas to be Drawn
			TCanvas *c_var = new TCanvas(("cvar_"+v.name).c_str(), ("cvar_"+v.name).c_str(),1200,1200);

			c_var->cd();

			TString legend_label = (is_bestfit_in)? "Data - (NueMC+Best Fit)": "Data - NueMC";
			legend_label += " w. Stat & Sys. Error";
			legend->AddEntry(data, legend_label, "lp" );


			//No need to adjust errors
			for(int ib=1; ib<v.n_bins+1; ib++){
				if(debug_message)std::cout<<" rescaled NueMC "<<ib<<" bin has "<<AllNue_MC->GetBinContent(ib)<< " with error "<<AllNue_MC->GetBinError(ib)<<" and the excess Error "<<data->GetBinError(ib)<<std::endl; 
//				data->SetBinError(ib, data_err);
			}

			data->Draw("E1");//black dot

			//set max&min, and scale MC evnets if isNorm is true;
			std::vector< double > mc_scale_factor(num_MC);
			for(size_t index = 0; index <num_MC; ++index){//the data here means the excess;


				mc_scale_factor[index] = data->Integral()/MC[index]->Integral();
				if(isNorm) MC[index]->Scale( mc_scale_factor[index] );
				//determine max
				double maxy_of_MC = MC[index]->GetBinContent(MC[index]->GetMaximumBin());
				if( maxy < maxy_of_MC ) maxy = maxy_of_MC;

				MC[index]->Draw("hist same");
			}
			if( maxy < maxy_of_data ) maxy = maxy_of_data;
			double miny_of_data = data->GetBinContent(data->GetMinimumBin());
			if(miny > miny_of_data ) miny = miny_of_data;


			//Configure labels yaxis;
			data->SetTitle(" ");

			//MC->GetXaxis()->SetTitle(v.unit.c_str());
			data->GetYaxis()->SetTitle((isNorm)? "Events [MC Rescaled]" : "Events");
			data->GetYaxis()->SetTitleOffset(1.5);
			data->SetMaximum(maxy*yaxis_factor);
			data->SetMinimum(miny);

			legend->SetBorderSize(0);
			legend->SetLineColor(kWhite);
			legend->SetFillStyle(0);
			legend->SetNColumns(2);

			legend->Draw();


			//add double axis, specially for R^3
			if((v.unit).find("R/500")!= std::string::npos){
				std::cout<<"Draw extrax axis"<<std::endl;
				double minv = v.plot_min;
				double maxv = v.plot_max;

				//						double mapped_maxv = pow(maxv,(1/3.0))*500;
				TF1 *temp_axis = new TF1("temp_axis","pow(x,3)",0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;
				//						TF1 *temp_axis = new TF1("temp_axis","pow(x,3.0)",0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;

				TGaxis *new_axis = new TGaxis(minv,maxy*yaxis_factor,maxv,maxy*yaxis_factor,"temp_axis", 403, "-");//This says where to draw axis (xmin,ymin,xmax,ymax,TF1,num_division =N1 + 100*N2 + 10000*N3 [by default, ROOT will take numbers not larger than these to optimize the axis] ,option"-");https://root.cern.ch/doc/master/classTGaxis.html
				new_axis->SetTitle("Radius [cm]");
				new_axis->SetLabelSize(0.03);
				new_axis->SetTitleSize(0.03);
				new_axis->SetTitleOffset(1.3);
				new_axis->Draw();
				//std::cout<<"min "<<minv<<" max "<<maxv<<" maxy "<<maxy<<std::endl;
			}
			c_var->Print(("vars/"+std::to_string(nv)+"_"+v.safe_unit+"_stage_"+std::to_string(stage)+".pdf").c_str(),"pdf");

			if(do_pair_plot){//data-NueMC vs each MC
				c_var->Clear();

				for(size_t index = 0; index <num_MC; ++index){
					maxy=yaxis_factor*std::max(data->GetBinContent(data->GetMaximumBin()),MC[index]->GetBinContent(MC[index]->GetMaximumBin()));

					data->Draw("E1");
					data->SetMaximum(maxy);
					data->SetMinimum(miny);

					MC[index]->Draw("hist same");

					//legend
					TLegend *le = new TLegend(0.13,0.75,0.42,0.88);
					le->AddEntry(MC[index], (MCfiles[index]->plot_name).c_str(),"lf");	
//					le->AddEntry(MC[index], (MCfiles[index]->plot_name+" w. Stat & Sys. Error").c_str(),"lf");	

					TString legend_label = (is_bestfit_in)? "Data - (NueMC+Best Fit)": "Data - NueMC";
					legend_label += " w. Stat & Sys. Error";
					le->AddEntry(data, legend_label,"lp");	
					le->SetLineColor(kWhite);
					le->SetFillStyle(0);
					le->SetNColumns(1);

					le->Draw();

					//get ks
					std::string ks = "KS: "+to_string_prec(MC[index]->KolmogorovTest(data), 3);

					//get chi_square
					double mychi=0;
					int ndof = v.n_bins - 1;

					for(int ib=1; ib<v.n_bins+1; ib++){
						double dav = data->GetBinContent(ib);//the excess
						double mcv = MC[index]->GetBinContent(ib);

						//data_err^2 = sigma_data^2+sigma_bkg^2, before substraction
						double mc_wgt_err = MC[index]->GetBinError(ib);//This is sum of sqrt(wgt); after scaling;
						double data_err = data->GetBinError(ib);
//					double mc_stat_err = sqrt(MC[index]->GetBinContent(ib)/mc_scale_factor[index]);
//						if(dav<1e-20){
//							ndof--;
//							continue;
//						}

//						double curchi = pow(dav-mcv,2)/(pow(mc_stat_err,2)+pow(mc_wgt_err,2));
						double curchi = pow(dav-mcv,2)/(pow(data_err,2)+pow(mc_wgt_err,2));
						if(debug_message) std::cout<<"Chi2 Bin "<<ib<<", data "<<dav<<", err "<<data_err<<",MC "<<mcv<<",err "<<mc_wgt_err<<",scale factor "<<mc_scale_factor[index]<<",chi2 "<<curchi<<std::endl;



//						double curchi = pow(dav-mcv,2)/(dav);
						mychi += curchi;
					}

					std::string chi_square = "#chi^{2}/n#it{DOF}: "+to_string_prec(mychi, 2) +"/"+to_string_prec(ndof)+", #chi^{2} P^{val}: "+to_string_prec(TMath::Prob(mychi,ndof),3);

					//draw ks, chi_square
					TLatex *estimators = new TLatex(0.34, 0.85,(ks+", "+chi_square).c_str());
					estimators->SetNDC();
					estimators->SetTextColor(kRed-7); 
					estimators->SetTextSize(0.03);
					if(draw_estimator) estimators->Draw("same");

					//add double axis, specially for R^3
					if((v.unit).find("R/500")!= std::string::npos){
					std::cout<<"Draw extrax axis"<<std::endl;
						double minv = v.plot_min;
						double maxv = v.plot_max;

//						double mapped_maxv = pow(maxv,(1/3.0))*500;
						TF1 *temp_axis = new TF1("temp_axis","pow(x,3)",0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;
//						TF1 *temp_axis = new TF1("temp_axis","pow(x,3.0)",0,maxv*500);//What to draw on the given empty axis: name, function, range of the function;

						TGaxis *new_axis = new TGaxis(minv,maxy,maxv,maxy,"temp_axis", 403, "-");//This says where to draw axis (xmin,ymin,xmax,ymax,TF1,num_division =N1 + 100*N2 + 10000*N3 [by default, ROOT will take numbers not larger than these to optimize the axis] ,option"-");https://root.cern.ch/doc/master/classTGaxis.html
						new_axis->SetTitle("Radius [cm]");
						new_axis->SetLabelSize(0.03);
						new_axis->SetTitleSize(0.03);
						new_axis->SetTitleOffset(1.3);
						new_axis->Draw();
						//std::cout<<"min "<<minv<<" max "<<maxv<<" maxy "<<maxy<<std::endl;
					}

					c_var->Print(("vars/"+std::to_string(nv)+"_"+v.safe_unit+"_"+MCfiles[index]->tag+"_stage_"+std::to_string(stage)+".pdf").c_str(),"pdf");
				}
				if(debug_message){

					std::vector< TH1*> pHists(MC);
					pHists.push_back(data);
					std::vector<bdt_file*> pFiles(MCfiles);
					pFiles.push_back(datafile);
					std::cout<<"------------------- Summary -------------------"<<std::endl;
					gadget_summary( v, pFiles, pHists);
//					int gap=9;
//					std::cout<<std::left<<std::setw(5)<<"Bin#";
//					std::cout<<std::left<<std::setw(gap)<<"Excess";
//					std::cout<<std::left<<std::setw(gap+8)<<"(Excess+MC)_E";
//
//					for(size_t index = 0; index <num_MC; ++index){
//						std::cout<<std::left<<std::setw(gap)<<MCfiles[index]->tag;
//						std::cout<<std::left<<std::setw(gap+6)<<MCfiles[index]->tag+"_statE";
//						std::cout<<std::left<<std::setw(gap+6)<<"scale_factor";
//					}
//					std::cout<<std::endl;
//
//					for(int ib=1; ib<v.n_bins+1; ib++){
////						if(data->GetBinContent(ib) < 10e-20) continue;
//						std::cout<<std::left<<std::setw(5)<<ib;
//						std::cout<<std::left<<std::setw(gap)<<data->GetBinContent(ib);
//						std::cout<<std::left<<std::setw(gap+8)<<data->GetBinError(ib);
//						for(size_t index = 0; index <num_MC; ++index){
//							std::cout<<std::left<<std::setw(gap)<<MC[index]->GetBinContent(ib);
//							std::cout<<std::left<<std::setw(gap+6)<<MC[index]->GetBinError(ib);
//							std::cout<<std::left<<std::setw(gap+6)<<mc_scale_factor[index];
//						}
//						std::cout<<std::endl;
//					}
				}

			}
			delete c_var;
			nv++;

		}

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
