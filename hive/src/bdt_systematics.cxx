#include <bdt_systematics.h>

/*
 * Some gadgets to keep things organized;
 */

TString gadget_nameroot(TString plot_type, TString dir, bdt_variable var){

	TString binning_text = to_string_prec(var.n_bins+var.plot_min+var.plot_max,0);
	TString file_name = dir+"/" + plot_type+"_"+var.safe_name +  binning_text + ".root";

	return file_name;
}


TString gadget_addindex(TString varname, int nth){

	TString nam = varname;
//	std::cout<<"before: "<<nam<<std::endl;	
	TString lab = "[" + to_string_prec(nth,0)+"]";
	
	nam.ReplaceAll("[0]", lab);
//	std::cout<<nam<<std::endl;	
	return nam;

}
/*
 * Constructor for struct bdt_sys, CHECK
 *
 */
bdt_sys::bdt_sys(int index, MVALoader XMLconfig)
	:
	tag				(XMLconfig.sys_tag[index]),
	dir           	(XMLconfig.sys_dir[index]),
	filename      	(XMLconfig.sys_filename[index]),
	treename      	(XMLconfig.sys_treename[index]),
	vars          	(XMLconfig.sys_vars[index]),
	vars_name     	(XMLconfig.sys_vars_name[index]),
	pot           	(XMLconfig.sys_pot[index]),
	throws        	(XMLconfig.sys_throws[index]),
	its_CV         	(XMLconfig.sys_its_CV[index]),
	its_multithrows	(XMLconfig.sys_its_multithrows[index]),
	its_OM	(XMLconfig.sys_its_OpticalModel[index])
	{
		num_vars = vars.size();
		hist.resize(num_vars);
		twodhist.resize(num_vars);
	
	};

/*
 * InitSys version 2, use branches;
 *
 */

void InitSys(std::vector<bdt_variable> vars, std::vector<std::string> precuts, std::vector<bdt_sys*> syss, double plot_pot, TString dir_root, TString dir_drawn){
	
	bool message = true;
	
	//make a string of precuts, CHECK, it will be more complicated when using BDT cuts;
	std::string temp_cut("("+precuts[0]+")");
	for(int index = 1; index < precuts.size(); ++index	){
		temp_cut += "&&("+precuts[index] + ")";
	}
	TString applycuts = temp_cut.c_str();

	//the histogram goes two ways: 1dhist or 2dhist
	
	bool its_1dhist = (vars.size()<2)? true: false;

	 bdt_variable* var1 = &vars[0];
	 bdt_variable* var2;
	 if(!its_1dhist) var2 = &vars[1];
	//Part I make histiograms;
	//work on each systematic


	if(message) std::cout<<"\n ---- Initialize systematics for "<<var1->unit<<std::endl;

	{//CHECK , need to update to more intelegent;
		TString check_covroot_name = dir_root+"/../cov_"+var1->safe_name + ".root";
		TFile *check_cov_root = (TFile*) TFile::Open(check_covroot_name,"READ"); //one var one hist root;
		if(check_cov_root != NULL ){ //oh we had the file! so .. 
			std::cout<<" Being lazy, no need to make cov matrix, because we had it already!"<<std::endl;
			check_cov_root->Close();
			goto ok_wedone;
		}
		delete check_cov_root;
	}

{
	TString outfile_name = gadget_nameroot("hist",dir_root, *var1);//, temps);

	TFile *hist_root = (TFile*) TFile::Open(outfile_name,"READ"); //one var one hist root;
	if(hist_root!=NULL ){ //quick check of existing files; 
		std::cout<<outfile_name<<" exist, no need to regenerate."<<std::endl;

		for( size_t index = 0; index < syss.size(); ++index){
			bdt_sys* temps = syss[index];
			int num_throws = temps->throws;
			bool its_multithrows = temps->its_multithrows;//multi weights or optical model

			for(int jndex = 0; jndex < temps->num_vars; ++jndex){//loop through weights

				TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
				TString Tdir_name = temps->tag+"_"+weights_nam;
				for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;

					TString hist_nam_label = (its_multithrows)? std::to_string(kndex+1) : "";
					TString hist_nam = weights_nam + hist_nam_label;
					TH1F* temp_th1 =  (TH1F*) hist_root->Get(Tdir_name+"/"+hist_nam);

					(temps->hist[jndex]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
					if(message) std::cout<<"\r\tLoad histogram "<<hist_nam;
				}//next throws
				if(message) std::cout<<std::endl;
			}//next systematic weights;
		}//next systematics
	}else{
		
		std::cout<<outfile_name<<" does not exist. Make it now."<<std::endl;

		hist_root = (TFile*) TFile::Open(outfile_name,"RECREATE"); //one var one hist root;

		for( size_t index = 0; index < syss.size(); ++index){
			bdt_sys* temps = syss[index];

			if(message) std::cout<<"\n"<<temps->tag<<" has "<<temps->num_vars<<" variables."<<std::endl;

			bool its_multithrows = temps->its_multithrows;//multi weights or optical model

			// type 1 single throw one file
			//STEP 1 Load files, TTree
			gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
			TFile *infile = (TFile*) TFile::Open(temps->dir+"/"+temps->filename,"READ"); 
			gSystem->RedirectOutput(0,0);//this let ROOT warns
			if(infile ==NULL){
				std::cout<<temps->dir+"/"+temps->filename<<" does not exist, please check! Program abort."<<std::endl;
				exit(EXIT_FAILURE);
			}
			TTree* temptree = (TTree*) infile->Get(temps->treename);//temptree is from the systematic weights file;

			//STEP 2 get the TH1F
			for(int jndex = 0; jndex < temps->num_vars; ++jndex){//loop through weights
				//subdirectory in root for each systematic weighs;

				TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
				TString Tdir_name = temps->tag+"_"+weights_nam;

				TDirectory *Tdir = hist_root->mkdir(Tdir_name);
				Tdir->cd();

				//STEP 2.2 prepare TH1F
				//			TFile *outfile = (TFile*) TFile::Open(outfile_name,"RECREATE"); 
				//			std::cout<<dir_root+"/"+outfile_name<<std::endl;
				const int num_throws = temps->throws;
				const int num_swvars = temps->num_vars;

				Long64_t nentries = temptree->GetEntries();
				if(message) std::cout<<" We have entries: "<<nentries<<std::endl;

				
				//STEP 2.1 make empty histgrams to be filled.
				std::vector<TH1F*> hist1d(num_throws);

				int count_hist = 0;
				for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
					TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];
					hist1d[kndex] = new TH1F(hist_name,hist_name, var1->n_bins, var1->plot_min, var1->plot_max);
					count_hist++;
				}
				if(message) std::cout<<"Creating "<<count_hist<<" histograms "<<std::endl;


				std::vector<TTreeFormula*> wgtforms(num_throws);
				
				for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
					TString temp_label = (kndex>0)? "["+std::to_string(kndex)+"]":"";
					TString temp_wgt = temps->vars[jndex] + temp_label;
					TString temp_wgtname = temps->vars_name[jndex]+std::to_string(kndex);

					wgtforms[kndex] = new TTreeFormula(temp_wgtname, temp_wgt+"*"+gadget_addindex(applycuts,kndex), temptree);//get the x-axis, with precuts;
					wgtforms[kndex]->GetNdata();
//					std::cout<< wgtforms[kndex]->EvalInstance()<<std::endl;
//					if(kndex== 10) sleep(10);
				}

				for(Long64_t kentry = 0; kentry< nentries; ++kentry){

					temptree->GetEntry(kentry);

					for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
						
						//CHECK, add continue statement for precuts;

						TTreeFormula varformula(var1->unit.c_str(), gadget_addindex(var1->mininame, kndex), temptree);//get the x-axis
						varformula.GetNdata();
						double temp_value = varformula.EvalInstance();
						double temp_weight = wgtforms[kndex]->EvalInstance()*plot_pot/temps->pot;

						hist1d[kndex]->Fill(temp_value,temp_weight);//jndex - nth sw; kndex - nth throw

						std::cout<<"\r\tLooping entries: "<<std::setw(5)<<kentry+1<<"/"<<std::setw(5)<<nentries;
						std::cout<<" "<<gadget_addindex(var1->mininame, kndex);
						std::cout<<"="<<std::setw(5)<<temp_value;
						std::cout<<" "<<temps->vars[jndex]<<"["<<kndex<<"]";
						std::cout<<"="<<std::setw(5)<<temp_weight;
						std::cout<<std::setw(4)<<kndex+1<<" throw sw name: "<<temps->vars_name[jndex];
						//					std::cout<<std::endl;
						//					if(kndex==10) exit(0);
					}//next throw;
				}

				std::cout<<"\n Finish filling "<<num_swvars<<" histograms"<<std::endl;

				//write hisograms to the diretory
				for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
					//configure th1
					TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];
					hist1d[kndex]->SetStats(0);
					hist1d[kndex]->GetXaxis()->SetTitle(var1->unit.c_str());
					hist1d[kndex]->GetYaxis()->SetTitle("Events");
				
					//save the th1 to the root file, class bdt_sys;
					Tdir->cd();
					hist1d[kndex]->Write(hist_name);
					//This is important for making covariance matrix using hist2cov();
					(temps->hist[jndex]).push_back( (TH1F*) hist1d[kndex]->Clone());// each bdt_sys now contains all hist;
					hist1d[kndex]->Delete();
				}

				std::cout<<"Save "<<num_throws;
				std::cout<<" histograms to the file: "<<outfile_name<<" dir:"<<Tdir_name<<std::endl;

			}//next systematic weight 
			infile->Close();
		}//next systematic files
	}
	std::cout<<"Finish making systematic histograms."<<std::endl;
	//
	//Part II: Overlay histograms & draw covaraicne matrix;
	//
	hist2cov( *var1, syss, dir_root, dir_drawn);
	hist_root->Close();
}
ok_wedone://CHECK
	std::cout<<" finish systematic matrix preparation"<<std::endl;
}


/*
 * Build covariance matrices based on histograms;
 */

void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss, TString dir_root, TString dir_drawn){
		
	bool message = true;
	
	if(message) std::cout<<"Making covariance matices."<<std::endl;
	
	TString unique_label = var.safe_name+ to_string_prec(var.n_bins+var.plot_min+var.plot_max,0)+"_";
	//
	//STEP 1 categorize hists; 
	std::vector<TString> tags;
	std::multimap< TString, bdt_sys* > sys2SWmap;
	std::map< TString, bdt_sys* > sys2CVmap;
	std::map< TString, bdt_sys* > sys2OMCVmap;

//	bool only_one_cv = true;
	for(size_t index = 0; index < syss.size(); ++index){//first classify what tags are available, then group them next;
		TString temptag = syss[index]->tag;
		tags.push_back(temptag);

		bdt_sys* temps = syss[index];
		if(message) std::cout<< "Looking at "<<temptag<<" ";	
		if(temps->its_CV){//save CV to map;
			sys2CVmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<" it is CV!";
		}else if(temps->its_OM){//Optical model has another CV
			sys2OMCVmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<" it is a OMCV!";

		} else{
			sys2SWmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<" it is Systematic weights!";
		}//save SW to map
		if(message)std::cout<<std::endl;
	}

	//STEP 1.1 clean up tags;
//	std::cout<<" size of tags "<<tags.size()<<std::endl;
	sort(tags.begin(), tags.end() );
	tags.erase( unique( tags.begin(), tags.end()), tags.end() );
//	std::cout<<" size of tags (after clean) "<<tags.size()<<std::endl;

	//STEP 2 pair up SW/CV to make covariance matrix;
	//prepare root file to store  covariance matrix;
	TString covroot_name = gadget_nameroot("cov",dir_root, var);//the root to save following covariance marices;
	TFile *cov_root = (TFile*) TFile::Open(covroot_name,"RECREATE"); //one var one hist root;
	TString finalcovroot_name = dir_root+"/../cov_"+var.safe_name + ".root";
	TFile *finalcov_root = (TFile*) TFile::Open(finalcovroot_name,"RECREATE"); //one var one hist root;


	int nb=var.n_bins;//cv_hist->GetNbinsX();//# of bins
	int nl=var.plot_min;//cv_hist->GetBinLowEdge(1);//lower bound
	int nh=var.plot_max;//cv_hist->GetBinLowEdge(nb+1);//upper bound
	
	TString finalcov_name = (var.safe_name).c_str();
	TH2D* finalcov =  new TH2D(finalcov_name, finalcov_name ,nb,nl,nh,nb,nl,nh);//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
	

	for(size_t index = 0; index < tags.size(); ++index){//each tag means one group of covariance matices;
		//prepare cov matrix here;
		//STEP 2.1 identify which bdt_sys syss[?] gives the CV, SW; and get the histogram;
		TString cur_tag = tags[index];
		//STEP 2.1.1 CV first
		bdt_sys* tempsCV = sys2CVmap.find(tags[index])->second;//get the sys for a given tag;


		TH1F* cv_hist = (TH1F*) (tempsCV->hist[0][0])->Clone();

		//STEP 2.1.2 Load SW;
		std::map< TString, std::vector<TH1F*> > sw_name2throws;//one sw - many_throws TH1F*
		//			std::vector< std::vector<TH1F*> > sw_hist;// each bdt_file->hist is a vector< vector<TH1F*> >
		//			std::vector< TString > sw_hist_name;
		for( std::_Rb_tree_iterator<std::pair<const TString, bdt_sys*> >  itr = sys2SWmap.begin();  itr!=sys2SWmap.end(); itr++){
			//go through each tag;
			if(itr->first == cur_tag){//find the tag, then check the corresponding bdt_sys;
				for( size_t jndex = 0; jndex < (itr->second)->hist.size(); ++jndex){
					//go through each SW(with many throws) of a bdt_sys
					TString temp_sw_name  = (itr->second)->tag+"_"+(itr->second)->vars_name[jndex];
					std::vector<TH1F*> temp_sw_hists = (itr->second)->hist[jndex] ;
					sw_name2throws.insert(std::make_pair( temp_sw_name, temp_sw_hists) );
				}//next sw (with nthrows)
			}//next systematics

		}//next bdt_sys (with various sw) under the same tag 

		//STEP 2.2 calculate covariance matrix and draw histograms;
		//note cv_hist is one histogram, sw_name2throws->second are different sw with N throws of histograms each;

		TCanvas *histCanvas = new TCanvas("histCanvas", "", 900, 400);
		TCanvas *covCanvas = new TCanvas("covCanvas","",900,400);


		int nby = 1.5*cv_hist->GetBinContent(cv_hist->GetMaximumBin());//histogram height

		for( std::_Rb_tree_iterator<std::pair<const TString, std::vector<TH1F*> > >  itr = sw_name2throws.begin();  itr!=sw_name2throws.end(); itr++){//loop over differenct covariance matrix
			
			//STEP 2.2.1 prepare covariance matrices and histograms;
			//work through the map< TString, std::vector<TH1F*> > sw_name2throws;
			TString temp_histname = itr->first;
			TString temp_covname = itr->first;//"<tag>_<sw> Covariance Matrix", name of cov matrix
			if(message) std::cout<<"Calculate covariance matrix for "<<itr->first<<std::endl;
			std::vector<TH1F*> temp_swhist = itr->second;
			int num_throws = temp_swhist.size();

			//for ovarlaying all histograms
			TH2D* all_hist = new TH2D(temp_histname, temp_histname, nb, nl,nh, nby, 0, nby);
			//for summing up cov matrice
			TH2D* covmatrices =  new TH2D(temp_covname+" Covariance Matrix", temp_covname+" Covaraince Matrix",nb,nl,nh,nb,nl,nh);//binning for xnbins,xmin,xmax,ybins,ymin.ymax;

			//to record the cov matrix;
			for(size_t jndex = 0; jndex < num_throws; ++jndex){//go through different throws

				//for overlay hsitograms
				for(int kndex = 1; kndex < nb+1; ++kndex){//get everything from a single histogram;
					all_hist->Fill(temp_swhist[jndex]->GetBinLowEdge(kndex) , temp_swhist[jndex]->GetBinContent(kndex) );
				}//next bin

				//for cov matrices
				TH2D* cov_temp =MakeCov("cov"+temp_covname+std::to_string(jndex), temp_swhist[jndex], cv_hist);
				cov_temp->Scale(1.0/(double) num_throws);
				covmatrices->Add(cov_temp);

//				std::cout<<"\r"<<jndex<<"th throws of histogram"<<std::endl;
			}//next throws
//			std::cout<<std::endl;
			

			//STEP 2.2.2 Draw histograms & covaraince matrix, save the covarianc matrix
			gStyle->SetPalette(kRainBow);//CHECK, this change colz colors;
			
			histCanvas->Clear();
			histCanvas->cd();
			all_hist->SetStats(false);
			all_hist->Draw("COLZ");
			cv_hist->Draw("L same");	
			cv_hist->SetLineColor(6);
			cv_hist->SetLineWidth(2);
			//Save the histogram


			histCanvas->SaveAs( dir_drawn + "/hist_"+unique_label+temp_histname+".pdf" ,"pdf");
			
			covCanvas->Clear();
			covCanvas->cd();

			 
			if(sys2OMCVmap.count(tags[index]) > 0){//its a Optical Model CV, propagate the matrix to a fractional covariance matrix; and add it to finalcov.
				bdt_sys* tempsOMCV = sys2OMCVmap.find(tags[index])->second;//get the sys for a given tag;

				if(message) std::cout<<"Making fractional covariance matrix for "<<tags[index]<<std::endl;
				TH1F* OMCV = (TH1F*) (tempsOMCV->hist[0][0])->Clone();
				TH2D* fracCov = MakeFracCov(temp_covname, covmatrices, cv_hist, OMCV);

				fracCov->SetStats(false);
				fracCov->Draw("COLZ");
				fracCov->SetTitle("Propagated "+temp_covname + " Covaraince Matrix");
				covCanvas->SaveAs( dir_drawn + "/cov_frac"+unique_label+temp_histname+".pdf" ,"pdf");
				if(message) std::cout<<"1,1 in the matrix "<< tags[index]<<" "<<fracCov->GetBinContent(1,1)<<std::endl;

				finalcov->Add(fracCov);

				cov_root->cd();
				fracCov->Write();

			}else{//normally, add the cov matrix to the final total matrix; 

				covmatrices->SetStats(false);
				covmatrices->GetXaxis()->SetTitle((var.unit).c_str());
				covmatrices->Draw("COLZ");
				covCanvas->SaveAs( dir_drawn + "/cov_"+unique_label+temp_histname+".pdf" ,"pdf");

				if(message) std::cout<<"1,1 in the matrix "<< tags[index]<<" "<<covmatrices->GetBinContent(1,1)<<std::endl;

				finalcov->Add(covmatrices);
				cov_root->cd();
				covmatrices->Write();
			}


		}//next sw covaraince matrix
	}//next tag;
	cov_root->Close();
	
	//set axis and name for the final covaraince matrix;
	finalcov->SetStats(false);

	finalcov->SetTitle(" Total Covariance Matrix");
	finalcov->GetXaxis()->SetTitle( (var.unit).c_str());

	TCanvas *c1 = new TCanvas("c1","",900,400);
	finalcov->Draw("COLZ");

	TString saved_root = gadget_nameroot("Total_cov",dir_root, var);//, temps);
	c1->SaveAs(saved_root,"pdf");
	if(message) std::cout<<"1,1 in the final matrix "<<finalcov->GetBinContent(1,1)<<std::endl;

	finalcov_root->cd();
	finalcov->Write();


	//turn the finalcov in the format of a matrix
	TMatrixD covMatrix(nb,nb); 
	TArrayD nums(nb*nb);

	for(int index = 0; index < nb; ++index){
		for(int jndex = 0; jndex < nb; ++jndex){
			nums[jndex*nb+index] = finalcov->GetBinContent(index+1,jndex+1);
		}
	}
	covMatrix.SetMatrixArray(nums.GetArray());
	covMatrix.Write((var.safe_name+"matrix").c_str());


	finalcov_root->Close();
	std::cout<<"Final total covariance matrix is saved at "<<finalcovroot_name<<std::endl;
}


TH2D* MakeCov(TString name,TH1F* hist, TH1F* cv){

	int nb=hist->GetNbinsX();
	int nl=hist->GetBinLowEdge(1);
	int nh=hist->GetBinLowEdge(nb+1);
//	cout<<"Bin config # bins: "<<nb<<" low "<<nl<<" high "<<nh<<" entries:"<<hist->Integral()<<endl;

//	cout<<"Total events in (MakeCov) "<<hist->Integral()<<endl;
	int counter = 0;
	if(nb!=cv->GetNbinsX() || nl !=cv->GetBinLowEdge(1)||nh!=cv->GetBinLowEdge(nb+1) ){
		std::cerr<<" Histograms don't match to the CV histogram. No Covariance matrix is calculated."<<std::endl;
		exit(EXIT_FAILURE);
	}
	TH2D* covmatrix =  new TH2D(name,name,nb,nl,nh,nb,nl,nh);//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
	for(int index = 1; index<nb+1; ++index){
		for(int jndex = 1; jndex<nb+1; ++jndex){
		double entry = (hist->GetBinContent(index)-cv->GetBinContent(index) )*(hist->GetBinContent(jndex)-cv->GetBinContent(jndex) );
//		covmatrix->Fill(index,jndex,entry);//Fill goes by values
		covmatrix->SetBinContent(index,jndex,entry);

//		cout<<"\r Fill in "<<name<<" "<<index<<","<<jndex<<" with "<<entry;
//		if(entry!=0) cout<<" # of good entry   --->"<<counter++;
		}
	}

	return covmatrix;
}

/*
 *
 * covariane matrix propagation, create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TString name,TH2D* cov_temp, TH1F* oldcv, TH1F* newcv){

	int nb=newcv->GetNbinsX();
	int nl=newcv->GetBinLowEdge(1);
	int nh=newcv->GetBinLowEdge(nb+1);
	
	TH2D* fractional_cov = (TH2D*) cov_temp->Clone("Fractional_Covaraince_Matrix_"+name);

//	fractional_cov->Reset("ICESM");//clear contents;

	for(int jndex = 1; jndex < nb+1; ++jndex){
		for(int kndex = 1; kndex < nb+1; ++kndex){
			double cvj = oldcv->GetBinContent(jndex);
			double cvk = oldcv->GetBinContent(kndex);
			double fjk = cov_temp->GetBinContent(jndex,kndex)/(cvj*cvk);

			double cvj_new = newcv->GetBinContent(jndex);
			double cvk_new = newcv->GetBinContent(kndex);
			fractional_cov->SetBinContent(jndex,kndex,fjk*cvj_new*cvj_new);
		}
	}

	return fractional_cov;
}

/*
 * Make Covaraince matrix from 2d histograms
 *
 */

TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv){

}

