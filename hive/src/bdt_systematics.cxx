#include <bdt_systematics.h>

int sys_env::checkEnv(){
	//check dirs, stage, hash
	try{
		if(fstage<0){ throw "Stage";}
		if(fcut_hash.Length() < 1){throw "Hash";}
		if(top_dir.Length()<1){ throw "Directory";}
		if(drawn_dir.Length()<1){ throw "Drawn directory";}
		if(root_dir.Length()<1){ throw "Root directory";}
	}

	catch (std::string errorMessage){
		std::cout<<"SysErr:"; 
		std::cout<<errorMessage;
		std::cout<<" is not set. Please check."<<std::endl;
		exit(EXIT_FAILURE);
	}
	return 0;
};

/*
 * Some gadgets to keep things organized;
 */
TString gadget_labelroot(TString plot_type, bdt_variable var, int stage, std::string hash){

//	TString binning_text = "_"+to_string_prec(var.n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);
	
	TString binning_text = to_string_prec((var.plot_max-var.plot_min)/var.int_n_bins,0);//to_string_prec(var.int_n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);

	//<type>_<var>_<bins, min, max>
	TString file_name = plot_type+"_"+var.safe_name + "_bgap"+ binning_text + "_s" + to_string_prec(stage,0) + "_"+hash;

	return file_name;
	//for variable binng that can be derived from this binnings, use these root files;
};

//Identical to the above;
TString sys_env::getFileName(TString plot_type, std::vector<bdt_variable> var){
	if(var.size()<2){
		TString bin_gap = to_string_prec((var[0].plot_max-var[0].plot_min)/var[0].int_n_bins,0);

		//<type>_<var[0]>_bgap<gap>_s<stage>_<hash>
		TString file_name = plot_type+"_"+var[0].safe_name + "_bgap"+ bin_gap + "_s" + to_string_prec(this->fstage,0) + "_"+this->fcut_hash;

		return file_name;
	}

};



TString gadget_updateindex(TString varname, int nth){

	TString nam = varname;
//	std::cout<<"before: "<<nam<<std::endl;	
	TString lab = "[" + to_string_prec(nth,0)+"]";
	
	nam.ReplaceAll("[0]", lab);
//	std::cout<<nam<<std::endl;	
	return nam;
}

//unsigned long  gadget_jenkins_hash(std::string key) {
//    size_t length = key.size();
//    size_t i = 0;
//    unsigned long hash = 0;
//    while (i != length) {
//        hash += key[i++];
//        hash += hash << 10;
//        hash ^= hash >> 6;
//    }
//    hash += hash << 3;
//    hash ^= hash >> 11;
//    hash += hash << 15;
//    return hash;
//}

/*
 * Constructor for struct bdt_sys, 
 *
 */
bdt_sys::bdt_sys(int index, MVALoader XMLconfig, int bdtfile_index, bdt_flow inflow)
	: bdt_file(bdtfile_index, XMLconfig, inflow){
		//members to be read off immediately;
		systag			=XMLconfig.sys_tag[index];
		tag				=XMLconfig.sys_tag[index] + "_" + bdt_file::tag;//double tags
		dir           	=XMLconfig.sys_dir[index];
		filename      	=XMLconfig.sys_filename[index];
		treename      	=XMLconfig.sys_treename[index];
		vars          	=XMLconfig.sys_vars[index];
		vars_name     	=XMLconfig.sys_vars_name[index];
		pot           	=XMLconfig.sys_pot[index];
		throws        	=XMLconfig.sys_throws[index];
		its_CV         	=XMLconfig.sys_its_CV[index];
		its_multithrows	=XMLconfig.sys_its_multithrows[index];
		its_OM			=XMLconfig.sys_its_OpticalModel[index];

		fullyloaded = true;
		//members that need to be derived slightly;
		num_vars = vars.size();
		hists.resize(num_vars);
		twodhists.resize(num_vars);
		histNames.resize(num_vars);
		twodhistNames.resize(num_vars);
		TdirNames.resize(num_vars);

		loads.assign(num_vars, false);
		for(int index = 0 ; index < num_vars; ++index){//loop through weights
			TdirNames[index] = tag+"_"+vars_name[index];

			for(int jndex = 0; jndex < throws; ++jndex){//save names;
				TString hist_name = vars_name[index];
				if(its_multithrows) hist_name += std::to_string(jndex);
				histNames[index].push_back(hist_name);
				twodhistNames[index].push_back(hist_name);
			}//next throw
		}//next weight

	};



bool bdt_sys::Load1dhist(TFile* cur_file){

	for(int index = 0; index< this->num_vars; ++index){//loop through weights;
		if(this->loads[index]) continue;
		this->loads[index] = true;
		for(int jndex = 0; jndex< this->throws; ++jndex){//loop through throws;
			TString nth_label = (this->its_multithrows)? std::to_string(jndex+1) : "";
			TString hist_name = this->vars_name[index]+nth_label;

			TH1F* temp_th1 =  (TH1F*) cur_file->Get(this->TdirNames[index]+"/"+hist_name);

			if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
				if(verbose) std::cout<<"\t"<<hist_name+" does not exist. Next"<<std::endl;
				this->loads[index] = false;
				(this->hists[jndex]).clear();
				break;
			}
			

			(this->hists[index]).push_back( (TH1F*) temp_th1->Clone() );// load it!
			if(verbose) std::cout<<"\r\tLoad Tdir: "<<this->TdirNames[index]<<" histogram: "<<hist_name;
		}//next throw
		if(verbose) std::cout<<std::endl;
		this->fullyloaded *= this->loads[index];
	}//next weight

	return (this->fullyloaded);
}

/*
 * Prepare 1d histograms here;
 * No POT weighting, nor variable binning (in this case, use equal binning then rebin later;
 *
 * Configuration check list:
 * - stage
 * - verbose
 * - 
 */
void bdt_sys::Make1dhist(bdt_variable *var){//TString histfilename, bdt_variable* var,  double plot_pot, std::vector<double> bdt_cuts){}

	
	if(verbose) std::cout<<"\nMaking systematic histograms for "<<this->tag<<" with "<<this->num_vars<<" variables at stage "<<this->fstage<<std::endl;

	//configure what type of histograms to make;
	int nbins =  var->int_n_bins;
	int lowbinedge = var->plot_min;
	int highbinedge = var->plot_max;

//	bool its_multithrows = this->its_multithrows;//multi weights or optical model
	const int num_throws = this->throws;
	const int num_swvars = this->num_vars;

	//STEP 1 Load files, prepare input TTree
	gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
	TFile *infile = (TFile*) TFile::Open(this->dir+"/"+this->filename,"READ"); 
	gSystem->RedirectOutput(0,0);//this let ROOT warns

	if(infile ==NULL){
		std::cout<<this->dir+"/"+this->filename<<" does not exist, please check! Program aborts."<<std::endl;
		exit(EXIT_FAILURE);
	}

	TTree* temptree = (TTree*) infile->Get(this->treename);//temptree is from the systematic weights file;

	//STEP 2 prepare TH1F
//	for(int index = this->start_with_weight; index < this->num_vars; ++index){//loop through weights}
	for(int index = 0 ; index < this->num_vars; ++index){//loop through weights

		if(this->loads[index]) continue;//well, this means no need make histograms;

		TFile* hist_root = (TFile*) TFile::Open(this->sys_env::root_dir + hist_prefix + ".root","UPDATE"); //Frequent Open and Close, hope this help prevent corrupted links;
		std::cout<<"Updating "<<sys_env::root_dir + sys_env::hist_prefix + ".root"<<std::endl;
		//subdirectory in root for each systematic weighs; start_with_weight labels which weight to start, i.e. index;
		//STEP 2.1
		//second check of directory, because there might be some leftover from previous run;
		

//		TString weights_nam = this->vars_name[index];//name of key for the histogram in root;
		TString Tdir_name = TdirNames[index];//this->tag+"_"+weights_nam;

		TDirectory* Tdir = hist_root->GetDirectory(Tdir_name);

//		if(Tdir !=NULL){ //TDirectory exists;
//			std::cout<<"\tDirectory "<<Tdir_name<<" exists; is it a good leftover?"<<std::endl; 
//
////			TString hist_name = (its_multithrows)? this->vars_name[index]+std::to_string(num_throws) :this->vars_name[index];
//
//			TH1F* check_temp_th1 =  (TH1F*) hist_root->Get(Tdir_name+"/"+this->histNames[index]);
//
//			if(check_temp_th1 !=NULL){//good we have the file, load it up;
//				std::cout<<"\tYes, proceed to next directory."<<std::endl;
//				hist_root->Close();
//				continue;
//			}else{//no we don't have the histogram; then remove the direcory
//				std::cout<<"\tNo, reproduce histograms in current directory "<<Tdir_name<<std::endl;
//				Tdir->Delete(Tdir_name + ";*");
//			}
//		} else{//No such directory, create it;
		if(Tdir==NULL){
			if(verbose) std::cout<<"\tCreate directory "<<Tdir_name<<std::endl;
			Tdir = hist_root->mkdir(Tdir_name);
		}else{
			if(verbose)  std::cout<<"\tSwitch to directory "<<Tdir_name<<std::endl;
		}	
//		}
		Tdir->cd();

		//STEP 2.2 prepare TH1F
		Long64_t nentries = temptree->GetEntries();
		if(verbose) std::cout<<" We have entries: "<<nentries<<"   ";

		//STEP 2.2.1 prepare empty histgrams to be filled, each histogram with corresponding weight[kndex] (Note, not POT reweighting).
		std::vector<TH1F*> hist1d(num_throws);//hist
		std::vector<TTreeFormula*> wgtforms(num_throws);//variable
		std::vector<TTreeFormula*> varformula(num_throws);//weight, but not with POT reweighting;

		int count_hist = 0;
		for(int jndex = 0; jndex <num_throws;++jndex){//loop throws of each weight;
			//histogram of a throw
//			TString hist_name = (its_multithrows)? this->vars_name[index]+std::to_string(jndex+1) :this->vars_name[index];

			hist1d[jndex] = new TH1F(this->histNames[index][jndex],this->histNames[index][jndex], nbins, lowbinedge, highbinedge);//Will be rebinned later for variable binning;
			count_hist++;

			//variables of a throw
			varformula[jndex] = new TTreeFormula(var->unit.c_str(), gadget_updateindex(var->mininame, jndex), temptree);//get the x-axis
			varformula[jndex]->GetNdata();

			//weights of a throw
//			TString temp_label = (jndex>0)? "["+std::to_string(jndex)+"]":"";
			TString temp_wgt = gadget_updateindex(this->vars[index], jndex);//this->vars[index] + temp_label;
			TString temp_wgtname = this->vars_name[index]+std::to_string(jndex);

			wgtforms[jndex] = new TTreeFormula(temp_wgtname, temp_wgt+"*"+this->getStageCutsIndex(this->fstage,fbdt_cuts, jndex), temptree);//get weight with precuts;
			wgtforms[jndex]->GetNdata();
//			std::cout<<temp_wgt+"*"+this->getStageCutsIndex(this->stage,fbdt_cuts, jndex)<<std::endl;
			//					std::cout<< wgtforms[jndex]->EvalInstance()<<std::endl;
		}//next throw

		if(verbose) std::cout<<"Creating "<<count_hist<<" histograms "<<std::endl;

		//STEP 2.2.2 fill histograms!
		for(Long64_t kentry = 0; kentry< nentries; ++kentry){

			temptree->GetEntry(kentry);

			for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;

//				TTreeFormula varformula(var->unit.c_str(), gadget_updateindex(var->mininame, kndex), temptree);//get the x-axis
//				varformula.GetNdata();
				double temp_xvalue = varformula[kndex]->EvalInstance();
				double temp_weight = wgtforms[kndex]->EvalInstance();//*plot_pot/this->pot; do this when making cov matrix
//				std::cout<<"CHECK! "<<kndex<<" throws "<<wgtforms[kndex]->EvalInstance()<<" "<<plot_pot<<" "<<this->pot<<std::endl;

				hist1d[kndex]->Fill(temp_xvalue,temp_weight);//index - nth sw; kndex - nth throw

				if(temp_xvalue*temp_weight>0 && verbose){
					std::cout<<"\r>>Looping entries: "<<std::setw(5)<<kentry+1<<"/"<<std::setw(5)<<nentries;
					if(debug_verbose){
						std::cout<<" "<<std::setw(10)<<gadget_updateindex(var->mininame, kndex);
						std::cout<<"="<<std::setw(10)<<temp_xvalue;
						std::cout<<" "<<std::setw(10)<<gadget_updateindex(this->vars[index],kndex);
						std::cout<<"="<<std::setw(10)<<temp_weight;
						//					std::cout<<this->getStageCutsIndex(this->fstage,fbdt_cuts, kndex);
						std::cout<<std::setw(4)<<kndex+1<<" throw of SW: "<<this->vars_name[index];
					}

				}
			}//next throw;
		}

		if(verbose) std::cout<<"\n Finish filling "<<num_swvars<<" systematic weights."<<std::endl;

		//write hisograms to the diretory
		for(int jndex = 0; jndex <num_throws;++jndex){//throw in each weight;
			//configure th1
//			TString hist_name = (its_multithrows)? this->vars_name[index]+std::to_string(jndex+1) :this->vars_name[index];
			hist1d[jndex]->SetStats(0);
			hist1d[jndex]->GetXaxis()->SetTitle(var->unit.c_str());
			hist1d[jndex]->GetYaxis()->SetTitle("Events");

			//save the th1 to the root file, class bdt_sys;
			Tdir->cd();
			hist1d[jndex]->Write(this->histNames[index][jndex]);
			//This is important for making covariance matrix using hist2cov();
//			(this->hist[index]).push_back( (TH1F*) hist1d[jndex]->Clone());// each bdt_sys now contains all hist;
//			std::cout<<"Push back hist "<<hist_name<<std::endl;

			hist1d[jndex]->Delete();
		}

		if(verbose){
			std::cout<<"Save "<<num_throws;
			std::cout<<" histograms to the dir: "<<Tdir_name<<std::endl;
		}

		hist_root->Close();
	}//next systematic weight 
	infile->Close();
			

//	std::cout<<"\t  Now, load histograms for tag "<<this->tag<<std::endl;
//	TFile* temp_hist_root = (TFile*) TFile::Open(histfilename,"UPDATE"); //Frequent Open and Close, hope this help prevent corrupted links;
//	for(int index = 0; index < this->num_vars; ++index){//loop through weights
//		TString weights_nam = this->vars_name[index];//name of key for the histogram in root;
//		TString Tdir_name = this->tag+"_"+weights_nam;
//
//		for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
//			TString hist_nam_label = (its_multithrows)? std::to_string(kndex+1) : "";
//			TString hist_nam = weights_nam + hist_nam_label;
//			TH1F* temp_th1 =  (TH1F*) temp_hist_root->Get(Tdir_name+"/"+hist_nam);
//			if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
//				std::cout<<"Broken hist?? No way! Come check "<<__FILE__<<__LINE__<<std::endl;
//				exit(EXIT_FAILURE);
//			} else{
//
//				(this->hist[index]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
//				(this->hist[index]).back()->SetDirectory(0);//prevent the histogram from being deleted when file is closed. https://root-forum.cern.ch/t/use-of-th1-clone/11234/3
//				if(verbose) std::cout<<"\r\tLoading histogram "<<hist_nam<<" from dir: "<<Tdir_name;
//			}
//		}//next throws
//		if(verbose) std::cout<<std::endl;
//	}//next systematic weights;
//	temp_hist_root->Close();

}

/*
 * Build covariance matrices based on histograms;
 * syss contains all the 1d histograms;
 */

void sys_env::hist2cov( bdt_variable var, bool rescale, bool smooth_matrix){//, TString dir_root, TString dir_drawn, double plot_pot, unsigned long hashdraw){}
	

	bool checkbins = true;//going to reject non-rebinnable runs;
	bool force_rebin = false;//test different smoothing; true - bingap 25MeV->50MeV

	if(!rescale && verbose) std::cout<<"\n Warning: raw histograms are used, no normalization"<<std::endl;
	//STEP 0 Configuration
	
	int stage = fstage;

	int nb = var.n_bins;//n_bins is the target binning;
	int nl = var.plot_min;
	int nh = var.plot_max;
	bool do_rebin = var.is_custombin;//this will be updated, if the sys histogram binning can be used.


	//Set binnings for output histograms/covariance marix;
	std::vector<double> cur_binning;
	if(do_rebin){//binnings are set in the xml
		cur_binning = var.edges;
		cur_binning.erase(cur_binning.begin());//remove the first element, then we get the binning;
	} else{//binnings need to be calcualted
		double bin_left = nl; 
		double bingap = (nh-nl)/nb;

		for(int jndex = 0; jndex < nb+1; jndex++){
			cur_binning.push_back(bin_left);
			bin_left+=bingap;
		}
	}

	if(verbose) std::cout<<"\n\nMaking covariance matices."<<std::endl;
	//tags are already prepared as: tag_collection, tag2SWmap, tag2SWmap;
	//STEP 1 prepare root file to store  covariance matrix and configure TH2D*;

	TFile *cov_root = (TFile*) TFile::Open(root_dir+cov_prefix+".root","RECREATE"); //Collection of resulting covariance matrices;

	TString final_covroot_name = top_dir + final_prefix + ".root";
	TFile *finalcov_root = (TFile*) TFile::Open(final_covroot_name,"RECREATE");//one final covariance matrix output for plotting;

	//final covariance matrix that sums up all systematics;
	TH2D* finalcov =  new TH2D((var.safe_name).c_str(), (var.safe_name).c_str() ,nb,&(cur_binning).front(),nb,&(cur_binning).front());//binning for xnbins,xmin,xmax,ybins,ymin.ymax;

	TCanvas *drawCanvas = new TCanvas("drawCanvas", "", 600, 400);//Cavans to draw;

	//STEP 2 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	for(size_t index = 0; index < tag_collection.size(); ++index){}//each tag means one group of covariance matices;
	for(auto cur_tag : tag_collection){//each tag means one group of covariance matices;
		
		if(verbose) std::cout<<"\n-->Working on tag "<<cur_tag<<std::endl;

		//STEP 2.1 Find the bdt_sys syss that gives the CV; then find and draw up SW on the flight;
		bdt_sys* tempsCV = tag2CVmap.find(cur_tag)->second;//bdt_sys that gives CV - tempsCV

		bool do_smooth = ( (var.name).compare(3,5,"EnuQE")==0 )&& smooth_matrix && ((tempsCV->systag).Contains("Unisim"));

//		if(do_smooth) std::cout<<tempsCV->systag<<" might need smoothing as long as initial binning >1, current bins "<<var.int_n_bins<<std::endl;

		//Got the CV! 
		TH1F* cv_hist = (TH1F*) (tempsCV->hists[0][0])->Clone();
		if(force_rebin) cv_hist->Rebin(2);

		if(checkbins){//check syss's hisogram bin edges; it is good, if all cur_binning edges are found.

			int binchecker = 0;//monitor if there are difference between var binning and histogram binning;
			for(double binedge : cur_binning){
				for(int bindex = 1; bindex < cv_hist->GetNbinsX()+1; ++bindex){
					if(debug_verbose)std::cout<<cv_hist->GetBinLowEdge(bindex)<<" vs "<<binedge;
					if(abs(cv_hist->GetBinLowEdge(bindex) - binedge) < 10e-10){//edge matches;
						binchecker++;
						if(debug_verbose) std::cout<<"Good match"<<std::endl;
						break;
					}
						if(debug_verbose) std::cout<<std::endl;
				}
			}
			checkbins = false;//only need to check once.

			if(binchecker == nb+1){//the left edge does not count;
				if(debug_verbose) std::cout<<"Matched bin edges "<<binchecker<<"/"<<nb+1<<std::endl;
				if(binchecker < cv_hist->GetNbinsX()){
					if(verbose) std::cout<<"Histograms are going be rebinned: "<<cv_hist->GetNbinsX()<<" -> "<<nb<<std::endl;
					do_rebin = true;
				}
			} else{
				std::cout<<"Binning not adjustable: matched bin edges "<<binchecker<<"/"<<nb+1<<std::endl;
				std::cout<<"Please check the input hist_*.root."<<std::endl; 
				exit(EXIT_FAILURE);
			}
		}

		//Update cv_hist, if the binning are not exactly what we want;
//		TH1F* smoothRef_hist = (TH1F*) cv_hist->Clone();
		if(do_rebin){ 
			cv_hist->Rebin(nb, cur_tag+"CV",&(cur_binning).front());//tags[index]+"CV" is just some name, can be anything;
			cv_hist = (TH1F*) gDirectory->Get(cur_tag+"CV");
		}

		if(rescale)cv_hist->Scale(out_POT/tempsCV->pot);//scale hist to data POT

//		//STEP 2.2 Load SW and draw histograms, covariance matrix only the flight;
		int nby = 1.5*cv_hist->GetBinContent(cv_hist->GetMaximumBin());//set histogram height
//		int total_throws = 0;
//		std::vector< bdt_sys*> to_plot_list;
//
		for( std::_Rb_tree_iterator<std::pair<const TString, bdt_sys*> >  itr = tag2SWmap.begin();  itr!=tag2SWmap.end(); itr++){//group bdt_sys SW under the same tag
//
			if(itr->first == cur_tag){//identify the related bdt_sys based on bdtfiles;
				std::cout<<"Identify a file for tag "<<cur_tag<<std::endl;
				total_throws += (itr->second)->throws;
				to_plot_list.push_back( itr->second);//get the bdt_sys that will be used;
				if(to_plot_list.size()>1 && (itr->second)->vars_name.size()>1){
					std::cout<<"Now the code does not work on multi bdt_sys with multi systematic variables.";
					std::cout<<"Need to upgrade if needed, or just modify the xml file for now."<<std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		TH2D* all_hist = Export_allhist(cv_hist, to_plot_list[0]);
		if(size_t jndex = 1; jndex < to_plot_list.size(); jndex ++){

		}
//		if(to_plot_list.size()> 1) std::cout<<to_plot_list.size()<<" bdt_sys's are added under the tag "<<cur_tag<<std::endl;


//		//STEP 2.2.1 prepare SW of throws;
//			bdt_sys* tempsSW = itr->second;//tempSW may contains  manySW x throws;
//
//		bdt_sys* tempsSW = to_plot_list[0];

//		//for Unisims, there will be two elements 
//		int numSW = (tempsSW->hist).size();//hist is a 2d vector w dimension numSW x throws;
		for(int jndex = 0; jndex < numSW; ++jndex){//go through different SW under same bdt_sys
//			//USE THIS LABEL!
			TString temp_sw_name  = cur_tag+"_"+tempsSW->vars_name[jndex];

			TH2D* all_hist = new TH2D(temp_sw_name, temp_sw_name, nb, &(cur_binning).front(), nby, 0, nby);//to store many throws SW
			TH2D* covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix", temp_sw_name+"_CovarainceMatrix",nb,&(cur_binning).front(), nb,&(cur_binning).front());//to store covariance matrix, equal width;
//			TH2D* ori_covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix_original", temp_sw_name+"_CovarainceMatrix_original",nb,&(cur_binning).front(), nb,&(cur_binning).front());//to store covariance matrix, equal width;

//			for(int kndex = 0; kndex < tempsSW->throws; kndex++){//go through different throws
//				TH1F* sw_hist = (TH1F*) (tempsSW->hists[jndex][kndex])->Clone("copy"+tempsSW->vars_name[jndex]);
//				if(force_rebin) sw_hist->Rebin(2);
////				TString hist_name = (tempsSW->its_multithrows)? tempsSW->vars_name[jndex]+std::to_string(kndex) :tempsSW->vars_name[jndex];
////				TH1F* sw_hist = new TH1F(hist_name+"copy",hist_name+"copy",nb,nl,nh);//(tempsCV->hists[0][0])->Clone();
////				int set_1stsw_bin = 1;
////				tempsSW->hists[jndex][kndex]->Rebin(2);
////				for(int bindex = 1; bindex < tempsSW->hists[jndex][kndex]->GetNbinsX()+1; ++bindex){
////					if(tempsSW->hists[jndex][kndex]->GetBinLowEdge(bindex+1)>nl||tempsSW->hists[jndex][kndex]->GetBinLowEdge(bindex-1)<nh){//select some bins in between;
////						sw_hist->SetBinContent(set_1stsw_bin++, tempsSW->hists[jndex][kndex]->GetBinContent(bindex));
//////						std::cout<<"CHECK "<<tempsSW->hists[jndex][kndex]->GetBinLowEdge(bindex)<<" "<<tempsSW->hists[jndex][kndex]->GetBinContent(bindex)<<std::endl;
////					}
////				}
//				///scale, smooth, rebin;
//
//				if(var.int_n_bins>1&& do_smooth ){ 
//					std::cout<<"\r Smoothing a sw histogram";
//					*sw_hist = SmoothSW(sw_hist, smoothRef_hist, omhistogram);
//				}
//
//				if(do_rebin){ 
//					TString rebin_label =  cur_tag+"_"+tempsSW->vars_name[jndex]+"sw"+to_string_prec(kndex,0);
//					sw_hist->Rebin(nb,rebin_label,&(cur_binning).front());//repeat label will get the previous histogram
//					sw_hist = (TH1F*) gDirectory->Get(rebin_label);
//				}
//
//				if(rescale)sw_hist->Scale(plot_pot/tempsSW->pot);
//
////				if(nb>1&& do_smooth ){ 
////					std::cout<<"\r Smoothing a sw histogram";
////					*sw_hist = SmoothSW(sw_hist, cv_hist, false);//omhistogram);
////				}
//
//				//STEP 3.2.2 make histograms & covriance matrix;
////				double temp_sw=0;//save first 14 bin contents;
////				double temp_cv=0;
////				for(int lndex = 1; lndex < nb+1; ++lndex){//fill sw_hist to all_hist bins by bins;
////					all_hist->Fill(sw_hist->GetBinLowEdge(lndex) , sw_hist->GetBinContent(lndex) );
////					if(lndex<14){
////				//	std::cout<<"\rGet bin up to "<<sw_hist->GetBinLowEdge(lndex+1);
////					temp_sw += sw_hist->GetBinContent(lndex);
////					temp_cv += cv_hist->GetBinContent(lndex);
////					}
////				}//next bin
//
//
////				temp_covcal+=pow(temp_sw-temp_cv,2)/(double) total_throws;
//
//				//Make cov matrices
//				TH2D* cov_temp = MakeCov("cov"+cur_tag+"_"+tempsSW->vars_name[jndex]+std::to_string(kndex), sw_hist, cv_hist);
//				//make propagated covariance matrix;
//
//				if(omhistogram && adjOMErr){//now substract the stat error, before scaling the covariance matrix;
////					std::cout<<"\rAdjust cov matrix by substituting the stat error ";
//					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
//						for(int mndex = 0; mndex < nb; ++mndex){//fill cv_hist to all_hist bins by bins;
//							double temp_covvalue = cov_temp->GetBinContent(lndex+1,mndex+1);
//							double temp_updatedvalue = temp_covvalue - uw_statErr[lndex]*uw_statErr[mndex] +(plot_pot/tempsCV->pot)*sqrt(cv_hist->GetBinContent(lndex+1)*cv_hist->GetBinContent(mndex+1));
//							cov_temp->SetBinContent(lndex+1,mndex+1, temp_updatedvalue);
//
//							if(lndex==mndex){
////							std::cout<<"Adjust Cov statError "<<temp_covvalue<<"->"<<temp_updatedvalue<<",  uwstatErr^2: "<<uw_statErr[lndex]*uw_statErr[mndex];
////							std::cout<<", cv_hist err^2:"<<sqrt(cv_hist->GetBinContent(lndex+1)*cv_hist->GetBinContent(mndex+1))<<std::endl;
//							}
//						}
//					}
//				}
//
//
//				double scale_factor = 1.0;
//				if(tempsSW->vars_name[jndex].Contains("KpProd") ) scale_factor = 1.5;//from Zarko's MaktrixMaker xml file, see /e898-data/data60/desktop_backup/tirpitz/scratch/zarko/analysis2012/MatrixMaker/xml/matrix_combined_analysis.xml & line 413 of MatrixMaker.cxx
//
//				cov_temp->Scale(pow(scale_factor,2)/(double) total_throws);
//
//				covmatrices->Add(cov_temp);
////				std::cout<<covmatrices->GetNbinsX()<<" vs. "<<cov_temp->GetNbinsX()<<std::endl;
////				exit(0);
//
//				if(kndex == tempsSW->throws-1 && to_plot_list.size()>1){//wait, one more throw! Reset
//					to_plot_list.erase(to_plot_list.begin());
//					kndex = -1;//will be 0 next round;	
//					tempsSW = to_plot_list[0];
//					std::cout<<"\nWait! More throws to append: "<<tempsSW->vars_name[0]<<std::endl;
//				}
//
//				//fill ori_covmatrices, no matter  smoothing or not;
//				if(omhistogram && false){
//					TH1F* originsw_hist = (TH1F*) (tempsSW->hists[jndex][kndex])->Clone("ori"+tempsSW->vars_name[jndex]);
//					///scale, smooth, rebin;
//					originsw_hist->Scale(plot_pot/tempsSW->pot);
//					if(var.is_custombin){ 
//						TString rebin_label =  "org"+cur_tag+"_"+tempsSW->vars_name[jndex]+"sw"+to_string_prec(kndex,0);
//						originsw_hist->Rebin(nb,rebin_label,&(cur_binning).front());//repeat label will get the previous histogram
//						originsw_hist = (TH1F*) gDirectory->Get(rebin_label);
//					}
//
//					//STEP 3.2.2 make histograms & covriance matrix;
//					for(int lndex = 1; lndex < nb+1; ++lndex){//fill originsw_hist to all_hist bins by bins;
//						all_hist->Fill(originsw_hist->GetBinLowEdge(lndex) , originsw_hist->GetBinContent(lndex) );
//					}//next bin
//
//					//Make cov matrices
//					TH2D* ori_cov_temp = MakeCov("ori_cov"+cur_tag+"_"+tempsSW->vars_name[jndex]+std::to_string(kndex), originsw_hist, cv_hist);
//
//					double scale_factor = 1.0;
//					if(tempsSW->vars_name[jndex].Contains("KpProd") ) scale_factor = 1.5;//from Zarko's MaktrixMaker.cxx
//					ori_cov_temp->Scale(scale_factor/(double) total_throws);
//
//					ori_covmatrices->Add(ori_cov_temp);
//				}
//
//			}//next throws;
//			std::cout<<std::endl;
//
////			if(false){//do one bin analysis here..
////				std::cout<<"\nOne-Bin Analysis "<< temp_covcal<<std::endl;
////			}
//
//
//
//			//STEP 3.2.3 all_hist & covmatrices are filled; now draw them;
//			gStyle->SetPalette(kRainBow);//this change colz colors;
//
//	gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
//			//Histograms;
//			drawCanvas->Clear();
//			drawCanvas->cd();
//			all_hist->SetStats(false);
//			all_hist->GetXaxis()->SetTitle(var.unit.c_str());
//			all_hist->GetYaxis()->SetTitle("Event Rate");
//			all_hist->Draw("COLZ");
//			cv_hist->Draw("L same");	
//			cv_hist->SetLineColor(6);
//			cv_hist->SetLineWidth(2);
//			drawCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("hist",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
//
//			//Covariance matrix;
//			drawCanvas->Clear();
//			drawCanvas->cd();
//			//check for OpticalModel
//
//			//make fractional covariance matrix, but this does not added to final covariance matrix;
//			TH2D* fracCov = MakeFracCov(temp_sw_name, covmatrices, cv_hist, false, ori_covmatrices);//omhistogram, ori_covmatrices);//omhistogram - false, then the last argument is useless.
//			fracCov->SetStats(false);
//			fracCov->Draw("COLZ");
//			fracCov->SetTitle("Fractional "+temp_sw_name + " Covaraince Matrix");
//			drawCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covFrac",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
//
//	  gSystem->RedirectOutput(0,0);//this let ROOT warns
//			if(verbose&& false){ 
//				std::cout<<"The Fractional matrix diagonal element "<< temp_sw_name<<": "<<std::endl;;
//				for(int bindex = 1; bindex<nb+1; ++bindex){
//					std::cout<<" "<<fracCov->GetBinContent(bindex,bindex)<<",";
//				}
//				std::cout<<std::endl;
//			}
//
//
//			drawCanvas->Clear();
//			drawCanvas->cd();
//
//			if(sys2OMCVmap.count(tags[index]) > 0){//its a Optical Model CV, propagate the matrix to a fractional covariance matrix; and add it to finalcov.
//				bdt_sys* tempsOMCV = (sys2OMCVmap.find(cur_tag))->second;//get the sys for a given tag;
//				double omcv_pot = tempsOMCV->pot;
//				if(verbose) std::cout<<"Propagating covariance matrix for "<<cur_tag<<std::endl;
//				TH1F* omcv_hist = (TH1F*) (tempsOMCV->hists[0][0])->Clone();
//				if(force_rebin) omcv_hist->Rebin(2);
//
//				//rebin, scale
//				if(do_rebin){ 
//					omcv_hist->Rebin(nb, temp_sw_name+"OMCV",&(cur_binning).front());
//					omcv_hist = (TH1F*) gDirectory->Get(temp_sw_name+"OMCV");
//				}
//
//				if(!adjOMErr){//capture statistical error for OpticalModel
//					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
//						uw_statErr[lndex] = omcv_hist->GetBinError(lndex+1);
//						std::cout<<"Pre-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<omcv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
//					}
//				};
//				if(rescale) omcv_hist->Scale(plot_pot/omcv_pot);
//
//
//				if(!adjOMErr){//capture statistical error for OpticalModel
//					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
//						uw_statErr[lndex] = omcv_hist->GetBinError(lndex+1);
//						std::cout<<"Post-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<omcv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
//					}
//				};
//
//				//make propagated covariance matrix;
//				TH2D* ProCov = PropagateCov(temp_sw_name, fracCov, omcv_hist);
//
//				if(!adjOMErr){//now substract the stat error, before scaling the covariance matrix;
////					std::cout<<"\rAdjust cov matrix by substituting the stat error ";
//					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
//						for(int mndex = 0; mndex < nb; ++mndex){//fill cv_hist to all_hist bins by bins;
//							double temp_covvalue = ProCov->GetBinContent(lndex+1,mndex+1);
//							double temp_updatedvalue = temp_covvalue - uw_statErr[lndex]*uw_statErr[mndex] + sqrt(omcv_hist->GetBinContent(lndex+1)*omcv_hist->GetBinContent(mndex+1));
//							ProCov->SetBinContent(lndex+1,mndex+1, temp_updatedvalue);
//
//							if(lndex==mndex){
//							std::cout<<"Adjust Cov statError "<<temp_covvalue<<"->"<<temp_updatedvalue<<",  uwstatErr^2: "<<uw_statErr[lndex]*uw_statErr[mndex];
//							std::cout<<", omcv_hist err^2:"<<sqrt(omcv_hist->GetBinContent(lndex+1)*omcv_hist->GetBinContent(mndex+1))<<std::endl;
//							}
//						}
//					}
//				}
//
//				ProCov->SetStats(false);
//				ProCov->Draw("COLZ");
//				ProCov->SetTitle("Propagated "+temp_sw_name + " Covaraince Matrix");
//
//				gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
//				drawCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covProp",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
//				gSystem->RedirectOutput(0,0);//this let ROOT warns
//
//				if(verbose){ 
//					std::cout<< temp_sw_name<<" matrix diagonal element: "<<std::endl;;
//					double temp_calcul = 0;
//					for(int bindex = 1; bindex<nb+1; ++bindex){
//						std::cout<<" "<<ProCov->GetBinContent(bindex,bindex)<<",";
//						if(ProCov->GetBinContent(bindex,bindex)>0&&bindex<14) temp_calcul+=ProCov->GetBinContent(bindex,bindex);
//					}
//					std::cout<<std::endl;
////					std::cout<<"Quadrature sum: "<<temp_calcul<<std::endl;
//				}
//
//				finalcov->Add(ProCov);
//
//				cov_root->cd();
//				fracCov->Write();
//				ProCov->Write();
//			} else{//no need extra processing;
//
//				covmatrices->SetStats(false);
//				covmatrices->GetXaxis()->SetTitle((var.unit).c_str());
//				covmatrices->Draw("COLZ");
//				gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
//				drawCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("cov",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
//				gSystem->RedirectOutput(0,0);//this let ROOT warns
//
//				if(verbose){ 
//					std::cout<<"The matrix diagonal element "<<temp_sw_name<<": "<<std::endl;
//					double temp_calcul = 0;
//					for(int bindex = 1; bindex<nb+1; ++bindex){
//						std::cout<<" "<<covmatrices->GetBinContent(bindex,bindex)<<",";
//						if(covmatrices->GetBinContent(bindex,bindex)>0&&bindex<14) temp_calcul+=covmatrices->GetBinContent(bindex,bindex);
//					}
//					std::cout<<std::endl;
////					std::cout<<"Quadrature sum: "<<temp_calcul<<std::endl;
//				}
//			
//				finalcov->Add(covmatrices);
//				cov_root->cd();
//				covmatrices->Write();
//			}
//			//one covmatrix is made;
//		}//next SW
//	}//next tag;
//
//	//STEP 4 draw out the final covariance matrix;
//
//	cov_root->Close();//basic covariance matrix is recorded, now prepare the one for plotting;
//	
//	//set axis and name for the final covaraince matrix;
//	finalcov->SetStats(false);
//
//	finalcov->SetTitle("Total Covariance Matrix for");
//	finalcov->GetXaxis()->SetTitle( (var.unit).c_str());
//
//	TCanvas *c1 = new TCanvas("c1","",600,400);
//	finalcov->Draw("COLZ");
//
//	c1->SaveAs(dir_drawn+"/"+gadget_labelroot("Totalcov", var, stage, hashdraw)+".pdf","pdf");
//
//	finalcov_root->cd();
//	finalcov->Write();
//
//
//	//turn the finalcov in the format of a matrix
//	TMatrixD covMatrix(nb,nb); 
//	TArrayD nums(nb*nb);
//
//	for(int index = 0; index < nb; ++index){
//		for(int jndex = 0; jndex < nb; ++jndex){
//			nums[jndex*nb+index] = finalcov->GetBinContent(index+1,jndex+1);
//		}
//	}
//	covMatrix.SetMatrixArray(nums.GetArray());
//	covMatrix.Write((var.safe_name+"matrix").c_str());
//
//
//	finalcov_root->Close();
//	std::cout<<"Final total covariance matrix is saved at >>>>   "<<final_covroot_name<<std::endl;
}


TH2D* MakeCov(TString name,TH1F* hist, TH1F* cv){

	int nb=hist->GetNbinsX();

	//	int nl=hist->GetBinLowEdge(1);
	//	int nh=hist->GetBinLowEdge(nb+1);
	//in case of customized binnings; use vector..
	std::vector<double> bins;
	for( int index = 1; index < nb+2; index++){

		bins.push_back(hist->GetBinLowEdge(index));//1st bin to n+1th bin
	}

	int counter = 0;
	TH2D* covmatrix =  new TH2D(name,name,nb, &(bins.front()),nb, &(bins.front()) );//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
	for(int index = 1; index<nb+1; ++index){
		for(int jndex = 1; jndex<nb+1; ++jndex){
			double entry = (hist->GetBinContent(index)-cv->GetBinContent(index) )*(hist->GetBinContent(jndex)-cv->GetBinContent(jndex) );
			covmatrix->SetBinContent(index,jndex,entry);

		}
	}

	return covmatrix;
}

/*
 * Create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TString name,TH2D* inputcov, TH1F* oldcv, bool ec_smoothing, TH2D* origin_cov){
	
	int nb=oldcv->GetNbinsX();
	int nl=oldcv->GetBinLowEdge(1);
	int nh=oldcv->GetBinLowEdge(nb+1);
	
	TH2D* fractional_cov = (TH2D*) inputcov->Clone("Fractional_Covaraince_Matrix_"+name);
//	fractional_cov->Reset("ICESM");//clear contents;

	for(int jndex = 1; jndex < nb+1; ++jndex){
		for(int kndex = 1; kndex < nb+1; ++kndex){
			double cvj = oldcv->GetBinContent(jndex);
			double cvk = oldcv->GetBinContent(kndex);
			double covjk = inputcov->GetBinContent(jndex,kndex);
			double fjk = covjk/(cvj*cvk);//not divide by nb;
//			double fjk = (ec_smoothing && jndex==kndex)? (covjk+cvj)/(cvj*cvk): covjk/(cvj*cvk);//not divide by nb;

//			if(ec_smoothing&&jndex==kndex && fjk > origin_cov->GetBinContent(jndex,kndex)/(cvj*cvk) ){
//				std::cout<<"\nAdjust matrix element to a smaller value: "<<jndex<<" "<<kndex<<": "<<fjk<<" to ";
//				fjk = origin_cov->GetBinContent(jndex,kndex)/(cvj*cvk);
//				std::cout<<fjk<<std::endl;
//			}

			fractional_cov->SetBinContent(jndex,kndex,fjk);
		}
	}

	return fractional_cov;
}

/*
 * Propagate matrix
 */
TH2D* PropagateCov(TString name,TH2D* frac_cov, TH1F* newcv){

	int nb=newcv->GetNbinsX();
	int nl=newcv->GetBinLowEdge(1);
	int nh=newcv->GetBinLowEdge(nb+1);
	
	TH2D* Pro_cov = (TH2D*) frac_cov->Clone("Propagated_Covaraince_Matrix_"+name);

	for(int jndex = 1; jndex < nb+1; ++jndex){
		for(int kndex = 1; kndex < nb+1; ++kndex){
			double fjk = frac_cov->GetBinContent(jndex,kndex);

			double cvj_new = newcv->GetBinContent(jndex);
			double cvk_new = newcv->GetBinContent(kndex);
			Pro_cov->SetBinContent(jndex,kndex,fjk*cvj_new*cvk_new);
		}
	}

	return Pro_cov;
}

/*
 * Smooth the matrix
 */
//http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
int global_index = 0;
TH1F SmoothSW(TH1F* sw, TH1F* cv, bool special_option){//fit each bin with polynomials

	using namespace std;
	//use the special option for the old (but effective) smoothing strategy
	
	bool print_ratio = false;
	bool ec_smoothing = false;

//	std::vector<double> sm_binning ={200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1250, 1300, 1500, 1900, 3000};
	std::vector<double> sm_binning ={200, 250, 300, 375, 475, 550, 600, 675, 750, 800, 950, 1100, 1150, 1250, 1300, 1500, 1700, 1900, 3000};

	TString cur_tag = to_string(sw->Integral())+"Unisims";

	if(special_option){
		cv->Rebin(sm_binning.size()-1, cur_tag+"smCV",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
		cv = (TH1F*) gDirectory->Get(cur_tag+"smCV");

		sw->Rebin(sm_binning.size()-1, cur_tag+"smSW",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
		sw = (TH1F*) gDirectory->Get(cur_tag+"smSW");
	}

	Int_t degree = 1;
	Int_t first_bin = 1;
	Int_t last_bin = sw->GetNbinsX();//1;
//	while(sw->GetBinLowEdge(last_bin)<1275){//dont smooth bins with EnuE>1899
//		last_bin++;
//	}
	Int_t Nbins=last_bin-first_bin+1;


	if(ec_smoothing){ //use en-chuan's 

		gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
		TCanvas *rCanvas = new TCanvas("rCanvas", "", 600, 400);//Cavans to draw;
		TF1* smoothFitter = new TF1("smoothFitter", "pol4", cv->GetBinLowEdge(1), cv->GetBinLowEdge(Nbins+1));
		TH1F* copysw = (TH1F*) sw->Clone();
		copysw->Divide(cv);
		//debugging
		if(print_ratio){
			copysw->Draw();
		}

		copysw->Fit("smoothFitter","+","",cv->GetBinLowEdge(1)-0.5,cv->GetBinLowEdge(Nbins+1)+0.5);

		sw->Reset();
		for (Int_t bin=1;bin<Nbins+1;bin++) {//Here modify the sw
			sw->SetBinContent(bin, cv->GetBinContent(bin) * smoothFitter->Eval(sw->GetBinCenter(bin)));
		}


		TH1F* copy_rsw = (TH1F*) sw->Clone();
		copy_rsw->Divide(cv);
		copy_rsw->SetLineColor(kRed);
		copy_rsw->Draw("same");

		if(print_ratio){  
			gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
			rCanvas->SaveAs(("/scratch/condor-tmp/klin/hellstroms_hive/hive/build/src/fullmcsystematics/drawn/ratio_OM_"+to_string_prec(global_index,0)+".pdf").c_str(),"pdf");
			gSystem->RedirectOutput(0,0);//this let ROOT warns
			rCanvas->Delete();
			global_index++;
		}

	}else{//use Zarko's
		//solve M a = y
		
		//proceed to smoothing

		cout << "Smoothing bins " << first_bin << " - " << last_bin << " (N=" << Nbins<< ", "<<cv->GetBinLowEdge(last_bin)<<") with pol "<<degree<<endl;
		for (Int_t bin=0;bin<Nbins;bin++) {
			//    cout <<cv->GetBinContent(bin+first_bin)<<"\t"<<sw->GetBinContent(bin+first_bin)<<endl;
		}
		Double_t M[degree+1][degree+1];
		Double_t y[degree+1];
		for (Int_t j=0;j<degree+1;j++) {
			y[j]=0.;
			for (Int_t k=0;k<degree+1;k++) {
				M[j][k]=0.;
			}
		}
		for (Int_t j=0;j<degree+1;j++) {
			for (Int_t i=0;i<Nbins;i++) {
				if ( cv->GetBinContent(i+first_bin)>0 )
					y[j] += pow(Double_t(i+1),Double_t(j))*(sw->GetBinContent(i+first_bin)/cv->GetBinContent(i+first_bin));
				else
					y[j] += pow(Double_t(i+1),Double_t(j));
			}
		}

		for (Int_t j=0;j<degree+1;j++) {
			for (Int_t k=0;k<degree+1;k++) {
				for (Int_t i=0;i<Nbins;i++) {
					M[j][k] += pow(Double_t(i+1),j)*pow(Double_t(i+1),k);
				}
			}
		}

		TVectorD va(degree+1);
		TMatrixD mM(degree+1,degree+1);
		TVectorD vy(degree+1);

		for (Int_t j=0;j<degree+1;j++) {
			vy[j]=y[j];
			for (Int_t k=0;k<degree+1;k++) {
				mM[j][k]=M[j][k];
			}
		}

		TDecompLU lu(mM);
		lu.Decompose();
		Bool_t ok;
		va=lu.Solve(vy,ok);

		if ( ok ) {
			//    cout <<"success"<<endl;
			for (Int_t i=0;i<Nbins;i++) {
				Double_t num=0.;
				for (Int_t j=0;j<degree+1;j++) {
					num += va[j]*pow(Double_t(i+1),Double_t(j)); 
				}
				num *= cv->GetBinContent(i+first_bin);
//				      cout <<"smoothed "<<i+first_bin<<"\t"<<num<<endl;
				sw->SetBinContent(i+first_bin,num);
			}
		}
		else cout <<"failed to solve"<<endl;
	}
		//calculte chi^2 on sw;
	
		double chi2 = 0;
		for (Int_t bin=1;bin<Nbins+1;bin++) {//Here modify the sw
			if(cv->GetBinContent(bin)>0){
			chi2+=pow(sw->GetBinContent(bin)-cv->GetBinContent(bin), 2)/cv->GetBinContent(bin);
			}
		}
		double AIC = chi2+ 2*degree+ (2*degree)*(degree+1)/(Nbins-degree-1);
		cout<<"Check smoothed result:\n Chi2: "<<chi2<<" AIC = "<<AIC<<" with degree "<<degree<<std::endl;
		//the minimum AIC gives the best fitting result;


	return *sw;


}

/*
 * Make Covaraince matrix from 2d histograms
 *
 */

TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv){

}

/*
 * Main code, initialize systematics;
 * Check if the histograms exist, if yes, then can save some works.
 */

void sys_env::InitSys(std::vector<bdt_variable> vars, std::vector<bdt_sys*> syss){
	bool do_cov = true; //true - generate covariance matrix everytime, no matter what; false -this boolean is not functioning
	bool do_hist = false; //true - generate hist everytime, no matter what; false -this boolean is not functioning

	bool skip_sysEvaluation = false;//true - do nothing;

//	bool verbose = this->verbose;
//	bool debug_verbose = this->debug_verbose;
	

	//the histogram goes two ways: 1dhist or 2dhist; 1dhist for now;
	bdt_variable* var1 = &vars[0];
	bdt_variable* var2;
	if(vars.size()>1) var2 = &vars[1];

	this->checkEnv();//if triggered, need to set things up outside (i.e. hive.cxx)
	
	//set the prefix names for plots and roots storing plots:
	//they go like <type>_<var>_bgap<gap>_s<stage>_<hash>
	this->hist_prefix = this->getFileName("hist",vars); 
	this->cov_prefix = this->getFileName("cov",vars); 
	this->final_prefix = this->getFileName("finalcov",vars); 

	//Before we start, check if 1) finalcov*.root exist,  2) hist*.root exist;

	std::cout<<"\n ---- Initialize systematics for "<<var1->unit;
	std::cout<<"\tSystematic files tag: s"<<this->fstage<<"_"<<fcut_hash<<std::endl;

	//check 1) for the finalcov*root;
	TString check_covroot_name = var1->covar_file.c_str();
	
	//check if it is good;
	TFile *check_cov_root = (TFile*) TFile::Open(check_covroot_name,"READ"); //one var one hist root;

	if(check_cov_root != NULL ){ //looks good

		if(verbose){
			std::cout<<"We found fractional covaraince matrices root: "<<check_covroot_name<<std::endl;
			std::cout<<"Skip making cov matrices."<<std::endl;
		}
		skip_sysEvaluation = true;
	}else{//no such root;
		if(verbose){ 
			std::cout<<"Did not found fractional covaraince matrices root: "<<check_covroot_name<<std::endl;
			std::cout<<"Will generate matices in "<<(this->top_dir+this->final_prefix + ".root")<<std::endl;
		}
	}
	check_cov_root->Close();
	delete check_cov_root;

	check_covroot_name = (this->top_dir+this->final_prefix + ".root");

	//check 2) for the 1d histgrams, if 1) pass;
	if(do_cov || !skip_sysEvaluation){//true - explore further, false - exit the function.

		TString histfilename = this->root_dir + this->hist_prefix + ".root";
		TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;

		if(verbose){
			if(do_cov) std::cout<<"(^ 3 ^) I am kidding, we still make cov matrices for "<<check_covroot_name<<std::endl;
			std::cout<<"Loading up "<<histfilename<<"\n"<<std::endl;
		}

		if(hist_root == NULL){
			hist_root = (TFile*) TFile::Open(histfilename,"RECREATE");
		}
		//Group bdt_sys, and load histograms; mark empty bdt_sys;
		bool all_good = true;
		for(auto cur_sys : syss){
			TString temptag = cur_sys->tag;

			//classify bdt_sys
			(this->tag_collection).push_back(temptag);
//			if(verbose) std::cout<< temptag<<" is "<<std::endl;

			if(cur_sys->its_CV){//save CV to map;
				tag2CVmap.insert(std::make_pair (temptag, cur_sys));
//				if(verbose) std::cout<<"CV!";

			} else{//save SW to map
				tag2SWmap.insert(std::make_pair (temptag, cur_sys));
//				if(verbose) std::cout<<"Systematic weight!";
			}

			cur_sys->catchupEnv(this);

			if(do_hist) cur_sys->Make1dhist(var1);

			if( !(cur_sys->Load1dhist(hist_root) )) all_good = false ;//mark!;
		}

		if(!all_good){//now need to Make histograms
			hist_root->Close();
			hist_root = (TFile*) TFile::Open(histfilename,"UPDATE");
			for(auto cur_sys : syss){
				if(!cur_sys->fullyloaded){
					cur_sys->Make1dhist(var1);
					cur_sys->Load1dhist(hist_root);
				}
			}
		}

		std::cout<<"Finish preparing systematic histograms, and now they are in\n>>>> "<<histfilename<<std::endl;

		//Part II: Overlay histograms & draw covaraicne matrix;

		//but first, clean up this->tag_collection;
		sort(this->tag_collection.begin(), this->tag_collection.end() );
		this->tag_collection.erase( unique( this->tag_collection.begin(), this->tag_collection.end()), this->tag_collection.end() );

		bool rescale_out_hist = true;
		bool smooth_unisims = true;
		hist2cov( *var1, rescale_out_hist, smooth_unisims);//tag_collection, tag2CVmap, tag2SWmap pointed to bdt_sys.
		hist_root->Close();
	}

///////////////////////
/////////OLD//////////
//////////////////////
//	TString readfile_method = "RECREATE";//create file by default, if exist, update the current file;
//	TString histfilename = dir_root+"/"+gadget_labelroot("hist", *var1, stage, hashlabel)+".root";
//
//	if(do_hist || !skip_sysEvaluation){//no root files for the covariance matix, so going to generate it.
//		std::cout<<"\nLooking into histograms generation"<<std::endl;
//		TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;
//
//		size_t index = 0;//index for syss, determine where we left over last time the code break;
//
//		if(hist_root!=NULL && !do_hist ){ //quick check of existing files; 
//			std::cout<<histfilename<<" exist! Now check what is inside."<<std::endl;
//
//			bool found_corrupted_hist = false;
//			while( index < syss.size()){
//				//CHECK, might need to update particular directory;
//
//				bdt_sys* temps = syss[index];
//				int num_throws = temps->throws;
//				bool its_multithrows = temps->its_multithrows;//multi weights or optical model
//
//				for(int jndex = 0; jndex < temps->num_vars; ++jndex){//loop through weights
//					TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
//					TString Tdir_name = temps->tag+"_"+weights_nam;
//
//					for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
//						TString hist_nam_label = (its_multithrows)? std::to_string(kndex+1) : "";
//						TString hist_nam = weights_nam + hist_nam_label;
//						TH1F* temp_th1 =  (TH1F*) hist_root->Get(Tdir_name+"/"+hist_nam);
//						if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
////							temps->start_with_weight = jndex;
//
////							readfile_method = "UPDATE";
//							hist_root->Close();
//							std::cout<<"\n>> Find a broken directory, when looking at "<<syss[index]->tag<<". Update it now."<<std::endl;
//
//							goto now_make_hist;//do partial generation, index determines which systematic file to start;
//						} else{
//
//							(temps->hists[jndex]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
//							if(verbose) std::cout<<"\r\tChecking dir: "<<Tdir_name<<" histogram: "<<hist_nam;
//						}
//
//					}//next throws
//					if(verbose) std::cout<<std::endl;
//				}//next systematic weights;
//				index++;
//			}//next systematics
//		} else{
//
//			std::cout<<histfilename<<" does not exist. Create it now."<<std::endl;
//
//now_make_hist://start ith syss (index), jth weight(start_with_weight), from the given index; 
////			hist_root = (TFile*) TFile::Open(histfilename,readfile_method); //one var one hist root;
//			
////			if(readfile_method == "UPDATE"){ 
////				hist_root->Recover();//when use update, its better to recover the file first.
////				hist_root->mkdir("TRASH");
////				hist_root->Delete("TRASH;*");//make a directory and delete it, then no recovery is needed.
////				hist_root->Write();
////				hist_root->Close();
////				hist_root = (TFile*) TFile::Open(histfilename,readfile_method); //one var one hist root;
////			}
////			std::cout<<"Working on "<<histfilename<<std::endl;
//
//			while( index < syss.size()){//Histogram making is here!! 
//				//configure Make1dhist process
//
//				syss[index]->Make1dhist(histfilename, var1, plot_pot, fbdt_cuts);
//				index++;
//			}//next systematic files
//		}
//
//		std::cout<<"Finish making systematic histograms, and save them to\n>>>> "<<histfilename<<std::endl;
//		//
//		//Part II: Overlay histograms & draw covaraicne matrix;
//		//
//		hist2cov( *var1, syss, dir_root, dir_drawn, plot_pot, hashlabel);
//		hist_root->Close();
//	}

	std::cout<<"Finish making covariance matrix, see \n>>>>  "<<check_covroot_name<<std::endl;
	std::cout<<"\n ----------------------------- \n"<<std::endl;
}
