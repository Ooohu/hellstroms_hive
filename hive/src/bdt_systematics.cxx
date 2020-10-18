#include <bdt_systematics.h>

int global_index = 0;

/*
 * Some gadgets to keep things organized;
 */
TH1D gadget_vrebin(TH1D* hist, std::vector<double > binning){

	if(true){
		double hb1 = hist->GetBinLowEdge(1);
		double hb2 = hist->GetBinLowEdge(hist->GetNbinsX()+1);
		if(binning.front() < hb1 ||binning.back() > hb2){
			std::cout<<"\nReject to rebin underflow bins: Hist is rebinned from ("<<hb1<<","<<hb2;
			std::cout<<") to ("<<binning.front() <<","<<binning.back() <<"). Please check bin edges in xml and hist*.root."<<std::endl;
			return *hist;
		}
	}
	hist->Rebin(binning.size() -1, ("tmp"+to_string_prec(global_index)).c_str() ,&(binning).front());//tags[index]+"CV" is just some name, can be anything;
	hist = (TH1D*) gDirectory->Get(("tmp"+to_string_prec(global_index++)).c_str());

	return *hist;
}

/*
 * Check bins of cv_hist, see if it is handlable;
 *0 - no problem; (1,2,4) - (1,2,4)
 *1 - Rebinable: cv_hist out of range (1,2,3,4) - (2,3)
 *2 - Rebinable: (1,2,4) - (1,4)
 *3 - BAD: underflow bins: output_binning range too large (1,2,4) - (0,1,4)
 *4 - BAD: bins edge not found in cv_hist (1,2,5) - (1,2,3,4)
 */
int gadget_BinMatcher(TH1D* cv_hist, std::vector< double > output_binning){

			bool debug_message = false;

			int binchecker = 0;//monitor if there are difference between var binning and histogram binning;
			int mis_count = 0;
			int edges = output_binning.size();
			int histnb = cv_hist->GetNbinsX();
			int bindex_starter = 1;
			for(double binedge : output_binning){//go through each bin in output_binning;
				for(int bindex = bindex_starter; bindex < histnb+2; ++bindex){//loop through lowegdes of cv_hist 
					if(debug_message)std::cout<<cv_hist->GetBinLowEdge(bindex)<<" vs "<<binedge;
					if(abs(cv_hist->GetBinLowEdge(bindex) - binedge) < 10e-10){//edge matches;
						binchecker++;
						bindex_starter = bindex+1;
						if(debug_message) std::cout<<"Good match"<<std::endl;
						break;
					}
					if(binchecker>1 && binchecker < edges) mis_count++;
					if(debug_message) std::cout<<std::endl;
				}
			}

			if(binchecker == edges ){ 
				if(edges - 1 == histnb ) return 0; //perfect
				if(mis_count < 1) return 1;
				std::cout<<mis_count<<std::endl;
				return 2; //rebinnable; checker is less than the NbinsX();
			}
			
			bool lb = output_binning.front() < cv_hist->GetBinLowEdge(1);
			bool hb = output_binning.back() > cv_hist->GetBinLowEdge(1+ histnb);
			if(lb||hb) return 3;//underflow bins; 

			return 4;
}

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



//generate proper name label from variable;
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
	
	//bdt_sys has boolean: loaded & fullyloaded 
	//it is loaded, unless NULL TH1D is found.
	//fullyloaded is false, if at least one NULL is found;
	this->fullyloaded = true;
	for(int index = 0; index< this->num_vars; ++index){//loop through weights;
//		if(this->loads[index]) continue;//Cant skip this, if file close, all hists are gone.
		this->loads[index] = true;

		(this->hists[index]).clear();//make sure it is empty before filling in;
		for(int jndex = 0; jndex< this->throws; ++jndex){//loop through throws;
			TString nth_label = (this->its_multithrows)? std::to_string(jndex) : "";
			TString hist_name = this->vars_name[index]+nth_label;

			TH1D* temp_th1 =  (TH1D*) cur_file->Get(this->TdirNames[index]+"/"+hist_name);

			if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
				if(verbose) std::cout<<"\t"<<this->tag+"_"+hist_name+" does not exist. Next";
				this->loads[index] = false;
				this->fullyloaded = false;
				break;
			}


			(this->hists[index]).push_back( (TH1D*) temp_th1->Clone() );// load it!
			if(verbose) std::cout<<"\r\tLoad Tdir: "<<this->TdirNames[index]<<" histogram: "<<hist_name;
		}//next throw
		if(verbose) std::cout<<std::endl;
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
void bdt_sys::Make1dhist(bdt_variable *var, TFile* out_root){//TString histfilename, bdt_variable* var,  double plot_pot, std::vector<double> bdt_cuts){}

	
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

	//STEP 2 prepare TH1D
//	for(int index = this->start_with_weight; index < this->num_vars; ++index){//loop through weights}
	for(int index = 0 ; index < this->num_vars; ++index){//loop through weights

		if(this->loads[index]) continue;//well, this means no need make histograms;
		out_root->cd();
//		TFile* out_root = (TFile*) TFile::Open(this->sys_env::root_dir + hist_prefix + ".root","UPDATE"); //Frequent Open and Close, hope this help prevent corrupted links;
//		std::cout<<"Updating "<<sys_env::root_dir + sys_env::hist_prefix + ".root"<<std::endl;
		//subdirectory in root for each systematic weighs; start_with_weight labels which weight to start, i.e. index;
		//STEP 2.1
		//second check of directory, because there might be some leftover from previous run;
		
		TString Tdir_name = TdirNames[index];//this->tag+"_"+weights_nam;

		TDirectory* Tdir = out_root->GetDirectory(Tdir_name);

		if(Tdir==NULL){
			if(verbose) std::cout<<"\tCreate directory "<<Tdir_name<<std::endl;
			Tdir = out_root->mkdir(Tdir_name);
		}else{
			if(verbose)  std::cout<<"\tSwitch to directory "<<Tdir_name<<std::endl;
		}	
//		}
		Tdir->cd();

		//STEP 2.2 prepare TH1D
		Long64_t nentries = temptree->GetEntries();
		if(verbose) std::cout<<" We have entries: "<<nentries<<"   ";

		//STEP 2.2.1 prepare empty histgrams to be filled, each histogram with corresponding weight[kndex] (Note, not POT reweighting).
		std::vector<TH1D*> hist1d(num_throws);//hist
		std::vector<TTreeFormula*> wgtforms(num_throws);//variable
		std::vector<TTreeFormula*> varformula(num_throws);//weight, but not with POT reweighting;

		int count_hist = 0;
		for(int jndex = 0; jndex <num_throws;++jndex){//loop throws of each weight;
			//histogram of a throw
//			TString hist_name = (its_multithrows)? this->vars_name[index]+std::to_string(jndex+1) :this->vars_name[index];

			hist1d[jndex] = new TH1D(this->histNames[index][jndex],this->histNames[index][jndex], nbins, lowbinedge, highbinedge);//Will be rebinned later for variable binning;
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
				double temp_weight = (wgtforms[kndex]->EvalInstance()>29)? 0 : wgtforms[kndex]->EvalInstance();//No large weight..
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

			hist1d[jndex]->Delete();
		}

		if(verbose){
			std::cout<<"Save "<<num_throws;
			std::cout<<" histograms to the dir: "<<Tdir_name<<std::endl;
		}

	}//next systematic weight 
	infile->Close();
			

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
	double nl = var.plot_min;
	double nh = var.plot_max;
	bool do_rebin = var.is_custombin;//this will be updated, if the sys histogram binning can be used.


	//Set binnings for output histograms/covariance marix;
	if(do_rebin){//binnings are set in the xml
		output_binning = var.edges;
		output_binning.erase(output_binning.begin());//remove the first element, then we get the binning;
	} else{//binnings need to be calcualted
		double bin_left = nl; 
		double bingap = (nh-nl)/nb;

		for(int jndex = 0; jndex < nb+1; jndex++){
			output_binning.push_back(bin_left);
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
	TH2D* finalcov =  new TH2D((var.safe_name).c_str(), (var.safe_name).c_str() ,nb,&(output_binning).front(),nb,&(output_binning).front());//binning for xnbins,xmin,xmax,ybins,ymin.ymax;

	TCanvas *drawCanvas = new TCanvas("drawCanvas", "", 600, 400);//Cavans to draw;

	//STEP 2 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	//	for(size_t index = 0; index < tag_collection.size(); ++index){}//each tag means one group of covariance matices;
	for(auto cur_tag : tag_collection){//each tag means one group of covariance matices;

		if(verbose) std::cout<<"\n-->Working on tag "<<cur_tag<<std::endl;

		//STEP 2.1 Find the bdt_sys syss that gives the CV; then find and draw up SW on the flight;
		bdt_sys* tempsCV = tag2CVmap.find(cur_tag)->second;//bdt_sys that gives CV - tempsCV

		bool do_smooth = ( (var.name).compare(3,5,"EnuQE")==0 )&& smooth_matrix && ((tempsCV->systag).Contains("Unisim"));

		//		if(do_smooth) std::cout<<tempsCV->systag<<" might need smoothing as long as initial binning >1, current bins "<<var.int_n_bins<<std::endl;

		//Got the CV! 
		TH1D* cv_hist = (TH1D*) (tempsCV->hists[0][0])->Clone();
		if(force_rebin) cv_hist->Rebin(2);

		if(checkbins){//check syss's hisogram bin edges; it is good, if all output_binning edges are found.

			checkbins = false;
			int matching_code = gadget_BinMatcher(cv_hist, output_binning);
			/*0 - no problem; (1,2,4) - (1,2,4)
			 *1 - Rebinable: cv_hist out of range (1,2,3,4) - (2,3)
			 *2 - Rebinable: (1,2,4) - (1,4)
			 *3 - BAD: underflow bins: output_binning range too large (1,2,4) - (0,1,4)
			 *4 - BAD: bins edge not found in cv_hist (1,2,5) - (1,2,3,4)
			 */
			if(matching_code > 2){	
				std::cout<<"BinMatcher Code "<<matching_code<<std::endl;
				std::cout<<"Please check the input the binning of hist_*.root."<<std::endl; 
				exit(EXIT_FAILURE);
			} else if( matching_code > 0){
				if(verbose) std::cout<<"Histograms are going be rebinned: "<<cv_hist->GetNbinsX()<<" -> "<<nb<<std::endl;
				do_rebin = true;
			}
		}

		//Update cv_hist, if the binning are not exactly what we want;
		if(rescale) cv_hist->Scale(out_POT/tempsCV->pot);//scale hist to data POT
		TH1D* smoothRef_hist = (TH1D*) cv_hist->Clone();
		if(do_rebin) *cv_hist = gadget_vrebin(cv_hist, output_binning);


		//		//STEP 2.2 Load SW and draw histograms, covariance matrix only the flight;

		//sort out bdt_sys SW under the same tag, breakdown SW, mainly for MCUnisims: Qtcorr&disc
		std::map< TH1D*, double > hist2tscalemap;//one set of hists for one Throw scaling;
		std::vector< TString > hists_names;//(SW1,SW2,...)
		std::vector< std::vector<TH1D*> > sw_hists_coll;//dimension: bdt_file x SW
		int file_counter = 0;//2 bdt_file means these two are added up as one matrix (independent, but just in the same plot);
		for( std::_Rb_tree_iterator<std::pair<const TString, bdt_sys*> >  itr = tag2SWmap.begin();  itr!=tag2SWmap.end(); itr++){
			if(itr->first == cur_tag){//identify the related bdt_sys under the same tag <bdt_sys+bdt_file>
				file_counter++;
				bdt_sys* temp_sys = itr->second;
				int tmpN = temp_sys->num_vars;

				double tmp_pscale = out_POT/temp_sys->pot;
				double tmp_tscale = (1.0/temp_sys->throws);
				for(int index = 0; index < tmpN; ++index){

				TString temp_varName = temp_sys->vars_name[index];
				if(verbose) std::cout<<"Got "<<temp_sys->tag+" "+temp_varName<<std::endl;
					std::vector<TH1D*> tmp_h = temp_sys->hists[index];
					if(file_counter<2){//others
						//exp:  hists1,hists2,...
						hists_names.push_back( temp_varName);
						sw_hists_coll.push_back(tmp_h);
					}else{//MCUnisims, append histograms to current SW;
						//exp:  hists1,hists2,
						//      hists3,hists4,
						//      second round should has no more SW than the first round.
						sw_hists_coll[index].insert(sw_hists_coll[index].end(), tmp_h.begin(), tmp_h.end());
					}
					//handle rescaling of sw_hist;
					for( auto cur_hist : tmp_h){
						if(debug_verbose) std::cout<<"\rScale "<<cur_tag<<temp_varName<<" "<<tmp_pscale<<" ";
						
						if(rescale){ 

							cur_hist->Scale(tmp_pscale);
						}
						hist2tscalemap.insert( std::make_pair(cur_hist, tmp_tscale));//scale with throws
					}
					if(debug_verbose) std::cout<<std::endl;
				}
			}
		}//finish sorting;

		int nby = 1.5*cv_hist->GetBinContent(cv_hist->GetMaximumBin());//set histogram height
		for(size_t index = 0; index < hists_names.size(); ++index){//loop over sets of independent histograms

			TString temp_sw_name = cur_tag +"_"+ hists_names[index];
			TH2D* all_hist = new TH2D(temp_sw_name, temp_sw_name, nb, &(output_binning).front(), nby, 0, nby);//to store many throws SW
			TH2D* covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix", temp_sw_name+"_CovarainceMatrix",nb,&(output_binning).front(), nb,&(output_binning).front());//to store covariance matrix, equal width;

			int hist_counter = 1;
			for(auto cur_hist : sw_hists_coll[index]){//ok, take out SW hist and add them to all_hists, covmatrices;

				if(verbose) std::cout<<std::flush<<"\r Processing hists "+temp_sw_name<<hist_counter++;
				
				if(force_rebin) cur_hist->Rebin(2);

				bool trad_bin = false;//for R<5 selection.
				if(do_smooth&& trad_bin){//smooth first
					if(verbose) std::cout<<"\r Traditionally Smoothing a sw histogram";
					*cur_hist = SmoothSW(cur_hist, smoothRef_hist, trad_bin);
				}

				if(do_rebin) *cur_hist = gadget_vrebin(cur_hist, output_binning);

				if(do_smooth && !trad_bin){//smooth first
					if(verbose) std::cout<<"\r Smoothing a sw histogram";
					*cur_hist = SmoothSW(cur_hist, cv_hist, trad_bin);
				}
				for(int jndex = 1; jndex < nb+1;  ++ jndex){
					all_hist->Fill(cur_hist->GetBinLowEdge(jndex) , cur_hist->GetBinContent(jndex));
					if(debug_verbose){
						std::cout<<"\r Fill "<<cur_hist->GetBinLowEdge(jndex)<<" with "<<cur_hist->GetBinContent(jndex)<<std::endl;
					}
				}
				if(debug_verbose)std::cout<<std::endl;

				//std::cout<<cur_hist->GetNbinsX()<<std::endl;
				TH2D* cov_temp = MakeCov(cur_hist, cv_hist);
				double extra_factor = 1;
				//Special! By doing so, we can blow upthe KpProd;
				if(temp_sw_name.Contains("KpProd")){ 
					if(verbose) std::cout<<" with an extra factor 1.5*1.5";
					extra_factor = 1.5*1.5;
				}

				covmatrices->Add(cov_temp, hist2tscalemap[cur_hist]*extra_factor);//throw rescaling;
			}
			if(verbose) std::cout<<std::endl;
			//finish handlig sw_hist;
			TH2D* fracCov = MakeFracCov( covmatrices, cv_hist);


			//Now we can draw all_hist, fracCov, Cov;
			gStyle->SetPalette(kRainBow);//this change colz colors;
			gSystem->RedirectOutput("/dev/null");
			//all_hist 
			drawCanvas->Clear();
			drawCanvas->cd();
			all_hist->SetStats(false);
			all_hist->GetXaxis()->SetTitle(var.unit.c_str());
			all_hist->GetYaxis()->SetTitle("Event Rate");
			all_hist->Draw("COLZ");
			cv_hist->Draw("L same");	
			cv_hist->SetLineColor(6);
			cv_hist->SetLineWidth(2);
			drawCanvas->SaveAs( drawn_dir + hist_prefix + "_"+temp_sw_name+".pdf" ,"pdf");

			//Covariance matrix;
			drawCanvas->Clear();
			covmatrices->SetStats(false);
			covmatrices->GetXaxis()->SetTitle(var.unit.c_str());
			covmatrices->GetYaxis()->SetTitle(var.unit.c_str());
			covmatrices->Draw("COLZ");
			covmatrices->SetTitle(temp_sw_name + " Covaraince Matrix");
			drawCanvas->SaveAs( drawn_dir + cov_prefix + "_"+temp_sw_name+".pdf" ,"pdf");

			cov_root->cd();
			covmatrices->Write();

			//Fractional Covariance matrix;
			drawCanvas->Clear();
			fracCov->SetStats(false);
			fracCov->GetXaxis()->SetTitle(var.unit.c_str());
			fracCov->GetYaxis()->SetTitle(var.unit.c_str());
			fracCov->Draw("COLZ");
			fracCov->SetTitle(temp_sw_name + " Fractional Covaraince Matrix");
			drawCanvas->SaveAs( drawn_dir + cov_prefix + "_"+temp_sw_name+"_frac.pdf" ,"pdf");
			gSystem->RedirectOutput(0,0);//re-activate the pirnt out.
			fracCov->Write();


			//Save CV,FracCov in TF & TMatrix forms to finalcov_root;
			finalcov_root->cd();
			fracCov->Write(temp_sw_name + "_FracMatrix_drawn");
			cv_hist->Write(temp_sw_name + "_CV_drawn");
			//TMatrix
			TMatrixD covMatrix(nb,nb); 
			TArrayD nums(nb*nb);
			TMatrixD cvhistM(nb,1); 
			TArrayD numsCV(nb);

			if(debug_verbose) std::cout<<"Diagonal element of the fractional matrix: ";
			for(int index = 0; index < nb; ++index){
				numsCV[index] = cv_hist->GetBinContent(index+1);
				for(int jndex = 0; jndex < nb; ++jndex){
					nums[jndex*nb+index] = fracCov->GetBinContent(index+1,jndex+1);
					if(debug_verbose) std::cout<<nums[jndex*nb+index]<<", ";
				}
			}
			if(debug_verbose) std::cout<<std::endl;
			covMatrix.SetMatrixArray(nums.GetArray());
			cvhistM.SetMatrixArray(numsCV.GetArray());
			covMatrix.Write(temp_sw_name + "_FracMatrix");
			cvhistM.Write(temp_sw_name + "_CV");

		}//next set of histograms;
	}//next tag
	if(verbose) std::cout<<"Saving matrices with hash "<<hist_prefix<<std::endl;
	finalcov_root->Close();
	cov_root->Close();
}

/*
 * Make one covariance matrix
 */
TH2D* MakeCov(TH1D* hist, TH1D* cv){

	int nb=hist->GetNbinsX();

	//	int nl=hist->GetBinLowEdge(1);
	//	int nh=hist->GetBinLowEdge(nb+1);
	//in case of customized binnings; use vector..
	std::vector<double> bins;
	for( int index = 1; index < nb+2; index++){

		bins.push_back(hist->GetBinLowEdge(index));//1st bin to n+1th bin
	}

//	int counter = 0;
//	delete gROOT->FindObject("tmp");
	TH2D* covmatrix = new TH2D(("tmp"+to_string_prec(global_index++,0)).c_str() ,"tmp",nb, &(bins.front()),nb, &(bins.front()) );//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
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
TH2D* MakeFracCov(TH2D* inputcov, TH1D* oldcv){
	
	int nb=oldcv->GetNbinsX();
	int nl=oldcv->GetBinLowEdge(1);
	int nh=oldcv->GetBinLowEdge(nb+1);
	

	TH2D* fractional_cov = (TH2D*) inputcov->Clone();//get the clone copy
//	fractional_cov->Reset("ICESM");//clear contents;

	for(int jndex = 1; jndex < nb+1; ++jndex){
		for(int kndex = 1; kndex < nb+1; ++kndex){
			double cvj = oldcv->GetBinContent(jndex);
			double cvk = oldcv->GetBinContent(kndex);
			double covjk = inputcov->GetBinContent(jndex,kndex);
			double fjk = covjk/(cvj*cvk);//not divide by nb;

			fractional_cov->SetBinContent(jndex,kndex,fjk);
		}
	}

	return fractional_cov;
}


/*
 * Smooth the matrix
 */
//http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
//int global_index = 0;
TH1D SmoothSW(TH1D* sw, TH1D* fix_cv, bool special_option){//fit each bin with polynomials

	using namespace std;
	//use the special option for the old (but effective) smoothing strategy
	
	bool message = false;
	bool print_ratio = false;
	bool ec_smoothing = false;

	std::vector<double> sm_binning ={200, 250, 300, 375, 475, 550, 600, 675, 750, 800, 950, 1100, 1150, 1250, 1300, 1500, 1700, 1900, 3000};

	TString cur_tag = to_string(sw->Integral())+"Unisims";
	TH1D* cv =(TH1D*) fix_cv->Clone();

	if(special_option){
		if(true){
			double hb1 = cv->GetBinLowEdge(1);
			double hb2 = cv->GetBinLowEdge(cv->GetNbinsX()+1);
			if(sm_binning.front() < hb1 ||sm_binning.back() > hb2){
				std::cout<<"\nNo traditional smoothing for underflow bins from ("<<hb1<<","<<hb2;
				std::cout<<") to ("<<sm_binning.front() <<","<<sm_binning.back() <<")."<<std::endl;
				std::cout<<"No smoothing."<<std::endl;
			return *sw;
			}
		}

		*cv = gadget_vrebin(cv, sm_binning);
		*sw = gadget_vrebin(sw, sm_binning);
		if(false){
			for(int index = 1;index < cv->GetNbinsX()+2; ++index){
				std::cout<<" "<<cv->GetBinLowEdge(index);
			}
			std::cout<<std::endl;
		}

//		cv->Rebin(sm_binning.size()-1, cur_tag+"smCV",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
//		cv = (TH1D*) gDirectory->Get(cur_tag+"smCV");
//
//		sw->Rebin(sm_binning.size()-1, cur_tag+"smSW",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
//		sw = (TH1D*) gDirectory->Get(cur_tag+"smSW");
	}

	Int_t degree = 2;
	Int_t first_bin = 1;
//	Int_t last_bin = sw->GetNbinsX()-1;//1;
//	while(sw->GetBinLowEdge(last_bin)<1275){//dont smooth bins with EnuE>1899
//		last_bin++;
//	}
	Int_t Nbins=sw->GetNbinsX()-1;//this controls how many bins to use;


	if(ec_smoothing){ //use en-chuan's 

		gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
		TCanvas *rCanvas = new TCanvas("rCanvas", "", 600, 400);//Cavans to draw;
		TF1* smoothFitter = new TF1("smoothFitter", "pol4", cv->GetBinLowEdge(1), cv->GetBinLowEdge(Nbins+1));
		TH1D* copysw = (TH1D*) sw->Clone();
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


		TH1D* copy_rsw = (TH1D*) sw->Clone();
		copy_rsw->Divide(cv);
		copy_rsw->SetLineColor(kRed);
		copy_rsw->Draw("same");

//		if(print_ratio){  
//			gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
//			rCanvas->SaveAs(("/scratch/condor-tmp/klin/hellstroms_hive/hive/build/src/fullmcsystematics/drawn/ratio_OM_"+to_string_prec(global_index,0)+".pdf").c_str(),"pdf");
//			gSystem->RedirectOutput(0,0);//this let ROOT warns
//			rCanvas->Delete();
//			global_index++;
//		}

	}else{//use Zarko's
		//solve M a = y
		
		//proceed to smoothing

		cout << "Smoothing bins "<< first_bin << "-" << Nbins << " (N=" << Nbins<< ", "<<cv->GetBinLowEdge(Nbins+1)<<") with pol "<<degree;
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
		if(message) std::cout<<std::endl;
		for (Int_t bin=1;bin<sw->GetNbinsX()+1;bin++) {//Here modify the sw
			if(cv->GetBinContent(bin)>0){
			chi2+=pow(sw->GetBinContent(bin)-cv->GetBinContent(bin), 2)/cv->GetBinContent(bin);
			if(message)	std::cout<<" bin "<<bin<<" sw "<<sw->GetBinContent(bin)<<" cv "<<cv->GetBinContent(bin)<<std::endl;
			}
		}
		double AIC = chi2+ 2*degree+ (2*degree)*(degree+1)/(Nbins-degree-1);
		cout<<" Chi2: "<<chi2<<" AIC = "<<AIC<<" with degree "<<degree<<std::endl;
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
	bool do_cov = false; //true - generate covariance matrix everytime, no matter what; false -this boolean is not functioning
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
		check_cov_root->Close();
	}else{//no such root;
		if(verbose){ 
			std::cout<<"Did not found fractional covaraince matrices root: "<<check_covroot_name<<std::endl;
			std::cout<<"Will generate matices in "<<(this->top_dir+this->final_prefix + ".root")<<std::endl;
		}
	}
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

		if(do_hist || hist_root == NULL){
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

			if(do_hist){ 
				cur_sys->Make1dhist(var1, hist_root);
				cur_sys->fullyloaded = true;
			}

			if( !(cur_sys->Load1dhist(hist_root) )) all_good = false ;//mark!;
		}

		if(!all_good){//now need to Make any histograms if false
			hist_root->Close();
			hist_root = (TFile*) TFile::Open(histfilename,"UPDATE");//reopen the file, but in UPDATE mode;
			for(auto cur_sys : syss){
				if(verbose) std::cout<<"Updating "<<histfilename<<std::endl;

				std::cout<<" Make up histogram "<<cur_sys->tag<<std::endl;
				if(!cur_sys->fullyloaded){ 
				cur_sys->Make1dhist(var1, hist_root);
					cur_sys->fullyloaded = true;
				}
				cur_sys->Load1dhist(hist_root);
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
		hist_root->Close();//if close, reference to hists are lost.
	}


	std::cout<<"Finish making covariance matrix, see \n>>>>  "<<check_covroot_name<<std::endl;
	std::cout<<"\n ----------------------------- \n"<<std::endl;
}
