#include <bdt_systematics.h>

/*
 * Some gadgets to keep things organized;
 */
TString gadget_labelroot(TString plot_type, bdt_variable var, int stage, unsigned long hash){

//	TString binning_text = "_"+to_string_prec(var.n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);
	
	TString binning_text = to_string_prec((var.plot_max-var.plot_min)/var.int_n_bins,0);//to_string_prec(var.int_n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);

	//<type>_<var>_<bins, min, max>
	TString file_name = plot_type+"_"+var.safe_name + "_bgap"+ binning_text + "_s" + to_string_prec(stage,0) + "_h"+to_string_prec(hash,0);

	return file_name;
	//for variable binng that can be derived from this binnings, use these root files;
}


TString gadget_addindex(TString varname, int nth){

	TString nam = varname;
//	std::cout<<"before: "<<nam<<std::endl;	
	TString lab = "[" + to_string_prec(nth,0)+"]";
	
	nam.ReplaceAll("[0]", lab);
//	std::cout<<nam<<std::endl;	
	return nam;
}

unsigned long  gadget_jenkins_hash(std::string key) {
    size_t length = key.size();
    size_t i = 0;
    unsigned long hash = 0;
    while (i != length) {
        hash += key[i++];
        hash += hash << 10;
        hash ^= hash >> 6;
    }
    hash += hash << 3;
    hash ^= hash >> 11;
    hash += hash << 15;
    return hash;
}

/*
 * Constructor for struct bdt_sys, 
 *
 */
bdt_sys::bdt_sys(int index, MVALoader XMLconfig, int bdtfile_index, bdt_flow inflow)
	:
	bdt_file(bdtfile_index, XMLconfig, inflow){
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

		
		//members that need to be derived slightly;
		num_vars = vars.size();
		hist.resize(num_vars);
		twodhist.resize(num_vars);
	
	};

/*
 * Main code, initialize systematics;
 * Check if the histograms exist, if yes, then can save some works.
 */

void InitSys(std::vector<bdt_variable> vars, std::vector<std::string> precuts, std::vector<bdt_sys*> syss, double plot_pot, std::vector<double> bdt_cuts, TString dir_root, TString dir_drawn){
	
	bool do_cov = true; //true - generate covariance matrix everytime, no matter what; false -this boolean is not functioning
	bool do_hist = true;//true - generate hist, no matter what; false -this boolean is not functioning

	bool message = true;
	
	
	//Configure stage, cuts 

	int stage = syss[0]->stage;

	//make a string of precuts, CHECK, it will be more complicated when using BDT cuts;
	std::string precut_for_hashing("("+precuts[0]+")");
	for(int index = 1; index < precuts.size(); ++index){
		precut_for_hashing += "&&("+precuts[index] + ")";
	}


	unsigned long hashlabel = gadget_jenkins_hash(precut_for_hashing);


	//the histogram goes two ways: 1dhist or 2dhist; 1dhist for now;
	bool its_1dhist = (vars.size()<2)? true: false;

	bdt_variable* var1 = &vars[0];
	bdt_variable* var2;
	if(!its_1dhist) var2 = &vars[1];

	
	//Before we start, check if 1) covaraicne matrix exist,  2) histograms exist;
	
	bool skip_sysEvaluation = false;//true - do nothing;

	std::cout<<"\n ---- Initialize systematics for "<<var1->unit;
	std::cout<<"\tSystematic files tag: s"<<stage<<"_h"<<hashlabel<<std::endl;

	//check 1) for the covariance matrix;
	TString check_covroot_name = var1->covar_file.c_str();//dir_root+"/../"+gadget_labelroot("finalcov", *var1, stage, hashlabel)+".root";
	TFile *check_cov_root = (TFile*) TFile::Open(check_covroot_name,"READ"); //one var one hist root;

	if(!do_cov && check_cov_root != NULL ){ //oh we had the file! so .. 

		TH1F* check_temp_th1 =  (TH1F*) check_cov_root->Get(var1->covar_name.c_str());
		if(check_temp_th1 != NULL){
			std::cout<<" Being lazy, no need to make cov matrix, because we had it already!"<<std::endl;
			skip_sysEvaluation = true;
		}
		check_cov_root->Close();
	}
	delete check_cov_root;


//	TString readfile_method = "RECREATE";//create file by default, if exist, update the current file;
	TString histfilename = dir_root+"/"+gadget_labelroot("hist", *var1, stage, hashlabel)+".root";

	if(do_hist || !skip_sysEvaluation){//no root files for the covariance matix, so going to generate it.
		std::cout<<"\nLooking into histograms generation"<<std::endl;
		TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;

		size_t index = 0;//index for syss, determine where we left over last time the code break;

		if(hist_root!=NULL && !do_hist ){ //quick check of existing files; 
			std::cout<<histfilename<<" exist! Now check what is inside."<<std::endl;

			bool found_corrupted_hist = false;
			while( index < syss.size()){
				//CHECK, might need to update particular directory;

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
						if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
//							temps->start_with_weight = jndex;

//							readfile_method = "UPDATE";
							hist_root->Close();
							std::cout<<"\n>> Find a broken directory, when looking at "<<syss[index]->tag<<". Update it now."<<std::endl;

							goto now_make_hist;//do partial generation, index determines which systematic file to start;
						} else{

							(temps->hist[jndex]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
							if(message) std::cout<<"\r\tChecking dir: "<<Tdir_name<<" histogram: "<<hist_nam;
						}

					}//next throws
					if(message) std::cout<<std::endl;
				}//next systematic weights;
				index++;
			}//next systematics
		} else{

			std::cout<<histfilename<<" does not exist. Create it now."<<std::endl;

now_make_hist://start ith syss (index), jth weight(start_with_weight), from the given index; 
//			hist_root = (TFile*) TFile::Open(histfilename,readfile_method); //one var one hist root;
			
//			if(readfile_method == "UPDATE"){ 
//				hist_root->Recover();//when use update, its better to recover the file first.
//				hist_root->mkdir("TRASH");
//				hist_root->Delete("TRASH;*");//make a directory and delete it, then no recovery is needed.
//				hist_root->Write();
//				hist_root->Close();
//				hist_root = (TFile*) TFile::Open(histfilename,readfile_method); //one var one hist root;
//			}
//			std::cout<<"Working on "<<histfilename<<std::endl;

			while( index < syss.size()){//Histogram making is here!! 
				Make1dhist(histfilename, var1, syss[index], plot_pot, bdt_cuts);
				index++;
			}//next systematic files
		}

		std::cout<<"Finish making systematic histograms, and save them to\n>>>> "<<histfilename<<std::endl;
		//
		//Part II: Overlay histograms & draw covaraicne matrix;
		//
		hist2cov( *var1, syss, dir_root, dir_drawn, plot_pot, hashlabel);
		hist_root->Close();
	}

	std::cout<<"Finish making covariance matrix, see \n>>>>  "<<check_covroot_name<<std::endl;
	std::cout<<"\n ----------------------------- \n"<<std::endl;
}


/*
 * Prepare 1d histograms here;
 * No POT weighting, nor variable binning (in this case, use equal binning then rebin later;
 */
void Make1dhist(TString histfilename, bdt_variable* var, bdt_sys* temps, double plot_pot, std::vector<double> bdt_cuts){

	bool message = true;
	
	std::cout<<"\nMaking systematic histograms for "<<temps->tag<<" with "<<temps->num_vars<<" variables."<<std::endl;	

	//configure what type of histograms to make;
	int nbins =  var->int_n_bins;
	int lowbinedge = var->plot_min;
	int highbinedge = var->plot_max;

	bool its_multithrows = temps->its_multithrows;//multi weights or optical model
	const int num_throws = temps->throws;
	const int num_swvars = temps->num_vars;

	//STEP 1 Load files, prepare input TTree
	gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
	TFile *infile = (TFile*) TFile::Open(temps->dir+"/"+temps->filename,"READ"); 
	gSystem->RedirectOutput(0,0);//this let ROOT warns

	if(infile ==NULL){
		std::cout<<temps->dir+"/"+temps->filename<<" does not exist, please check! Program aborts."<<std::endl;
		exit(EXIT_FAILURE);
	}

	TTree* temptree = (TTree*) infile->Get(temps->treename);//temptree is from the systematic weights file;

	//STEP 2 prepare TH1F
//	for(int jndex = temps->start_with_weight; jndex < temps->num_vars; ++jndex){//loop through weights}
	for(int jndex = 0 ; jndex < temps->num_vars; ++jndex){//loop through weights
		TFile* hist_root = (TFile*) TFile::Open(histfilename,"UPDATE"); //Frequent Open and Close, hope this help prevent corrupted links;
		//subdirectory in root for each systematic weighs; start_with_weight labels which weight to start, i.e. jndex;
		//STEP 2.1
		//second check of directory, because there might be some leftover from previous run;
		

		TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
		TString Tdir_name = temps->tag+"_"+weights_nam;

		TDirectory* Tdir = hist_root->GetDirectory(Tdir_name);

		if(Tdir !=NULL){ 
			std::cout<<"\tDirectory "<<Tdir_name<<" exists; is it a good leftover?"<<std::endl; 

			TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(num_throws) :temps->vars_name[jndex];

			TH1F* check_temp_th1 =  (TH1F*) hist_root->Get(Tdir_name+"/"+hist_name);

			if(check_temp_th1 !=NULL){//good we have the file, load it up;
				std::cout<<"\tYes, proceed to next directory."<<std::endl;
				hist_root->Close();
				continue;
			}else{//no we don't have the histogram; then remove the direcory
				std::cout<<"\tNo, reproduce histograms in current directory "<<Tdir_name<<std::endl;
				Tdir->Delete(Tdir_name + ";*");
			}

		} else{//No such directory, create it;
			std::cout<<"\tCreate directory "<<Tdir_name<<std::endl;
			Tdir = hist_root->mkdir(Tdir_name);
		}
		Tdir->cd();

		//STEP 2.2 prepare TH1F
		Long64_t nentries = temptree->GetEntries();
		if(message) std::cout<<" We have entries: "<<nentries<<"   ";

		//STEP 2.2.1 prepare empty histgrams to be filled, each histogram with corresponding weight[kndex] (Note, not POT reweighting).
		std::vector<TH1F*> hist1d(num_throws);//hist
		std::vector<TTreeFormula*> wgtforms(num_throws);//variable
		std::vector<TTreeFormula*> varformula(num_throws);//weight, but not with POT reweighting;

		int count_hist = 0;
		for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
			//histogram of a throw
			TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];

			hist1d[kndex] = new TH1F(hist_name,hist_name, nbins, lowbinedge, highbinedge);//Will be rebinned later for variable binning;
			count_hist++;


			//variables of a throw
			varformula[kndex] = new TTreeFormula(var->unit.c_str(), gadget_addindex(var->mininame, kndex), temptree);//get the x-axis
			varformula[kndex]->GetNdata();

			//weights of a throw
			TString temp_label = (kndex>0)? "["+std::to_string(kndex)+"]":"";
			TString temp_wgt = temps->vars[jndex] + temp_label;
			TString temp_wgtname = temps->vars_name[jndex]+std::to_string(kndex);

			wgtforms[kndex] = new TTreeFormula(temp_wgtname, temp_wgt+"*"+temps->getStageCutsPlus(temps->stage,bdt_cuts, kndex), temptree);//get weight with precuts;
			wgtforms[kndex]->GetNdata();
//			std::cout<<temp_wgt+"*"+temps->getStageCutsPlus(temps->stage,bdt_cuts, kndex)<<std::endl;
			//					std::cout<< wgtforms[kndex]->EvalInstance()<<std::endl;
		}//next throw

		if(message) std::cout<<"Creating "<<count_hist<<" histograms "<<std::endl;

		//STEP 2.2.2 fill histograms!
		for(Long64_t kentry = 0; kentry< nentries; ++kentry){

			temptree->GetEntry(kentry);

			for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;

//				TTreeFormula varformula(var->unit.c_str(), gadget_addindex(var->mininame, kndex), temptree);//get the x-axis
//				varformula.GetNdata();
				double temp_xvalue = varformula[kndex]->EvalInstance();
				double temp_weight = wgtforms[kndex]->EvalInstance();//*plot_pot/temps->pot; do this when making cov matrix
//				std::cout<<"CHECK! "<<kndex<<" throws "<<wgtforms[kndex]->EvalInstance()<<" "<<plot_pot<<" "<<temps->pot<<std::endl;

				hist1d[kndex]->Fill(temp_xvalue,temp_weight);//jndex - nth sw; kndex - nth throw

				if(temp_xvalue*temp_weight>0){
				std::cout<<"\r\tLooping entries: "<<std::setw(5)<<kentry+1<<"/"<<std::setw(5)<<nentries;
//				std::cout<<" "<<std::setw(10)<<gadget_addindex(var->mininame, kndex);
				std::cout<<"="<<std::setw(10)<<temp_xvalue;
//				std::cout<<" "<<std::setw(10)<<temps->vars[jndex]<<"["<<kndex<<"]";
				std::cout<<"="<<std::setw(10)<<temp_weight;
//				std::cout<<temps->getStageCutsPlus(temps->stage,bdt_cuts, kndex);
				std::cout<<std::setw(4)<<kndex+1<<" throw sw name: "<<temps->vars_name[jndex];

				}
			}//next throw;
		}

		std::cout<<"\n Finish filling "<<num_swvars<<" systematic weights."<<std::endl;

		//write hisograms to the diretory
		for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
			//configure th1
			TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];
			hist1d[kndex]->SetStats(0);
			hist1d[kndex]->GetXaxis()->SetTitle(var->unit.c_str());
			hist1d[kndex]->GetYaxis()->SetTitle("Events");

			//save the th1 to the root file, class bdt_sys;
			Tdir->cd();
			hist1d[kndex]->Write(hist_name);
			//This is important for making covariance matrix using hist2cov();
//			(temps->hist[jndex]).push_back( (TH1F*) hist1d[kndex]->Clone());// each bdt_sys now contains all hist;
//			std::cout<<"Push back hist "<<hist_name<<std::endl;

			hist1d[kndex]->Delete();
		}

		std::cout<<"Save "<<num_throws;
		std::cout<<" histograms to the dir: "<<Tdir_name<<std::endl;

		hist_root->Close();
	}//next systematic weight 
	infile->Close();
			

	std::cout<<"\t  Now, load histograms for tag "<<temps->tag<<std::endl;
	TFile* temp_hist_root = (TFile*) TFile::Open(histfilename,"UPDATE"); //Frequent Open and Close, hope this help prevent corrupted links;
	for(int jndex = 0; jndex < temps->num_vars; ++jndex){//loop through weights
		TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
		TString Tdir_name = temps->tag+"_"+weights_nam;

		for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
			TString hist_nam_label = (its_multithrows)? std::to_string(kndex+1) : "";
			TString hist_nam = weights_nam + hist_nam_label;
			TH1F* temp_th1 =  (TH1F*) temp_hist_root->Get(Tdir_name+"/"+hist_nam);
			if(temp_th1 == NULL){//oops, some histograms are missing, label it and redo the previous one and the rest;
				std::cout<<"Broken hist?? No way! Come check "<<__FILE__<<__LINE__<<std::endl;
				exit(EXIT_FAILURE);
			} else{

				(temps->hist[jndex]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
				(temps->hist[jndex]).back()->SetDirectory(0);//prevent the histogram from being deleted when file is closed. https://root-forum.cern.ch/t/use-of-th1-clone/11234/3
				if(message) std::cout<<"\r\tLoading histogram "<<hist_nam<<" from dir: "<<Tdir_name;
			}
		}//next throws
		if(message) std::cout<<std::endl;
	}//next systematic weights;
	temp_hist_root->Close();

}

/*
 * Build covariance matrices based on histograms;
 * syss contains all the 1d histograms;
 */

void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss, TString dir_root, TString dir_drawn, double plot_pot, unsigned long hashdraw){
		
	bool message = true;
	bool smooth_matrix = true;
	bool checkbins = true;

	bool rescale = true;

	bool force_rebin = false;
	bool adjOMErr = true;

	if(!rescale) std::cout<<"\n Warning: raw histograms are used, no normalization"<<std::endl;
	//STEP 0 Configuration
	int stage = syss[0]->stage;

	int nb = var.n_bins;//n_bins is the target binning;
	int nl = var.plot_min;
	int nh = var.plot_max;
	bool do_rebin = var.is_custombin;//this will be updated, if the sys histogram binning can be used.

	std::vector<double> uw_statErr(nb,0);//unwated statistical error;

	//Set binnings for output histograms/covariance marix;
	std::vector<double> cur_binning;
	if(do_rebin){
		cur_binning = var.edges;
		cur_binning.erase(cur_binning.begin());//remove the first element, then we get the binning;
	} else{
		double bin_left = nl; 
		double bingap = (nh-nl)/nb;

		for(int jndex = 0; jndex < nb+1; jndex++){
			cur_binning.push_back(bin_left);
			bin_left+=bingap;
		}
	}

	if(message) std::cout<<"\n\nMaking covariance matices."<<std::endl;
	
	//STEP 1 group systematic hists for evaluating covariance matrices;
	std::vector<TString> tags;//tags = systags[?]+"_"+bdtftags[?];
	std::vector<TString> systags;//systematic type; not used yet.
	std::vector<TString> bdtftags;//data, dirt, pi0misd, etc. Not used yet.
	std::multimap< TString, bdt_sys* > sys2SWmap;
	std::map< TString, bdt_sys* > sys2CVmap;
	std::map< TString, bdt_sys* > sys2OMCVmap;

	for(size_t index = 0; index < syss.size(); ++index){//first classify what tags are available, then group them next;
		bdt_sys* temps = syss[index];

		TString temptag = temps->tag;

		tags	.push_back(temptag);//for individual systematic x bdtfiles; will repeat for one of each: sw, CV, and maybe also OMCV;
		systags	.push_back(temps->systag);//for multi-cov; will contain duplicated
		bdtftags.push_back(temps->bdt_file::tag);//for multi-cov; will contain duplicated

		if(message) std::cout<< "Looking at "<<temptag;
		if(temps->its_CV){//save CV to map;
			sys2CVmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<", CV!";
		}else if(temps->its_OM){//Optical model has another CV
			sys2OMCVmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<", OMCV!";
		} else{
			sys2SWmap.insert(std::make_pair (temptag, temps));
			if(message)std::cout<<", Systematic weights!";
		}//save SW to map
		if(message)std::cout<<std::endl;
	}

	//STEP 1.1 clean up tags;
	sort(tags.begin(), tags.end() );
	tags.erase( unique( tags.begin(), tags.end()), tags.end() );


	//STEP 2 pair up SW/CV to make covariance matrix;
	
	//STEP 2.1 prepare root file to store  covariance matrix and configure TH2D*;

	TFile *cov_root = (TFile*) TFile::Open(dir_root+"/"+gadget_labelroot("cov",var, stage, hashdraw)+".root","RECREATE"); //Collection of resulting covariance matrices;

	TString final_covroot_name = dir_root+"/../"+gadget_labelroot("finalcov", var, stage, hashdraw)+".root";
	TFile *finalcov_root = (TFile*) TFile::Open(final_covroot_name,"RECREATE");//one final covariance matrix output for plotting;

	
	//final covariance matrix that sums up all systematics;
	TH2D* finalcov =  new TH2D((var.safe_name).c_str(), (var.safe_name).c_str() ,nb,&(cur_binning).front(),nb,&(cur_binning).front());//binning for xnbins,xmin,xmax,ybins,ymin.ymax;

	TCanvas *histCanvas = new TCanvas("histCanvas", "", 600, 400);//Cavans to draw;
	TCanvas *covCanvas = new TCanvas("covCanvas","",600,400);

	//STEP 3 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	for(size_t index = 0; index < tags.size(); ++index){//each tag means one group of covariance matices;
		TString cur_tag = tags[index];
		std::cout<<"\n-->Working on tag "<<cur_tag<<std::endl;

		//STEP 3.1 Find the bdt_sys syss that gives the CV; then do SW, OMCV on the flight;
		bdt_sys* tempsCV = sys2CVmap.find(cur_tag)->second;//bdt_sys that gives CV - tempsCV

		bool omhistogram = (tempsCV->systag).Contains("OpticalModel");
		bool do_smooth = ( (var.name).compare(3,5,"EnuQE")==0 )&& smooth_matrix && ((tempsCV->systag).Contains("Unisim"));//||omhistogram);

//		std::cout<<"\nWorking on "<<tempsCV->systag<<std::endl;
		if(do_smooth) std::cout<<tempsCV->systag<<" might need smoothing as long as initial binning >1, current bins "<<var.int_n_bins<<std::endl;

		//Got the CV! 
		TH1F* cv_hist = (TH1F*) (tempsCV->hist[0][0])->Clone();
		if(force_rebin) cv_hist->Rebin(2);

		if(checkbins){//check syss's hisogram bin edges; see if it matches var's binning;
			//case 1 - non-rebinable
			//case 2 - rebinable, but need to merge "dynamically"
			//case 3 - rebinable, perfect;
			int binchecker = 0;//monitor if there are difference between var binning and histogram binning;
			int initial_bindex = 1;
			for(double binedge : cur_binning){
				//exp: (20,50,70) , the target binning;
				//test1: (20,30,51), cv_hist binning;
				//test2: (20,25,30,50,70)
				//test3: (20,50,70)
				// As long as we can find all bins from the hist, thne it is good.
				for(int bindex = initial_bindex; bindex < cv_hist->GetNbinsX()+1; ++bindex){
//					std::cout<<cv_hist->GetBinLowEdge(bindex)<<" vs "<<binedge;
					if(abs(cv_hist->GetBinLowEdge(bindex) - binedge) < 10e-20){//edge matches;
						binchecker++;
						bindex = cv_hist->GetNbinsX()+1;
//								std::cout<<"Good match"<<std::endl;
						continue;
					}
//							std::cout<<std::endl;
				}
				//initial_bindex++;
			}
			checkbins = false;//only need to check once.

			if(binchecker == nb+1){//the left edge does not count;
				std::cout<<"Binning is adjustable: nbins from "<<binchecker<<" to "<<nb+1<<std::endl;
				do_rebin = true;
			} else{
				std::cout<<"Binning not adjustable: nbins from "<<binchecker<<" to "<<nb+1<<std::endl;
				std::cout<<"Please check the input hist_*.root."<<std::endl; 
				exit(EXIT_FAILURE);
			}
		}

		//Update cv_hist, if the binning are not exactly what we want;
//		TH1F* cv_hist = new TH1F(tempsCV->vars_name[0]+"copy",tempsCV->vars_name[0]+"copy",nb,nl,nh);//(tempsCV->hist[0][0])->Clone();
//		int set_1st_bin = 1;
//		tempsCV->hist[0][0]->Rebin(2);
//		for(int bindex = 1; bindex < tempsCV->hist[0][0]->GetNbinsX()+1; ++bindex){
//			if(tempsCV->hist[0][0]->GetBinLowEdge(bindex+1)>nl||tempsCV->hist[0][0]->GetBinLowEdge(bindex-1)<nh){//select some bins in between;
//				cv_hist->SetBinContent(set_1st_bin++,tempsCV->hist[0][0]->GetBinContent(bindex));
//			}
//		}

		TH1F* smoothRef_hist = (TH1F*) cv_hist->Clone();
		if(do_rebin){ 
			cv_hist->Rebin(nb, cur_tag+"CV",&(cur_binning).front());//tags[index]+"CV" is just some name, can be anything;
			cv_hist = (TH1F*) gDirectory->Get(cur_tag+"CV");
		}
//		if(omhistogram){//capture statistical error for OpticalModel
//			for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
//				uw_statErr[lndex] = cv_hist->GetBinError(lndex+1);
//				std::cout<<"Pre-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<cv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
//			}
//		};
		if(omhistogram && adjOMErr){//capture statistical error for OpticalModel
			for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
				uw_statErr[lndex] = cv_hist->GetBinError(lndex+1);
				std::cout<<"Pre-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<cv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
			}
		};

		if(rescale)cv_hist->Scale(plot_pot/tempsCV->pot);//scale hist to data POT

		if(omhistogram && adjOMErr){//capture statistical error for OpticalModel
			for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
				uw_statErr[lndex] = cv_hist->GetBinError(lndex+1);
				std::cout<<"Post-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<cv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
			}
		};


		//STEP 3.2 Load SW and draw histograms, covariance matrix only the flight;
		int nby = 1.5*cv_hist->GetBinContent(cv_hist->GetMaximumBin());//set histogram height
		int total_throws = 0;
		std::vector< bdt_sys*> to_plot_list;

		for( std::_Rb_tree_iterator<std::pair<const TString, bdt_sys*> >  itr = sys2SWmap.begin();  itr!=sys2SWmap.end(); itr++){//group bdt_sys SW under the same tag

			if(itr->first == cur_tag){//identify the related bdt_sys based on bdtfiles;
				std::cout<<"Identify a file for tag "<<cur_tag<<std::endl;
				total_throws += (itr->second)->throws;
				to_plot_list.push_back( itr->second);//get the bdt_file;
				if(to_plot_list.size()>1 && (itr->second)->vars_name.size()>1){
					std::cout<<"Now the code does not work on multi bdt_sys with multi systematic variables.";
					std::cout<<"Need to upgrade if needed, or just modify the xml file for now."<<std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}
		if(to_plot_list.size()> 1) std::cout<<to_plot_list.size()<<" bdt_sys's are added under the tag "<<cur_tag<<std::endl;


		//STEP 3.2.1 prepare SW of throws;
		//				bdt_sys* tempsSW = itr->second;//tempSW may contains  manySW x throws;

		bdt_sys* tempsSW = to_plot_list[0];
		//
		//for Unisims, there will be two elements 
		int numSW = (tempsSW->hist).size();//hist is a 2d vector w dimension numSW x throws;
		for(int jndex = 0; jndex < numSW; ++jndex){//go through different SW under same bdt_sys
			//USE THIS LABEL!
			TString temp_sw_name  = cur_tag+"_"+tempsSW->vars_name[jndex];

			TH2D* all_hist = new TH2D(temp_sw_name, temp_sw_name, nb, &(cur_binning).front(), nby, 0, nby);//to store many throws SW
			TH2D* covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix", temp_sw_name+"_CovarainceMatrix",nb,&(cur_binning).front(), nb,&(cur_binning).front());//to store covariance matrix, equal width;
			TH2D* ori_covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix_original", temp_sw_name+"_CovarainceMatrix_original",nb,&(cur_binning).front(), nb,&(cur_binning).front());//to store covariance matrix, equal width;

//			double temp_covcal = 0;
			for(int kndex = 0; kndex < tempsSW->throws; kndex++){//go through different throws
				TH1F* sw_hist = (TH1F*) (tempsSW->hist[jndex][kndex])->Clone("copy"+tempsSW->vars_name[jndex]);
				if(force_rebin) sw_hist->Rebin(2);
//				TString hist_name = (tempsSW->its_multithrows)? tempsSW->vars_name[jndex]+std::to_string(kndex) :tempsSW->vars_name[jndex];
//				TH1F* sw_hist = new TH1F(hist_name+"copy",hist_name+"copy",nb,nl,nh);//(tempsCV->hist[0][0])->Clone();
//				int set_1stsw_bin = 1;
//				tempsSW->hist[jndex][kndex]->Rebin(2);
//				for(int bindex = 1; bindex < tempsSW->hist[jndex][kndex]->GetNbinsX()+1; ++bindex){
//					if(tempsSW->hist[jndex][kndex]->GetBinLowEdge(bindex+1)>nl||tempsSW->hist[jndex][kndex]->GetBinLowEdge(bindex-1)<nh){//select some bins in between;
//						sw_hist->SetBinContent(set_1stsw_bin++, tempsSW->hist[jndex][kndex]->GetBinContent(bindex));
////						std::cout<<"CHECK "<<tempsSW->hist[jndex][kndex]->GetBinLowEdge(bindex)<<" "<<tempsSW->hist[jndex][kndex]->GetBinContent(bindex)<<std::endl;
//					}
//				}
				///scale, smooth, rebin;

				if(var.int_n_bins>1&& do_smooth ){ 
					std::cout<<"\r Smoothing a sw histogram";
					*sw_hist = SmoothSW(sw_hist, smoothRef_hist, omhistogram);
				}

				if(do_rebin){ 
					TString rebin_label =  cur_tag+"_"+tempsSW->vars_name[jndex]+"sw"+to_string_prec(kndex,0);
					sw_hist->Rebin(nb,rebin_label,&(cur_binning).front());//repeat label will get the previous histogram
					sw_hist = (TH1F*) gDirectory->Get(rebin_label);
				}

				if(rescale)sw_hist->Scale(plot_pot/tempsSW->pot);

//				if(nb>1&& do_smooth ){ 
//					std::cout<<"\r Smoothing a sw histogram";
//					*sw_hist = SmoothSW(sw_hist, cv_hist, false);//omhistogram);
//				}

				//STEP 3.2.2 make histograms & covriance matrix;
//				double temp_sw=0;//save first 14 bin contents;
//				double temp_cv=0;
//				for(int lndex = 1; lndex < nb+1; ++lndex){//fill sw_hist to all_hist bins by bins;
//					all_hist->Fill(sw_hist->GetBinLowEdge(lndex) , sw_hist->GetBinContent(lndex) );
//					if(lndex<14){
//				//	std::cout<<"\rGet bin up to "<<sw_hist->GetBinLowEdge(lndex+1);
//					temp_sw += sw_hist->GetBinContent(lndex);
//					temp_cv += cv_hist->GetBinContent(lndex);
//					}
//				}//next bin


//				temp_covcal+=pow(temp_sw-temp_cv,2)/(double) total_throws;

				//Make cov matrices
				TH2D* cov_temp = MakeCov("cov"+cur_tag+"_"+tempsSW->vars_name[jndex]+std::to_string(kndex), sw_hist, cv_hist);
				//make propagated covariance matrix;

				if(omhistogram && adjOMErr){//now substract the stat error, before scaling the covariance matrix;
//					std::cout<<"\rAdjust cov matrix by substituting the stat error ";
					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
						for(int mndex = 0; mndex < nb; ++mndex){//fill cv_hist to all_hist bins by bins;
							double temp_covvalue = cov_temp->GetBinContent(lndex+1,mndex+1);
							double temp_updatedvalue = temp_covvalue - uw_statErr[lndex]*uw_statErr[mndex] +(plot_pot/tempsCV->pot)*sqrt(cv_hist->GetBinContent(lndex+1)*cv_hist->GetBinContent(mndex+1));
							cov_temp->SetBinContent(lndex+1,mndex+1, temp_updatedvalue);

							if(lndex==mndex){
//							std::cout<<"Adjust Cov statError "<<temp_covvalue<<"->"<<temp_updatedvalue<<",  uwstatErr^2: "<<uw_statErr[lndex]*uw_statErr[mndex];
//							std::cout<<", cv_hist err^2:"<<sqrt(cv_hist->GetBinContent(lndex+1)*cv_hist->GetBinContent(mndex+1))<<std::endl;
							}
						}
					}
				}


				double scale_factor = 1.0;
				if(tempsSW->vars_name[jndex].Contains("KpProd") ) scale_factor = 1.5;//from Zarko's MaktrixMaker xml file, see /e898-data/data60/desktop_backup/tirpitz/scratch/zarko/analysis2012/MatrixMaker/xml/matrix_combined_analysis.xml & line 413 of MatrixMaker.cxx

				cov_temp->Scale(pow(scale_factor,2)/(double) total_throws);

				covmatrices->Add(cov_temp);
//				std::cout<<covmatrices->GetNbinsX()<<" vs. "<<cov_temp->GetNbinsX()<<std::endl;
//				exit(0);

				if(kndex == tempsSW->throws-1 && to_plot_list.size()>1){//wait, one more throw! Reset
					to_plot_list.erase(to_plot_list.begin());
					kndex = -1;//will be 0 next round;	
					tempsSW = to_plot_list[0];
					std::cout<<"\nWait! More throws to append: "<<tempsSW->vars_name[0]<<std::endl;
				}

				//fill ori_covmatrices, no matter  smoothing or not;
				if(omhistogram && false){
					TH1F* originsw_hist = (TH1F*) (tempsSW->hist[jndex][kndex])->Clone("ori"+tempsSW->vars_name[jndex]);
					///scale, smooth, rebin;
					originsw_hist->Scale(plot_pot/tempsSW->pot);
					if(var.is_custombin){ 
						TString rebin_label =  "org"+cur_tag+"_"+tempsSW->vars_name[jndex]+"sw"+to_string_prec(kndex,0);
						originsw_hist->Rebin(nb,rebin_label,&(cur_binning).front());//repeat label will get the previous histogram
						originsw_hist = (TH1F*) gDirectory->Get(rebin_label);
					}

					//STEP 3.2.2 make histograms & covriance matrix;
					for(int lndex = 1; lndex < nb+1; ++lndex){//fill originsw_hist to all_hist bins by bins;
						all_hist->Fill(originsw_hist->GetBinLowEdge(lndex) , originsw_hist->GetBinContent(lndex) );
					}//next bin

					//Make cov matrices
					TH2D* ori_cov_temp = MakeCov("ori_cov"+cur_tag+"_"+tempsSW->vars_name[jndex]+std::to_string(kndex), originsw_hist, cv_hist);

					double scale_factor = 1.0;
					if(tempsSW->vars_name[jndex].Contains("KpProd") ) scale_factor = 1.5;//from Zarko's MaktrixMaker.cxx
					ori_cov_temp->Scale(scale_factor/(double) total_throws);

					ori_covmatrices->Add(ori_cov_temp);
				}

			}//next throws;
			std::cout<<std::endl;

//			if(false){//do one bin analysis here..
//				std::cout<<"\nOne-Bin Analysis "<< temp_covcal<<std::endl;
//			}



			//STEP 3.2.3 all_hist & covmatrices are filled; now draw them;
			gStyle->SetPalette(kRainBow);//this change colz colors;

	gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
			//Histograms;
			histCanvas->Clear();
			histCanvas->cd();
			all_hist->SetStats(false);
			all_hist->GetXaxis()->SetTitle(var.unit.c_str());
			all_hist->GetYaxis()->SetTitle("Event Rate");
			all_hist->Draw("COLZ");
			cv_hist->Draw("L same");	
			cv_hist->SetLineColor(6);
			cv_hist->SetLineWidth(2);
			histCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("hist",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");

			//Covariance matrix;
			covCanvas->Clear();
			covCanvas->cd();
			//check for OpticalModel

			//make fractional covariance matrix, but this does not added to final covariance matrix;
			TH2D* fracCov = MakeFracCov(temp_sw_name, covmatrices, cv_hist, false, ori_covmatrices);//omhistogram, ori_covmatrices);//omhistogram - false, then the last argument is useless.
			fracCov->SetStats(false);
			fracCov->Draw("COLZ");
			fracCov->SetTitle("Fractional "+temp_sw_name + " Covaraince Matrix");
			covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covFrac",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");

	  gSystem->RedirectOutput(0,0);//this let ROOT warns
			if(message&& false){ 
				std::cout<<"The Fractional matrix diagonal element "<< temp_sw_name<<": "<<std::endl;;
				for(int bindex = 1; bindex<nb+1; ++bindex){
					std::cout<<" "<<fracCov->GetBinContent(bindex,bindex)<<",";
				}
				std::cout<<std::endl;
			}


			covCanvas->Clear();
			covCanvas->cd();

			if(sys2OMCVmap.count(tags[index]) > 0){//its a Optical Model CV, propagate the matrix to a fractional covariance matrix; and add it to finalcov.
				bdt_sys* tempsOMCV = (sys2OMCVmap.find(cur_tag))->second;//get the sys for a given tag;
				double omcv_pot = tempsOMCV->pot;
				if(message) std::cout<<"Propagating covariance matrix for "<<cur_tag<<std::endl;
				TH1F* omcv_hist = (TH1F*) (tempsOMCV->hist[0][0])->Clone();
				if(force_rebin) omcv_hist->Rebin(2);

				//rebin, scale
				if(do_rebin){ 
					omcv_hist->Rebin(nb, temp_sw_name+"OMCV",&(cur_binning).front());
					omcv_hist = (TH1F*) gDirectory->Get(temp_sw_name+"OMCV");
				}

				if(!adjOMErr){//capture statistical error for OpticalModel
					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
						uw_statErr[lndex] = omcv_hist->GetBinError(lndex+1);
						std::cout<<"Pre-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<omcv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
					}
				};
				if(rescale) omcv_hist->Scale(plot_pot/omcv_pot);


				if(!adjOMErr){//capture statistical error for OpticalModel
					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
						uw_statErr[lndex] = omcv_hist->GetBinError(lndex+1);
						std::cout<<"Post-scaled Get Err^2:"<<pow(uw_statErr[lndex],2)<<" from "<<omcv_hist->GetBinContent(lndex+1)<<" events"<<std::endl;
					}
				};

				//make propagated covariance matrix;
				TH2D* ProCov = PropagateCov(temp_sw_name, fracCov, omcv_hist);

				if(!adjOMErr){//now substract the stat error, before scaling the covariance matrix;
//					std::cout<<"\rAdjust cov matrix by substituting the stat error ";
					for(int lndex = 0; lndex < nb; ++lndex){//fill cv_hist to all_hist bins by bins;
						for(int mndex = 0; mndex < nb; ++mndex){//fill cv_hist to all_hist bins by bins;
							double temp_covvalue = ProCov->GetBinContent(lndex+1,mndex+1);
							double temp_updatedvalue = temp_covvalue - uw_statErr[lndex]*uw_statErr[mndex] + sqrt(omcv_hist->GetBinContent(lndex+1)*omcv_hist->GetBinContent(mndex+1));
							ProCov->SetBinContent(lndex+1,mndex+1, temp_updatedvalue);

							if(lndex==mndex){
							std::cout<<"Adjust Cov statError "<<temp_covvalue<<"->"<<temp_updatedvalue<<",  uwstatErr^2: "<<uw_statErr[lndex]*uw_statErr[mndex];
							std::cout<<", omcv_hist err^2:"<<sqrt(omcv_hist->GetBinContent(lndex+1)*omcv_hist->GetBinContent(mndex+1))<<std::endl;
							}
						}
					}
				}

				ProCov->SetStats(false);
				ProCov->Draw("COLZ");
				ProCov->SetTitle("Propagated "+temp_sw_name + " Covaraince Matrix");

				gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
				covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covProp",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
				gSystem->RedirectOutput(0,0);//this let ROOT warns

				if(message){ 
					std::cout<< temp_sw_name<<" matrix diagonal element: "<<std::endl;;
					double temp_calcul = 0;
					for(int bindex = 1; bindex<nb+1; ++bindex){
						std::cout<<" "<<ProCov->GetBinContent(bindex,bindex)<<",";
						if(ProCov->GetBinContent(bindex,bindex)>0&&bindex<14) temp_calcul+=ProCov->GetBinContent(bindex,bindex);
					}
					std::cout<<std::endl;
//					std::cout<<"Quadrature sum: "<<temp_calcul<<std::endl;
				}

				finalcov->Add(ProCov);

				cov_root->cd();
				fracCov->Write();
				ProCov->Write();
			} else{//no need extra processing;

				covmatrices->SetStats(false);
				covmatrices->GetXaxis()->SetTitle((var.unit).c_str());
				covmatrices->Draw("COLZ");
				gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
				covCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("cov",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");
				gSystem->RedirectOutput(0,0);//this let ROOT warns

				if(message){ 
					std::cout<<"The matrix diagonal element "<<temp_sw_name<<": "<<std::endl;
					double temp_calcul = 0;
					for(int bindex = 1; bindex<nb+1; ++bindex){
						std::cout<<" "<<covmatrices->GetBinContent(bindex,bindex)<<",";
						if(covmatrices->GetBinContent(bindex,bindex)>0&&bindex<14) temp_calcul+=covmatrices->GetBinContent(bindex,bindex);
					}
					std::cout<<std::endl;
//					std::cout<<"Quadrature sum: "<<temp_calcul<<std::endl;
				}
			
				finalcov->Add(covmatrices);
				cov_root->cd();
				covmatrices->Write();
			}
			//one covmatrix is made;
		}//next SW
	}//next tag;

	//STEP 4 draw out the final covariance matrix;

	cov_root->Close();//basic covariance matrix is recorded, now prepare the one for plotting;
	
	//set axis and name for the final covaraince matrix;
	finalcov->SetStats(false);

	finalcov->SetTitle("Total Covariance Matrix for");
	finalcov->GetXaxis()->SetTitle( (var.unit).c_str());

	TCanvas *c1 = new TCanvas("c1","",600,400);
	finalcov->Draw("COLZ");

	c1->SaveAs(dir_drawn+"/"+gadget_labelroot("Totalcov", var, stage, hashdraw)+".pdf","pdf");

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
	std::cout<<"Final total covariance matrix is saved at >>>>   "<<final_covroot_name<<std::endl;
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
TH1F SmoothSW(TH1F* sw, TH1F* cv, bool ec_smoothing){//fit each bin with polynomials

	using namespace std;
	bool print_ratio = false;
	ec_smoothing = false;

//	std::vector<double> sm_binning ={200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1250, 1300, 1500, 1900, 3000};
	std::vector<double> sm_binning ={200, 250, 300, 375, 475, 550, 600, 675, 750, 800, 950, 1100, 1150, 1250, 1300, 1500, 1700, 1900, 3000};

	TString cur_tag = to_string(sw->Integral())+"Unisims";

//	cv->Rebin(sm_binning.size()-1, cur_tag+"smCV",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
//	cv = (TH1F*) gDirectory->Get(cur_tag+"smCV");
//
//	sw->Rebin(sm_binning.size()-1, cur_tag+"smSW",&(sm_binning).front());//tags[index]+"CV" is just some name, can be anything;
//	sw = (TH1F*) gDirectory->Get(cur_tag+"smSW");

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

