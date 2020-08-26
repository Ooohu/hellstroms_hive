#include <bdt_systematics.h>

/*
 * Some gadgets to keep things organized;
 */
TString gadget_labelroot(TString plot_type, bdt_variable var, int stage, unsigned long hash){

//	TString binning_text = "_"+to_string_prec(var.n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);
	
	TString binning_text = to_string_prec(var.int_n_bins,0)+"_"+to_string_prec(var.plot_min,0)+"_"+to_string_prec(var.plot_max,0);

	//<type>_<var>_<bins, min, max>
	TString file_name = plot_type+"_"+var.safe_name + "_"+ binning_text + "_s" + to_string_prec(stage,0) + "_h"+to_string_prec(hash,0);

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
 *
 */

void InitSys(std::vector<bdt_variable> vars, std::vector<std::string> precuts, std::vector<bdt_sys*> syss, double plot_pot, std::vector<double> bdt_cuts, TString dir_root, TString dir_drawn){
	
	bool do_cov = true; //true - generate covariance matrix everytime, no matter what; false -this boolean is not functioning
	bool do_hist = false;//true - generate hist, no matter what; false -this boolean is not functioning

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

	if(message) std::cout<<"\n ---- Initialize systematics for "<<var1->unit<<std::endl;

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


	TString readfile_method = "RECREATE";//create file by default, if exist, update the current file;
	TString histfilename = dir_root+"/"+gadget_labelroot("hist", *var1, stage, hashlabel)+".root";

	if(do_hist || !skip_sysEvaluation){//no root fiels for the covariance matix, so going to generate it.
		TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;

		size_t index = 0;//index for syss, determine where we left over last time the code break;

		if(hist_root!=NULL ){ //quick check of existing files; 
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
							temps->start_with_weight = jndex;
//							rmdir(temps->tag+"_"+weights_nam[jndex]);

							readfile_method = "UPDATE";

							hist_root->Close();

							std::cout<<"\n>> Find a broken directory, reproduce histograms.starting from "<<syss[index]->tag<<std::endl;

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
			hist_root = (TFile*) TFile::Open(histfilename,readfile_method); //one var one hist root;

			std::cout<<"Working on "<<histfilename<<std::endl;

			while( index < syss.size()){//Histogram making is here!! 
				Make1dhist(hist_root, var1, syss[index], plot_pot, bdt_cuts);
				index++;
			}//next systematic files
		}

		std::cout<<"Finish making systematic histograms, and save them to >>>> "<<histfilename<<std::endl;
		//
		//Part II: Overlay histograms & draw covaraicne matrix;
		//
		hist2cov( *var1, syss, dir_root, dir_drawn, plot_pot, hashlabel);
		hist_root->Close();
	}

	std::cout<<"Finish making covariance matrix, see >>>>  "<<check_covroot_name<<std::endl;
	std::cout<<"\n ----------------------------- \n"<<std::endl;
}


/*
 * Prepare 1d histograms here;
 * No POT weighting, nor variable binning (in this case, use equal binning then rebin later;
 */
void Make1dhist(TFile* hist_root, bdt_variable* var, bdt_sys* temps, double plot_pot, std::vector<double> bdt_cuts){

	bool message = true;
	
	std::cout<<"Making systematic histograms for "<<temps->tag<<" with "<<temps->num_vars<<" variables."<<std::endl;	


	//configure what type of histograms to make;
	bool its_multithrows = temps->its_multithrows;//multi weights or optical model


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
	bool prev_bad_dir = true;//previous directory must be corrupted;
	for(int jndex = temps->start_with_weight; jndex < temps->num_vars; ++jndex){//loop through weights
		//subdirectory in root for each systematic weighs; start_with_weight labels which weight to start, i.e. jndex;
		
		//STEP 2.1
		//second check of directory, because there might be some leftover from previous run;
		
		const int num_throws = temps->throws;
		const int num_swvars = temps->num_vars;

		TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
		TString Tdir_name = temps->tag+"_"+weights_nam;

		TDirectory* Tdir = hist_root->GetDirectory(Tdir_name);

		if(Tdir !=NULL){ 
			std::cout<<Tdir_name<<" exists, and it is properly a good leftover;"<<std::endl; 

			TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(num_throws) :temps->vars_name[jndex];

			TH1F* check_temp_th1 =  (TH1F*) gDirectory->Get(hist_name);

			if(check_temp_th1 !=NULL){//good we have the file,
				hist_root->cd();
				std::cout<<"\t  Proceed to next directory."<<std::endl;
				continue;
			}else{//no we don't have the histogram; then remove the direcory
				std::cout<<"\t Reproduce histograms in current directory "<<Tdir_name<<std::endl;
				Tdir->Delete(Tdir_name + ";*");//might not functioning correctly :( CHECK
			}

		} else{//No such directory, create it;
			Tdir = hist_root->mkdir(Tdir_name);
		}
		Tdir->cd();

		//STEP 2.2 prepare TH1F
		Long64_t nentries = temptree->GetEntries();
		if(message) std::cout<<" We have entries: "<<nentries<<std::endl;

		//STEP 2.2.1 prepare empty histgrams to be filled, each histogram with corresponding weight[kndex] (Note, not POT reweighting).
		std::vector<TH1F*> hist1d(num_throws);//hist
		std::vector<TTreeFormula*> wgtforms(num_throws);//variable
		std::vector<TTreeFormula*> varformula(num_throws);//weight, but not with POT reweighting;

		int count_hist = 0;
		for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
			//histogram of a throw
			TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];

				hist1d[kndex] = new TH1F(hist_name,hist_name, var->int_n_bins, var->plot_min, var->plot_max);//Will be rebinned later for variable binning;
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

				//					std::cout<<std::endl;
				//					if(kndex==10) exit(0);
				}
			}//next throw;
		}

		std::cout<<"\n Finish filling "<<num_swvars<<" histograms"<<std::endl;

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
			(temps->hist[jndex]).push_back( (TH1F*) hist1d[kndex]->Clone());// each bdt_sys now contains all hist;
			hist1d[kndex]->Delete();
		}

		std::cout<<"Save "<<num_throws;
		std::cout<<" histograms to the dir: "<<Tdir_name<<std::endl;

	}//next systematic weight 
	infile->Close();

}

/*
 * Build covariance matrices based on histograms;
 */

void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss, TString dir_root, TString dir_drawn, double plot_pot, unsigned long hashdraw){
		
	bool message = true;
	bool smooth_matrix = true;

	int stage = syss[0]->stage;
	
	if(message) std::cout<<"\tMaking covariance matices."<<std::endl;
	
	//STEP 1 group systematic hists for evaluating covariance matrices;
	std::vector<TString> tags;
	std::vector<TString> systags;
	std::vector<TString> bdtftags;
	std::multimap< TString, bdt_sys* > sys2SWmap;
	std::map< TString, bdt_sys* > sys2CVmap;
	std::map< TString, bdt_sys* > sys2OMCVmap;

//	bool only_one_cv = true;
	for(size_t index = 0; index < syss.size(); ++index){//first classify what tags are available, then group them next;

		bdt_sys* temps = syss[index];

		TString temptag = temps->tag;

		tags	.push_back(temptag);//for individual systematic x bdtfiles; will repeat for one of each: sw, CV, and maybe also OMCV;
		systags	.push_back(temps->systag);//for multi-cov; will contain duplicated
		bdtftags.push_back(temps->bdt_file::tag);//for multi-cov; will contain duplicated

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
	sort(tags.begin(), tags.end() );
	tags.erase( unique( tags.begin(), tags.end()), tags.end() );


	//STEP 2 pair up SW/CV to make covariance matrix;
	
	//STEP 2.1 prepare root file to store  covariance matrix and configure TH2D*;
	TString covroot_name = dir_root+"/"+gadget_labelroot("cov",var, stage, hashdraw)+".root";//the root to save following covariance marices;
	TFile *cov_root = (TFile*) TFile::Open(covroot_name,"RECREATE"); //Collection of individual covariance matrices;

	TString final_covroot_name = dir_root+"/../"+gadget_labelroot("finalcov", var, stage, hashdraw)+".root";
	TFile *finalcov_root = (TFile*) TFile::Open(final_covroot_name,"RECREATE"); //one var one hist root;


	int nb = var.n_bins;//n_bins also for rebinned binning for variable binning;
	int nl = var.plot_min;
	int nh = var.plot_max;
	
	TString finalcov_name = (var.safe_name).c_str();

	TH2D* finalcov =  new TH2D(finalcov_name, finalcov_name ,nb,nl,nh,nb,nl,nh);//binning for xnbins,xmin,xmax,ybins,ymin.ymax;

	TCanvas *histCanvas = new TCanvas("histCanvas", "", 600, 400);//Cavans to draw;
	TCanvas *covCanvas = new TCanvas("covCanvas","",600,400);

	//STEP 3 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	for(size_t index = 0; index < tags.size(); ++index){//each tag means one group of covariance matices;
		TString cur_tag = tags[index];
		std::cout<<"Working on tag "<<cur_tag<<std::endl;

		//Set binnings;
		std::vector<double> cur_binning;

		if(var.is_custombin){
			cur_binning = var.edges;
			cur_binning.erase(cur_binning.begin());//remove the first element;
		} else{
			double bin_left = nl; 
			double step = (nh-nl)/nb;

			//			std::cout<<"CHECK, binning:";
			for(int jndex = 0; jndex < nb+1; jndex++){
				cur_binning.push_back(bin_left);
				//				std::cout<<bin_left<<" "<<std::endl;
				bin_left+=step;
			}

		}

		//STEP 3.1 Find the bdt_sys syss that gives the CV; then do SW, OMCV on the flight;
		bdt_sys* tempsCV = sys2CVmap.find(cur_tag)->second;//bdt_sys that gives CV - tempsCV
		
		bool special_smoothing = (tempsCV->systag).Contains("OpticalModel");
		bool do_smooth = smooth_matrix && ((tempsCV->systag).Contains("Unisim")||special_smoothing);
		std::cout<<"\nWorking on "<<tempsCV->systag<<std::endl;
		if(do_smooth) std::cout<<tempsCV->systag<<" might need smoothing (if nbins>1)"<<std::endl;

		//Got the CV! Now Rebin and Scale with POT;
		TH1F* cv_hist = (TH1F*) (tempsCV->hist[0][0])->Clone();

		cv_hist->Scale(plot_pot/tempsCV->pot);//scale hist to data POT
		TH1F* smoothRef_hist = (TH1F*) cv_hist->Clone();
		if(var.is_custombin){ 
			cv_hist->Rebin(nb, cur_tag+"CV",&(cur_binning).front());//tags[index]+"CV" is just some name, can be anything;
			cv_hist = (TH1F*) gDirectory->Get(cur_tag+"CV");
		}


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

		//for Unisims, there will be two elements 

		int numSW = (tempsSW->hist).size();//hist is a 2d vector w dimension numSW x throws;
		for(int jndex = 0; jndex < numSW; ++jndex){//go through different SW under same bdt_sys
			//USE THIS LABEL!
			TString temp_sw_name  = cur_tag+"_"+tempsSW->vars_name[jndex];

			TH2D* all_hist = new TH2D(temp_sw_name, temp_sw_name, nb, &(cur_binning).front(), nby, 0, nby);//to store many throws SW
			TH2D* covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix", temp_sw_name+"_CovarainceMatrix",nb,nl,nh,nb,nl,nh);//to store covariance matrix, equal width;
			TH2D* ori_covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix_original", temp_sw_name+"_CovarainceMatrix_original",nb,nl,nh,nb,nl,nh);//to store covariance matrix, equal width;

			for(int kndex = 0; kndex < tempsSW->throws; kndex++){//go through different throws
				TH1F* sw_hist = (TH1F*) (tempsSW->hist[jndex][kndex])->Clone("copy"+tempsSW->vars_name[jndex]);
				///scale, smooth, rebin;
				sw_hist->Scale(plot_pot/tempsSW->pot);
				if(var.int_n_bins>1&& do_smooth ){ 
					std::cout<<"\r Smoothing a sw histogram";
					SmoothSW(sw_hist, smoothRef_hist, special_smoothing);
				}
				if(var.is_custombin){ 
					TString rebin_label =  cur_tag+"_"+tempsSW->vars_name[jndex]+"sw"+to_string_prec(kndex,0);
					sw_hist->Rebin(nb,rebin_label,&(cur_binning).front());//repeat label will get the previous histogram
					sw_hist = (TH1F*) gDirectory->Get(rebin_label);
				}

				//STEP 3.2.2 make histograms & covriance matrix;
				for(int lndex = 1; lndex < nb+1; ++lndex){//fill sw_hist to all_hist bins by bins;
					all_hist->Fill(sw_hist->GetBinLowEdge(lndex) , sw_hist->GetBinContent(lndex) );
				}//next bin

				//Make cov matrices
				TH2D* cov_temp = MakeCov("cov"+cur_tag+"_"+tempsSW->vars_name[jndex]+std::to_string(kndex), sw_hist, cv_hist);

				double scale_factor = 1.0;
				if(tempsSW->vars_name[jndex].Contains("KpProd") ) scale_factor = 1.5;//from Zarko's MaktrixMaker.cxx
				cov_temp->Scale(scale_factor/(double) total_throws);

				covmatrices->Add(cov_temp);

				if(kndex == tempsSW->throws-1 && to_plot_list.size()>1){//wait, one more throw! Reset
					to_plot_list.erase(to_plot_list.begin());
					kndex = -1;//will be 0 next round;	
					tempsSW = to_plot_list[0];
					std::cout<<"\nWait! More throws to append: "<<tempsSW->vars_name[0]<<std::endl;
				}

				//fill ori_covmatrices, no matter  smoothing or not;
				if(special_smoothing){
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

			//STEP 3.2.3 all_hist & covmatrices are filled; now draw them;
			gStyle->SetPalette(kRainBow);//this change colz colors;

			//Histograms;
			histCanvas->Clear();
			histCanvas->cd();
			all_hist->SetStats(false);
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
			TH2D* fracCov = MakeFracCov(temp_sw_name, covmatrices, cv_hist, special_smoothing, ori_covmatrices);//special_smoothing - false, then the last argument is useless.
			fracCov->SetStats(false);
			fracCov->Draw("COLZ");
			fracCov->SetTitle("Fractional "+temp_sw_name + " Covaraince Matrix");
			covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covFrac",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");

			if(message){ 
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

				if(message) std::cout<<"\nPropagating covariance matrix for "<<cur_tag<<std::endl;
				TH1F* omcv_hist = (TH1F*) (tempsOMCV->hist[0][0])->Clone();

				//rebin, scale
				if(var.is_custombin){ 
					omcv_hist->Rebin(nb, temp_sw_name+"OMCV",&(cur_binning).front());
					omcv_hist = (TH1F*) gDirectory->Get(temp_sw_name+"OMCV");
				}
				omcv_hist->Scale(plot_pot/tempsCV->pot);

				//make propagated covariance matrix;
				TH2D* ProCov = PropagateCov(temp_sw_name, fracCov, omcv_hist);
				ProCov->SetStats(false);
				ProCov->Draw("COLZ");
				ProCov->SetTitle("Propagated "+temp_sw_name + " Covaraince Matrix");
				covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covProp",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");

				if(message){ 
					std::cout<<"The matrix diagonal element "<< temp_sw_name<<": "<<std::endl;;
					for(int bindex = 1; bindex<nb+1; ++bindex){
						std::cout<<" "<<ProCov->GetBinContent(bindex,bindex)<<",";
					}
					std::cout<<std::endl;
				}

				finalcov->Add(ProCov);

				cov_root->cd();
				fracCov->Write();
				ProCov->Write();
			} else{//no need extra processing;

				covmatrices->SetStats(false);
				covmatrices->GetXaxis()->SetTitle((var.unit).c_str());
				covmatrices->Draw("COLZ");
				covCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("cov",var, stage, hashdraw)+"_"+temp_sw_name+".pdf" ,"pdf");

				if(message){ 
					std::cout<<"The matrix diagonal element "<<temp_sw_name<<": "<<std::endl;
					for(int bindex = 1; bindex<nb+1; ++bindex){
						std::cout<<" "<<covmatrices->GetBinContent(bindex,bindex)<<",";
					}
					std::cout<<std::endl;
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
			double fjk = (ec_smoothing && jndex==kndex)? (covjk+cvj)/(cvj*cvk): covjk/(cvj*cvk);//not divide by nb;

			if(ec_smoothing&&jndex==kndex && fjk > origin_cov->GetBinContent(jndex,kndex)/(cvj*cvk) ){
				std::cout<<"\nAdjust matrix element to a smaller value: "<<jndex<<" "<<kndex<<": "<<fjk<<" to ";
				fjk = origin_cov->GetBinContent(jndex,kndex)/(cvj*cvk);
				std::cout<<fjk<<std::endl;
			}

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
void SmoothSW(TH1F* sw, TH1F* cv, bool ec_smoothing){//fit each bin with polynomials
{
  using namespace std;
  bool print_ratio = false;


  Int_t degree = 2;
  Int_t first_bin = 1;
  Int_t last_bin = cv->GetNbinsX() ;
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
	  gSystem->RedirectOutput(0,0);//this let ROOT warns

	  sw->Reset();
	  for (Int_t bin=1;bin<Nbins+1;bin++) {//Here modify the sw
		  sw->SetBinContent(bin, cv->GetBinContent(bin) * smoothFitter->Eval(sw->GetBinCenter(bin)));
	  }

		
	  TH1F* copy_rsw = (TH1F*) sw->Clone();
	  copy_rsw->Divide(cv);
	  copy_rsw->SetLineColor(kRed);
	  copy_rsw->Draw("same");

	  if(print_ratio){  
		rCanvas->SaveAs(("/scratch/condor-tmp/klin/hellstroms_hive/hive/build/src/fullmcsystematics/drawn/ratio_OM_"+to_string_prec(global_index,0)+".pdf").c_str(),"pdf");
		  rCanvas->Delete();
		  global_index++;
	  }

  }else{//use Zarko's
  //solve M a = y

//  cout << "Smoothing bins " << first_bin << " - " << last_bin << " (N=" << Nbins<< ") with pol "<<degree<<endl;
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
//      cout <<"smoothed "<<i+first_bin<<"\t"<<num<<endl;
      sw->SetBinContent(i+first_bin,num);
    }
  }
  else cout <<"failed to solve"<<endl;
}
}

}

/*
 * Make Covaraince matrix from 2d histograms
 *
 */

TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv){

}

