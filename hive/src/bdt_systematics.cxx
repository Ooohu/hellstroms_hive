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
 * Constructor for struct bdt_sys, CHECK
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

	if(message) std::cout<<"\n ---- Initialize systematics for "<<var1->unit<<std::endl;

	//check 1) for the covariance matrix;
	TString check_covroot_name = dir_root+"/../"+gadget_labelroot("finalcov", *var1, stage, hashlabel)+".root";
	TFile *check_cov_root = (TFile*) TFile::Open(check_covroot_name,"READ"); //one var one hist root;

	if(check_cov_root != NULL ){ //oh we had the file! so .. 
		std::cout<<" Being lazy, no need to make cov matrix, because we had it already!"<<std::endl;
		check_cov_root->Close();

		skip_sysEvaluation = true;
	}
	delete check_cov_root;


	TString readfile_method = "RECREATE";//create file by default, if exist, update the current file;
	TString histfilename = dir_root+"/"+gadget_labelroot("hist", *var1, stage, hashlabel)+".root";

	if(do_hist || !skip_sysEvaluation){//no root fiels for the covariance matix, so going to generate it.
		TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;

		size_t index = 0;//index for syss, determine where we left over last time the code break;

		if(!do_cov && hist_root!=NULL ){ //quick check of existing files; 
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
//							if(jndex > 0){//well, at least one good histogram;
								temps->start_with_weight = jndex;
								rmdir(temps->tag+"_"+weights_nam[jndex]);
//
//								(temps->hist[jndex-1]).clear();

								readfile_method = "UPDATE";
//							}//if jndex == 0; it means just redo everything; readfile_method = "RECREATE";

//							(temps->hist[jndex]).clear();//clear the current histogram, which might be corrupted;
							hist_root->Close();

							std::cout<<"\n>> Find a broken directory, reproduce histograms.starting from "<<syss[index]->tag<<std::endl;

							goto now_make_hist;//do partial generation, index determines which systematic file to start;
						} else{

							(temps->hist[jndex]).push_back( (TH1F*) temp_th1->Clone() );// each bdt_sys now contains all hist;
							//						std::cout<<"Throw "<<kndex<<" "<<temp_th1->Integral()<<std::endl;
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
		std::cout<<temps->dir+"/"+temps->filename<<" does not exist, please check! Program abort."<<std::endl;
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

		TDirectory* temp_dir = (TDirectory*) hist_root->Get(Tdir_name);

		if(temp_dir !=NULL){ 
			std::cout<<Tdir_name<<" exists, and it is properly a good leftover; proceed to next directory."<<std::endl;
			continue;
		}

		TDirectory *Tdir = hist_root->mkdir(Tdir_name);
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

//			if(var->is_custombin){
//				hist1d[kndex] = new TH1F(hist_name,hist_name, (var->edges).size() -1, &(var->edges).front() );
//			} else{
				hist1d[kndex] = new TH1F(hist_name,hist_name, var->int_n_bins, var->plot_min, var->plot_max);//Will be rebinned later for variable binning;
//			}
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

	TCanvas *histCanvas = new TCanvas("histCanvas", "", 900, 400);//Cavans to draw;
	TCanvas *covCanvas = new TCanvas("covCanvas","",900,400);

	//STEP 3 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	for(size_t index = 0; index < tags.size(); ++index){//each tag means one group of covariance matices;
		TString cur_tag = tags[index];
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

		bool do_smooth =  (tempsCV->systag).Contains("Unisim")||(tempsCV->systag).Contains("OpticalModel");
		if(do_smooth) std::cout<<tempsCV->systag<<" might need smoothing (if nbins>1)"<<std::endl;

		//Got the CV! Now Rebin and Scale with POT;
		TH1F* cv_hist = (TH1F*) (tempsCV->hist[0][0])->Clone();

		if(var.is_custombin){ 
			cv_hist->Rebin(nb, cur_tag+"CV",&(cur_binning).front());//tags[index]+"CV" is just some name, can be anything;
			cv_hist = (TH1F*) gDirectory->Get(cur_tag+"CV");
		}

		cv_hist->Scale(plot_pot/tempsCV->pot);//scale hist to data POT

		//STEP 3.2 Load SW and draw histograms, covariance matrix only the flight;
		int nby = 1.5*cv_hist->GetBinContent(cv_hist->GetMaximumBin());//set histogram height

		for( std::_Rb_tree_iterator<std::pair<const TString, bdt_sys*> >  itr = sys2SWmap.begin();  itr!=sys2SWmap.end(); itr++){//loop over bdt_sys for SW;

			if(itr->first == cur_tag){//find the corresponding SW, many throws TH1F* of a associated bdt_sys;
				std::cout<<"Working on tag "<<cur_tag<<std::endl;

				//STEP 3.2.1 prepare SW of throws;
				bdt_sys* tempsSW = itr->second;//tempSW may contains  manySW x throws;

				TH2D* all_hist = new TH2D(cur_tag, cur_tag, nb, &(cur_binning).front(), nby, 0, nby);//to store many throws SW
				TH2D* covmatrices =  new TH2D(cur_tag+"_CovarianceMatrix", cur_tag+"_CovarainceMatrix",nb,nl,nh,nb,nl,nh);//to store covariance matrix, equal width;

				int numSW = (tempsSW->hist).size();//hist is a 2d vector w dimension numSW x throws;

				for(int jndex = 0; jndex < numSW; ++jndex){//go through different SW under same bdt_sys
					TString temp_sw_name  = cur_tag+"_"+tempsSW->vars_name[jndex];

					for(int kndex = 0; kndex < tempsSW->throws; ++kndex){//go through different throws
						TH1F* sw_hist = (TH1F*) (tempsSW->hist[jndex][kndex])->Clone();

						///scale, smooth, rebin;
						sw_hist->Scale(plot_pot/tempsSW->pot);
						if(nb>1&& do_smooth ){ 
							std::cout<<"\r Smoothing a sw histogram";
							SmoothSW(sw_hist, cv_hist);
						}
						if(var.is_custombin){ 
							sw_hist->Rebin(nb, cur_tag+"sw",&(cur_binning).front());
							sw_hist = (TH1F*) gDirectory->Get(cur_tag+"sw");
						}

						//STEP 3.2.2 make histograms & covriance matrix;
						for(int lndex = 1; lndex < nb+1; ++lndex){//fill sw_hist to all_hist bins by bins;
							all_hist->Fill(sw_hist->GetBinLowEdge(lndex) , sw_hist->GetBinContent(lndex) );
						}//next bin

						//Make cov matrices
						TH2D* cov_temp = MakeCov("cov"+cur_tag+std::to_string(kndex), sw_hist, cv_hist);
						cov_temp->Scale(1.0/(double) tempsSW->throws);
						covmatrices->Add(cov_temp);
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


					if(sys2OMCVmap.count(tags[index]) > 0){//its a Optical Model CV, propagate the matrix to a fractional covariance matrix; and add it to finalcov.
						bdt_sys* tempsOMCV = (sys2OMCVmap.find(cur_tag))->second;//get the sys for a given tag;

						if(message) std::cout<<"Making fractional covariance matrix for "<<cur_tag<<std::endl;
						TH1F* omcv_hist = (TH1F*) (tempsOMCV->hist[0][0])->Clone();

						//rebin, scale
						if(var.is_custombin){ 
							omcv_hist->Rebin(nb, cur_tag+"OMCV",&(cur_binning).front());
							omcv_hist = (TH1F*) gDirectory->Get(cur_tag+"OMCV");
						}
						omcv_hist->Scale(plot_pot/tempsCV->pot);

						//make fractional covariance matrix, but this does not added to final covariance matrix;
						TH2D* fracCov = MakeFracCov(temp_sw_name, covmatrices, cv_hist);
						fracCov->SetStats(false);
						fracCov->Draw("COLZ");
						fracCov->SetTitle("Fractional "+cur_tag + " Covaraince Matrix");
						covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covFrac",var, stage, hashdraw)+"_"+cur_tag+".pdf" ,"pdf");


						//make propagated covariance matrix;
						covCanvas->Clear();
						covCanvas->cd();
						TH2D* ProCov = PropagateCov(temp_sw_name, fracCov, omcv_hist);
						ProCov->SetStats(false);
						ProCov->Draw("COLZ");
						ProCov->SetTitle("Propagated "+cur_tag + " Covaraince Matrix");
						covCanvas->SaveAs( dir_drawn +  "/"+gadget_labelroot("covProp",var, stage, hashdraw)+"_"+cur_tag+".pdf" ,"pdf");

						if(message){ 
							std::cout<<"The matrix diagonal element "<< cur_tag<<": "<<std::endl;;
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
						covCanvas->SaveAs( dir_drawn + "/"+gadget_labelroot("cov",var, stage, hashdraw)+"_"+cur_tag+".pdf" ,"pdf");

						if(message){ 
							std::cout<<"The matrix diagonal element "<<cur_tag<<": "<<std::endl;
							for(int bindex = 1; bindex<nb+1; ++bindex){
								std::cout<<" "<<covmatrices->GetBinContent(bindex,bindex)<<",";
							}
							std::cout<<std::endl;
						}

						finalcov->Add(covmatrices);
						cov_root->cd();
						covmatrices->Write();
					}
				}//next SW
			}//end of if statement; this finishes making plots, bc we got the corresponding bdt_sys
		}//next bdt_sys;
	}//next tag;

	//STEP 4 draw out the final covariance matrix;

	cov_root->Close();//basic covariance matrix is recorded, now prepare the one for plotting;
	
	//set axis and name for the final covaraince matrix;
	finalcov->SetStats(false);

	finalcov->SetTitle("Total Covariance Matrix");
	finalcov->GetXaxis()->SetTitle( (var.unit).c_str());

	TCanvas *c1 = new TCanvas("c1","",900,400);
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

		bins.push_back(hist->GetBinLowEdge(index));
	}
	//	cout<<"Bin config # bins: "<<nb<<" low "<<nl<<" high "<<nh<<" entries:"<<hist->Integral()<<endl;

	//	cout<<"Total events in (MakeCov) "<<hist->Integral()<<endl;
	int counter = 0;
	if(nb!=cv->GetNbinsX() || hist->GetBinLowEdge(1) !=cv->GetBinLowEdge(1)||hist->GetBinLowEdge(nb+1)!=cv->GetBinLowEdge(nb+1) ){
		std::cerr<<" Histograms don't match to the CV histogram. No Covariance matrix is calculated."<<std::endl;
		exit(EXIT_FAILURE);
	}
	TH2D* covmatrix =  new TH2D(name,name,nb, &(bins.front()),nb, &(bins.front()) );//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
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
 * Create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TString name,TH2D* cov_temp, TH1F* oldcv){

	int nb=oldcv->GetNbinsX();
	int nl=oldcv->GetBinLowEdge(1);
	int nh=oldcv->GetBinLowEdge(nb+1);
	
	TH2D* fractional_cov = (TH2D*) cov_temp->Clone("Fractional_Covaraince_Matrix_"+name);

//	fractional_cov->Reset("ICESM");//clear contents;

	for(int jndex = 1; jndex < nb+1; ++jndex){
		for(int kndex = 1; kndex < nb+1; ++kndex){
			double cvj = oldcv->GetBinContent(jndex);
			double cvk = oldcv->GetBinContent(kndex);
			double fjk = cov_temp->GetBinContent(jndex,kndex)/(cvj*cvk);

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
void SmoothSW(TH1F* sw, TH1F* cv){//fit each bin with polynomials

	Int_t nb=sw->GetNbinsX();
	Int_t degree = 5;

//	std::vector< Double_t >					vectorA(degree+1, 0.0);
//	std::vector< std::vector<Double_t> >	matrixA(degree+1, 
//							(std::vector<Double_t>)(degree+1, 0.0));

	TVectorD vectorA(degree+1);
	TMatrixD matrixA(degree+1, degree+1);

	for (Int_t j=0;j<degree+1;j++) {
		vectorA[j] = 0.;//initialize the vector;
//		std::cout<<"Initialize vector A "<<vectorA[j]<<std::endl;
		for (Int_t k=0;k<degree+1;k++) {
			matrixA[j][k] = 0.;//initialize the matrix;

			for (Int_t bini=0; bini<nb; bini++) {
				if(k==0){// do it once for each j;
					Double_t temp_cv = cv->GetBinContent(bini+1);
					if(temp_cv < 10e-20){//in case of divide by 0;
						temp_cv = sw->GetBinContent(bini+1);
					}
					if(temp_cv > 10e-20){
						//calculate vectorA
//						std::cout<<"Check vector A "<<vectorA[j]<<std::endl;
						vectorA[j] += pow( bini+1, j) * sw->GetBinContent(bini+1)/temp_cv;
//						std::cout<<"Vector A, before calculation "<<pow( bini+1, j)<<" "<<sw->GetBinContent(bini+1)<<" "<<temp_cv<<" "<<vectorA[j]<<std::endl;
					}
				}
				//calculate matrixA
				matrixA[j][k] += pow(bini+1,j)* pow(bini+1,k);
			}
		}
	}

	TVectorD vectorB(degree+1);
	TDecompLU lu(matrixA);
	lu.Decompose();
	Bool_t ok;

	vectorB = lu.Solve(vectorA,ok);
	if(ok){
		for (Int_t bini = 0; bini<nb; bini++){
			Double_t num = 0.;
			for(Int_t j=0; j<degree+1; j++){
//			std::cout<<" vector A"<<vectorA[j]<<std::endl;
				num+=vectorB[j]*pow(bini+1, j)*cv->GetBinContent(bini+1);
			}
//			std::cout<<"Smooth "<<bini+1<<"\t"<<sw->GetBinContent(bini+1)<<" to "<<num<<std::endl;
			sw->SetBinContent(bini+1 , num);
		}
	} else{
//		std::cout<<"Cant smooth this one"<<std::endl;
	}
	//step 1 take ratio;

	//Modify ratio

	//Recover SW

}

/*
 * Make Covaraince matrix from 2d histograms
 *
 */

TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv){

}

