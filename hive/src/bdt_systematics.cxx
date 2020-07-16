#include <bdt_systematics.h>

/*
 * Some gadgets to keep things organized;
 */

TString gadget_nameroot(TString plot_type, TString dir, bdt_variable var){

	TString binning_text = to_string_prec(var.n_bins+var.plot_min+var.plot_max,0);
	TString file_name = dir+"/" + plot_type+"_"+var.safe_name +  binning_text + ".root";

	return file_name;
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
	its_multifiles	(XMLconfig.sys_its_multifiles[index])
	{
		num_vars = vars.size();
		hist.resize(num_vars);
		twodhist.resize(num_vars);
	
	};

//TH1F GetATH1(TTree &temptree, bdt_variable var, TString hist_nam, TString modified_wgt){
//
//	TH1F th1(TH1F(hist_nam, hist_nam, var.n_bins, var.plot_min, var.plot_max));
//	temptree.Draw(var.mininame+">>"+hist_nam, modified_wgt,"goff");//This draw caus memory leak;
//
//	return th1;
//}

/*
 * Initialize Systematics of given variable 1dhistogram of weights by preparing histograns;
 */
void InitSys(bdt_variable var, std::vector<bdt_sys*> syss, double plot_pot, TString dir_root, TString dir_drawn){
	
	bool message = true;
	
	//Part I make histiograms;
	//work on each systematic

	if(message) std::cout<<"\n ---- Initialize systematics for "<<var.unit<<std::endl;

	TString outfile_name = gadget_nameroot("hist",dir_root, var);//, temps);


	TFile *hist_root = (TFile*) TFile::Open(outfile_name,"READ"); //one var one hist root;
	if(hist_root!=NULL && false){ //quick check of existing files;
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

			bool its_CV = temps->its_CV;//use Totalweight, one throw; seems not useful, CHECK
			bool its_multithrows = temps->its_multithrows;//multi weights or optical model
			bool its_multifiles = temps->its_multifiles;//optical model has multifiles;

			// type 1 single throw one file
			//STEP 1 Load files, TTree
			gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
			TFile *infile = (TFile*) TFile::Open(temps->dir+"/"+temps->filename,"READ"); 
			gSystem->RedirectOutput(0,0);//this let ROOT warns
			TTree* temptree = (TTree*) infile->Get(temps->treename);


			//STEP 2 get the TH1F
			for(int jndex = 0; jndex < temps->num_vars; ++jndex){//loop through weights
				//subdirectory in root for each systematic weighs;

				TString weights_nam = temps->vars_name[jndex];//name of key for the histogram in root;
				TString Tdir_name = temps->tag+"_"+weights_nam;

				TDirectory *Tdir = hist_root->mkdir(Tdir_name);
				Tdir->cd();
				//			TString binning_text = to_string_prec(var.n_bins+var.plot_max+var.plot_max,0);
				//			TString outfile_name = dir_root+"/"+var.safe_name + binning_text  + temps->tag + "_" + weights_nam +".root";
				//			TString outfile_name = gadget_nameroot("hist",dir_root, var);//, temps);

				//CHECK, add skipping for previously done root;


				//STEP 2.2 prepare TH1F
				//			TFile *outfile = (TFile*) TFile::Open(outfile_name,"RECREATE"); 
				//			std::cout<<dir_root+"/"+outfile_name<<std::endl;
				int num_throws = temps->throws;
				TDirectory::AddDirectory(0);
				for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;

					//if(kndex%5==1){//reopen the input root file doesnt help;
					//	std::cout<<"Take a break"<<std::endl;
					//	infile->Close();	
					//	infile = (TFile*) TFile::Open(temps->dir+"/"+temps->filename,"READ"); //one var one hist root;
					//	temptree =(TTree*) infile->Get(temps->treename); 
					//	hist_root->cd();
					//}
					//if(kndex%5==1){//reopen the output root file doesnt help;
					//	std::cout<<"Take a break"<<std::endl;
					//	hist_root->Close();	
					//	hist_root = (TFile*) TFile::Open(outfile_name,"RECREATE"); //one var one hist root;

					//	hist_root->cd();
					//}
					TString hist_nam_label = (its_multithrows)? std::to_string(kndex+1) : "";
					TString hist_nam = weights_nam + hist_nam_label;
					TString modified_wgt = (its_multithrows)? temps->vars[jndex] + "["+std::to_string(kndex) +"]": temps->vars[jndex];

//					TH1F th1 = GetATH1(*temptree,var, hist_nam, modified_wgt);


//					TCanvas *tempC = new TCanvas("tempC", "", 900, 400);
//					TPad *tempP = new TPad("tempP", "tempP", 0, 0, 1, 1.0);
//					tempP->cd();
//					tempC->cd();
					//why no popup window? ,because I didn't call root

//					TH1F* th1 = new TH1F(hist_nam, hist_nam, var.n_bins, var.plot_min, var.plot_max);
//					temptree->Draw(var.mininame+">>"+hist_nam, modified_wgt,"goff");//This draw cause memory leak;
					
					temptree->Draw(var.name+">>"+hist_nam+var.binning, modified_wgt,"goff");//This draw cause memory leak; Whenever a histogram is drawn, the memory will be consumed and be gone. Once the memory is out, then the code will be killed - fail.

					TH1F* th1= (TH1F*)gDirectory->Get(hist_nam); 
//					th1->AddDirectory(0);
					if(kndex==0){ 
						gSystem->RedirectOutput("/dev/null");
						th1->Sumw2();//the error will be [sum of sqrt(weights)] // CHECK
						gSystem->RedirectOutput(0,0);
					}

					th1->Scale(plot_pot/temps->pot);

					//configure th1
					th1->SetStats(0);
					th1->GetXaxis()->SetTitle(var.unit.c_str());
					th1->GetYaxis()->SetTitle("Events");
				
					th1->SetDirectory(0);
					//save the th1 to the root file, class bdt_sys;
					th1->Write(hist_nam);

					(temps->hist[jndex]).push_back( (TH1F*) th1->Clone() );// each bdt_sys now contains all hist;

					if(message){
						std::cout<<"\r\tPreparing "<<std::setw(9)<<weights_nam<<", throws: "<<kndex+1<<"/"<<num_throws;
						//save them according to types:
					}
					th1->Reset();	
//					th1->Delete();
//					std::cout<<th1<<std::endl;

//					th1->~TH1();
//					gDirectory->DeleteAll();
//					gDirectory->Delete(hist_nam);//useless
//					delete th1;

//					std::cout<<"Going to free"<<std::endl;
//					sleep(4);
//					tempP->Close();
//					tempC->Close();
//					std::cout<<"free"<<std::endl;
//					sleep(4);

				//	delete hist_nam_label;
				//	delete hist_nam;
//				//	delete modified_wgt;
//					delete th1;
//					free(th1);//nothing is free.. LOL
					gDirectory->Remove(th1);
				}//next throw;
				std::cout<<"\n Save "<<num_throws;
				std::cout<<" histograms to the file: "<<outfile_name<<std::endl;
//				sleep(10);
			Tdir->DeleteAll();
			delete Tdir;
			}//next systematic weight 

			// type 3 single throw multi files 

			infile->Close();
		}//next systematic files
	}	
	std::cout<<"Finish making systematic histograms."<<std::endl;
	sleep(10);
	//
	//
	//Part II: Overlay histograms & draw covaraicne matrix;
	//
	//
	//CHECK, need to load syss[?]->hist manually, if not made above;
	hist2cov( var, syss, dir_root, dir_drawn);
	hist_root->Close();
}

/*
 * InitSys version 2, use branches;
 *
 */

void InitSys2(bdt_variable var, std::vector<bdt_sys*> syss, double plot_pot, TString dir_root, TString dir_drawn){
	
	bool message = true;
	
	//Part I make histiograms;
	//work on each systematic

	if(message) std::cout<<"\n ---- Initialize systematics for "<<var.unit<<std::endl;

	TString outfile_name = gadget_nameroot("hist",dir_root, var);//, temps);


	TFile *hist_root = (TFile*) TFile::Open(outfile_name,"READ"); //one var one hist root;
	if(hist_root!=NULL){ //quick check of existing files;
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

			bool its_CV = temps->its_CV;//use Totalweight, one throw; seems not useful, CHECK
			bool its_multithrows = temps->its_multithrows;//multi weights or optical model
			bool its_multifiles = temps->its_multifiles;//optical model has multifiles;

			// type 1 single throw one file
			//STEP 1 Load files, TTree
			gSystem->RedirectOutput("/dev/null");//no warning, shut up! thanks https://root-forum.cern.ch/t/suppress-all-root-info-warning-error-output/30766
			TFile *infile = (TFile*) TFile::Open(temps->dir+"/"+temps->filename,"READ"); 
			gSystem->RedirectOutput(0,0);//this let ROOT warns
			TTree* temptree = (TTree*) infile->Get(temps->treename);


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
				const int num_mininames = (var.mininame).size();

				Long64_t nentries = temptree->GetEntries();
				if(message) std::cout<<" We have entries: "<<nentries<<std::endl;

				//Get branches
				temptree->SetMakeClass(1);

				temptree->SetBranchStatus("*",0);

				std::vector< Float_t * > swvars(num_swvars); 
				std::vector< Float_t * > varvars(num_mininames); 
				//To load Float_t []: https://arduino.stackexchange.com/questions/3774/how-can-i-declare-an-array-of-variable-size-globally
				for(size_t kndex = 0; kndex < (temps->vars).size(); ++ kndex){//get sw branches
					swvars[kndex] = (Float_t*) malloc(nentries*sizeof(Float_t*)*num_throws);
					temptree->SetBranchStatus(temps->vars[kndex],1);
					temptree->SetBranchAddress(temps->vars[kndex], swvars[kndex]);
				}

				for(size_t kndex = 0; kndex < (var.mininame).size(); ++ kndex){//get sw branches
					varvars[kndex] = (Float_t*) malloc(nentries*sizeof(Float_t*)*num_throws);
					temptree->SetBranchStatus(var.mininame[kndex],1);
					temptree->SetBranchAddress(var.mininame[kndex], varvars[kndex]);
				}

				std::vector<TH1F*> hist1d(num_throws);
				std::vector<TH2F*> hist2d(num_throws); //if num_mininames  == 2;

				int count_hist = 0;
				for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
					TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];
					hist1d[kndex] = new TH1F(hist_name,hist_name, var.n_bins, var.plot_min, var.plot_max);
					count_hist++;
				}
				if(message) std::cout<<"Creating "<<count_hist<<" histograms "<<std::endl;

				for(Long64_t kentry = 0; kentry< nentries; ++kentry){

					temptree->GetEntry(kentry);

//					if(kentry>2) exit(0);
					for(int kndex = 0; kndex < num_throws; ++kndex){//loop through throws;
					hist1d[kndex]->Fill(varvars[0][0],swvars[jndex][kndex]*plot_pot/temps->pot);//jndex - nth sw; kndex - nth throw
					std::cout<<"\r\tLooping entries: "<<std::setw(5)<<kentry+1<<"/"<<std::setw(5)<<nentries;
					std::cout<<std::setw(5)<<kndex+1<<" throw sw name: "<<temps->vars_name[jndex];
					}//next throw;
				}

				std::cout<<"\n Finish filling "<<num_swvars<<" histograms"<<std::endl;

				//write hisograms to the diretory
				for(int kndex = 0; kndex <num_throws;++kndex){//throw in each weight;
					//configure th1
					TString hist_name = (its_multithrows)? temps->vars_name[jndex]+std::to_string(kndex+1) :temps->vars_name[jndex];
					hist1d[kndex]->SetStats(0);
					hist1d[kndex]->GetXaxis()->SetTitle(var.unit.c_str());
					hist1d[kndex]->GetYaxis()->SetTitle("Events");
				
					//save the th1 to the root file, class bdt_sys;
					Tdir->cd();
					hist1d[kndex]->Write(hist_name);
					(temps->hist[jndex]).push_back( (TH1F*) hist1d[kndex]->Clone());// each bdt_sys now contains all hist;
					hist1d[kndex]->Delete();
				}

				std::cout<<"Save "<<num_throws;
				std::cout<<" histograms to the file: "<<outfile_name<<" dir:"<<Tdir_name<<std::endl;
				//				sleep(10);
				//Tdir->DeleteAll();
				//delete Tdir;
				//free memory
				for(size_t kndex = 0; kndex < (temps->vars).size(); ++ kndex){//get sw branches
					free(swvars[kndex]); 
				}

				for(size_t kndex = 0; kndex < (var.mininame).size(); ++ kndex){//get sw branches
					free(varvars[kndex]);
				}
			}//next systematic weight 
			infile->Close();
		}//next systematic files
	}
	std::cout<<"Finish making systematic histograms."<<std::endl;
	//
	//Part II: Overlay histograms & draw covaraicne matrix;
	//
	hist2cov( var, syss, dir_root, dir_drawn);
	hist_root->Close();
}
/*
 * Build covariance matrices based on histograms;
 */

void hist2cov( bdt_variable var, std::vector<bdt_sys*> syss, TString dir_root, TString dir_drawn){
		
	bool message = true;
	
	if(message) std::cout<<"Making covariance matices."<<std::endl;
	
	//
	//STEP 1 categorize hists; 
	std::vector<TString> tags;
	std::multimap< TString, bdt_sys* > sys2SWmap;
	std::map< TString, bdt_sys* > sys2CVmap;

//	bool only_one_cv = true;
	for(size_t index = 0; index < syss.size(); ++index){//first classify what tags are available, then group them next;
		TString temptag = syss[index]->tag;
		tags.push_back(temptag);

		bdt_sys* temps = syss[index];
		
		if(temps->its_CV){//save CV to map;
			sys2CVmap.insert(std::make_pair (temptag, temps));

//			if(only_one_cv){ 
//				only_one_cv = false;
//			} else{
//				std::cerr<<" Error: did you put 2 CVs under the same tag?"<<std::endl; 
//				exit(EXIT_FAILURE);
//			}

		} else{
			sys2SWmap.insert(std::make_pair (temptag, temps));
		}//save SW to map
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

		//std::cout<<__LINE__<<" size of bdt_sys->hist "<< tempsCV->hist.size()<<std::endl;
		//std::cout<<__LINE__<<" size of bdt_sys->hist[0] "<< tempsCV->hist[0].size()<<std::endl;
		//std::cout<<__LINE__<<" address of bdt_sys->hist[0][0] "<< tempsCV->hist[0][0]<<std::endl;
		//std::cout<<__LINE__<<" content of bdt_sys->hist[0][0] "<< tempsCV->hist[0][0]->GetBinContent(1)<<std::endl;

		TH1F* cv_hist = (TH1F*) (tempsCV->hist[0][0])->Clone();
		std::cout<<__LINE__<<std::endl;

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

		for( std::_Rb_tree_iterator<std::pair<const TString, std::vector<TH1F*> > >  itr = sw_name2throws.begin();  itr!=sw_name2throws.end(); itr++){
			
		std::cout<<__LINE__<<std::endl;
			//STEP 2.2.1 prepare covariance matrices and histograms;
			//work through the map< TString, std::vector<TH1F*> > sw_name2throws;
			TString temp_histname = itr->first;
			TString temp_covname = itr->first;//"<tag>_<sw> Covariance Matrix"
			std::vector<TH1F*> temp_swhist = itr->second;
			int num_throws = temp_swhist.size();

			//for ovarlaying all histograms
			TH2D* all_hist = new TH2D(temp_histname, temp_histname, nb, nl,nh, nby, 0, nby);
		std::cout<<__LINE__<<std::endl;
			//for summing up cov matrice
			TH2D* covmatrices =  new TH2D(temp_covname, temp_covname ,nb,nl,nh,nb,nl,nh);//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
		std::cout<<__LINE__<<std::endl;

			//to record the cov matrix;
		std::cout<<__LINE__<<std::endl;
			for(size_t jndex = 0; jndex < num_throws; ++jndex){//go through different throws

				//for overlay hsitograms
				for(int kndex = 1; kndex < nb+1; ++kndex){//get everything from a single histogram;
					all_hist->Fill(temp_swhist[jndex]->GetBinLowEdge(kndex) , temp_swhist[jndex]->GetBinContent(kndex) );
				}//next bin

				//for cov matrices
				TH2D* cov_temp =MakeCov("cov"+temp_covname+std::to_string(jndex), temp_swhist[jndex], cv_hist);
				cov_temp->Scale(1.0/(double) num_throws);
				covmatrices->Add(cov_temp);

			}//next throws
			
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
			histCanvas->SaveAs( dir_drawn + "/hist_"+temp_histname+".pdf" ,"pdf");
			
			covCanvas->Clear();
			covCanvas->cd();
			covmatrices->SetStats(false);
			covmatrices->Draw("COLZ");
			covCanvas->SaveAs( dir_drawn + "/cov_"+temp_histname+".pdf" ,"pdf");
			
			finalcov->Add(covmatrices);

			cov_root->cd();
			covmatrices->Write();

			if(false){//do matrix propagation;

			}
		}//next sw
	}//next tag;
	cov_root->Close();
	
	//set axis and name for the final covaraince matrix;
	finalcov->SetStats(false);
	finalcov->SetTitle(" Total Covariance Matrix");
	finalcov->GetXaxis()->SetTitle( (var.unit).c_str());

	TCanvas *c1 = new TCanvas("c1","",900,400);
	finalcov->Draw("COLZ");
	c1->SaveAs(dir_drawn+"/Total_cov_"+var.safe_name+".pdf","pdf");

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
 * Make Covaraince matrix from 2d histograms
 *
 */

TH2D* Make2DCov(TString name,TH2D* hist, TH2D* cv){

}

