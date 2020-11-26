#include <bdt_systematics.h>

int global_index = 0;

/*
 * Some gadgets to keep things organized;
 */

//rebin hitogram to the given binning;
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
 *3 - BAD: underflow bins: cur_binning range too large (1,2,4) - (0,1,4)
 *4 - BAD: bins edge not found in cv_hist (1,2,5) - (1,2,3,4)
 */
int gadget_BinMatcher(TH1D* cv_hist, std::vector< double > cur_binning){

			bool debug_message = false;

			int binchecker = 0;//monitor if there are difference between var binning and histogram binning;
			int mis_count = 0;
			int edges = cur_binning.size();
			int histnb = cv_hist->GetNbinsX();
			int bindex_starter = 1;
			for(double binedge : cur_binning){//go through each bin in cur_binning;
					if(debug_message)std::cout<<binedge<<" vs ";
				for(int bindex = bindex_starter; bindex < histnb+2; ++bindex){//loop through lowegdes of cv_hist 
					if(debug_message)std::cout<<cv_hist->GetBinLowEdge(bindex)<<",";
					if(abs(cv_hist->GetBinLowEdge(bindex) - binedge) < 10e-10){//edge matches;
						binchecker++;
						bindex_starter = bindex+1;
						if(debug_message) std::cout<<"Good match"<<std::endl;
						break;
					}
					if(binchecker>1 && binchecker < edges) mis_count++;
					//if(debug_message) std::cout<<std::endl;
				}
			}

			if(binchecker == edges ){ 
				if(edges - 1 == histnb ) return 0; //perfect
				if(mis_count < 1) return 1;
				std::cout<<mis_count<<std::endl;
				return 2; //rebinnable; checker is less than the NbinsX();
			}
			
			bool lb = cur_binning.front() < cv_hist->GetBinLowEdge(1);
			bool hb = cur_binning.back() > cv_hist->GetBinLowEdge(1+ histnb);
			if(lb||hb) return 3;//underflow bins; 

			return 4;
}


/*
 * return a covariance matrix whose diagonal elements are 
 * the systematic errors;
 * basically, it is (total_fractional matrices_ii)*(MChist_i)
 */
TMatrixD gadget_PrepareMatrix(std::vector<bdt_sys*> syss, TFile* matrix_root , TH1* MChist, double outPOT, std::string bdt_tag){
	//Extract Error matrix from the fractionla matrix and adjust optical Model;
	//Only the diagonal element matters.
	
	bool message = true;
	bool adjustOM = false;

	int nb = MChist->GetNbinsX();
	TMatrixD total_fm(nb,nb);
	std::cout<<"\n Getting fracitonal matrices: ";
	for(size_t index = 0; index < syss.size(); ++index){//loop over systematics

		if(bdt_tag.compare(syss[index]->bdt_file::tag) != 0 ) continue;
		if(syss[index]->its_CV) continue;//work with non-CV systematic weights, bdt_sys only provides names to retrieve contents;
		bdt_sys* cur_sys = syss[index];

		std::vector<std::stringstream> summary(nb+1);//print out buffer
		std::vector<bool> summary_headers(nb+1, true);
		
		for(size_t jndex = 0; jndex < (cur_sys->vars).size(); ++jndex){//loop over spcific type of SW 

			TString temp_tag = cur_sys->tag+"_"+cur_sys->vars_name[jndex];

			TMatrixD temp_matrix(nb,nb); 
			TMatrixD temp_CV(nb,1); 

			TMatrixD* read_matrix = (TMatrixD*) matrix_root->Get(temp_tag + "_FracMatrix");
			TMatrixD* read_CV = (TMatrixD*) matrix_root->Get(temp_tag + "_CV");

			if(read_matrix == NULL){//it is skipped because it has been merged to another matrices.
				std::cout<<"\tSkip loading matrix for "<<temp_tag+ "_FracMatrix"<<std::endl;
				continue;
			}

			if(read_CV->GetNrows() != nb){
				std::cout<<"Pls check input matrix: cov. has dimension "<<read_CV->GetNrows()<<" vs histogram "<<nb<<std::endl;
				exit(EXIT_FAILURE);
			}


			TH1D* temp_CVhist  = (TH1D*) matrix_root->Get(temp_tag + "_CV_drawn");
			if(true){//check bins;
				std::vector< double > MCbinning;
				for(int kndex = 1; kndex < nb+2; kndex ++){
					MCbinning.push_back(MChist->GetBinLowEdge(kndex) );
				}

				int matching_code = gadget_BinMatcher(temp_CVhist, MCbinning);
				/*0 - no problem; (1,2,4) - (1,2,4)
				 *1 - Rebinable: temp_CVhist out of range (1,2,3,4) - (2,3)
				 *2 - Rebinable: (1,2,4) - (1,4)
				 *3 - BAD: underflow bins: xoutput_binning range too large (1,2,4) - (0,1,4)
				 *4 - BAD: bins edge not found in temp_CVhist (1,2,5) - (1,2,3,4)
				 */
				switch (matching_code){
					case 0:
						temp_matrix = *read_matrix;
						temp_CV = *read_CV;
						std::cout<<"\rCovariance matrix binnings are all good.";
						break;
					case 1:
						{//take out temp_matrix & temp_CV to match MChist binning;
							std::cout<<"Covariance matrix is partially useful, pick bin edges started from ";
							int initial_bin = 0;
							for(int lndex = 1; lndex < temp_CVhist->GetNbinsX()+1; lndex++){
								if(abs(temp_CVhist->GetBinLowEdge(lndex) - MCbinning.front()) < 10e-5){//got it 
									initial_bin = lndex - 1;
									std::cout<<"bin "<<lndex<<std::endl;
									break;
								} 
							}
							for(int lndex = 0; lndex < nb; lndex ++){
								temp_CV(lndex,0) = (*read_CV)(initial_bin + lndex ,0);
								for(int mndex = 0; mndex < nb; mndex ++){
								temp_matrix(lndex,mndex) = (*read_matrix)(initial_bin+lndex,initial_bin+mndex);
								}
							}
						}
						break;
					default://>1;
						std::cout<<"BinMatcher Code "<< matching_code<<std::endl;
						std::cout<<"Binnings from matrices:(" <<temp_CVhist->GetBinLowEdge(1)<<",";
						std::cout<<temp_CVhist->GetBinLowEdge(1+ temp_CVhist->GetNbinsX())<<") vs MC (";
						std::cout<<"("<<MChist->GetBinLowEdge(1)<<","<<MChist->GetBinLowEdge(nb+1);
						std::cout<<"). Dont match, need to generate the matrices."<<std::endl;
						exit(EXIT_FAILURE);
						break;	
				}
			}

			if(cur_sys->its_OM && adjustOM ){//modify the stat part of the matrix for Optical Model
				std::cout<<"\tAdjust Optical Model statistical error in the fractional matrix"<<std::endl;
				for(int kndex =  0; kndex < nb; ++kndex){
					double temp_mii = temp_matrix(kndex,kndex);//[kndex*nb+kndex];
					double origin_cvi = temp_CV(kndex,0)/outPOT*cur_sys->pot;//which is OM weights POT
					double MCi = MChist->GetBinContent(kndex+1);
					std::cout<<"\tThe ("<<kndex<<","<<kndex<<") element: "<<temp_mii<<" -> ";

					temp_matrix(kndex,kndex) = ((temp_mii*pow(origin_cvi,2)-origin_cvi)*pow(MCi/origin_cvi,2)+MCi)/pow(MCi,2);
					std::cout<<temp_matrix(kndex,kndex)<<std::endl;
				}//next bin
			}
//			temp_matrix.SetMatrixArray(nums.GetArray());
			total_fm+= temp_matrix;
			for(int kndex = 0; kndex < nb; ++kndex){
				double temp_error = sqrt(temp_matrix(kndex,kndex))*MChist->GetBinContent(kndex+1);
				if(summary_headers[kndex+1]){
					if(kndex == 0) summary[kndex]<<"Bin ";
					summary[kndex+1]<<std::setw(4)<<std::left<<kndex+1;
					summary_headers[kndex+1] = false;
				}
				summary[kndex+1]<<std::setw(9)<<temp_error<<"("<<std::setw(11)<<std::left<<temp_matrix(kndex,kndex)<<") ";
			}//next bin
				summary[0]<<std::setw(23)<<std::left<<cur_sys->vars_name[jndex]+"_Err (fm) ";
		}//next SW

		std::cout<<" "<<cur_sys->systag<<std::endl;
		std::streambuf* orig_buf = std::cout.rdbuf();//save the original buffer,
		for(int kndex = 0; kndex < summary.size(); ++kndex){
			std::cout<<summary[kndex].rdbuf()<<"\n";
		}
		std::cout.rdbuf(orig_buf);//swap back after finish;
	}//next systematic

	for(int index = 0; index < nb; ++index){
		std::cout<<"Total Fractional Error Matrix ("<<index<<","<<index<<"): "<<(total_fm)(index,index);
		std::cout<<" CV:"<<MChist->GetBinContent(index+1)<<" Frac Err:"<<sqrt((total_fm)(index,index))*MChist->GetBinContent(index+1)<<std::endl;
	}
	std::cout<<std::endl;
	return total_fm;

}


/*
 * This function separate the matrix cov into shape, mixed, and norm three components.
 * output: vector<TH2D> {shape, mixed, norm}
 */
std::vector<TMatrixD> gadget_SeparateMatrix(TMatrixD* frac_cov, TH1* hist, TString label){
	
	bool message = true;
	bool db_message = false;

	std::cout<<"Handle and Separate the fractional covaraince matrix."<<std::endl;
	int nb = hist->GetNbinsX();
	double nl = hist->GetBinLowEdge(1);
	double nh = hist->GetBinLowEdge(nb+1);

	double total_event = hist->Integral();
	std::cout<<"Total event "<<total_event<<std::endl;

	TH2D* shape = new TH2D("shape", "shape", nb, nl, nh, nb, nl, nh);
	TH2D* mixed = new TH2D("mixed", "mixed", nb, nl, nh, nb, nl, nh);
	TH2D* norm = new TH2D("norm", "norm", nb, nl, nh, nb, nl, nh);

	TMatrixD shapeM(nb,nb);
	TMatrixD mixedM(nb,nb);
	TMatrixD normM(nb,nb);

	std::vector<std::stringstream> pshapeM(nb);
	std::vector<std::stringstream> pmixedM(nb);
	std::vector<std::stringstream> pnormM(nb);
	for(int index=0; index< nb;index++){
		for(int jndex=0; jndex< nb;jndex++){
			double bini = hist->GetBinContent(index+1);
			double binj = hist->GetBinContent(jndex+1);
			double fm_ij = bini*binj*(*frac_cov)(index,jndex);//-hist->GetBinError(index+1)*hist->GetBinError(jndex+1);//CHECK, statErr removed

			double msum_ik = 0;
			double msum_kj = 0;
			double msum_kl = 0;
			for(int kndex=0; kndex< nb;kndex++){
				msum_ik += hist->GetBinContent(index+1)*hist->GetBinContent(kndex+1)*(*frac_cov)(index,kndex);
				msum_kj += hist->GetBinContent(kndex+1)*hist->GetBinContent(jndex+1)*(*frac_cov)(kndex,jndex);
				for(int lndex=0; lndex< nb;lndex++){

					msum_kl += hist->GetBinContent(kndex+1)*hist->GetBinContent(lndex+1)*(*frac_cov)(kndex,lndex);
				}
			}
			double norm_ij = bini*binj/pow(total_event,2)*msum_kl;
			double mixed_ij = binj/total_event*msum_ik + bini/total_event*msum_kj - 2*norm_ij;

			double shape_ij = fm_ij - mixed_ij - norm_ij;
			if(db_message){
				std::cout<<"("<<index+1<<","<<jndex+1<<") fm_ij:"<<fm_ij;
				std::cout<<" bini:"<<bini<<" binj:"<<binj;
				std::cout<<" msum_ik:"<<msum_ik;
				std::cout<<" msum_kj:"<<msum_kj;
				std::cout<<" msum_kl:"<<msum_kl;
				std::cout<<" shape_ij:"<<shape_ij;
				std::cout<<" mixed_ij:"<<mixed_ij;
				std::cout<<" norm_ij:"<<norm_ij<<std::endl;;
			}

			shape->SetBinContent(index+1,jndex+1, shape_ij);
			mixed->SetBinContent(index+1,jndex+1, mixed_ij);
			norm->SetBinContent (index+1,jndex+1, norm_ij);
			shapeM(index,jndex) = shape_ij;
			mixedM(index,jndex) =  mixed_ij;
			normM(index,jndex) =  norm_ij;

			pshapeM[index]<<std::setw(12)<<shape_ij<<" ";
			pmixedM[index]<<std::setw(12)<<mixed_ij<<" ";
			pnormM[index]<<std::setw(12)<<norm_ij<<" ";
		}
	}
	if(message){
		std::cout<<"Note: (1,1) is the top left "<<std::endl;
		std::cout<<"Shape only"<<std::endl;
		for(int printdex = 0; printdex < nb; ++printdex){
		std::cout<<pshapeM[printdex].rdbuf()<<std::endl;
		}
		std::cout<<"Mixed only"<<std::endl;
		for(int printdex = 0; printdex < nb; ++printdex){
			std::cout<<pmixedM[printdex].rdbuf()<<std::endl;
		}
		std::cout<<"Norm only"<<std::endl;
		for(int printdex = 0; printdex < nb; ++printdex){
			std::cout<<pnormM[printdex].rdbuf()<<std::endl;
		}

	}
	
	if(false){//print out?
	TCanvas *canvas_out = new TCanvas("tmp_3matrices","tmp_3matrices",1800,1600);
	shape->SetStats(false);
	shape->Draw("COLZ");
	canvas_out->SaveAs((label+"m_shape.pdf"),"pdf");
	canvas_out->Clear();

	mixed->SetStats(false);
	mixed->Draw("COLZ");
	canvas_out->SaveAs((label+"m_mixed.pdf"),"pdf");
	canvas_out->Clear();

	norm->SetStats(false);
	norm->Draw("COLZ");
	canvas_out->SaveAs((label+"m_norm.pdf"),"pdf");
	canvas_out->Delete();
	}

	if(true){	
		TMatrixD test(nb,nb);
		TVectorD eigenv(nb);
		test = shapeM+mixedM+normM;
		test.EigenVectors(eigenv);

		std::cout<<"Check eigen values:";
		for(int index=0; index< nb;index++){
			std::cout<<eigenv(index)<<" ";
		}
		std::cout<<"\n"<<std::endl;
	}

	return {shapeM, mixedM, normM};
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
TString sys_env::getFileName(TString plot_type, std::vector<bdt_variable> vars){
	if(vars.size()<2){
		TString bin_gap = to_string_prec((vars[0].plot_max-vars[0].plot_min)/vars[0].int_n_bins,0);

		//<type>_<vars[0]>_bgap<gap>_s<stage>_<hash>
		TString file_name = plot_type+"_"+vars[0].safe_name + "_bgap"+ bin_gap + "_s" + to_string_prec(this->fstage,0) + "_"+this->fcut_hash;

		return file_name;
	}else if (vars.size()<3){
		TString bin_gap = to_string_prec((vars[0].plot_max-vars[0].plot_min)/vars[0].int_n_bins,0) + "_" + to_string_prec((vars[1].plot_max-vars[1].plot_min)/vars[1].int_n_bins,1);;

		//<type>_<vars[0]>_bgap<gap>_s<stage>_<hash>
		TString file_name = plot_type+"_"+vars[0].safe_name + "_"+ vars[1].safe_name+ "_bgap"+ bin_gap + "_s" + to_string_prec(this->fstage,0) + "_"+this->fcut_hash;

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
//		twodhists.resize(num_vars);
		histNames.resize(num_vars);
//		twodhistNames.resize(num_vars);
		TdirNames.resize(num_vars);

		histloaded.assign(num_vars, false);
		for(int index = 0 ; index < num_vars; ++index){//loop through weights
			TdirNames[index] = tag+"_"+vars_name[index];

			for(int jndex = 0; jndex < throws; ++jndex){//save names;
				TString hist_name = vars_name[index];
				if(its_multithrows) hist_name += std::to_string(jndex);
				histNames[index].push_back(hist_name);
//				twodhistNames[index].push_back(hist_name);
			}//next throw
		}//next weight

	};



// Check if TH1D exist; return true when all required TH1D are found.
//record - true -> save the TH1D into a vector;
bool bdt_sys::Load1dhist(TFile* cur_file, bool savehist){
	
	//bdt_sys has boolean: loaded & fullyloaded 
	//it is loaded, unless NULL TH1D is found.
	//fullyloaded is false, if at least one NULL is found;
	this->fullyloaded = true;
	for(int index = 0; index< this->num_vars; ++index){//loop through weights;
//		if(this->histloaded[index]) continue;//Cant skip this, if file close, all hists are gone.
		this->histloaded[index] = true;

		(this->hists[index]).clear();//make sure it is empty before filling in;
		for(int jndex = 0; jndex< this->throws; ++jndex){//loop through throws;
			TString nth_label = (this->its_multithrows)? std::to_string(jndex) : "";
			TString hist_name = this->vars_name[index]+nth_label;

			TH1D* temp_th1 =  (TH1D*) cur_file->Get(this->TdirNames[index]+"/"+hist_name);

			if(temp_th1 == NULL){//oops, find a histogram is missing; label this bdt_sys as unloaded with histograms;
				if(verbose) std::cout<<"\t"<<this->TdirNames[index]+"/"+hist_name<<" does not exist. Next"<<std::endl;
				this->histloaded[index] = false;
				this->fullyloaded = false;
				break;
			}

			if(savehist){ 
				(this->hists[index]).push_back( (TH1D*) temp_th1->Clone() );// load it!
				if(verbose) std::cout<<"\r\tLoad Tdir: "<<this->TdirNames[index]<<" histogram: "<<hist_name;
			}
		}//next throw
		if(verbose&&savehist) std::cout<<std::endl;
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
void bdt_sys::Make1dhist(std::vector<bdt_variable> vars, TFile* out_root, bool twovars){

	bdt_variable var = vars[0];

	if(verbose) std::cout<<"\nMaking systematic histograms for "<<this->tag<<" with "<<this->num_vars<<" systematic variable at stage "<<this->fstage<<std::endl;

	//configure bin shifting parameter for 2 variables;
	//1. prelong hist
	//2. edit label
	int twod_histext = 1;
	//	TString twodlabel = "";
	if(twovars){
		twod_histext = vars[1].int_n_bins;
		//		twodlabel = "2d"+vars[0].safe_name+"_"+vars[1].safe_name;
	}

	//configure what type of histograms to make;
	int nbins =  var.int_n_bins;
	int lowbinedge = var.plot_min;
	int highbinedge = var.plot_max;

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
	for(int index = 0 ; index < this->num_vars; ++index){//loop through weights

		if(this->histloaded[index]) continue;//well, this means no need make histograms;
		out_root->cd();
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

		Tdir->cd();

		//STEP 2.2 prepare TH1D
		Long64_t nentries = temptree->GetEntries();
		if(verbose) std::cout<<" We have entries: "<<nentries<<"   ";

		//STEP 2.2.1 prepare empty histgrams to be filled, each histogram with corresponding weight[kndex] (Note, not POT reweighting).
		std::vector<TH1D*> hist1d(num_throws);//hist
		std::vector<TTreeFormula*> wgtforms(num_throws);//variable
		std::vector<std::vector<TTreeFormula*>> varformula(twod_histext, std::vector<TTreeFormula*>(num_throws));//weight, but not with POT reweighting;
		std::vector<std::vector<TTreeFormula*>> varcondition(twod_histext, std::vector<TTreeFormula*>(num_throws));//condition {type1, type2{throw1,throw2,...}, ...};

		for(int kndex = 0; kndex < twod_histext; ++kndex){//loop over conditions to prepare TTreeFormula, i.e. varx[?]
			int count_hist = 0;
			for(int lndex = 0; lndex <num_throws;++lndex){//loop throws of each weight;
				//histogram of a throw

				if(kndex==0){
					hist1d[lndex] = new TH1D(this->histNames[index][lndex],this->histNames[index][lndex], nbins*(twod_histext), lowbinedge, highbinedge+(highbinedge-lowbinedge)*(twod_histext-1));//Will be rebinned later for variable binning;
				}
				count_hist++;

				//variables of a throw
				TString cur_formula = var.name;
				TString cur_condition = "1";
				if(twovars){
					double condition_step = (vars[1].plot_max-vars[1].plot_min)/(vars[1].int_n_bins*1.0);
					double con_min = vars[1].plot_min+condition_step*kndex;
					double con_max = vars[1].plot_min+condition_step*(kndex+1);
					cur_condition = "("+vars[1].name+">"+to_string_prec(con_min,5)+"&&"+vars[1].name+"<"+to_string_prec(con_max,5)+")";

					TString con_adj = "+"+to_string_prec(kndex,0)+"*"+to_string_prec(var.plot_max-var.plot_min,5);
					cur_formula = "("+var.name+con_adj+")";
				}
				//as a convention, kndex for different constrains, and lndex for different throws;
				varformula[kndex][lndex] = new TTreeFormula(var.unit.c_str(), gadget_updateindex(cur_formula, lndex), temptree);//get the x-axis

				varcondition[kndex][lndex] = new TTreeFormula("tmp_condition", gadget_updateindex(cur_condition, lndex), temptree);
				varcondition[kndex][lndex]->GetNdata();

				if(lndex==0) std::cout<<"\nVariable is "<<gadget_updateindex(cur_formula, lndex)<<" with condition "<<gadget_updateindex(cur_condition, lndex);

				varformula[kndex][lndex]->GetNdata();

				//weights of a throw
				//			TString temp_label = (lndex>0)? "["+std::to_string(lndex)+"]":"";
				TString temp_wgt = gadget_updateindex(this->vars[index], lndex);//this->vars[index] + temp_label;
				TString temp_wgtname = this->vars_name[index]+std::to_string(lndex);

				wgtforms[lndex] = new TTreeFormula(temp_wgtname, temp_wgt+"*"+this->getStageCutsIndex(this->fstage,fbdt_cuts, lndex), temptree);//get weight with precuts;
				wgtforms[lndex]->GetNdata();
				//			std::cout<<temp_wgt+"*"+this->getStageCutsIndex(this->stage,fbdt_cuts, lndex)<<std::endl;
				//					std::cout<< wgtforms[lndex]->EvalInstance()<<std::endl;
			}//next throw
		}//next condition
		std::cout<<std::endl;

		//			if(verbose) std::cout<<"Creating "<<count_hist<<" histograms "<<std::endl;

		//STEP 2.2.2 fill histograms!

		for(Long64_t lentry = 0; lentry< nentries; ++lentry){

			temptree->GetEntry(lentry);

			for(int mndex = 0; mndex < num_throws; ++mndex){//loop through throws;
				for(int kndex = 0; kndex < twod_histext; ++kndex){

					double temp_xvalue = varformula[kndex][mndex]->EvalInstance();
					double temp_weight = (wgtforms[mndex]->EvalInstance()>29)? 0 : wgtforms[mndex]->EvalInstance();//No large weight..
					bool temp_condition = varcondition[kndex][mndex]->EvalInstance();
					//					std::cout<<"CHECK! "<<mndex<<" throws "<<wgtforms[mndex]->EvalInstance()<<" "<<this->pot<<std::endl;
					if(temp_weight>0 && verbose){
						std::cout<<"\r>>Looping entries: "<<std::setw(5)<<lentry+1<<"/"<<std::setw(5)<<nentries;
						if(debug_verbose && temp_condition>0){//print out
							std::cout<<" "<<std::setw(10)<<gadget_updateindex(var.name, mndex);
							std::cout<<"="<<std::setw(10)<<temp_xvalue;
							std::cout<<" "<<std::setw(10)<<gadget_updateindex(this->vars[index],mndex);
							std::cout<<"="<<std::setw(10)<<temp_weight;
							std::cout<<std::setw(4)<<mndex+1<<" throw of SW: "<<this->vars_name[index];
							std::cout<<" in? "<<temp_condition;
						}
						//std::cout<<std::endl;
					}

					if(!temp_condition) continue;
					hist1d[mndex]->Fill(temp_xvalue,temp_weight);//index - nth sw; mndex - nth throw
				}

			}//next throw;
		}//next entry

		if(verbose) std::cout<<"\n Finish filling "<<num_swvars<<" systematic weights."<<std::endl;

		//write hisograms to the diretory
		for(int jndex = 0; jndex <num_throws;++jndex){//throw in each weight;
			//configure th1
			//			TString hist_name = (its_multithrows)? this->vars_name[index]+std::to_string(jndex+1) :this->vars_name[index];
			hist1d[jndex]->SetStats(0);
			hist1d[jndex]->GetXaxis()->SetTitle(var.unit.c_str());
			hist1d[jndex]->GetYaxis()->SetTitle("Events");

			//save the th1 to the root file, class bdt_sys;
			Tdir->cd();
			hist1d[jndex]->Write(this->histNames[index][jndex]);

			//				if(twod_histext-1 == kndex) hist1d[jndex]->Delete();//remove hist1d when leaving the condition loop
		}

		if(verbose){
			std::cout<<"Save "<<num_throws;
			std::cout<<" histograms to the dir: "<<Tdir_name<<std::endl;
		}
	}//next systematic weight 
	infile->Close();

}

/*
 * covariance matrix from 2d histograms;
 *
 */
void sys_env::hist2d2cov( std::vector<bdt_variable> vars, bool rescale, bool smooth_matrix){

	bool checkbins = true;//going to reject non-rebinnable runs;
	bool force_rebin = false;//test different smoothing; true - bingap 25MeV->50MeV

//	TString twodlabel = "2d"+vars[0].safe_name+"_"+vars[1].safe_name;

	if(!rescale && verbose) std::cout<<"\n Warning: raw histograms are used, no normalization"<<std::endl;
	//STEP 0 Configuration

	int stage = fstage;
	int xnb = vars[0].n_bins;//n_bins is the target binning;
	bool xdo_rebin = vars[0].is_custombin;//this will be updated, if the sys histogram binning can be used.

	std::vector<double> xoutput_binning = gadget_CalBinning(vars);

	int nb = xoutput_binning.size()-1;//totalbins on xaxis;
	int ynb = (vars.size()<2)? 1 : vars[1].n_bins;//n_bins is the target binning;

	//the following are the same for one var or two vars;
	if(verbose) std::cout<<"\n\nMaking covariance matices."<<std::endl;
	//tags are already prepared as: tag_collection, tag2SWmap, tag2SWmap;
	//STEP 1 prepare root file to store  covariance matrix and configure TH2D*;

	TFile *cov_root = (TFile*) TFile::Open(root_dir+cov_prefix+".root","RECREATE"); //Collection of resulting covariance matrices;

	TString final_covroot_name = top_dir + final_prefix +".root";
	TFile *finalcov_root = (TFile*) TFile::Open(final_covroot_name,"RECREATE");//one final covariance matrix output for plotting;

	TCanvas *drawhistCanvas = new TCanvas("drawhistCanvas", "", 1200, 400);//Cavans to draw;
	TCanvas *drawCanvas = new TCanvas("drawCanvas", "", 800, 400);//Cavans to draw;

	//STEP 2 Get histograms, rebin (if is_custombin is true) and rescale with POT; then make covariance matrices;
	//	for(size_t index = 0; index < tag_collection.size(); ++index){}//each tag means one group of covariance matices;
	for(auto cur_tag : tag_collection){//each tag means one group of covariance matices;

		if(verbose) std::cout<<"\n-->Working on tag "<<cur_tag<<std::endl;

		//STEP 2.1 Find the bdt_sys syss that gives the CV; then find and draw up SW on the flight;
		bdt_sys* tempsCV = tag2CVmap.find(cur_tag)->second;//bdt_sys that gives CV - tempsCV


		//Got the CV! 
		TH1D* cv_hist = (TH1D*) (tempsCV->hists[0][0])->Clone();
		if(force_rebin) cv_hist->Rebin(2);

		if(checkbins){//check syss's hisogram bin edges; it is good, if all xoutput_binning edges are found.

			checkbins = false;
			int matching_code = gadget_BinMatcher(cv_hist, xoutput_binning);
			/*0 - no problem; (1,2,4) - (1,2,4)
			 *1 - Rebinable: cv_hist out of range (1,2,3,4) - (2,3)
			 *2 - Rebinable: (1,2,4) - (1,4)
			 *3 - BAD: underflow bins: xoutput_binning range too large (1,2,4) - (0,1,4)
			 *4 - BAD: bins edge not found in cv_hist (1,2,5) - (1,2,3,4)
			 */
			if(matching_code > 2){	
				std::cout<<"BinMatcher Code "<<matching_code<<" plothist bins "<<xoutput_binning.size()-1<<" cv bins "<<cv_hist->GetNbinsX()<<std::endl;
				std::cout<<"Please check the input the binning of hist_*.root."<<std::endl; 
				exit(EXIT_FAILURE);
			} else if( matching_code > 0){
				if(verbose) std::cout<<"Histograms are going be rebinned: "<<cv_hist->GetNbinsX()<<" -> "<<nb-1<<std::endl;
				xdo_rebin = true;
			}
		}

		//Update cv_hist, if the binning are not exactly what we want;
		if(rescale){ 
			cv_hist->Scale(out_POT/tempsCV->pot);//scale hist to data POT
			if(debug_verbose) std::cout<<"\rScale CV "<<cur_tag<<" POT: "<<out_POT/tempsCV->pot<<std::endl;
		}
		TH1D* smoothRef_hist = (TH1D*) cv_hist->Clone();
		if(xdo_rebin) *cv_hist = gadget_vrebin(cv_hist, xoutput_binning);


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
						if(debug_verbose) std::cout<<"\rScale "<<cur_tag<<temp_varName<<" POT: "<<tmp_pscale<<" throws:"<<tmp_tscale;
						
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

			TString temp_sw_name = cur_tag +"_"+ hists_names[index];//+twodlabel;

			bool do_smooth = ( (vars[0].name).compare(3,5,"EnuQE")==0 )&& smooth_matrix && ((tempsCV->systag).Contains("Unisim"));
//		if(cur_tag.Contains("Optical")) do_smooth = true;

			int plotbins = xnb*ynb;
			TH2D* all_hist = new TH2D(temp_sw_name, temp_sw_name, plotbins, &(xoutput_binning).front(), nby, 0, nby);//to store many throws SW
			TH2D* covmatrices =  new TH2D(temp_sw_name+"_CovarianceMatrix", temp_sw_name+"_CovarainceMatrix",plotbins,&(xoutput_binning).front(), plotbins,&(xoutput_binning).front());//to store covariance matrix, equal width;
			TH2D* cormatrices =  new TH2D(temp_sw_name+"_CorrelationMatrix", temp_sw_name+"_CorrelationMatrix",plotbins,&(xoutput_binning).front(), plotbins,&(xoutput_binning).front());//to store covariance matrix, equal width;

			int hist_counter = 1;
			for(auto cur_hist : sw_hists_coll[index]){//ok, take out SW hist and add them to all_hists, covmatrices;

				if(verbose) std::cout<<std::flush<<"\r Processing hists "+temp_sw_name<<hist_counter++;
				
				if(force_rebin) cur_hist->Rebin(2);

				bool trad_bin = false;//for R<5 selection.
				if(do_smooth&& trad_bin){//smooth first
					if(verbose) std::cout<<"\r Traditionally Smoothing a sw histogram";
					*cur_hist = SmoothSW(cur_hist, smoothRef_hist, trad_bin);
				}

				if(xdo_rebin) *cur_hist = gadget_vrebin(cur_hist, xoutput_binning);

				if(do_smooth && !trad_bin){//smooth first
					if(verbose) std::cout<<"\r Smoothing a sw histogram";
					*cur_hist = SmoothSW(cur_hist, cv_hist, trad_bin);
				}
				for(int jndex = 1; jndex < nb+1;  ++ jndex){
					all_hist->Fill(cur_hist->GetBinLowEdge(jndex) , cur_hist->GetBinContent(jndex));
					if(debug_verbose){
						std::cout<<"\r Fill bin "<<cur_hist->GetBinLowEdge(jndex)<<" with content "<<cur_hist->GetBinContent(jndex);
					}
				}
//				if(debug_verbose)std::cout<<std::endl;

				//std::cout<<cur_hist->GetNbinsX()<<std::endl;
				TH2D* cov_temp = Make2VarCov(cur_hist, cv_hist, xnb);
				TH2D* cor_temp = MakeCor(cov_temp, cur_hist, cv_hist);
				double extra_factor = 1;
				//Special! By doing so, we can blow upthe KpProd;
				if(temp_sw_name.Contains("KpProd")){ 
					if(verbose) std::cout<<" with an extra factor 1.5*1.5";
					extra_factor = 1.5*1.5;
				}

				covmatrices->Add(cov_temp, hist2tscalemap[cur_hist]*extra_factor);//throw rescaling;
				cormatrices->Add(cor_temp, hist2tscalemap[cur_hist]);//throw rescaling;


			}
			if(verbose) std::cout<<std::endl;
			//finish handlig sw_hist;
			TH2D* fracCov = MakeFracCov( covmatrices, cv_hist);


			//Now we can draw all_hist, fracCov, Cov;
			gStyle->SetPalette(kRainBow);//this change colz colors;
			gSystem->RedirectOutput("/dev/null");
			//all_hist 
			drawhistCanvas->Clear();
			drawhistCanvas->cd();
			all_hist->SetStats(false);
			all_hist->GetXaxis()->SetTitle(vars[0].unit.c_str());
			all_hist->GetYaxis()->SetTitle("Event Rate");
			all_hist->Draw("COLZ");
			cv_hist->Draw("L same");	
			cv_hist->SetLineColor(6);
			cv_hist->SetLineWidth(2);
			drawhistCanvas->SaveAs( drawn_dir + hist_prefix + "_"+temp_sw_name+".pdf" ,"pdf");

			//Covariance matrix;
			drawCanvas->Clear();
			drawCanvas->cd();
			covmatrices->SetStats(false);
			covmatrices->GetXaxis()->SetTitle(vars[0].unit.c_str());
			covmatrices->GetYaxis()->SetTitle(vars[0].unit.c_str());
			covmatrices->Draw("COLZ");
			covmatrices->SetTitle(temp_sw_name + " Covaraince Matrix");
			drawCanvas->SaveAs( drawn_dir + cov_prefix + "_"+temp_sw_name+".pdf" ,"pdf");

			cov_root->cd();
			covmatrices->Write();

			//Correlation matrix;
			drawCanvas->Clear();
			cormatrices->SetStats(false);
			cormatrices->GetXaxis()->SetTitle(vars[0].unit.c_str());
			cormatrices->GetYaxis()->SetTitle(vars[0].unit.c_str());
			cormatrices->Draw("COLZ");
			cormatrices->SetTitle(temp_sw_name + " Correlation Matrix");
			drawCanvas->SaveAs( drawn_dir + cov_prefix + "_"+temp_sw_name+"_cor.pdf" ,"pdf");

//			cov_root->cd();
			cormatrices->Write();


			//Fractional Covariance matrix;
			drawCanvas->Clear();
			fracCov->SetStats(false);
			fracCov->GetXaxis()->SetTitle(vars[0].unit.c_str());
			fracCov->GetYaxis()->SetTitle(vars[0].unit.c_str());
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
					if(debug_verbose&&index==jndex) std::cout<<nums[jndex*nb+index]<<", ";
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
 * Make one covariance matrix for two variables
 */
TH2D* Make2VarCov(TH1D* hist, TH1D* cv, int one_set_bins){

	bool db_message = false;
	int totalbins=hist->GetNbinsX();
	int x_ini = 1;
	int y_ini = 1;
	//	int nl=hist->GetBinLowEdge(1);
	//	int nh=hist->GetBinLowEdge(nb+1);
	//in case of customized binnings; use vector..
	std::vector<double> bins;
	for( int index = 1; index < totalbins+2; index++){

		bins.push_back(hist->GetBinLowEdge(index));//1st bin to n+1th bin
	}

//	int counter = 0;
//	delete gROOT->FindObject("tmp");
	TH2D* covmatrix = new TH2D(("tmp"+to_string_prec(global_index++,0)).c_str() ,"tmp",totalbins, &(bins.front()),totalbins, &(bins.front()) );//binning for xnbins,xmin,xmax,ybins,ymin.ymax;
	for(int index = x_ini; index<one_set_bins+1; ++index){
//		if(db_message)	std::cout<<"Make2VarCov check: Onesetbins :"<<one_set_bins<<" total bins :"<<totalbins<<std::endl;
		for(int jndex = y_ini; jndex<one_set_bins+1; ++jndex){
			double entry = (hist->GetBinContent(index)-cv->GetBinContent(index) )*(hist->GetBinContent(jndex)-cv->GetBinContent(jndex) );
			covmatrix->SetBinContent(index,jndex,entry);
			if(db_message && index ==4 &&jndex ==4){ 
			std::cout<<"CHECK ";
			std::cout<<index<<" cv "<<cv->GetBinContent(index) << " sw "<<hist->GetBinContent(index);
			std::cout<<" "<<jndex<<" cv "<<cv->GetBinContent(jndex) << " sw "<<hist->GetBinContent(jndex);
			std::cout<<" cov:"<<entry<<std::endl;
			}

		}
		if( index == totalbins) break;
		if(index == one_set_bins){//update index for 2nd blocks;
			x_ini = one_set_bins+1;
			y_ini = one_set_bins+1;
			one_set_bins*=2;
		}
	}
	return covmatrix;
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

			if(index == jndex && jndex ==4){
			std::cout<<"CHECK"<<std::endl;
			std::cout<<index<<" cv "<<cv->GetBinContent(index) << " sw "<<hist->GetBinContent(index);
			std::cout<<" "<<jndex<<" cv "<<cv->GetBinContent(jndex) << " sw "<<hist->GetBinContent(jndex);
			std::cout<<" cov:"<<entry<<std::endl;
			}

		}
	}
	return covmatrix;
}

/*
 * Make one correlation matrix
 */
TH2D* MakeCor(TH2D* twodhist, TH1D* hist, TH1D* cv){

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
		double entry = twodhist->GetBinContent(index,jndex)/sqrt(twodhist->GetBinContent(index,index)*twodhist->GetBinContent(jndex,jndex));
//			double entry = twodhist->GetBinContent(index,jndex)/sqrt((pow(hist->GetBinError(index),2)+pow(cv->GetBinError(index),2))*(pow(hist->GetBinError(jndex),2)+pow(cv->GetBinError(jndex),2)));

			covmatrix->SetBinContent(index,jndex,entry);

		}
	}
	return covmatrix;
}



/*
 * Create a fractional covariance matrix;
 */
TH2D* MakeFracCov(TH2D* inputcov, TH1D* oldcv){
	
	bool p_message = false;
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
			double fjk = (cvj==0 || cvk ==0 )? 0 : covjk/(cvj*cvk);//not divide by nb;
			if(kndex == jndex && jndex ==4 && p_message){
			std::cout<<jndex<<" value "<<cvj;
			std::cout<<" "<<kndex<<" value "<<cvk;
			std::cout<<" cov:"<<covjk<<std::endl;
			}

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
	
	bool message = true;
	bool print_ratio = false;
	bool ec_smoothing = false;
	Int_t degree = 2;

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

		if(message) cout << "Smoothing bins "<< first_bin << "-" << Nbins << " (N=" << Nbins<< ", "<<cv->GetBinLowEdge(Nbins+1)<<") with pol "<<degree;
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
		if(message) cout<<" Chi2: "<<chi2<<" AIC = "<<AIC<<" with degree "<<degree<<std::endl;
		//the minimum AIC gives the best fitting result;
	return *sw;
}


/*
 * Main code, initialize systematics;
 * Check if the histograms exist, if yes, then can save some works.
 */

int sys_env::InitSys(std::vector<bdt_variable> vars, std::vector<bdt_sys*> syss){
	bool skip_sysEvaluation = false;//true - exit the function right away, 
	if(skip_sysEvaluation) return 0;

	bool do_cov = false; //exit the function early only when fractional covar matrix exists & and do_cov is false;
	bool do_hist = false; //true - generate hist everytime, no matter what; false -this boolean is not functioning


	//the histogram goes two ways: 1dhist - 1 variable, or 2dhist (sidebyside 1dhist) - 2 variables;
	bool do_2dhist = false;
	bdt_variable varx = vars[0];

	//settings for making cov. matrices;
	bool rescale_out_hist = true;
	bool smooth_unisims = false;

	this->checkEnv();//if triggered, need to set things up outside (i.e. hive.cxx)
	
	//set the prefix names for plots and roots storing plots:
	//they go like <type>_<var>_bgap<gap>_s<stage>_<hash>
	//names will change as do_2hist change.
	this->hist_prefix  = this->getFileName("hist",vars); 
	this->cov_prefix   = this->getFileName("cov",vars); 
	this->final_prefix = this->getFileName("finalcov",vars); 

	std::cout<<"\n ---- Initialize systematics for "<<varx.unit;
	if(vars.size()>1){ 
		do_2dhist = true;	
		std::cout<<", "<<vars[1].unit;
	}
	std::cout<<"\tSystematic files selection tag: s"<<this->fstage<<"_"<<fcut_hash<<std::endl;

	//I. Before we start, check if 1) finalcov*.root exist,  2) hist*.root exist;
	//Check 1) for the finalcov*root;
	TString xml_covroot_name = varx.covar_file.c_str();
	TString out_covroot_name = (this->top_dir+this->final_prefix + ".root");
	
	//check if it is good;
	TFile *check_cov_root = (TFile*) TFile::Open(xml_covroot_name,"READ"); //one var one hist root;

	if(verbose) std::cout<<"xml suggests covaraince matrices root: \n\t"<<xml_covroot_name<<std::endl;
	if(check_cov_root != NULL){//The job was previously done;

		if(verbose) std::cout<<"\tFound root, so skip making cov matrices."<<std::endl;
		check_cov_root->Close();

		if(!do_cov){ 
			return 0;//exit the function, not generating hists/cov. matrices;
		} else std::cout<<"(^ 3 ^) I am kidding, we still generate them in \n\t"<<out_covroot_name<<std::endl;
		
	} else if(verbose) std::cout<<"\tDid not found root. Generate them in \n\t"<<out_covroot_name<<std::endl;

	//at this point, it is going to generate covariance matrices, 
	//but still know if we need to generate 1dhistograms;
	//Check 2) if hist*.root exist

	TString histfilename = this->root_dir + this->hist_prefix + ".root";
	TFile *hist_root = (TFile*) TFile::Open(histfilename,"READ");//one var one hist root;

	if(verbose) std::cout<<"Examinating systematic 1dhistograms in: \n\t"<<histfilename<<std::endl;

	if(do_hist || hist_root == NULL){//persume no 1dhist in either cases, and recreate a new root file;
		if(hist_root!=NULL) hist_root->Close();
		hist_root = (TFile*) TFile::Open(histfilename,"RECREATE");
	} else{//ok, we have a file, just read it to see whats inside
		hist_root = (TFile*) TFile::Open(histfilename,"READ");
	}

	//II. Pair up CV&SW of all bdt_sys and mark empty bdt_sys;
	//1. check & mark
	//2. make 1dhist if bdt_sys is marked
	//3. load 1dhists;
	bool all_good = true;
	if(verbose) std::cout<<"Make maps for central value root & multisims weights root: "<<std::endl;
	for(auto cur_sys : syss){//the basic stuffs;
		cur_sys->catchupEnv(this);//set up bdt_sys environment according to the marco setting in bdt_env;
		for(int index = 0; index < cur_sys->num_vars; index ++){
			cur_sys->TdirNames[index];//+=twodlabel ;//add to vars_name[index], TdirNames[index];
			cur_sys->vars_name[index];//add to vars_name[index], TdirNames[index];
		}
		TString temptag = cur_sys->tag;

		//classify bdt_sys
		(this->tag_collection).push_back(temptag);
		if(debug_verbose) std::cout<< temptag<<" is ";

		if(cur_sys->its_CV){//save CV to map;
			tag2CVmap.insert(std::make_pair (temptag, cur_sys));
			if(debug_verbose) std::cout<<"CV!"<<std::endl;

		} else{//save SW to map
			tag2SWmap.insert(std::make_pair (temptag, cur_sys));
			if(debug_verbose) std::cout<<"multisim weight!"<<std::endl;
		}

		//mark empty Tdirectory in the 1dhist root file;
		if(!(cur_sys->Load1dhist(hist_root,false) )){ //just check dont load; all good only if every bdt_sys have all 1dhists;
			all_good = false;
			if(verbose) std::cout<<" Some histograms are missing"<<std::endl;
		}
	}
	hist_root->Close();
	
	//now we know whether or not there are empty bdt_sys (no 1dhist in the root);
	//III. take action on 1dhist;
	//some or no 1dhists exist - make 1dhist ("UPDATE");
	//all 1dhists exist - no making, just load 1dhists("READ");

	if(!all_good){//now need to Make some or all histograms; open and close root for one set of bdt_sys;
		if(verbose) std::cout<<"Updating "<<histfilename<<std::endl;
		for(auto cur_sys : syss){
			if(debug_verbose) std::cout<<" Check 1dhists of the systematics "<<cur_sys->tag<<std::endl;

			if(!cur_sys->fullyloaded){ 
				hist_root = (TFile*) TFile::Open(histfilename,"UPDATE");//reopen the file, but in UPDATE mode;
				cur_sys->Make1dhist(vars, hist_root, do_2dhist);
				cur_sys->fullyloaded = true;
				hist_root->Write();
				hist_root->Close();
			}
		}
	}

	//load
	hist_root = (TFile*) TFile::Open(histfilename,"READ");
	for(auto cur_sys : syss){//load everything;
		cur_sys->Load1dhist(hist_root,true);
	}
	std::cout<<"Finish preparing systematic histograms, and now they are in\n>>>> "<<histfilename<<std::endl;

	//IV: Overlay histograms & draw covaraicne matrix;
	//but first, clean up this->tag_collection;
	sort(this->tag_collection.begin(), this->tag_collection.end() );
	this->tag_collection.erase( unique( this->tag_collection.begin(), this->tag_collection.end()), this->tag_collection.end() );

	hist2d2cov( vars, rescale_out_hist, smooth_unisims);
	hist_root->Close();//if close, reference to hists are lost.

	std::cout<<"Finish making covariance matrices, see \n\t"<<out_covroot_name<<std::endl;
	std::cout<<"\n ----------------------------- \n"<<std::endl;
}
