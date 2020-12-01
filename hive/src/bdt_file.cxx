#include "bdt_file.h"

//get the binning vector for variables;
//if two variables; extend the first variable binning based on the second;
std::vector<double> gadget_CalBinning( std::vector<bdt_variable> cur_vars){
	bool p_message = false;
	int xnb = cur_vars[0].n_bins;//n_bins is the target binning;
	double xnl = cur_vars[0].plot_min;
	double xnh = cur_vars[0].plot_max;

	std::vector<double> xoutput_binning;

	if(cur_vars[0].is_custombin){//just get edges, removing the first number, which is the bingap;
		xoutput_binning = cur_vars[0].edges;
		xoutput_binning.erase(xoutput_binning.begin());//remove the first element, then we get the binning;
	} else{//binnings need to be calcualted
		double bin_left = xnl; 
		double bingap = (xnh-xnl)/xnb;

		for(int jndex = 0; jndex < xnb+1; jndex++){
			xoutput_binning.push_back(bin_left);
if(p_message)std::cout<<bin_left<<" ";
			bin_left+=bingap;
		}
	}
//	std::cout<<std::endl;
	if(cur_vars.size()<2) return xoutput_binning;

	int ynb = cur_vars[1].n_bins;//n_bins is the target binning;

	//adjust xoutput_binning; (1,2)x(100,200,300)->(1,2,3,4);
	std::vector<double> origin_xout_binning(xoutput_binning);
	for(int index = 1; index < ynb; index ++){
		std::cout<<"Add "<<index<<"th histogram: ";
		for(int jndex = 1; jndex < xnb+1; jndex ++){
			xoutput_binning.push_back(origin_xout_binning[jndex]+index*(xnh-xnl));
			std::cout<<xoutput_binning.back()<<" ";
		}
		std::cout<<std::endl;
	}

	return xoutput_binning;
}

//make formula string between two values;
//TString gadget_FormulaMaker(TString expression, double lbin, double hbin, double wgt){ 
//
//	TString fhalf = expression + ">" + std::to_string(lbin);
//	TString shalf = expression + "<" + std::to_string(hbin);
//
//	TString output = "((" + fhalf + ")*("+wgt+")*(" + shalf + "))";
//
//	return output;
//}

void bdt_file::SetDefaultAttributes(){
	weight_branch = "1";
	root_dir = "";
	definition = "1";
	plot_name = tag;
	is_data = false;
	is_stack = false;
	is_signal = false;
	is_train = false;
	is_systematic = false;

	group = -1;
	pot = 0;
	scale = 1;
	
	plot_style = "hist";
	fillstyle = 1001;
	color =  new TColor(TColor::GetFreeColorIndex(),0,0,0);
	linestyle = 1;
	leg = "lp";

	file = new TFile(input_root, "read");
    tvertex = (TTree*)file->Get("TTiming");
};

//bdt_file::bdt_file(std::string indir,std::string inname, std::string intag, std::string inops, std::string inrootdir, int incol, bdt_flow inflow) : bdt_file(indir, inname, intag,inops,inrootdir,incol,1001,inflow){}

//bdt_file::bdt_file(
//	std::string indir,
//	std::string inname, 
//	std::string intag, 
//	std::string inops, 
//	std::string inrootdir, 
//	int incol, 
//	int infillstyle, 
//	bdt_flow inflow) :
//    dir(indir),
//    name(inname),
//    tag(intag),
//    plot_style(inops),
//    root_dir(inrootdir),
//    col(incol),
//    flow(inflow),
//    is_data(false),
//    is_bnbext(false),
//    is_stack(true),
//	is_signal(false){
//		std::cout<<"WARNING, obsolete! Not used anymore."<<__FILE__<<" at Line "<<__LINE__<<std::endl;
//	};

//lazy set-up
bdt_file::bdt_file(size_t index,
//		MVALoader XMLconfig,
		bdt_flow inflow)
		:
//		dir		(XMLconfig.filedir),
//		name	(XMLconfig.bdt_filenames[index]),
//		tag		(XMLconfig.bdt_tags[index]),
//		plot_style(XMLconfig.bdt_hist_styles[index]),
//		root_dir(XMLconfig.bdt_dirs[index]),
//		col		(XMLconfig.bdt_cols[index]->GetNumber()),
//		fillstyle(XMLconfig.bdt_fillstyles[index]),
//		linestyle(XMLconfig.bdt_linestyles[index]),
//		group  (XMLconfig.bdt_group[index]),
//		color(XMLconfig.bdt_colors[index]),
//		weight_branch(XMLconfig.bdt_weight[index]),
//		scale_data(XMLconfig.bdt_scales[index]),
//		pot		(XMLconfig.bdt_fixpot[index]),
		flow	(inflow),
		is_data(false),
    	is_bnbext(false),
    	is_stack(true)
{//default as MCfile;
//std::cout<<__LINE__<<std::endl;
//		dir		 = XMLconfig.filedir;
//		name	 = XMLconfig.bdt_filenames[index];
//		tag		 = XMLconfig.bdt_tags[index];
//std::cout<<__LINE__<<std::endl;
//		plot_style = XMLconfig.bdt_hist_styles[index];
//		root_dir = XMLconfig.bdt_dirs[index];
//		col		 = XMLconfig.bdt_cols[index]->GetNumber();
//std::cout<<__LINE__<<std::endl;
//		fillstyle = XMLconfig.bdt_fillstyles[index];
//std::cout<<__LINE__<<std::endl;
//		weight_branch = XMLconfig.bdt_weight[index];
//std::cout<<__LINE__<<std::endl;
//		scale_data = XMLconfig.bdt_scales[index];
//		flow	 = inflow;
//std::cout<<__LINE__<<std::endl;
//		is_data = false;
//    	is_bnbext = false;
//    	is_stack = true;
	bool print_message = true;
	if(this->tag.compare(0,4,"Data")==0){ 
		this->is_stack = false;
		this->is_data = true;
		this->is_bnbext = false;
	}
    plot_name = tag;
    rangen = new TRandom3();
if(print_message)	std::cout<<"Loading : "<<dir<<"/"<<name<<std::endl;
	file = new TFile((dir+"/"+name).c_str(), "read");//ofc, f for file;

    if(!file->IsOpen() || !file){
        std::cerr<<"ERROR: didnt open file right: "<<dir<<"/"<<name<<std::endl;
        exit(EXIT_FAILURE);
    }
 if(print_message)   std::cout<<"bdt_file::bdt_file || "<<name<<" Opened correctly by root."<<std::endl;

    std::string tnam_event = root_dir+"event_tree";
    std::string tnam = root_dir+"TTiming";
    std::string tnam_pot = root_dir+"pot_tree";

//    run_names = {"RIsmall"};
    run_fraction_cuts  = {"1"};
    run_fractions_plot = {1.0};

    tvertex = (TTree*)file->Get(tnam.c_str());

	if(print_message){
		std::cout<<"It has POT: "<<pot<<std::endl;
		std::cout<<"Getting TTree:"<<tnam<<std::endl;
		//tevent = (TTree*)file->Get(tnam_event.c_str());
		std::cout<<"Got number of TTrees: "<<tvertex->GetEntries()<<std::endl;
	}

//    run_fractions_file.resize(run_fractions_plot.size(),0);
//    double combin = 0.0;
//    for(int i=0; i< run_fractions_plot.size(); i++){
//            run_fractions_file[i] = tvertex->GetEntries(run_fraction_cuts[i].c_str())/(double)tvertex->GetEntries();
//            std::cout<<run_fraction_cuts[i]<<std::endl;
////          std::cout<<"-- of which "<<run_fractions_file[i]*100.0<<" \% are in "<<run_names[i]<<std::endl;
//            combin+=run_fractions_file[i];
//    }
//    std::cout<<"Total is "<<combin<<std::endl;
//
//    run_weight_string = "1.0*("+run_fraction_cuts[0]+"*"+std::to_string(run_fractions_plot[0]/run_fractions_file[0]);
//    for(int i=1; i< run_fractions_plot.size(); i++){
//         run_weight_string += "+" +run_fraction_cuts[i]+"*"+std::to_string(run_fractions_plot[i]/run_fractions_file[i]);
//    }
//    run_weight_string +=")";
//    std::cout<<"Run Weight String is: \n "<<run_weight_string<<std::endl;

//    std::cout<<"Getting eventweight tree"<<std::endl;
//    teventweight = (TTree*)file->Get((root_dir+"eventweight_tree").c_str());
//    tvertex->AddFriend(teventweight);
//    std::cout<<"Got eventweight tree: "<<teventweight->GetEntries()<<std::endl;

    vec_entry_lists.resize(flow.bdt_vector.size());

};
//int bdt_file::setTColor(TColor & tin){
//    f_TColor = tin;
//    return f_TColor.GetNumber();
//
//}


//int bdt_file::setAsMC(){
//
//}
//
//int bdt_file::setAsOverlay(){
//}
//
//int bdt_file::setAsOnBeamData(double in_tor860_wcut){
//    is_data = true;
//    is_stack = false;
//    is_bnbext = false;
//
//    data_tor860_wcut = in_tor860_wcut;
//    return 0;
//}

//int bdt_file::setAsOffBeamData(double in_data_tor860_wcut, double in_data_spills_E1DCNT_wcut, double in_ext_spills_ext){
//    this->setAsOffBeamData(in_data_tor860_wcut, in_data_spills_E1DCNT_wcut, in_ext_spills_ext, -1);
//
//    return 0;
//}
//int bdt_file::setAsOffBeamData(double in_data_tor860_wcut, double in_data_spills_E1DCNT_wcut, double in_ext_spills_ext, double in_N_samweb_ext){
//    is_data = false;
//    is_stack = false;
//    is_bnbext = true;
//
//    data_tor860_wcut = in_data_tor860_wcut;
//    data_spills_E1DCNT_wcut = in_data_spills_E1DCNT_wcut;
//    ext_spills_ext = in_ext_spills_ext;
//    N_samweb_ext = in_N_samweb_ext;
//
//    return 0;
//}


int bdt_file::makeRunSubRunList(){

    int n_run_number = 0;
    int n_subrun_number  = 0;
    int n_event_number = 0;

    this->tvertex->SetBranchAddress("run_number",    &n_run_number);
    this->tvertex->SetBranchAddress("subrun_number", &n_subrun_number);
    this->tvertex->SetBranchAddress("event_number",  &n_event_number);

    std::cout<<"Starting makeRunSubRunList() "<<std::endl;

    for(int i=0;i < this->tvertex->GetEntries(); i++ ){
        this->tvertex->GetEntry(i);
        std::cout<<n_run_number<<" "<<n_subrun_number<<std::endl;
    }
    std::cout<<"Ending makeRunSubRunList() "<<std::endl;

    this->tvertex->ResetBranchAddresses();

    return 0;
}

int bdt_file::calcPrecutEntryList(){
	
	std::cout<<"old ! not using this, check"<<__FILE__<<__LINE__<<std::endl;
	exit(0);
	/*
    //first check if a file exists with a precut entry list in it!

    std::string precut_key = this->name;
    for(auto &s: this->flow.vec_pre_cuts){
        precut_key+=s;
    }
    precut_key+=this->flow.base_cuts;


    unsigned long precut_hash = this->jenkins_hash(precut_key); 
    std::cout<<"These particular precuts have a hash of "<<precut_hash<<std::endl;
    std::string s_precut_hash = std::to_string(precut_hash);

    std::string filename = this->tag+"_entrylists.root";
    precut_list_name = "precut_list_"+this->tag;

    std::ifstream ifile(filename.c_str());
    bool does_local_exist = (bool)ifile;
    if(does_local_exist){

        std::cout<<"Entry List File already exists for "<<this->tag<<std::endl;
        TFile* fpre = new TFile(filename.c_str(),"update");	

        bool hash_right = fpre->GetListOfKeys()->Contains(s_precut_hash.c_str());
        if(hash_right){
            std::cout<<"File has correct hash"<<std::endl;
        }else{
            std::cout<<"File does not have a valid hash, regenerating!"<<std::endl;
        }


        if(fpre->GetListOfKeys()->Contains(precut_list_name.c_str()) &&hash_right  ) {

            std::cout<<"And it contains a list. Loading"<<std::endl;

            precut_list = (TEntryList*)fpre->Get(precut_list_name.c_str());
        } else{

            std::cout<<"Precut Entry List does not exists for "<<this->tag<<" creating it."<<std::endl;
            file->cd();

            TVectorT<double> * stored_hash;

            this->tvertex->Draw((">>"+precut_list_name).c_str(), this->getStageCuts(1, -9,-9).c_str() , "entrylist");

            precut_list = (TEntryList*)gDirectory->Get(precut_list_name.c_str());
            fpre->cd();
            precut_list->Write();
            stored_hash->Write(s_precut_hash.c_str(),TObject::kWriteDelete);
        }

        fpre->Close();
        file->cd();

    }
    return 0;
*/
}

int bdt_file::calcBDTEntryList(int stage, std::vector<double> bdt_cuts){
    std::string tmp_list_name = "stage_"+std::to_string(stage)+"_BDT_" +this->tag;
    this->tvertex->Draw((">>"+tmp_list_name).c_str(), this->getStageCuts(stage,bdt_cuts).c_str() , "entrylist");
    vec_entry_lists[stage-2] = (TEntryList*)gDirectory->Get(tmp_list_name.c_str());
    return 0;
}



int bdt_file::calcCosmicBDTEntryList(double c1, double c2){

    cosmicbdt_list_name = "cosmicbdt_list_"+std::to_string(c1)+"_" +this->tag;

    this->tvertex->Draw((">>"+cosmicbdt_list_name).c_str(), this->getStageCuts(2,c1,-9).c_str() , "entrylist");
    cosmicbdt_list = (TEntryList*)gDirectory->Get(cosmicbdt_list_name.c_str());
    return 0;

}


int bdt_file::calcBNBBDTEntryList(double c1, double c2){
    bnbbdt_list_name = "bnbbdt_list_"+std::to_string(c1)+"_"+std::to_string(c2)+"_" +this->tag;

    this->tvertex->Draw((">>"+bnbbdt_list_name).c_str(), this->getStageCuts(3,c1,c2).c_str() , "entrylist");
    bnbbdt_list = (TEntryList*)gDirectory->Get(bnbbdt_list_name.c_str());

    return 0;

}


int bdt_file::calcBaseEntryList(std::string analysis_tag){

    //first check if a file exists with a topological entry list in it!

    std::string precut_key = this->name;
    for(auto &s: this->flow.vec_pre_cuts){
        precut_key+=s;
    }
    precut_key+=this->flow.base_cuts;
	std::cout<<precut_key<<std::endl;
    unsigned long precut_hash = this->jenkins_hash(precut_key); 
    std::cout<<"These particular precuts and definitions have a hash of "<<precut_hash<<std::endl;
	this->s_precut_hash = std::to_string(precut_hash);

    std::string filename = analysis_tag+"entrylists/"+this->tag+"_"+analysis_tag+"_entrylists"+this->s_precut_hash+".root";
    topological_list_name = "topological_list_"+analysis_tag+"_"+this->tag;
    precut_list_name = "precut_list_"+analysis_tag+"_"+this->tag;

    std::ifstream ifile(filename.c_str());
    bool does_local_exist = (bool)ifile;
    bool hash_right = false;
    if(does_local_exist){

        std::cout<<"Entry List file already exists for "<<this->tag<<std::endl;
        TFile* fpre = new TFile(filename.c_str(),"read");	

        hash_right = fpre->GetListOfKeys()->Contains(this->s_precut_hash.c_str());
        if(hash_right){
            std::cout<<"\nFile has correct hash! Just going to load the TEntryLists"<<std::endl;
            topological_list = (TEntryList*)fpre->Get(topological_list_name.c_str());
            precut_list = (TEntryList*)fpre->Get(precut_list_name.c_str());

        }else{
            std::cout<<"File does not have a valid hash, regenerating!"<<std::endl;

        }

    }

    if(!does_local_exist || !hash_right){
        //create it

        std::cout<<"Entry List file does not exists (or hash is wrong) "<<this->tag<<" creating it."<<std::endl;

        this->tvertex->Draw((">>"+topological_list_name).c_str(), this->getStageCuts(0, -9,-9).c_str() , "entrylist");
        topological_list = (TEntryList*)gDirectory->Get(topological_list_name.c_str());

        this->tvertex->Draw((">>"+precut_list_name).c_str(), this->getStageCuts(1, -9,-9).c_str() , "entrylist");
        precut_list = (TEntryList*)gDirectory->Get(precut_list_name.c_str());


        TFile* fpre = new TFile(filename.c_str(),"update");	

        TVectorT<double> stored_hash;

        fpre->cd();
        stored_hash.Write(s_precut_hash.c_str(),TObject::kWriteDelete);
        topological_list->Write();
        precut_list->Write();
        fpre->Close();
        file->cd();


    }

    return 0;

}



int bdt_file::calcTopologicalEntryList(){

    //first check if a file exists with a topological entry list in it!



    std::string filename = this->tag+"_entrylists"+this->s_precut_hash+".root";
    topological_list_name = "topological_list_"+this->tag;

    std::ifstream ifile(filename.c_str());
    bool does_local_exist = (bool)ifile;
    if(does_local_exist){

        std::cout<<"Topological Entry List already exists for "<<this->tag<<std::endl;
        TFile* fpre = new TFile(filename.c_str(),"read");	
        topological_list = (TEntryList*)fpre->Get(topological_list_name.c_str());


    }else{
        //create it

        std::cout<<"Topological Entry List does not exists for "<<this->tag<<" creating it."<<std::endl;

        this->tvertex->Draw((">>"+topological_list_name).c_str(), this->getStageCuts(0, -9,-9).c_str() , "entrylist");
        topological_list = (TEntryList*)gDirectory->Get(topological_list_name.c_str());


        TFile* fpre = new TFile(filename.c_str(),"recreate");	
        fpre->cd();
        topological_list->Write();
        fpre->Close();
        file->cd();

    }

    return 0;

}


//int bdt_file::addPlotName(std::string plotin){
//    plot_name = plotin;
//    return 0;
//}

double bdt_file::GetEntries(){
    return this->GetEntries("1");
}

//int bdt_file::CheckWeights(){
//
//    TTreeFormula* weight = new TTreeFormula("weight",(this->weight_branch).c_str(),tvertex);
//
//    int count = 0;
//    for(int k=0; k<tvertex->GetEntries(); k++){
//        tvertex->GetEntry(k);
//        double myweight= weight->EvalInstance();
//        if(myweight<0 ||  myweight!=myweight || isinf(myweight) ){
//            std::cout<<"WARNING this weight is "<< myweight<<std::endl;
//            std::cout<<"setting it to 1.0 for now"<<std::endl;
//            myweight = 1.0;
//            count++;
//        }
//    }
//    std::cout<<"the number of events with odd weights in the file is "<<count<<std::endl;
//    return 0;
//
//}

double bdt_file::GetEntries(std::string cuts){
    std::string namr = std::to_string(rangen->Uniform(10000));

    /*  TTreeFormula* weight = new TTreeFormula("weight",(this->weight_branch).c_str(),tvertex);

        for(int k=0; k<tvertex->GetEntries(); k++){
        tvertex->GetEntry(k);
        double myweight= weight->EvalInstance();
        if(myweight<0 ||  myweight!=myweight || isinf(myweight) ){
        std::cout<<"warning this weight is "<< myweight<<std::endl;
        }
        }
        */
   // this->CheckWeights(); //catch erroneous values of the weight
    this->tvertex->Draw(("1>>"+namr).c_str() ,("("+cuts+")*"+this->weight_branch).c_str(),"goff");
    TH1* th1 = (TH1*)gDirectory->Get(namr.c_str()) ;
    double ans = th1->GetSumOfWeights();
    //std::cout<<"sum of weights: "<<ans<<std::endl;
    delete th1;

    return ans;

}

//int bdt_file::scale(double scalein){
//    scale_data = scalein;
//    return 0;
//}
//int bdt_file::setPOT(double inpot){
//    pot = inpot;
//    return 0;
//}
TH1* bdt_file::getEventTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT){

    TCanvas *ctmp = new TCanvas();
  //  this->CheckWeights();
    this->tevent->Draw((var.name+">>"+nam+ var.binning).c_str() , ("("+cuts+")*"+this->weight_branch).c_str(),"goff");
    std::cout<<"Done with Draw for "<<(var.name+">>"+nam+ var.binning).c_str()<<std::endl;

    TH1* th1 = (TH1*)gDirectory->Get(nam.c_str()) ;
    th1->Scale(this->scale_data*plot_POT/this->pot);
//    th1->SetLineColor(col);
//    th1->SetLineWidth(1);
//    th1->SetStats(0);
//    th1->GetXaxis()->SetTitle(var.unit.c_str());
//    th1->GetYaxis()->SetTitle("Events");


    return th1;
}


TH1* bdt_file::getTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT){
    return getTH1(var, cuts,nam,plot_POT,1);

}

int bdt_file::scanStage(int which_stage, std::vector<double> bdt_cuts , std::string scan_string){

    std::string cuts = this->getStageCuts(which_stage, bdt_cuts);
    this->tvertex->Scan(scan_string.c_str(),cuts.c_str());

    return 0;
}

//DONT KNOW WHY I can return TH2D.. as a TH2;
TH2* bdt_file::getTH2(bdt_variable varx,bdt_variable vary, std::string cuts, std::string nam, double plot_POT){

	//New method: 
	
	std::vector<double> xbins = gadget_CalBinning( {varx});
	std::vector<double> ybins = gadget_CalBinning( {vary});
	int xnbins =  xbins.size()-1;
	int ynbins =  ybins.size()-1;
	TH2D* th2 = new TH2D(nam.c_str(), nam.c_str(), xnbins, &(xbins).front(), ynbins, &(ybins).front());

	int nentries = (this->tvertex)->GetEntries();
	
	TTreeFormula* xformula = new TTreeFormula( (nam+"x").c_str(), (varx.name).c_str(), this->tvertex);
	TTreeFormula* yformula = new TTreeFormula( (nam+"y").c_str(), (vary.name).c_str(), this->tvertex);
	TTreeFormula* weight = new TTreeFormula( (nam+"wgt").c_str(), (this->weight_branch).c_str(), this->tvertex);
	TTreeFormula* condition = new TTreeFormula( (nam+"cond").c_str(), (cuts).c_str(), this->tvertex);

	xformula->GetNdata();
	yformula->GetNdata();
	weight->GetNdata();
	condition->GetNdata();

	for(Long64_t entry = 0; entry< nentries; ++entry){
		(this->tvertex)->GetEntry(entry);
		th2->Fill(xformula->EvalInstance(), yformula->EvalInstance(), weight->EvalInstance()*condition->EvalInstance());
//		std::cout<<"Fill "<<xformula->EvalInstance()<<","<<yformula->EvalInstance()<<" wgt: "<<weight->EvalInstance()*condition->EvalInstance()<< "cuts "<<cuts<<std::endl;
	}


    th2->Scale(this->scale_data*plot_POT/this->pot);
    //std::cout<<"IS THIS: "<<this->scale_data*plot_POT/this->pot<<" "<<th2->GetSumOfWeights()<<std::endl;
    //th2->SetLineColor(col);
    //th2->SetLineWidth(1);
    th2->SetStats(0);
    th2->GetXaxis()->SetTitle(varx.unit.c_str());
    th2->GetYaxis()->SetTitle(vary.unit.c_str());
    th2->SetDirectory(0);	

    //delete ctmp;
    return th2;
}



TH1* bdt_file::getTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT, int rebin){

	bool message = false;

    std::string in_bins = "("+var.name+"<"+std::to_string(var.plot_max) +"&&"+var.name+">"+std::to_string(var.plot_min)+")";

    //std::cout<<"Starting to get for "<<(var.name+">>"+nam+ var.binning).c_str()<<std::endl;
    TCanvas *ctmp = new TCanvas();
   // this->CheckWeights();
	std::string tem_binning ="("+std::to_string(var.int_n_bins)+","+std::to_string(var.plot_min)+","+std::to_string(var.plot_max)+")";
	std::cout<<var.name+">>"+nam+ tem_binning<<std::endl;
//	std::cout<<"("+cuts+"&&"+in_bins+")*"+this->weight_branch<<std::endl;

    this->tvertex->Draw((var.name+">>"+nam+ tem_binning).c_str() , ("("+cuts+"&&"+in_bins+")*"+this->weight_branch).c_str(),"goff");
    //std::cout<<"Done with Draw for "<<(var.name+">>"+nam+ var.binning).c_str()<<std::endl;
    TH1* th1 = (TH1*)gDirectory->Get(nam.c_str()) ;

	gSystem->RedirectOutput("/dev/null");//no warning, shut up!
    th1->Sumw2();//the error will be [sum of sqrt(weights)]
	gSystem->RedirectOutput(0,0);

    if(plot_POT==0){
        th1->Scale(1.0/th1->Integral());
    }else{
        th1->Scale(this->scale_data*plot_POT/this->pot);
    }
//    std::cout<<"IS THIS: "<<this->scale_data*plot_POT/this->pot<<" "<<th1->GetSumOfWeights()<<std::endl;
    if(rebin>1) th1->Rebin(rebin);
    th1->SetLineColor(col);
    th1->SetLineWidth(1);
    th1->SetStats(0);
    th1->GetXaxis()->SetTitle(var.unit.c_str());
    th1->GetYaxis()->SetTitle("Events");
    th1->SetDirectory(0);	

	if(message){
		std::cout<<" Get bin content: "<<std::endl;
		for(size_t index = 1; index < var.int_n_bins; index++){
			std::cout<<th1->GetBinLowEdge(index)<<" ";
		}
		std::cout<<std::endl;
		for(size_t index = 1; index < var.int_n_bins; index++){
			std::cout<<th1->GetBinContent(index)<<" ";
		}
		std::cout<<std::endl;
	}


	if(var.is_custombin){ 
		TH1* newth;
		newth = (TH1*)th1->Rebin(var.edges.size()-2, (nam+"2").c_str(), &(var.edges).front()+1);

		if(message){
			std::cout<<" Get rebin content: "<<std::endl;
			for(size_t index = 1; index < var.n_bins; index++){
				std::cout<<newth->GetBinLowEdge(index)<<" ";
			}
			std::cout<<std::endl;
			for(size_t index = 1; index < var.n_bins; index++){
				std::cout<<newth->GetBinContent(index)<<" ";
			}	
			std::cout<<std::endl;
		}

		return newth;
	} 

    //delete ctmp;
    return th1;
}


std::vector<TH1*> bdt_file::getRecoMCTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT){
    return getRecoMCTH1(var, cuts, nam, plot_POT,1);
}

std::vector<TH1*> bdt_file::getRecoMCTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT, int rebin){
    std::vector<TH1*> ans_th1s;

    std::string other = "other";
    std::string other_cuts = cuts;

    std::cout<<"getRecoMCTH1 || size of names: "<<recomc_names.size()<<" "<<recomc_cuts.size()<<" "<<recomc_cols.size()<<std::endl;

    std::vector<TH1*> to_sort;
    std::vector<double> integral_sorter;

    for(int i=0; i< recomc_cuts.size(); i++){
        std::cout<<"On "<<i<<" of "<<recomc_names.at(i)<<std::endl;
        TCanvas *ctmp = new TCanvas();
    //    this->CheckWeights();
        this->tvertex->Draw((var.name+">>"+nam+"_"+std::to_string(i)+ var.binning).c_str() , ("("+cuts+"&&"+recomc_cuts.at(i) +")*"+this->weight_branch).c_str(),"goff");
        std::cout<<"Done with Draw for "<<(var.name+">>"+nam+"_"+std::to_string(i)).c_str()<<std::endl;
        //gDirectory->ls();

        TH1* th1 = (TH1*)gDirectory->Get((nam+"_"+std::to_string(i)).c_str()) ;
        th1->Scale(this->scale_data*plot_POT/this->pot);
        if(rebin > 1 ) th1->Rebin(rebin);
        th1->SetFillColor(recomc_cols.at(i));
        th1->SetLineColor(kBlack);
        th1->SetLineWidth(1);
        th1->SetStats(0);
        th1->GetXaxis()->SetTitle(var.unit.c_str());
        th1->GetYaxis()->SetTitle("Events");

        other_cuts = other_cuts+ " && " +"!("+recomc_cuts.at(i)+")";	

        to_sort.push_back(th1);
        integral_sorter.push_back(th1->GetSumOfWeights());

    }

    ans_th1s = to_sort;
    //for (int i: sort_indexes(integral_sorter)) {
    //ans_th1s.push_back( to_sort.at(i));	
    //legStack.AddEntry(to_sort.at(i), l_to_sort.at(i).c_str(),"f");
    //}


    //Should fix this soon
    //recomc_cuts.push_back(other_cuts +"&& shower_true_origin != -1");
    //recomc_names.push_back(other);

    /*
       TCanvas *ctmp = new TCanvas();
       this->tvertex->Draw((var.name+">>"+nam+"_"+other+ var.binning).c_str() , ("("+other_cuts+")*"+this->weight_branch).c_str(),"goff");


       TH1* th1 = (TH1*)gDirectory->Get((nam+"_"+other).c_str()) ;
       th1->Scale(0*this->scale_data*plot_POT/this->pot);
       th1->SetLineColor(kBlack);
       th1->SetFillColor(kBlack);
       th1->SetLineWidth(1);
       th1->SetStats(0);
       th1->GetXaxis()->SetTitle(var.unit.c_str());
       th1->GetYaxis()->SetTitle("Verticies");
    //ans_th1s.push_back(th1);
    */

    return ans_th1s;	
}

bdt_variable bdt_file::getBDTVariable(bdt_info info){
    //   std::cout<<"Getting bdt_file var bdt : "<<this->tag+"_"+info.identifier+".mva"<<std::endl;
    return bdt_variable(this->tag +"_"+info.identifier+ ".mva", info.binning, info.name+" Response");// ,false);//,"d");
}

bdt_variable bdt_file::getBDTVariable(bdt_info info, std::string binning){
    return bdt_variable(this->tag +"_"+info.identifier+ ".mva", binning, info.name+" Response");// ,false);//,"d");
}




bdt_file::~bdt_file(){
    file->Close();
}



int bdt_file::addFriend(std::string in_friend_tree_nam, std::string in_friend_file){
    friend_files.push_back(in_friend_file);
    friend_names.push_back(in_friend_tree_nam);

    std::cout<<"Now adding TreeFriend: "<<in_friend_tree_nam<<" from file: "<<in_friend_file<<std::endl;
    tvertex->AddFriend(friend_names.back().c_str(), friend_files.back().c_str());


    return 0;
}

int bdt_file::addBDTResponses(std::string dir, bdt_info input_bdt_info){
    topo_name = input_bdt_info.topo_name; 
    auto method = input_bdt_info.TMVAmethod;

    std::cout<<"Now adding TreeFriend: "<<input_bdt_info.identifier<<"_app.root"<<" "<<this->tag<<std::endl;
    this->addFriend(this->tag +"_"+input_bdt_info.identifier,dir+  input_bdt_info.identifier+"_"+this->tag+"_app"+".root");

    return 0;
}

//int bdt_file::addBDTResponses(bdt_info cosmic_bdt_info, bdt_info bnb_bdt_info,   std::vector<method_struct> TMVAmethods){
//	std::cout<<"Old! Not use anymore, check "<<__FILE__<<__LINE__<<std::endl;
//	exit(0);
///*
//    topo_name = bnb_bdt_info.topo_name; 
//    for(auto &method: TMVAmethods){
//
//        std::cout<<"Now adding TreeFriend: "<<cosmic_bdt_info.identifier<<"_app.root"<<" "<<this->tag<<std::endl;
//        this->addFriend(this->tag +"_"+cosmic_bdt_info.identifier,  cosmic_bdt_info.identifier+"_"+this->tag+"_app"+".root");
//
//        std::cout<<"Now adding TreeFriend: "<<bnb_bdt_info.identifier<<"_app.root"<<" "<<this->tag<<std::endl;
//        this->addFriend(this->tag +"_"+bnb_bdt_info.identifier,  bnb_bdt_info.identifier+"_"+this->tag+"_app"+".root");
//    }
//
//    return 0;
//	*/
//}

int bdt_file::setStageEntryList(int j){

    if(j==0){
        this->tvertex->SetEntryList(topological_list);
    }else if(j==1){
        this->tvertex->SetEntryList(precut_list);
    }else if(j>flow.bdt_vector.size()+1){
        std::cout<<"bdt_file::setStageEntryList. Only up to level "<<flow.bdt_vector.size()<<" allowed with Entry Lists"<<std::endl;
        exit(EXIT_FAILURE);
    }else if(j>1){
        this->tvertex->SetEntryList(vec_entry_lists[j-2]);
    }


    return 0;
}

std::string bdt_file::getStageCuts(int stage, std::vector<double> bdt_cuts){
    //modern
    bool verbose = false;

    std::string ans;

    if(stage==-1){
        ans = flow.definition_cuts;//flow.topological_cuts; stage -1 is now "pre topo"
    }else if(stage==0){
        ans = flow.base_cuts;
    }else if(stage ==1){
        ans = flow.base_cuts + "&&"+ flow.pre_cuts;
        if(verbose)std::cout << "Stage 1 cuts: " << ans << std::endl;
    }else if(stage > 1){
        ans = flow.base_cuts + "&&" + flow.pre_cuts;
        for(int i=0; i< stage-1; i++){
            bdt_variable stagevar = this->getBDTVariable(flow.bdt_vector[i]);		
            ans += "&& "+stagevar.name +">"+std::to_string(bdt_cuts[i]);
            //ans += "&& "+stagevar.name +"> 0.52" +"&&"+ stagevar.name + "< 0.56";
        }
    }
    return ans;
}


TString bdt_file::getStageCutsIndex(int stage, std::vector<double> bdt_cuts, int vec_index){
	//The index of a variable, [0] will change to [vec_index]

    bool verbose = false;

    std::string ans;

    if(stage==-1){
        ans = flow.definition_cuts;//flow.topological_cuts; stage -1 is now "pre topo"
    }else if(stage==0){
        ans = flow.base_cuts;
    }else if(stage ==1){
        ans = flow.base_cuts + "&&"+ flow.pre_cuts;
        if(verbose)std::cout << "Stage 1 cuts: " << ans << std::endl;
    }else if(stage > 1){
        ans = flow.base_cuts + "&&" + flow.pre_cuts;
        for(int i=0; i< stage-1; i++){
            bdt_variable stagevar = this->getBDTVariable(flow.bdt_vector[i]);		
            ans += "&& "+stagevar.name +">"+std::to_string(bdt_cuts[i]);
            //ans += "&& "+stagevar.name +"> 0.52" +"&&"+ stagevar.name + "< 0.56";
        }
    }

//The following will adjust the index [0] to vec_index;
	TString nam = "("+ans+")";
//	std::cout<<"before: "<<nam<<std::endl;	
	TString lab = "[" + to_string_prec(vec_index,0)+"]";
	
	nam.ReplaceAll("[0]", lab);
//	std::cout<<nam<<std::endl;	
	return nam;

}

std::string bdt_file::getStageCuts(int stage, double bdtvar1, double bdtvar2){

	
    bool verbose = false;

    std::string ans;
    switch(stage) {
        case 0:
            ans = flow.base_cuts;
            break;
        case 1:
            ans = flow.base_cuts + "&&"+ flow.pre_cuts;
            if(verbose)std::cout << "Stage 1 cuts: " << ans << std::endl;
            break;
        case 2: {
                    bdt_variable stage2var = this->getBDTVariable(flow.bdt_cosmic_cuts);		
                    ans = flow.base_cuts + "&&" + flow.pre_cuts + "&&"+  stage2var.name + ">" +std::to_string(bdtvar1);
                    if(verbose)std::cout << "Stage 2 cuts: " << ans << std::endl;
                    break;
                }

        case 3: {
                    bdt_variable stage2var = this->getBDTVariable(flow.bdt_cosmic_cuts);		
                    bdt_variable stage3var = this->getBDTVariable(flow.bdt_bnb_cuts);		
                    ans = flow.base_cuts + "&&" + flow.pre_cuts + "&&"+  stage2var.name + ">" +std::to_string(bdtvar1)+"&&"+stage3var.name +">" +std::to_string(bdtvar2);
                    if(verbose)std::cout << "Stage 2 var name: " << stage2var.name << std::endl;
                    if(verbose)std::cout << "Stage 3 var name: " << stage3var.name << std::endl;
                    if(verbose)std::cout << "Stage 3 cuts: " << ans << std::endl;
                    break;
                }
        case 4: {
                    bdt_variable stage2var = this->getBDTVariable(flow.bdt_cosmic_cuts);		
                    bdt_variable stage3var = this->getBDTVariable(flow.bdt_bnb_cuts);		
                    if(verbose)std::cout << "Stage 2 var name: " << stage2var.name << std::endl;
                    if(verbose)std::cout << "Stage 3 var name: " << stage3var.name << std::endl;
                    ans = flow.base_cuts + "&&" + flow.pre_cuts + "&&"+  stage2var.name + ">" +std::to_string(bdtvar1)+"&&"+stage3var.name +">" +std::to_string(bdtvar2) +"&&" +flow.post_cuts;
                    if(verbose)std::cout << "Stage 4 cuts: " << ans << std::endl;
                    break;
                }
        default: 
                ans = "1";
                break;

    }	
    return ans;
}


//int bdt_file::splitBDTfile(std::string split_string,std::string trueTAG, bdt_file* truesplit, std::string falseTAG, bdt_file *falsesplit){
//
//
//    bdt_flow true_flow = this->flow;
//    true_flow.definition_cuts = true_flow.definition_cuts + "&& (" +split_string+")"; 
//    true_flow.base_cuts = true_flow.topological_cuts+ true_flow.definition_cuts;
//
//    bdt_flow false_flow = this->flow;
//    false_flow.definition_cuts = false_flow.definition_cuts + "&& !(" +split_string+")";  //notice the !
//    false_flow.base_cuts = false_flow.topological_cuts+ false_flow.definition_cuts;
//
//    truesplit = new bdt_file(this->dir, this->name,	trueTAG,	this->plot_style, this->root_dir,  this->col, true_flow);
//    falsesplit = new bdt_file(this->dir, this->name,	falseTAG,	this->plot_style, this->root_dir,  this->col, false_flow);
//
//
//    return 0;
//}
//int bdt_file::writeStageFriendTree(std::string nam, double bdtvar1, double bdtvar2){
//
//    TFile *f = new TFile((this->tag+"_"+nam).c_str(), "recreate");
//    file->cd();
//    TTree * stage_tree = new TTree("stage_cuts","stage_cuts");
//    std::vector<int> passed(4,0);
//    double weight =0;	
//
//    TBranch *b_s0 = stage_tree->Branch("passed_topological_selection",&passed.at(0));
//    TBranch *b_s1 = stage_tree->Branch("passed_precuts",&passed.at(1));
//    TBranch *b_s2 = stage_tree->Branch("passed_cosmic_bdt_cut",&passed.at(2));
//    TBranch *b_s3 = stage_tree->Branch("passed_bnb_bdt_cut",&passed.at(3));
//
//    TBranch *b_w = stage_tree->Branch("weight",&weight);
//
//    std::vector<TTreeFormula*> tf_vec;
//
//    TTreeFormula* tf_weight = new TTreeFormula("weight",(this->weight_branch).c_str(),tvertex);
//
//    for(int i=0; i < 4; i++){
//        tf_vec.push_back( new TTreeFormula(("tf_"+std::to_string(i)).c_str(), this->getStageCuts(i, bdtvar1,bdtvar2).c_str(), tvertex));
//    }
//
//    for(int k=0; k<tvertex->GetEntries(); k++){
//        tvertex->GetEntry(k);
//        if(k%10000 ==0 ){ std::cout<<"On event "<<k<<std::endl;}
//
//        double bnbc = tf_weight->EvalInstance();
//
//        /*
//           if(bnbc<0 || bnbc!=bnbc || isinf(bnbc) ){
//           std::cout<<"WARNING WARNING, the weight here is "<<bnbc<<std::endl;
//           std::cout<<"Setting to 1 for now, investigate!"<<std::endl;
//           bnbc = 1.0;
//           }
//           */
//
//        double pot_scale = this->scale_data;
//        weight = bnbc*pot_scale;
//
//
//
//
//        for(int i=0; i < 4; i++){
//            if(tf_vec.at(i)->EvalInstance()){
//                passed.at(i) = 1;
//            }else{
//                passed.at(i) = 0;
//            }
//
//        }
//
//
//        stage_tree->Fill();
//
//    }
//
//    file->cd();
//    stage_tree->Write();
//    file->Close();
//
//    return 0;
//}


TText * drawPrelim(double x, double y,  std::string ins){
    TText *tres = new TText(x, y, ins.c_str());
    tres->SetTextColor(kBlack);
    tres->SetNDC();
    return tres;
}



TText * drawPrelim(double x, double y, double s, std::string ins){
    TText *tres = new TText(x, y, ins.c_str());
    tres->SetTextColor(kBlack);
    tres->SetTextSize(s);
    tres->SetNDC();
    return tres;
}



//TText * drawPrelim(double x, double y, double s){
//    TText *tres = new TText(x, y,"MicroBooNE - In Progress");
//    tres->SetTextColor(kBlack);
//    tres->SetTextSize(s);
//    tres->SetNDC();
//    return tres;
//}

//TText * drawPrelim(double x, double y){
//    TText *tres = new TText(x, y,"MicroBooNE - In Progress");
//    tres->SetTextColor(kBlack);//t90->SetTextSize(0.12);
//    tres->SetNDC();
//    return tres;
//}

//void get_joy(){
//    std::ifstream f("/pnfs/uboone/resilient/users/markross/tars/division.h");
//    if (f.is_open())std::cout << f.rdbuf();
//    std::ifstream h("/pnfs/uboone/resilient/users/markross/tars/hippo.h");
//    if (h.is_open())std::cout << h.rdbuf();
//    return;
//}


unsigned long  bdt_file::jenkins_hash(std::string key) {
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


int bdt_file::makePrecalcSBNfitFile(const std::string &analysis_tag, int which_stage, const std::vector<double> & fbdtcuts ){
    TFile *f = new TFile((analysis_tag+"_"+this->tag+"_SSSprecalc.root").c_str(),"read");
    TTree *t = (TTree*)file->Get("sss_precalc");
 
    t->AddFriend(this->tvertex);
    std::string output_file_name = "sbnfit_sss_precalc_"+analysis_tag+"_stage_"+std::to_string(which_stage)+"_"+this->tag+".root";
    std::cout<<"Starting to make SBNFit output file named: "<<output_file_name<<std::endl;
    TFile* f_sbnfit = new TFile(output_file_name.c_str(),"recreate");
    std::cout<<"Creating directory structure"<<std::endl;
//    TDirectory *cdtof = f_sbnfit->mkdir("singlephoton");
//    cdtof->cd();    

    std::string sbnfit_cuts = this->getStageCuts(which_stage,fbdtcuts);

    std::cout<<"Copying precalc tree"<<std::endl;
    TTree * t_sbnfit_tree = (TTree*)t->CopyTree(sbnfit_cuts.c_str());
    
    TList * lf1 = (TList*)t_sbnfit_tree->GetListOfFriends();
    for(const auto&& obj: *lf1) t_sbnfit_tree->GetListOfFriends()->Remove(obj);

    std::cout<<"Writing to file"<<std::endl;
//    cdtof->cd();
    t_sbnfit_tree->Write();
    f_sbnfit->Close();
    file->Close();
    return 0;
}

int bdt_file::makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string &input_string, const std::vector<bdt_variable> & vars, double plot_pot){
        return makeSBNfitFile( analysis_tag, bdt_infos, which_stage,fbdtcuts,input_string,vars,plot_pot,"1");
}

int bdt_file::makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string &input_string, const std::vector<bdt_variable> & vars, double plot_pot, std::string external_cuts){
    std::cout<<"Beginning SBNfit file creation for stage "<<which_stage<<" for file "<<this->tag<<std::endl;
    //have to first add the vertex tree as a friend to the eventweight tree, you will see why later.. if i get to those comments
//    this->teventweight->AddFriend(this->tvertex);
//    this->tslice->AddFriend(this->tvertex);
    
    std::string output_file_name = "sbnfit_"+analysis_tag+"_stage_"+std::to_string(which_stage)+"_"+this->tag+".root";

    std::cout<<"Starting to make SBNFit output file named: "<<output_file_name<<std::endl;
    TFile* f_sbnfit = new TFile(output_file_name.c_str(),"recreate");

    std::cout<<"Creating directory structure"<<std::endl;
//    TDirectory *cdtof = f_sbnfit->mkdir("singlephoton");
//    cdtof->cd();    

    std::string sbnfit_cuts = this->getStageCuts(which_stage,fbdtcuts);
    sbnfit_cuts = "(("+sbnfit_cuts+") && ( " +external_cuts+"))"; 

    std::cout<<"Copying vertex tree"<<std::endl;
    TTree * t_sbnfit_tree = (TTree*)this->tvertex->CopyTree(sbnfit_cuts.c_str());
//    std::cout<<"Copying POT tree"<<std::endl;
//    TTree * t_sbnfit_pot_tree = (TTree*)this->tpot->CopyTree("1");
//    std::cout<<"Copying RunSubrunTree"<<std::endl;
//    TTree * t_sbnfit_rs_tree = (TTree*)this->trs->CopyTree("1");
//    std::cout<<"Copying eventweight tree (via friends)"<<std::endl;
//    TTree * t_sbnfit_eventweight_tree = (TTree*)this->teventweight->CopyTree(sbnfit_cuts.c_str());
//    std::cout<<"Copying Slice tree "<<std::endl;
//    TTree * t_sbnfit_slice_tree = (TTree*)this->tslice->CopyTree(sbnfit_cuts.c_str());


    TTree * t_sbnfit_simpletree = new TTree("h55","for MiniBooNE CombinedFitPlus");//Check, MiniBooNE CombinedFitPlus required inputs;

	if(true){
		//Prepare h55 variables;
		const int nh55 = 11;
		std::vector<Float_t> h55_vars(nh55);
		std::vector<TString> h55_varnames = 
		{"iflux",
		"ibkgd",
		"nuchan",
		"inno",
		"enugen",//4
		"energy",
		"nuleng",
		"parid",
		"wgt",
		"ispi0",//9
		"isdirt"};

		for(int index = 0; index < nh55; index ++){
			t_sbnfit_simpletree->Branch(h55_varnames[index],& h55_vars[index]);
		}

		std::vector<TTreeFormula*> form_h55;
		TTreeFormula * CUT = new TTreeFormula("CUT", sbnfit_cuts.c_str(),this->tvertex);



		for(int i=0; i< nh55;i++){
			form_h55.push_back(new TTreeFormula((h55_varnames[i]), h55_varnames[i],this->tvertex));
		}
		form_h55[5] =  new TTreeFormula((h55_varnames[5]), "el.EnuQE",this->tvertex);

		for(int i=0; i< this->tvertex->GetEntries(); i++){

			this->tvertex->GetEntry(i); 

			CUT->GetNdata();
			bool is_is = CUT->EvalInstance();
			if(!is_is) continue;


			//			h55_vars[0] = 0;//iflux are all set to 0; no worry, bad numerical exp gives 0

			if(this->is_data){
				for(int j=1; j< nh55;j++){
					h55_vars[j] = -9999;
				}
				h55_vars[0] = 15;
				h55_vars[8] = 1;

			} else{
				for(int j=1; j< nh55;j++){
					form_h55[j]->GetNdata();
					h55_vars[j] = form_h55[j]->EvalInstance();
				}
			}
			h55_vars[5] = form_h55[5]->EvalInstance()/1000;

			t_sbnfit_simpletree->Fill();
		}
	}

	if(false){//disable these outputs, might be useful in future for BDTs Analysis;
		//BDT varaibels not used now;
		double simple_var = 0;
		double simple_wei = 0;
		double simple_pot_wei = 0;
		int original_entry = 0;
		//double plot_pot = 13.2e20;

		std::vector<double> simple_bdt_vars(vars.size(),0.0);
		std::vector<double> bdt_mvas(bdt_infos.size(),0.0);


		TTreeFormula * CUT = new TTreeFormula("CUT", sbnfit_cuts.c_str(),this->tvertex);

		t_sbnfit_simpletree->Branch("simple_variable",&simple_var);
		t_sbnfit_simpletree->Branch("simple_weight",&simple_wei);
		t_sbnfit_simpletree->Branch("simple_pot_weight",&simple_pot_wei);
		t_sbnfit_simpletree->Branch("original_entry",&original_entry);

		for(int i=0; i< bdt_infos.size(); i++){
			std::string nam = "simple_"+bdt_infos[i].identifier+"_mva";
			t_sbnfit_simpletree->Branch(nam.c_str(), &(bdt_mvas[i]));
		}

		for(int i=0; i< vars.size(); i++){
			std::string tnam = "simple_bdt_var_"+std::to_string(i);
			t_sbnfit_simpletree->Branch(tnam.c_str(),&(simple_bdt_vars[i]));
		}

		TTreeFormula* weight = new TTreeFormula("weight_formula ",this->weight_branch.c_str(),this->tvertex);
		TTreeFormula* var = new TTreeFormula("var_formula ",input_string.c_str(),this->tvertex);

		std::vector<TTreeFormula*> form_vec;
		std::vector<TTreeFormula*> form_vec_vars;

		for(int i=0; i< bdt_infos.size();i++){
			std::string nam = this->tag+"_"+bdt_infos[i].identifier+".mva";
			form_vec.push_back(new TTreeFormula((bdt_infos[i].identifier+"_mva_formula").c_str(), nam.c_str(),this->tvertex));
		}

		for(int i=0; i< vars.size();i++){
			form_vec_vars.push_back(new TTreeFormula((vars[i].safe_unit).c_str(), vars[i].name.c_str(),this->tvertex));
		}


		std::string var_string = input_string;
		if(var_string == "") var_string = "reco_vertex_size";
		std::cout<<"Starting to make a simpletree with variable "<<var_string<<std::endl;
		for(int i=0; i< this->tvertex->GetEntries(); i++){
			this->tvertex->GetEntry(i); 

			CUT->GetNdata();
			bool is_is = CUT->EvalInstance();

			if(!is_is) continue;

			weight->GetNdata();
			var->GetNdata();
			simple_wei = weight->EvalInstance();

			/*
			   if(simple_wei<0 || simple_wei!=simple_wei || isinf(simple_wei) ){
			   std::cout<<"WARNING WARNING, the weight here is "<<simple_wei<<std::endl;
			   std::cout<<"Setting to 1 for now, investigate!"<<std::endl;
			   simple_wei = 1.0;
			   }
			   */

			simple_var = var->EvalInstance();
			simple_pot_wei = simple_wei*this->scale_data*plot_pot/this->pot;
			original_entry = i;

			for(int j=0; j< bdt_infos.size();j++){
				form_vec[j]->GetNdata();
				bdt_mvas[j] = form_vec[j]->EvalInstance();
			}

			for(int j=0; j< vars.size();j++){
				form_vec_vars[j]->GetNdata();
				simple_bdt_vars[j] = form_vec_vars[j]->EvalInstance();
			}
			t_sbnfit_simpletree->Fill();
		}


		TList * lf1 = (TList*)t_sbnfit_tree->GetListOfFriends();
		for(const auto&& obj: *lf1) t_sbnfit_tree->GetListOfFriends()->Remove(obj);
	}
//    TList * lf2 = (TList*)t_sbnfit_eventweight_tree->GetListOfFriends();
//    for(const auto&& obj: *lf2) t_sbnfit_eventweight_tree->GetListOfFriends()->Remove(obj);


    std::cout<<"Writing to file"<<std::endl;
//    cdtof->cd();
    t_sbnfit_tree->Write();
//    t_sbnfit_pot_tree->Write();
//    t_sbnfit_rs_tree->Write();
//    t_sbnfit_eventweight_tree->Write(); 
//    t_sbnfit_slice_tree->Write();
    t_sbnfit_simpletree->Write();
//    weight->Write();
//    var->Write();
//    CHECK, disable some branches;

    TVectorD POT_value(1);
    POT_value[0] = this->pot;
//    POT_value.Write("POT_value");

    f_sbnfit->Close();
    std::cout<<"Done!"<<std::endl;

    return 0;
}


