#ifndef BDT_FILE_H
#define BDT_FILE_H

#include <vector>
#include <string>
#include <set>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <numeric>
/******** Our includes *****/

#include  "bdt_flow.h"
#include  "bdt_var.h"
#include  "bdt_info.h"
#include  "method_struct.h"
/******** Root includes *****/

#include "TTreeFormula.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMVA/Types.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TFriendElement.h"
#include "TText.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TEntryList.h"
#include "TSystem.h"

//#include "load_mva_param.h"

template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out <<std::fixed<< std::setprecision(n) << a_value;
    return out.str();
}

template<typename T>
bool marks_compare_vec_nonsense(std::vector<T>& v1, std::vector<T>& v2)
{
        std::sort(v1.begin(), v1.end());
            std::sort(v2.begin(), v2.end());
                return v1 == v2;
}


template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}


TText * drawPrelim(double x, double y,double s, std::string in);
TText * drawPrelim(double x, double y, std::string in);


class bdt_file{
    public:
		
		//essential members
		TString input_root;	
		std::string tag;

		//members with default value
		std::string weight_branch;
		std::string root_dir;//TDiretory
		std::string definition;
        bool is_data;
		bool is_signal;
		bool is_train; //determined in <training> block
		bool is_systematic;	//determined in <systematic> block
		int group;
        double pot;//determined in <fixPOT> section>
		double scale;

		// - plotting special
		std::string plot_name;
		std::string plot_style;
        bool is_stack;
		bool bdt_on_top;
        TColor* color;//this is also the fill color
        int fillstyle;
		int linecolor;
        int linestyle;
        std::string leg;//legend style; l - line; p - Marker; f - box
		
		//derived members
        TFile *file;
        TTree *tvertex;
		std::string s_precut_hash;


        //Optional members
        std::vector<std::string> friend_files;
        std::vector<std::string> friend_names;
        int numberofevents_raw;//before POT scaling?

		//items done in hive.cxx
        TEntryList * topological_list;
        TEntryList * precut_list;

		//function
		//constructor;
		bdt_file(std::string ininput_root, std::string in_tag)
		:
		input_root(ininput_root.c_str()),
		tag(in_tag){};

        ~bdt_file();


        void calcBaseEntryList(std::string);

		//copy constructor forbdt_sys; THIS CAUSE LINKING PROBLEM...
//		bdt_file( const bdt_file& copy);

		void SetDefaultAttributes();
		
		
		//do these in hive.cxx;
		std::vector< std::string > stage_cuts;
		//out this, use new variable above
        bdt_flow flow;





		//CHECK to do;
		void outputRoot(int stage);

        int addFriend(std::string in_friend_tree_nam, std::string in_friend_file);

        TH1* getTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT);
        TH2* getTH2(bdt_variable varx, bdt_variable vary, std::string cuts, std::string nam, double plot_POT);

        unsigned long jenkins_hash(std::string key);





//        int makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string & inpu);
        std::string dir;//file location
        std::string name;
//        std::string tag;
//        std::string plot_name;
//        std::string plot_ops;
//        std::string root_dir;//TDirectory

//        std::string weight_branch;




        TRandom3* rangen;
        std::string topo_name;

        //This is slightly deprecisated
//        std::string friend_tree_file;
//        std::string friend_tree_name;


        //These are now passed into bdt_recomc as they should be really.
        //Still used here (in get recomcbdts) but are filled from bdt_reco not in constructor.
        std::vector<std::string> recomc_cuts;
        std::vector<std::string> recomc_names;
        std::vector<int> recomc_cols;


	int col;
        int rebin;	//combing bins, maybe not use;

        bool is_bnbext;//bkg data; CHECK, might not use
        int numberofevents;


        //copy tvertex into topovertex, but with topological cut.
        TTree *topovertex;

        TTree *tevent;
        TTree *tpot;
        TTree *trs;
        TTree *teventweight;
        TTree *tslice;

        std::string cosmicbdt_list_name;
        TEntryList * cosmicbdt_list;
        std::string bnbbdt_list_name;
        TEntryList * bnbbdt_list;

        std::vector<TEntryList*> vec_entry_lists;


        //Run management stuff, there is 5 Runs R1,R2,R3a , R3b,R4
        std::vector<double> run_fractions_plot; //fractions to plot
        std::vector<double> run_fractions_file; //fractions in file
        std::vector<std::string> run_names;
        std::vector<std::string> run_fraction_cuts;
        std::string run_weight_string;

        //a function that splits a BDT file based on string and !string
//        int splitBDTfile(std::string split_string,std::string trueTAG, bdt_file* truesplit, std::string falseTAG, bdt_file *falsesplit);


        int setStageEntryList(int j);
        int setStageEntryList(int j, double, double);
        int calcBDTEntryList(int stage, std::vector<double> bdt_cuts);

        int scanStage(int which_stage, std::vector<double> bdt_cuts , std::string scan_string);


        double data_tor860_wcut;
        double data_spills_E1DCNT_wcut;
        double ext_spills_ext;
        double N_samweb_ext;

//        int setAsMC();
//        int setAsOverlay();
//        int setAsOnBeamData(double in_tor860_wcut);
//        int setAsOffBeamData(double in_data_tor860_wcut, double in_data_spills_E1DCNT_wcut, double in_ext_spills_ext, double N_samweb_ext);
//        int setAsOffBeamData(double in_data_tor860_wcut, double in_data_spills_E1DCNT_wcut, double in_ext_spills_ext);

//      int calcPOT();

//      int calcPOT(std::vector<std::string> run_names, std::vector<std::string> run_cuts, std::vector<double> run_fractions);

        int makeRunSubRunList();

        double scale_data;

        bdt_variable getBDTVariable(bdt_info info);
        bdt_variable getBDTVariable(bdt_info info, std::string bin);
        //legacy code, and damned lazy too
        //bdt_variable getBDTVariable(std::string cut);
		//new
		bdt_file(size_t index,
//			MVALoader XMLconfig,
			bdt_flow inflow);

        bdt_file(std::string indir,
		std::string inname, 
		std::string intag, 
		std::string inops, 
		std::string inrootdir,  
		int incol, 
		bdt_flow inflow);	
        
//		bdt_file(std::string indir,
//		std::string inname, 
//		std::string intag, 
//		std::string inops, 
//		std::string inrootdir, 
//		int incol, 
//		int fillstyle,
//		bdt_flow inflow);	

        //legacy code OBSOLETE
        //bdt_file(std::string indir,std::string inname, std::string intag, std::string inops, std::string inrootdir, std::string infriend, std::string infriendtree, int incol, bool indata);	


//        int scale(double scalein);
//        int setPOT(double inpot);

        TH1* getEventTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT);
//        int CheckWeights();
     
        double GetEntries(std::string cuts);
        double GetEntries();
        TH1* getTH1(std::string invar, std::string cuts, std::string nam, double plot_POT, int rebin);
        TH1* getTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT, int rebin);

        std::vector<TH1*> getRecoMCTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT);
        std::vector<TH1*> getRecoMCTH1(bdt_variable var, std::string cuts, std::string nam, double plot_POT, int rebin);

        int addBDTResponses(std::string dir, bdt_info input_bdt_info);


        int makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string & inpu);
        int makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string & inpu, const std::vector<bdt_variable> &vars ,double plot_pot);
        int makeSBNfitFile(const std::string &analysis_tag, const std::vector<bdt_info>& bdt_infos, int which_stage, const std::vector<double> & fbdtcuts, const std::string & inpu, const std::vector<bdt_variable> &vars ,double plot_pot,std::string external_cuts);


    int makePrecalcSBNfitFile(const std::string &analysis_tag, int which_stage, const std::vector<double> & fbdtcuts );



        std::string getStageCuts(int stage, double bdtvar1, double bdtvar2);
        std::string getStageCuts(int stage, std::vector<double> bdt_cuts);
		TString getStageCutsIndex(int stage, std::vector<double> bdt_cuts, int vec_index);


//        int writeStageFriendTree(std::string nam,double,double);
//        int addPlotName(std::string plotin);
//        int setTColor(TColor &);
};
//void get_joy();

//get binning vector for variables;
std::vector<double> gadget_CalBinning( std::vector<bdt_variable> cur_vars);
#endif
