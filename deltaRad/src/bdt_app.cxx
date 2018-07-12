#include "bdt_app.h"


bdt_app_tree_struct::bdt_app_tree_struct(std::string const & tree_name, bool const extra) {
	tree = new TTree(tree_name.c_str(), "");
	if(extra) {
	//	tree->Branch("is_delta_rad", &is_delta_rad, "is_delta_rad/I");
	//	tree->Branch("true_nu_vtx_fid_contained", &true_nu_vtx_fid_contained, "true_nu_vtx_fid_contained/I");
	//	tree->Branch("reco_asso_tracks", &reco_asso_tracks, "reco_asso_tracks/I");
	}
	tree->Branch("mva", &mva, "mva/D");
	//tree->Branch("summed_associated_reco_shower_energy", &summed_associated_reco_shower_energy, "summed_associated_reco_shower_energy/D");
}

bdt_app_tree_struct::~bdt_app_tree_struct() {
	delete tree;
}



void bdt_app_update_formula(std::vector<TTreeFormula*> & tfv, std::vector<float *> & rvv) {

	for(size_t i = 0; i < tfv.size(); ++i) {
		 //This line must be here. see https://sft.its.cern.ch/jira/browse/ROOT-7465
		 tfv.at(i)->GetNdata();

		*rvv.at(i) = tfv.at(i)->EvalInstance();
	}
	  
}

void bdt_app_update(void_vec const & tvv, std::vector<float *> & rvv) {

	for(size_t i = 0; i < tvv.size(); ++i) {
		if(tvv.get_type(i) == typeid(int).name()) {
			*rvv.at(i) = *tvv.at<int>(i);
		}
		else if(tvv.get_type(i) == typeid(double).name()) {
			*rvv.at(i) = *tvv.at<double>(i);
		}
	}
}

void bdt_topological_update(std::vector<double*> *topv){


}
		

int bdt_app_tree(std::string identifier, TTree * tree, bdt_flow flow, std::string otree_name, std::vector<bdt_variable> vars, std::vector<method_struct> const & methods) {

	TMVA::Reader * reader = new TMVA::Reader("!Color:!Silent");

	void_vec tree_var_v;
	std::vector<float *> reader_var_v;
	std::vector<TTreeFormula*> tree_formulas_v;

	for(bdt_variable &p : vars) {
		/*
		if(p.type == "d") {

			double * dp = new double(-1);
			tree_var_v.push_back(dp);
			tree->SetBranchAddress(p.name.c_str(), dp);
		}
		else if(p.type == "i") {
			int * ip = new int(-1);
			tree_var_v.push_back(ip);
			tree->SetBranchAddress(p.name.c_str(), ip);
		}
		else {
			std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ <<"\n"
				<< "ERROR: invalid type: " << p.type << "\n";
		}
		*/	
       	        tree_formulas_v.push_back(new TTreeFormula(p.safe_name.c_str(), p.name.c_str() ,tree));
		reader_var_v.push_back(new float(-1));
		reader->AddVariable(p.name.c_str(), reader_var_v.back());
	}

	/*
	int shr1_index = 0;
	int shr2_index = 0;
	float f1 =0;
	float f2=0;

	int oshr1_index = 0;
	int oshr2_index = 0;
	std::vector<int> *of =0;
        TBranch *bvpx = 0;

	TTreeFormula *t1 = new TTreeFormula("i1","most_energetic_shower_index",tree);
	TTreeFormula *t2 = new TTreeFormula("i2","second_most_energetic_shower_index",tree);

	TTreeFormula *t3 = new TTreeFormula("i3","true_shower_parent_pdg[most_energetic_shower_index]",tree);
	TTreeFormula *t4 = new TTreeFormula("i4","true_shower_parent_pdg[second_most_energetic_shower_index]",tree);

	std::cout<<"GetNdata: "<<t3->GetNdata()<<" "<<t4->GetNdata()<<std::endl;
	
	tree->SetBranchAddress("most_energetic_shower_index",&oshr1_index);	
	tree->SetBranchAddress("second_most_energetic_shower_index",&oshr2_index);	
	tree->SetBranchAddress("true_shower_parent_pdg",&of,&bvpx);
	*/

	//TTreeFormula * tf = new TTreeFormula("tf", cut.c_str(), tree);
	TTreeFormula * tf_topological = new TTreeFormula("tf_top", flow.topological_cuts.c_str(), tree);
	TTreeFormula * tf_definition = new TTreeFormula("tf_def", flow.definition_cuts.c_str(), tree);
        std::cout << "Toplogical cuts: " << flow.topological_cuts.c_str() << std::endl;
        std::cout << "Definition cuts: " << flow.definition_cuts.c_str() << std::endl;

	for(method_struct const & method : methods) {
		reader->BookMVA(method.str.c_str(), ("BDTxmls_"+identifier+"/weights/"+identifier+"_"+method.str+".weights.xml").c_str());
		bdt_app_tree_struct ts(otree_name, false);
	
        int passed = 0;
        int N = tree->GetEntries();
        for(int i = 0; i < N; ++i) {
            if(i%50000==0){std::cout<<i<<"/"<<N<<std::endl;}
            tree->GetEntry(i);
            //bdt_app_update(tree_var_v, reader_var_v);
            ts.mva = -999;
            if(i%25000==0){
                std::cout<<i<<"/"<<N<<std::endl;
            }
	
	
	   /* std::cout<<"TTreeFormula:\t"<<t1->EvalInstance()<<" "<<t2->EvalInstance()<<" "<<t3->EvalInstance()<<" "<<t4->EvalInstance()<<std::endl;	
	    std::cout<<"DIRECT:\t\t"<<oshr1_index<<" "<<oshr2_index<<" "<<of->at(oshr1_index)<<" "<<of->at(oshr2_index)<<std::endl;	
	
	    if(t3->EvalInstance() == 111 && t4->EvalInstance() == 111){
		std::cout<<"Should Pass: check next line"<<std::endl;
	    }
	   */

	    //These two lines must be here. or FAIL. see https://sft.its.cern.ch/jira/browse/ROOT-7465
	    tf_topological->GetNdata();
	    tf_definition->GetNdata();

            if(tf_topological->EvalInstance()) {
            	bdt_app_update_formula(tree_formulas_v, reader_var_v);
            
                if(tf_definition->EvalInstance()) {
		            //std::cout<<"Pass Defin"<<std::endl;
                    ts.mva = reader->EvaluateMVA(method.str.c_str());
                    passed++;
                }
            }
            ts.tree->Fill();
        }
	    ts.tree->Write();
        std::cout << "Number of events passed: " << passed << std::endl;
    }
        /*
		for(int i = 0; i < tree->GetEntries(); ++i) {
			tree->GetEntry(i);
			bdt_app_update_formula(tree_var_v, reader_var_v);
			ts.mva = -999;
			if(tf->EvalInstance()) ts.mva = reader->EvaluateMVA(method.str.c_str());
			ts.tree->Fill();

		}
        */
	

	//delete reader;
	//delete tf;
	//for(float * f : reader_var_v) delete f;

}


int bdt_app(bdt_info info, std::vector<bdt_file*> files, std::vector<bdt_variable> vars, std::vector<method_struct> & method){

	std::string identifier = info.identifier;
	
	TFile * app_ofile = TFile::Open((identifier+"_app"+".root").c_str(), "update");
	for(size_t i = 0; i < files.size(); ++i) {
		std::cout<<"On file: "<<files.at(i)->tag<<std::endl;
		std::string bdt_response_friend_tree_name = files.at(i)->tag+"_"+info.identifier;
		bdt_app_tree(identifier, files.at(i)->tvertex, files.at(i)->flow, bdt_response_friend_tree_name , vars, method);
	}
	app_ofile->Close();
	delete app_ofile;



	return 0;
}
