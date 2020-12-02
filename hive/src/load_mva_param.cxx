#include "load_mva_param.h"
#include "method_struct.h"
		
bool gadget_boolreader( const char* intext){

	if( intext !=NULL){
		if (strcmp(intext ,"yes")==0||strcmp(intext,"true")==0){
			return true;
		}
	}
	return false;

};

//MVALoader::MVALoader(std::string xmlname): MVALoader(xmlname, 1) {}//activate verbosity

void MVALoader::LoadEnv(){

    //Setup TiXml documents
    doc = new TiXmlDocument(whichxml.c_str() );
    bool loadOkay = doc->LoadFile();

    if(loadOkay){
        std::cout<<"MVALoader::LoadEnv() || Loaded "<<whichxml<<std::endl;
    }else{
        std::cerr<<"ERROR: MVALoader::LoadEnv() || Failed to load "<<whichxml<<std::endl;
        std::cerr<<"ERROR: MVALoader::LoadEnv() || You probably just forgot to add a --xml my_analysis.xml , or maybe you missed a </>"<<std::endl;
        exit(EXIT_FAILURE);
    }

	//Verbosity & Input Directory
    TiXmlElement *pVerbo;
    pVerbo = doc->FirstChildElement("verbosity");
    if(pVerbo)
    {
        verbosity =std::stod(pVerbo->GetText());
    }

    TiXmlElement *pFileDir;
    pFileDir = doc->FirstChildElement("inputdir");
    if(pFileDir)
    {
        inputdir =std::string(pFileDir->GetText());
    }else{
        inputdir = "./";
    }

	//ALIAS
    if(verbosity>0) std::cout<<"\n #### <<Alias>> ####"<<std::endl;
    TiXmlElement *pAlias;
    pAlias = doc->FirstChildElement("alias");
    int n_alias = 0;
    while(pAlias)
    {
        std::string key = pAlias->Attribute("key");
        std::string val = pAlias->Attribute("value");
        aliasMap[key] = val; 
        std::cout<<"XML-Alias "<<key<<" --> "<<val<<std::endl;
        n_alias++;
        pAlias = pAlias->NextSiblingElement("alias");
    }


	//TOPOLOGY
    std::cout<<"\n#### <<Topology>> ####"<<std::endl;

    TiXmlElement *pTopoCut;
    pTopoCut = doc->FirstChildElement("topology");

    if(pTopoCut){
        const char* pTopo_tag = pTopoCut->Attribute("tag");
        const char* pTopo_def = pTopoCut->Attribute("def");
        const char* pTopo_name = pTopoCut->Attribute("name");
		try{
			if(pTopo_tag = NULL) throw "tag";
			if(pTopo_def = NULL) throw "def";
			if(pTopo_name = NULL) throw "name";
		}
		catch( std::string topo_Error){
			std::cerr<<"ERROR: MVALoader::LoadEnv() || ";
			std::cerr<<" attribute of "<<topo_Error<<" is missing, please check your xml"<<std::endl;
			exit(EXIT_FAILURE);
		}
		analysis_tag = pTopo_tag;
		topo_name = pTopo_name;

		stage_cuts.push_back(pTopo_def);
		stage_cuts.push_back("");//leave a space for precus;

		if(verbosity > 0) std::cout<<"Analysis tag: "<<analysis_tag<<" with "<<topo_name<<" definition:"<<pTopo_def<<std::endl;

        const char* pTopo_bdtcuts = pTopoCut->Attribute("bdtcut");
		if(pTopo_bdtcuts!=NULL){ 
			std::vector<std::string > pTopo_bdtcuts_v = gadget_Tokenizer<std::string> (pTopo_bdtcuts);
			stage_cuts.insert( stage_cuts.end(), pTopo_bdtcuts_v.begin(), pTopo_bdtcuts_v.end() );
			if(verbosity > 0){
				std::cout<<"BDT cuts are placed at: "<<pTopo_bdtcuts_v[0];
				for(int lndex = 1; lndex < pTopo_bdtcuts_v.size(); lndex++){
					std::cout<<","<<pTopo_bdtcuts_v[lndex];
				}
				std::cout<<std::endl;
			}
		}
//        pTopoCut = pTopoCut->NextSiblingElement("topology");
    }

    std::cout<<"\n#### <<Precuts>> ####"<<std::endl;

    TiXmlElement *pPreCut;
    pPreCut = doc->FirstChildElement("precut");
	std::string cut_total = "(1)";

    while(pPreCut )
    {
		const char* cut_def =  pPreCut->Attribute("def");
		if(cut_def !=NULL){
			pre_cuts.push_back( this->AliasParse( cut_def));
			cut_total += "&&("+(std::string)cut_def+")";
			const char* cut_name = pPreCut->Attribute("name");
			if(verbosity > 0){ 
				if(cut_name!=NULL)std::cout<<cut_name;
				std::cout<<" is defined as "<<cut_def<<std::endl;
			}
		}
		
		stage_cuts[1] = cut_total;

        pPreCut = pPreCut->NextSiblingElement("precut");
    }
};


void MVALoader::LoadBDTfiles(){

    std::cout<<"\n#### <<BDT Files>> ####"<<std::endl;
    //Grab the first element. Note very little error checking here! make sure they exist.
    TiXmlElement *pBDTfile; 
    pBDTfile = doc->FirstChildElement("bdtfile");

    if(!pBDTfile){
        std::cerr<<"ERROR: MVALoader::LoadBDTfiles() || "<<whichxml<<" contains no BDT files! "<<std::endl;
        exit(EXIT_FAILURE);
    }
	
	while(pBDTfile){
        const char* pBDTf_inputr = pBDTfile->Attribute("filename");
        const char* pBDTf_tag = pBDTfile->Attribute("tag");
		try{
			if(pBDTf_inputr = NULL) throw "filename";
			if(pBDTf_tag = NULL) throw "tag";
		}
		catch( std::string topo_Error){
			std::cerr<<"ERROR: MVALoader::LoadBDTfiles() || ";
			std::cerr<<" attribute of "<<topo_Error<<" is missing, please check your xml"<<std::endl;
			exit(EXIT_FAILURE);
		}
		std::string tmp_input = inputdir + "/" + (std::string) pBDTf_inputr;
		//Build the bdt_file class;
		bdt_file tmp_BDTf( tmp_input, (std::string) pBDTf_tag);
		tmp_BDTf.SetDefaultAttributes();

		const char* pBDTf_wgt = pBDTfile->Attribute("wgt");
		if(pBDTf_wgt!=NULL) tmp_BDTf.weight_branch = (std::string) pBDTf_wgt;

		const char* pBDTf_Tdir = pBDTfile->Attribute("Tdir");
		if(pBDTf_Tdir!=NULL) tmp_BDTf.root_dir = (std::string) pBDTf_Tdir;

		const char* pBDTf_isdata = pBDTfile->Attribute("isdata");
		if(pBDTf_isdata!=NULL) tmp_BDTf.is_data = gadget_boolreader(pBDTf_isdata);

		const char* pBDTf_issignal = pBDTfile->Attribute("issignal");
		if(pBDTf_issignal!=NULL) tmp_BDTf.is_signal = gadget_boolreader(pBDTf_issignal);

		const char* pBDTf_group = pBDTfile->Attribute("group");
		if(pBDTf_group!=NULL) tmp_BDTf.group = (int) std::stod(pBDTf_group);

		const char* pBDTf_scale = pBDTfile->Attribute("scale");
		if(pBDTf_scale!=NULL) tmp_BDTf.scale = (double) std::stod(pBDTf_scale);

		//go to plotstyle block
		TiXmlElement *pPlotstyle = pBDTfile->FirstChildElement("plotstyle");
		if(pPlotstyle){
			const char* pPlots_isstack = pPlotstyle->Attribute("isstack");
			if(pPlots_isstack!=NULL) tmp_BDTf.is_stack = gadget_boolreader(pPlots_isstack);

			const char* pPlots_color = pPlotstyle->Attribute("color");
			if(pPlots_color!=NULL){ 
				std::vector< double > v_col = gadget_Tokenizer< double > (pPlots_color);
				tmp_BDTf.color = new TColor(TColor::GetFreeColorIndex(),v_col[0],v_col[1],v_col[2]) ;
			}

			const char* pPlots_histstyle = pPlotstyle->Attribute("histstyle");
			if(pPlots_histstyle!=NULL) tmp_BDTf.plot_style = (std::string) pPlots_histstyle;

			const char* pPlots_fillstyle = pPlotstyle->Attribute("fillstyle");
			if(pPlots_fillstyle!=NULL) tmp_BDTf.fillstyle = (int) std::stod(pPlots_fillstyle);

			const char* pPlots_linecolor = pPlotstyle->Attribute("linecolor");
			if(pPlots_linecolor!=NULL) tmp_BDTf.linecolor = (int) std::stod(pPlots_linecolor);

			const char* pPlots_linestyle = pPlotstyle->Attribute("linestyle");
			if(pPlots_linestyle!=NULL) tmp_BDTf.linestyle = (int) std::stod(pPlots_linestyle);

			const char* pPlots_legend = pPlotstyle->Attribute("legend");
			if(pPlots_legend!=NULL) tmp_BDTf.leg = (std::string) pPlots_legend;

			const char* pPlots_on_top = pPlotstyle->Attribute("on_top");
			if(pPlots_on_top!=NULL) tmp_BDTf.bdt_on_top = gadget_boolreader(pPlots_on_top);

			const char* pPlots_plot_name = pPlotstyle->Attribute("plot_name");
			if(pPlots_plot_name!=NULL) tmp_BDTf.plot_name = (std::string) pPlots_plot_name;
		}


		//go to training block
      TiXmlElement *pTrain = pBDTfile->FirstChildElement("training");
      while(pTrain){
          tmp_BDTf.is_train = true;
          TiXmlElement *pCut = pTrain->FirstChildElement("cut");
			while(pCut){
				const char* text =  pCut->GetText();
				if(text != NULL) tmp_BDTf.definition += "&&(" + this->AliasParse(text) + ")";
				pCut = pCut->NextSiblingElement("cut");
			}

          pTrain = pTrain->NextSiblingElement("training");
      }//multiple <training> blocks will be concatenated;

		//go to definition block
		TiXmlElement *pDefinition = pBDTfile->FirstChildElement("definition");
		while(pDefinition){
			TiXmlElement *pCut = pDefinition->FirstChildElement("cut");
			while(pCut){
				const char* text =  pCut->GetText();
				if(text != NULL) tmp_BDTf.definition += "&&(" + this->AliasParse(text) + ")";
				pCut = pCut->NextSiblingElement("cut");
			}
			pDefinition = pDefinition->NextSiblingElement("definition");
		}//multiple <definition> blocks will be concatenated;

		//go to fixPOT block
		TiXmlElement *pfixPOT = pBDTfile->FirstChildElement("fixPOT");
		if(pfixPOT){
			std::string fixpot = std::string(pfixPOT->GetText());
			tmp_BDTf.pot = (double) stod(fixpot);
		}
		

		//go to systematic block build new object based on current bdtfiles;
		TiXmlElement *pSystematics = pBDTfile->FirstChildElement("systematics");
		if(pSystematics){
			isSys += pSystematics;//one on, all on;
          tmp_BDTf.is_systematic = true;

			bdt_sys tmp_BDTsys( tmp_BDTf);
			const char* pSys_systag = pSystematics->Attribute("systag");
			if(pSys_systag!=NULL) tmp_BDTsys.systag = (std::string) pSys_systag;

			const char* pSys_var = pSystematics->Attribute("var");
			if(pSys_var!=NULL) tmp_BDTsys.vars =gadget_Tokenizer<TString>( this->AliasParse(pSys_var) );

			const char* pSys_varname = pSystematics->Attribute("varname");
			if(pSys_varname!=NULL) tmp_BDTsys.vars_name =gadget_Tokenizer<TString>( this->AliasParse(pSys_varname) );

			const char* pSys_throws = pSystematics->Attribute("throws");
			if(pSys_throws!=NULL) tmp_BDTsys.throws = (int) std::stod(pSys_throws);

			const char* pSys_isCV = pSystematics->Attribute("isCV");
			if(pSys_isCV!=NULL) tmp_BDTsys.is_CV = gadget_boolreader(pSys_isCV);

			const char* pSys_isOM = pSystematics->Attribute("isOpticalModel");
			if(pSys_isOM!=NULL) tmp_BDTsys.is_OM = gadget_boolreader(pSys_isOM);
			
			all_systematics.push_back( tmp_BDTsys);
		}	

		all_files.push_back( tmp_BDTf);
        pBDTfile = pBDTfile->NextSiblingElement("bdtfile");
	}
}
	
//this feeds the bdt_info;
void MVALoader::LoadBDTSettings(){


    std::cout<<"\n#### <<BDT Training Setting>> ####"<<std::endl;
	TiXmlElement *pMVA; 
	
	//    //Grab the first element. Note very little error checking here! make sure they exist.
    pMVA = doc->FirstChildElement("mva");
	
	    if(!pMVA && verbosity > 0) std::cerr<<"Warnnig: MVALoader::LoadBDTSettings() || XMl contains no mva's! No BDT is used."<<whichxml<<std::endl;

    TMVA::Types  type_instance = TMVA::Types::Instance();
    int n_bdt = 0;
    while(pMVA )
    {
        if( gadget_boolreader(pMVA->Attribute("use")) ){
            const char* mva_type = pMVA->Attribute("type");
            const char* bdt_tag = pMVA->Attribute("tag");
            const char* bdt_name = pMVA->Attribute("name");
            const char* bdt_binning = pMVA->Attribute("binning");
			try{
				if(mva_type = NULL) throw "type";
				if(bdt_tag = NULL) throw "tag";
				if(bdt_name = NULL) throw "name";
				if(bdt_binning = NULL) throw "binning";
			}
			catch( std::string pMVA_Error){
				std::cerr<<"ERROR: MVALoader::LoadBDTSettings() || ";
				std::cerr<<" attribute of "<<pMVA_Error<<" is missing, please check your xml"<<std::endl;
				exit(EXIT_FAILURE);
			}
			bdt_info tmp_bdtInfo();

            TMVA::Types::EMVA tmva_type = type_instance.GetMethodType(mva_type.c_str());
			tmp_bdtInfo.type = tmva_type;

			tmp_bdtInfo.identifier = analysis_tag + (std::string) bdt_tag;
			tmp_bdtInfo.name = (std::string) bdt_name;
			tmp_bdtInfo.binning = (std::string) bdt_binning;

			//options are concatenated with : as the connector

           if(verbosity>0)  std::cout<<"\n["<<n_bdt<<"] BDT TAG: "<<bdt_tag<<", name: "<<bdt_name<<", binning: "<<bdt_binning<<std::endl;

            TiXmlElement *pMethod = pMVA->FirstChildElement("method");

				
            while(pMethod ){

				if( gadget_boolreader(pMethod->Attribute("use")) ){
                    const char* method_type = pMethod->Attribute("type");		
					if(method_type != NULL) tmp_bdtInfo.str = method_type;

                    std::string param_string = "!H:!V";
                    std::string xgboost_string = "";
                    TiXmlElement *pParam = pMethod->FirstChildElement("param");
                    while(pParam){
                        param_string += ":"+ std::string(pParam->GetText()) ;
						xgboost_string += "," + std::string(pParam->GetText());
                        pParam = pParam->NextSiblingElement("param");
                    }//-->end param loop

					if(method_type=="XGBoost"){
						//Loop over all parameters, splitting by "=" sign and saving parameters into a vector of pairs of strings.
						//x = 1 , y = 2-> pair(x,1) (y,2)
						if(verbosity > 0) std::cout<<"\tReading XGBoost config \n\t";

						std::replace(xgboost_string.begin(), xgboost_string.end(), "=",",");
						std::vector<std::string> vec_params = gadget_Tokenizer(std::string) (xgboost_string);

						for(int nndex = 0; nndex<vec_params/2;nndex ++){
							std::pair<std::string,std::string> pairs = std::make_pair(vec_params[2*nndex],vec_params[2*nndex+1]);
							(tmp_bdtInfo.xg_config).push_back(pairs);
							if(verbosity>1) std::cout<<vec_params[2*nndex]<<" = "<<vec_params[2*nndex+1]<<", ";
						}
						if(verbosity>1) std::cout<<std:endl;
					}
					//choose what bdt_files are used  in this training;
					TiXmlElement *pMVAfile = pMVA->FirstChildElement("file");
					while(pMVAfile){
						const char* tmp_bkg_train_tag = pMVAfile->Attribute("bkg_train_tag");
						const char* tmp_bkg_test_tag = pMVAfile->Attribute("bkg_test_tag");
						const char* tmp_bkg_test_cut = pMVAfile->Attribute("bkg_test_cut");
						const char* tmp_sig_train_tag = pMVAfile->Attribute("sig_train_tag");
						const char* tmp_sig_test_tag = pMVAfile->Attribute("sig_test_tag");
						const char* tmp_sig_test_cut = pMVAfile->Attribute("sig_test_cut");

						if(tmp_bkg_train_tag!=NULL)		tmp_bdtInfo.bkg_train_tag = tmp_bkg_train_tag;
						if(tmp_bkg_test_tag !=NULL)		tmp_bdtInfo.bkg_test_tag  = tmp_bkg_test_tag;
						if(tmp_bkg_test_cut !=NULL)		tmp_bdtInfo.bkg_test_cut  = tmp_bkg_test_cut;
						if(tmp_sig_train_tag!=NULL)		tmp_bdtInfo.sig_train_tag = tmp_sig_train_tag;
						if(tmp_sig_test_tag !=NULL)		tmp_bdtInfo.sig_test_tag  = tmp_sig_test_tag;
						if(tmp_sig_test_cut !=NULL)		tmp_bdtInfo.sig_test_cut  = tmp_sig_test_cut;
						pMVAfile = pMVA->NextSiblingElement("file");
					}//end file


            }
//            while(pMethod ){
//                if((std::string)pMethod->Attribute("use") == "yes"){
//
//
//                    std::vector<std::pair<std::string,std::string>> xg_config;
//
//                    method_struct temp_struct = {tmva_type , method_type, param_string};
//                    temp_struct.bdt_tag = bdt_tag;
//                    temp_struct.bdt_name = bdt_name;
//                    temp_struct.bdt_binning = bdt_binning;
//                    temp_struct.precuts = precuts;
//                    temp_struct.xg_config = xg_config;
//
//                    vec_methods.push_back(temp_struct);		
//
//                    if(verbosity > 0) std::cout<<" MVALoader::MVALoader || Loading a method: "<<mva_type<<"::"<<method_type<<" with params: "<<param_string<<std::endl;
//
//
//                }
//                pMethod = pMethod->NextSiblingElement("method");
//
//            }//-->end method loop	
//
//
//            TiXmlElement *pMVAfile = pMVA->FirstChildElement("file");
//            while(pMVAfile){
//                //filename, foldername, traning_cut, training_fraction
//
//                //Old style
//                //vec_methods.back().filename = pMVAfile->Attribute("filename"); 
//                //vec_methods.back().foldername = pMVAfile->Attribute("foldername"); 
//                //vec_methods.back().training_cut = this->AliasParse(pMVAfile->Attribute("training_cut")); 
//                //vec_methods.back().training_fraction= strtof(pMVAfile->Attribute("training_fraction"),NULL); 
//
//                vec_methods.back().bkg_train_tag = pMVAfile->Attribute("bkg_train_tag");
//                vec_methods.back().bkg_test_tag = pMVAfile->Attribute("bkg_test_tag");
//                vec_methods.back().bkg_test_cut = pMVAfile->Attribute("bkg_test_cut");
//                vec_methods.back().sig_train_tag = pMVAfile->Attribute("sig_train_tag");
//                vec_methods.back().sig_test_tag = pMVAfile->Attribute("sig_test_tag");
//                vec_methods.back().sig_test_cut = pMVAfile->Attribute("sig_test_cut");
//
//                pMVAfile = pMVA->NextSiblingElement("file");
//            }//end file
//
//
//            TiXmlElement *pMVAscan = pMVA->FirstChildElement("scan");
//            if(!pMVAscan){
//                vec_methods.back().scan_max = -99; 
//                vec_methods.back().scan_min = -99; 
//                vec_methods.back().scan_steps = -99; 
//            }else{
//                while(pMVAscan){
//                    vec_methods.back().scan_max = strtof(pMVAscan->Attribute("max"),NULL); 
//                    vec_methods.back().scan_min = strtof(pMVAscan->Attribute("min"),NULL); 
//                    vec_methods.back().scan_steps = strtof(pMVAscan->Attribute("steps"),NULL); 
//                    std::cout<<"\rScan params "<<vec_methods.back().scan_steps<<" "<<vec_methods.back().scan_min<<" "<<vec_methods.back().scan_max<<std::endl;
//                    pMVAscan = pMVA->NextSiblingElement("scan");
//                }//end scan
//            }
//
//            vec_methods.back().topological_definition = topo_def;
//            vec_methods.back().topological_name = topo_name;
//
        }
        n_bdt++;
        pMVA = pMVA->NextSiblingElement("mva");

    }//--> end mva loop



};




void MVALoader::LoadVariables(){};


//MVALoader::MVALoader(std::string xmlname, int Verbose_in) :whichxml(xmlname) {

//
//    std::cout<<"\n####################### Variables ########################################"<<std::endl;
//     //first lets see if there is a covariance general file (GLOBAL)
//    TiXmlElement *pCovar = doc.FirstChildElement("covar");
//
//    bool has_global_covar = false;
//    std::string global_covar_dir;
//    std::string global_covar_name;
//    std::string global_leg_name;
//    while(pCovar){
//        has_global_covar = true;
//
//        const char* var_covar_dir = pCovar->Attribute("dir");
//        const char* var_covar_name = pCovar->Attribute("name");
//        const char* var_leg_name = pCovar->Attribute("plotname");
//        if (var_covar_dir==NULL || var_covar_name==NULL){
//            has_global_covar = false;
//        }else{
//            global_covar_dir = var_covar_dir;
//            global_covar_name = var_covar_name;
//            global_leg_name = var_leg_name;
//            std::cout<<"Loading a GLOBAL covariance matrix direectory"<<std::endl;
//        }
//    
//        pCovar = pCovar->NextSiblingElement("covar");
//
//    }
//
//
//
//	//bdt variables 
//    TiXmlElement *pVar = doc.FirstChildElement("var");
//    std::vector<bdt_variable> bdt_all_vars;
////    std::vector<bdt_variable> bdt_train_vars;
////    std::vector<bdt_variable> bdt_spec_vars;
//
//    int n_var = 0;
//    while(pVar){
//		std::vector< std::string> read_labels = {"def", "binning", "unit"}; 
//
//		std::vector< std::string> pVar_input_values(read_labels.size()); 
//		for(int jndex  = 0; jndex < read_labels.size(); jndex++){
//
//			pVar_input_values[jndex] = pVar->Attribute( read_labels[jndex].c_str() );
//
//			if(read_labels[jndex] == "def" ) pVar_input_values[jndex] = this->AliasParse(pVar_input_values[jndex]);
//		}
//
//        bdt_variable tvar( pVar_input_values[0], pVar_input_values[1], pVar_input_values[2]);
//
//		tvar.id = n_var;
//		std::cout<<"["<<n_var<<"] "<< tvar.name<<" ("<<tvar.unit<<") in "<<tvar.binning<<std::endl;
//
//		tvar.is_logplot = gadget_boolreader( pVar->Attribute("logplot"));
//		if(tvar.is_logplot) std::cout<<"\tLogplot!"<<std::endl;
//
//		tvar.is_custombin = gadget_boolreader( pVar->Attribute("custombin"));
//		if(tvar.is_custombin) std::cout<<"\tVariable binning!"<<std::endl;
//		
//		
//		//The following const char* is for Attribute that can be NULL;
//		const char* t_covar_files =  pVar->Attribute("covarfile");//can be empty;
//		if( t_covar_files !=NULL && strlen(t_covar_files) > 0){
//			tvar.covar_file = t_covar_files;
//			tvar.has_covar = true;
//            std::cout<<"\tLink fractional covaraince matrices from "<<tvar.covar_file<<std::endl;
//		}
//
//        const char* t_cats = pVar->Attribute("group");
//		if( t_cats !=NULL && strlen(t_cats) > 0 ){ 
//			std::string t_cats_str (t_cats);
//			std::vector<std::string> t_cats_v = gadget_Tokenizer<std::string> ( t_cats_str);
//			std::cout<<"\tGroup "<<t_cats<<std::endl;
//
//			tvar.cats.clear();
//			for(int jndex = 0; jndex < t_cats_v.size(); jndex++){
//				(tvar.cats).push_back( std::stod(t_cats_v[jndex]));
//			}
//		}
//
//        const char* t_pheight = pVar->Attribute("plotheight");
//		if(t_pheight!=NULL && strlen(t_pheight) > 0 && std::stod(t_pheight) > 0){ 
//			std::cout<<"\tSet plotheight: "<<t_pheight<<std::endl;
//			tvar.plot_height =(double) std::stod(t_pheight);
//		}
//
//
//        const char* t_train = pVar->Attribute("training");
//        std::vector<int> var_train_int;
//		if(t_train!=NULL && strlen(t_train) > 0){
//
//			std::string t_train_str (t_train);
//			std::vector<std::string> t_train_v = gadget_Tokenizer<std::string>( t_train_str);
//
//			for(int jndex = 0; jndex < t_train_v.size(); jndex++){
//				var_train_int.push_back( (int) std::stod( t_train_v[jndex]));
//			}
//			//for (auto && c : t_train_v) {
//			//	var_train_int.push_back((int)c - '0');
//			//}
//			tvar.is_train = true;
//			std::cout<<" -- Variable training string is "<<t_train_str<<std::endl;
//		}
//
//        //Loop over vec_methods, fill in bdt_train_vars;
//        for(int p=0; p< vec_methods.size(); p++){
//
////            bool is_train = false;
//            for(int k=0; k< var_train_int.size(); k++){
//                if(var_train_int[k]> vec_methods.size()){
//                    std::cout<<"ERROR! BDT variable: "<<tvar.name <<" has been assigned to train on BDT "<<var_train_int[k]<<" But only "<<vec_methods.size()<<" has been defined!"<<std::endl;
//                    exit(EXIT_FAILURE);
//                }
////                if(p==var_train_int[k]){
////                    is_train = true;
////                    break;
////                }
//            }
//            if(tvar.is_train){
//                vec_methods[p].bdt_train_vars.push_back(tvar);
//                std::cout<<" -- so adding "<<tvar.name<<" as training to method "<<vec_methods[p].bdt_name<<std::endl;
//            }
//        }
//
//		//finish up binning info;
//		tvar.load_bininfo();
//
//		//add tvar;
//        bdt_all_vars.push_back(tvar);
//
//        n_var++;
//        pVar = pVar->NextSiblingElement("var");
//    }//Next variable;
//
//
//    for(auto &v: vec_methods){
//        v.bdt_all_vars = bdt_all_vars;
//    }
//
//
//    std::cout<<"\n####################### Efficiency  ########################################"<<std::endl;
//
//    TiXmlElement *pEff = doc.FirstChildElement("efficiency");
//
//    while(pEff){
// 
//        const char* t_denom_stage = pEff->Attribute("denom_stage");
//            if(t_denom_stage==NULL){
//                std::cout<<"Setting default denominator stage of -1"<<std::endl;
//                v_eff_denom_stage.push_back(-1);
//            }else{
//                v_eff_denom_stage.push_back(std::stoi(t_denom_stage,NULL));
//                std::cout<<"Setting Denom Stage of "<<v_eff_denom_stage.back()<<std::endl;
//            }
//
//        const char* t_denom_cut = pEff->Attribute("denom_cut");
//            if(t_denom_cut==NULL){
//                std::cout<<"Setting default denominator cut of `1` "<<std::endl;
//                v_eff_denom_cut.push_back("1");
//            }else{
//                v_eff_denom_cut.push_back(t_denom_cut);
//                std::cout<<"Setting Denom Cut of "<<t_denom_cut<<std::endl;
//            }
//  
//        const char* t_numer_stage = pEff->Attribute("numer_stage");
//            if(t_numer_stage==NULL){
//                std::cout<<"Setting default numerinator stage of -1"<<std::endl;
//                v_eff_numer_stage.push_back(-1);
//            }else{
//                v_eff_numer_stage.push_back(std::stoi(t_numer_stage,NULL));
//                std::cout<<"Setting Numer Stage of "<<v_eff_numer_stage.back()<<std::endl;
//            }
//
//        const char* t_numer_cut = pEff->Attribute("numer_cut");
//            if(t_numer_cut==NULL){
//                std::cout<<"Setting default numerinator cut of `1` "<<std::endl;
//                v_eff_numer_cut.push_back("1");
//            }else{
//                v_eff_numer_cut.push_back(t_numer_cut);
//                std::cout<<"Setting Numer Cut of "<<t_numer_cut<<std::endl;
//            }
//
//
//        pEff = pEff->NextSiblingElement("efficiency");
//
//    }
//
//
////    std::cout<<"\n ####################### RECO-MC Matching ########################################"<<std::endl;
////
////
////
////    TiXmlElement *pRecoMC = doc.FirstChildElement("recomc");
////    while(pRecoMC){
////
////
////        TiXmlElement *pDef = pRecoMC->FirstChildElement("def");
////        while(pDef){
////
////            const char* t_recomc_name = pDef->Attribute("name");
////            if(t_recomc_name==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || recomc has no `name` attribute! "<<std::endl; exit(EXIT_FAILURE);}
////            recomc_names.push_back(t_recomc_name);
////
////            const char* t_col = pDef->Attribute("col");
////            if(t_col==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || recomc has no `col` attribute! Come and set a color! "<<std::endl; exit(EXIT_FAILURE);}
////            std::string s_col = t_col;
////            s_col.erase(std::remove(s_col.begin(), s_col.end(), '('), s_col.end());
////            s_col.erase(std::remove(s_col.begin(), s_col.end(), ')'), s_col.end());
////            std::vector<double> v_col; 
////
////            size_t pos = 0;
////            std::string delim = ",";
////            std::string token;
////            while ((pos = s_col.find(delim)) != std::string::npos) {
////                token = s_col.substr(0, pos);
////                v_col.push_back(std::stod(token));
////                s_col.erase(0, pos + delim.length());
////            }
////            v_col.push_back(std::stod(s_col));
////            for(int cc=0; cc< v_col.size();cc++) {
////                if(v_col[cc]>1) v_col[cc] = v_col[cc]/255.0;
////            }
////            recomc_cols.push_back(new TColor(TColor::GetFreeColorIndex(),v_col[0],v_col[1],v_col[2]) );
////
////
////            std::string unpar =  pDef->GetText();
////            std::string parsed = this->AliasParse(unpar);
////            recomc_defs.push_back(parsed);
////            pDef = pDef->NextSiblingElement("def");
////
////            std::cout<<"RecoMC cut "<<recomc_names.back()<<" "<<recomc_defs.back()<<" Color: ";recomc_cols.back()->Print();std::cout<<std::endl;
////
////        }
////        pRecoMC = pRecoMC->NextSiblingElement("recomc");
////    }//-->end definition
//	
//	//systematics
//    std::cout<<"####################### Systematics  ########################################"<<std::endl;
//    TiXmlElement *pSys = doc.FirstChildElement("systematics");
//
//    n_sys = 0;
//    while(pSys){
//		sys_tag.push_back(pSys->Attribute("tag"));
//		sys_dir.push_back(pSys->Attribute("dir"));
//		sys_filename.push_back(pSys->Attribute("file"));
//		sys_treename.push_back(pSys->Attribute("tree"));
//
//		//numbers
//		sys_pot.push_back(atof( pSys->Attribute("pot") ));
//		sys_throws.push_back(atof( pSys->Attribute("throws") ));
//	
//		//boolean 
//		sys_its_CV.push_back(gadget_boolreader( pSys->Attribute("isCV")));
//		sys_its_multithrows.push_back(gadget_boolreader( pSys->Attribute("isMultithrows")));
//		sys_its_OpticalModel.push_back(gadget_boolreader(pSys->Attribute("OpticalModel")));
//
////tokenize the following var, varnam, bdtfiles
////
//
//        std::string bdtfiles_need_alias = pSys->Attribute("bdtfiles");
//        std::string aliased_bdtfiles = this->AliasParse(bdtfiles_need_alias); 
//		sys_for_files.push_back(gadget_Tokenizer<TString>( aliased_bdtfiles));
//
//		sys_vars.push_back(gadget_Tokenizer<TString>( this->AliasParse(pSys->Attribute("var"))));
//		sys_vars_name.push_back(gadget_Tokenizer<TString>( this->AliasParse(pSys->Attribute("varnam"))));
//
//		if(true){//print the systemati files information out;
//			std::cout<<"File directory: "<<sys_dir[n_sys];
//			std::cout<<"/"<<sys_filename[n_sys]<<std::endl;
//
//			std::cout<<"For "<<sys_for_files[n_sys].size()<<" bdtfiles: ";
//			for(size_t index = 0; index < sys_for_files[n_sys].size(); ++index){
//				std::cout<<sys_for_files[n_sys][index]<<" ";
//			}
//
//			std::cout<<"\nVariables: ";	
//			for(size_t index = 0; index < sys_vars[n_sys].size(); ++index){
//				std::cout<<sys_vars[n_sys][index]<<" ";
//			}
//			std::cout<<"\t Names: ";
//			for(size_t index = 0; index < sys_vars_name[n_sys].size(); ++index){
//				std::cout<<sys_vars_name[n_sys][index]<<" ";
//			}
//			std::cout<<std::endl;
//
//			std::cout<<std::setw(18)<<std::left<<"Tag "<<std::setw(9)<<"TTree ";
//			std::cout<<std::setw(12)<<"POT "<<std::setw(6)<<"throws ";
//			std::cout<<std::setw(5)<<"CV "<<std::setw(8)<<"Nthrows ";
//			std::cout<<std::setw(8)<<"Nfiles "<<std::endl;
//
//			std::cout<<std::setw(18)<<sys_tag[n_sys];
//			std::cout<<std::setw(9)<<sys_treename[n_sys];
//			std::cout<<std::setw(12)<<sys_pot[n_sys];
//			std::cout<<std::setw(6)<<sys_throws[n_sys];
//
//
//			std::vector< std::string > t_of_f = {" Yes ", " No "};
//			int toff_index = 1;
//			toff_index = (sys_its_CV[n_sys])? 0: 1;
//			std::cout<<std::setw(5)<<t_of_f[toff_index];
//
//			toff_index = (sys_its_multithrows[n_sys])? 0: 1;
//			std::cout<<std::setw(8)<<t_of_f[toff_index];
//
//			toff_index = (sys_its_OpticalModel[n_sys])? 0: 1;
//			std::cout<<std::setw(8)<<t_of_f[toff_index]<<"\n"<<std::endl;
//			
//		}
//        n_sys++;
//        pSys = pSys->NextSiblingElement("systematics");
//    }//end of systematics
//	std::cout<<std::endl;
//};



//std::vector<method_struct> MVALoader::GetMethods(){
//    return vec_methods; 
//}



std::string MVALoader::AliasParse(std::string in){
    std::string delim = "#";
    std::string compressed = in;

    size_t pos = 0;
    size_t cur = 0;
    size_t n_del = std::count(compressed.begin(), compressed.end(), '#');

    if(n_del%2!=0){
        std::cerr<<"ERROR! AliasParse has found an ODD number of deliminters # "<<n_del<<" in the string "<<in<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::string key;
    while ((pos = compressed.find(delim,0)) != std::string::npos) {

        size_t start = pos+1;
        size_t end = compressed.find(delim,start);

        key = compressed.substr(start, end-start);
        //std::cout<<start<<" "<<end<<" "<<key<<std::endl;

        if(aliasMap.count(key)==0){
            std::cerr<<"ERROR! AliasParse has found a key "<<key<<" that is not in the map. Check your spelling"<<std::endl;
            exit(EXIT_FAILURE);
        }

        compressed.replace(start-1,key.size()+2,aliasMap[key]);
        //std::cout<<compressed<<std::endl; 
    }

    return compressed;
}
