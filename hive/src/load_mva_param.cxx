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

void MVALoader::LoadEnv(){};
void MVALoader::LoadBDTfiles(){};
void MVALoader::LoadVariables(){};
void MVALoader::LoadSystematics(){};


//MVALoader::MVALoader(std::string xmlname, int Verbose_in) :whichxml(xmlname) {
//
//    verbosity = Verbose_in;
//
//
//    //Setup TiXml documents
//    TiXmlDocument doc( xmlname.c_str() );
//    bool loadOkay = doc.LoadFile();
//
//    if(loadOkay){
//        if(verbosity > 0) std::cout<<"MVALoader::MVALoader || Loaded "<<whichxml<<std::endl;
//    }else{
//        std::cerr<<"ERROR: MVALoader::MVALoader || Failed to load "<<whichxml<<std::endl;
//        std::cerr<<"ERROR: MVALoader::MVALoader || You probably just forgot to add a --xml my_analysis.xml , or maybe you missed a </>"<<std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//    TiXmlHandle hDoc(&doc);
//
//    std::cout<<"\n########################### Alias ###########################"<<std::endl;
//    TiXmlElement *pAlias;
//    pAlias = doc.FirstChildElement("alias");
//    int n_alias = 0;
//    while(pAlias)
//    {
//        std::string key = pAlias->Attribute("key");
//        std::string val = pAlias->Attribute("value");
//        aliasMap[key] = val; 
//        std::cout<<"XML-Alias "<<key<<" --> "<<val<<std::endl;
//        n_alias++;
//        pAlias = pAlias->NextSiblingElement("alias");
//    }
//
//
//
//
//    std::cout<<"\n########################### Topology ###########################"<<std::endl;
//
//    TiXmlElement *pTopoCut;
//    pTopoCut = doc.FirstChildElement("topology");
//    std::string topo_def;
//    std::string topo_name;
//    while(pTopoCut )
//    {
//
//        topo_def = pTopoCut->Attribute("def");
//        topo_name = pTopoCut->Attribute("name");
//        analysis_tag = pTopoCut->Attribute("tag");
//
//        const char* t_cut = pTopoCut->Attribute("bdtcut");
//        if(t_cut==NULL){ 
//
//        }else{
//            std::string s_cuts = t_cut;
//            s_cuts.erase(std::remove(s_cuts.begin(), s_cuts.end(), '('), s_cuts.end());
//            s_cuts.erase(std::remove(s_cuts.begin(), s_cuts.end(), ')'), s_cuts.end());
//
//            size_t pos = 0;
//            std::string delim = ",";
//            std::string token;
//            while ((pos = s_cuts.find(delim)) != std::string::npos) {
//                token = s_cuts.substr(0, pos);
//                bdt_cuts.push_back(std::stod(token));
//                s_cuts.erase(0, pos + delim.length());
//            }
//            bdt_cuts.push_back(std::stod(s_cuts));
//        }
//
//
//        std::cout<<"Loading Topology "<<topo_name<<" with definition "<<topo_def<<std::endl;
//        pTopoCut = pTopoCut->NextSiblingElement("topology");
//    }
//
//    TiXmlElement *pFileDir;
//    pFileDir = doc.FirstChildElement("inputdir");
//    if(pFileDir)
//    {
//        inputdir =std::string(pFileDir->GetText());
//        pFileDir = pFileDir->NextSiblingElement("inputdir");
//    }else{
//        inputdir = "./";
//    }
//
//
//
//
//
//    TiXmlElement *pPreCut;
//    pPreCut = doc.FirstChildElement("precut");
//
//    std::vector<std::string> precuts;
//
//    std::cout<<"\n########################### Precuts ###########################"<<std::endl;
//    while(pPreCut )
//    {
//
//        std::string cut_def_unparsed = pPreCut->Attribute("def");
//        std::string cut_def = this->AliasParse(cut_def_unparsed); 
//
//
//        std::string cut_name = pPreCut->Attribute("name");
//
//        std::cout<<"Loading Precut number "<<precuts.size()<<" "<<cut_name<<std::endl;
//        std::cout<<"--- Define: "<<cut_def<<std::endl;
//        precuts.push_back(cut_def);
//        pPreCut = pPreCut->NextSiblingElement("precut");
//
//    }
//
//    std::cout<<"\n######################## BDT's to Train on ########################################"<<std::endl;
//
//
//    TiXmlElement *pMVA; 
//
//
//    //Grab the first element. Note very little error checking here! make sure they exist.
//    pMVA = doc.FirstChildElement("mva");
//    if(!pMVA) {
//        std::cerr<<"Warnnig: MVALoader::MVALoader || XMl contains no mva's! No BDT is used."<<whichxml<<std::endl;
////        std::cerr<<"ERROR: MVALoader::MVALoader || XMl contains no mva's! "<<whichxml<<std::endl;
////        exit(EXIT_FAILURE);
//    }
//
//    TMVA::Types  type_instance = TMVA::Types::Instance();
//    int n_bdt = 0;
//    while(pMVA )
//    {
//        if( (std::string)pMVA->Attribute("use") == "yes"){
//            std::string mva_type = pMVA->Attribute("type");
//            //for each type, find all methods to be used
//
//            std::string bdt_tag = pMVA->Attribute("tag");
//            std::string bdt_name = pMVA->Attribute("name");
//            std::string bdt_binning = pMVA->Attribute("binning");
//
//            std::cout<<"\n["<<n_bdt<<"] BDT TAG: "<<bdt_tag<<", name: "<<bdt_name<<", binning: "<<bdt_binning<<std::endl;
//
//            //use TMVA instance to get the right EMVA type
//            TMVA::Types::EMVA tmva_type = type_instance.GetMethodType(mva_type.c_str());
//
//            TiXmlElement *pMethod = pMVA->FirstChildElement("method");
//            while(pMethod ){
//                if((std::string)pMethod->Attribute("use") == "yes"){
//
//                    std::string method_type = pMethod->Attribute("type");		
//
//
//                    std::vector<std::string> vec_params;	
//                    std::string param_string = "!H:!V";
//
//                    TiXmlElement *pParam = pMethod->FirstChildElement("param");
//                    while(pParam){
//                        vec_params.push_back( std::string(pParam->GetText()) );
////                        std::cout<<vec_params.back()<<", ";
//                        pParam = pParam->NextSiblingElement("param");
//                    }//-->end param loop
////					std::cout<<std::endl;
//
//                    for(std::string p: vec_params){
//                        param_string = param_string + ":" +p;
//                    }	
//
//
//                    std::vector<std::pair<std::string,std::string>> xg_config;
//                    if(method_type=="XGBoost"){
//                        //Loop over all parameters, splitting by "=" sign and saving parameters into a vector of pairs of strings.
//                            std::cout<<"\tReading XGBoost config \n\t";
//
//                        for(auto &p: vec_params){
//                            size_t pos = 0;
//                            std::string delim = "=";
//                            std::string firstone;
//                            while((pos = p.find(delim)) != std::string::npos) {
//                                firstone = p.substr(0, pos);
//                                p.erase(0, pos + delim.length());
//                            }
//                            std::string secondone = p;
//                            std::pair<std::string,std::string> pairs = std::make_pair(firstone,secondone);
//                            xg_config.push_back(pairs);
//                            std::cout<<firstone<<" = "<<secondone<<", ";
//                        }
//						std::cout<<std::endl;
//                    }
//
//
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
//        }
//        n_bdt++;
//        pMVA = pMVA->NextSiblingElement("mva");
//
//    }//--> end mva loop
//
//
////    std::cout<<"\n#######################  BDT_Files  ########################################"<<std::endl;
////    std::cout<<"(Print out will be below when constructing bdt_file class)"<<std::endl;
//
//    TiXmlElement *pBDTfile; 
//
//
//    //Grab the first element. Note very little error checking here! make sure they exist.
//    pBDTfile = doc.FirstChildElement("bdtfile");
//    if(!pBDTfile) {
//        std::cerr<<"ERROR: MVALoader::MVALoader || XMl contains no BDT_files! "<<whichxml<<std::endl;
//        exit(EXIT_FAILURE);
//    }
//    n_bdt_files = 0;
//
//	//essential information: filename, tag;
//	//with default: systematic, definition;
//    while(pBDTfile)
//    {
//
//        const char* t_tag = pBDTfile->Attribute("tag");
//        if(t_tag==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has no `tag` attribute! "<<std::endl; exit(EXIT_FAILURE);}
//        bdt_tags.push_back(t_tag);
//
//        const char* t_filename = pBDTfile->Attribute("filename");
//        if(t_filename==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has no `filename` attribute! "<<std::endl; exit(EXIT_FAILURE);}
//        bdt_filenames.push_back(t_filename);
//
//        const char* t_hist_styles = pBDTfile->Attribute("hist_style");
//        if(t_hist_styles==NULL){ bdt_hist_styles.push_back("hist");}else{
//            bdt_hist_styles.push_back(t_hist_styles);
//        }
//
//        const char* t_fillstyles = pBDTfile->Attribute("fillstyle");
//		if(t_fillstyles==NULL){ 
//			bdt_fillstyles.push_back(1001);
//		}else{
//			bdt_fillstyles.push_back((int)std::stoi(t_fillstyles,nullptr,10));
//		}
//
//        const char* t_linecol = pBDTfile->Attribute("linecol");
//		if(t_linecol==NULL){ 
//			bdt_linecols.push_back(1);
//		}else{
//			bdt_linecols.push_back((int)std::stoi(t_linecol,nullptr,10));
//		}
//
//        const char* t_linestyles = pBDTfile->Attribute("linestyle");
//		if(t_linestyles==NULL){ 
//			bdt_linestyles.push_back(1);
//		}else{
//			bdt_linestyles.push_back((int)std::stoi(t_linestyles,nullptr,10));
//		}
//
//        const char* t_dirs = pBDTfile->Attribute("dirs");
//        if(t_dirs==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has no `dirs` attribute! "<<std::endl; exit(EXIT_FAILURE);}
//        bdt_dirs.push_back(t_dirs);
//
//        const char* t_col = pBDTfile->Attribute("col");
//        if(t_col==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has no `col` attribute! "<<std::endl; exit(EXIT_FAILURE);}
//        std::string s_col = t_col;
//        s_col.erase(std::remove(s_col.begin(), s_col.end(), '('), s_col.end());
//        s_col.erase(std::remove(s_col.begin(), s_col.end(), ')'), s_col.end());
//        std::vector<double> v_col; 
//
//        size_t pos = 0;
//        std::string delim = ",";
//        std::string token;
//        while ((pos = s_col.find(delim)) != std::string::npos) {
//            token = s_col.substr(0, pos);
//            v_col.push_back(std::stod(token));
//            s_col.erase(0, pos + delim.length());
//        }
//        v_col.push_back(std::stod(s_col));
//        //for(auto &d:v_col)std::cout<<d<<std::endl;
//        for(int cc=0; cc< v_col.size();cc++) {
//            if(v_col[cc]>1) v_col[cc] = v_col[cc]/255.0;
//        }
//        bdt_cols.push_back(new TColor(TColor::GetFreeColorIndex(),v_col[0],v_col[1],v_col[2]) );
//
//        const char* t_scales = pBDTfile->Attribute("scale");
//        if(t_scales==NULL){
//            bdt_scales.push_back(1.0);
//        }else{
//            bdt_scales.push_back(atof(t_scales));
//        }
//
//        const char* t_weight = pBDTfile->Attribute("weight");
//        if(t_weight==NULL){
//            bdt_weight.push_back("1");
//        }else{
//            bdt_weight.push_back(t_weight);
//        }
//
//
//        const char* t_signals = pBDTfile->Attribute("signal");
//        if(t_signals==NULL){
//            bdt_is_signal.push_back(false);
//        }else{
//            std::string sig = t_signals;
//            if(sig=="true"){
//                bdt_is_signal.push_back(true);
//            }else{
//                bdt_is_signal.push_back(false);
//
//            }
//        }
//
//        const char* t_ptop = pBDTfile->Attribute("plot_on_top");
//        if(t_ptop==NULL){
//            bdt_on_top.push_back(false);
//        }else{
//            std::string sig = t_ptop;
//            if(sig=="true"){
//                bdt_on_top.push_back(true);
//            }else{
//                bdt_on_top.push_back(false);
//
//            }
//        }
//
//        const char* t_valid = pBDTfile->Attribute("validate");
//        if(t_valid==NULL){
//            bdt_is_validate_file.push_back(false);
//        }else{
//            std::string sig = t_valid;
//            if(sig=="true"){
//                bdt_is_validate_file.push_back(true);
//            }else{
//                bdt_is_validate_file.push_back(false);
//            }
//        }
//
//
//        const char* t_plotname = pBDTfile->Attribute("plot_name");
//        if(t_plotname==NULL){
//            bdt_plotnames.push_back(bdt_tags.back());
//        }else{
//            bdt_plotnames.push_back(t_plotname);
//        }
//
//		const char* t_groupname = pBDTfile->Attribute("group");//Keng
//        if(t_groupname==NULL){
//            bdt_group.push_back(-1);
//        }else{
//            bdt_group.push_back(atoi(t_groupname));
//        }
//
//
//
//        TiXmlElement *pDefinition = pBDTfile->FirstChildElement("definition");
//        std::vector<std::string> this_denom; 
//        while(pDefinition){
//            TiXmlElement *pCut = pDefinition->FirstChildElement("cut");
//            while(pCut){
//                std::string unpar =  pCut->GetText();
//                std::string parsed = this->AliasParse(unpar);
//                this_denom.push_back(parsed);
//                pCut = pCut->NextSiblingElement("cut");
//            }
//            pDefinition = pDefinition->NextSiblingElement("definition");
//        }//-->end definition
//        bdt_definitions.push_back(this_denom);
//
//        //next lets check if its the Signal Training
//        TiXmlElement *pTrain = pBDTfile->FirstChildElement("training");
//        std::vector<std::string> this_tcut; 
//        bool is_train = false;
//        while(pTrain){
//            is_train = true;
//            TiXmlElement *pTCut = pTrain->FirstChildElement("cut");
//            while(pTCut){
//                std::string unpar =  pTCut->GetText();
//                std::string parsed = this->AliasParse(unpar);
//
//                this_tcut.push_back(parsed);
//                pTCut = pTCut->NextSiblingElement("cut");
//            }
//
//            pTrain = pTrain->NextSiblingElement("training");
//        }//-->end training cuts           
//        bdt_training_cuts.push_back(this_tcut);
//        if(is_train){
//            bdt_is_training_signal.push_back(true);
//        }else{
//            bdt_is_training_signal.push_back(false);
//        }
//
//		//Add fixpot for sample comparison, CHECK
//		TiXmlElement *pfixPOT = pBDTfile->FirstChildElement("fixPOT");
//		if(pfixPOT){
//			std::string fixpot = std::string(pfixPOT->GetText());
//			bdt_fixpot.push_back(stod(fixpot));
//		}else{
//			bdt_fixpot.push_back(0);
//		}
//
//        //So, some book-keeping if its data!
//        bool is_data = false;
//        TiXmlElement *pData = pBDTfile->FirstChildElement("data");
//
//        bdt_is_onbeam_data.push_back(false);
//        bdt_is_offbeam_data.push_back(false);
//        bdt_onbeam_pot.push_back(-999);
//        bdt_offbeam_spills.push_back(-999);
//        bdt_onbeam_spills.push_back(-999);
//
//
//        while(pData){//<data > section
//
//            const char* t_use = pData->Attribute("use");
//            if(t_use=="no"){break;}
//
//            const char* t_type = pData->Attribute("type");
//            if(t_type==NULL){std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has been designated data, but not given a `type` attribute! wither OnBeam or OffBeam!! "<<std::endl; exit(EXIT_FAILURE);}
//            std::string s_type = t_type;
//            if(s_type == "OnBeam"){
//                bdt_is_onbeam_data.back() = true;
//                is_data = true;
//            }else if(s_type == "OffBeam"){
//                bdt_is_offbeam_data.back() = true;
//                is_data = true;
//            }else{
//                std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has been designated data, but its `type` attribute is neither OnBeam or OffBeam!! "<<t_type<<" "<<std::endl; exit(EXIT_FAILURE);
//            }
//            bdt_onbeam_pot.back() = 0;
//            bdt_offbeam_spills.back() = 0;
//            bdt_onbeam_spills.back() = 0;
//
//            TiXmlElement *pTor = pData->FirstChildElement("tor860_wcut");
//            while(pTor){//tor860_wcut section
//                std::string tor = std::string(pTor->GetText());
//                bdt_onbeam_pot.back() = stod(tor);    
//                pTor = pTor->NextSiblingElement("tor860_wcut");
//            }
//
//            TiXmlElement *pSpills = pData->FirstChildElement("E1DCNT_wcut");
//            while(pSpills){
//                std::string tor = std::string(pSpills->GetText());
//                bdt_onbeam_spills.back() = stod(tor);    
//                pSpills = pSpills->NextSiblingElement("E1DCNT_wcut");
//            }
//            if(bdt_onbeam_spills.back()==0 && bdt_is_offbeam_data.back()==true){
//                std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has been designated "<<t_type<<" data, but no ON beam Spills has been given for normalization! Use Zarko's tool "<<std::endl; exit(EXIT_FAILURE);
//            }
//
//
//            TiXmlElement *pExt = pData->FirstChildElement("EXT");
//            while(pExt){
//                std::string tor = std::string(pExt->GetText());
//                bdt_offbeam_spills.back() = stod(tor);    
//                pExt = pExt->NextSiblingElement("EXT");
//            }
//            if(bdt_offbeam_spills.back()==0 && bdt_is_offbeam_data.back()==true){
//                std::cerr<<"ERROR: MVALoader::MVALoader || bdt_file has been designated "<<t_type<<" data, but no OFF beam Spills has been given for normalization! Use Zarko's tool "<<std::endl; exit(EXIT_FAILURE);
//            }
//
//            pData = pData->NextSiblingElement("data");
//        }//-->end data
//
//
//        pBDTfile = pBDTfile->NextSiblingElement("bdtfile");
//        n_bdt_files++;
//
//    }//--> end bdt_file
//
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
