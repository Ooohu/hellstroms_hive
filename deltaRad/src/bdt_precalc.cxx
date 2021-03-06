#include "bdt_precalc.h"
double bragg_fnc(double *x, double *par) {
	return par[0]*pow(x[0],par[1]);
}



int addPreFriends(bdt_file* filein,std::string which){
	std::string filename = filein->tag+"_friends.root";

	std::ifstream ifile(filename.c_str());
	bool does_local_exist = (bool)ifile;
	if(does_local_exist){
		filein->addFriend( which+"_info"  , filename);
	}else{
		filein->addFriend( which+"_info"  , filein->dir+"friends/"+filename);
	}

	return 0;
}



int bdt_precalc::gen(std::string which ){

}



int bdt_precalc::genTrackInfo(){
	TFile *friend_file_out = new TFile(friend_file_out_name.c_str(),"update");
	std::cout<<"bdt_precalc::genTrackInfo() || Starting. "<<std::endl;

	file->tvertex->ResetBranchAddresses();
	//Some constants, not chosen for any huge reason
	double max_realistic_dEdx = 25;
	double min_realistic_hits = 5;
	double min_dEdx = 0.01;

	friend_file_out->cd();
	TTree * friend_tree = new TTree("track_info","track_info");


	TFile *fileconv = new TFile((bnbcorrection_dir+"/bnbcorrection/proton_conversion.root").c_str(), "read");
	TGraph * gconv = (TGraph*)fileconv->Get("Graph");

	friend_file_out->cd();


	//Some branches for some basic info
	//	double asso_r
	int reco_asso_tracks = 0;
	int reco_asso_showers = 0;
	//std::vector<double> longest_asso_track_displacement = 0;	
	//Branches for the dEdx calo only!
	std::vector<std::vector<double>> *vtrk_dEdx = 0;
	std::vector<std::vector<double>> *vtrk_resrange = 0;

	std::vector<std::vector<double>> *vtrk_X = 0;
	std::vector<std::vector<double>> *vtrk_Y = 0;
	std::vector<std::vector<double>> *vtrk_Z = 0;
	//std::vector<double>* true_track_pdg=0;

	TBranch *bdEdx = 0;
	TBranch *bresrange = 0;
	TBranch *btrackpdg = 0;

	TBranch *bX =0;
	TBranch *bY =0;
	TBranch *bZ =0;

	file->tvertex->SetBranchAddress("reco_asso_tracks", &reco_asso_tracks);
	file->tvertex->SetBranchAddress("reco_asso_showers", &reco_asso_showers);
	file->tvertex->SetBranchAddress("reco_track_calo_dEdx",&vtrk_dEdx,&bdEdx);
	file->tvertex->SetBranchAddress("reco_track_calo_resrange",&vtrk_resrange,&bresrange);
	//file->tvertex->SetBranchAddress("true_track_pdg",&true_track_pdg,&btrackpdg);

	file->tvertex->SetBranchAddress("reco_track_X",&vtrk_X,&bX);
	file->tvertex->SetBranchAddress("reco_track_Y",&vtrk_Y,&bY);
	file->tvertex->SetBranchAddress("reco_track_Z",&vtrk_Z,&bZ);


	double reco_track_startx[100], reco_track_starty[100], reco_track_startz[100];
	double reco_track_endx[100], reco_track_endy[100], reco_track_endz[100];
	file->tvertex->SetBranchAddress("reco_track_startx", &reco_track_startx);
	file->tvertex->SetBranchAddress("reco_track_starty", &reco_track_starty);
	file->tvertex->SetBranchAddress("reco_track_startz", &reco_track_startz);
	file->tvertex->SetBranchAddress("reco_track_endx", &reco_track_endx);
	file->tvertex->SetBranchAddress("reco_track_endy", &reco_track_endy);
	file->tvertex->SetBranchAddress("reco_track_endz", &reco_track_endz);


	// New branches for FRIEND TREEEE 
	double v_bragg_parA=0;
	double v_bragg_parD=0;
	int v_reco_asso_tracks = 0;

	std::vector<double> v_principal_0;	
	std::vector<double> v_principal_1;	
	std::vector<double> v_principal_2;	
	std::vector<double> v_straightness_chi2;	
	std::vector<double> v_num_spacepoints;	
	std::vector<double> v_straightness_sumangle;	

	std::vector<double> v_mean_dEdx;	
	std::vector<double> v_start_mean_dEdx;	
	std::vector<double> v_end_mean_dEdx;	
	std::vector<double> v_track_range;
	std::vector<double> v_track_PIDA;
	std::vector<double> v_track_braggA;
	std::vector<double> v_track_good_calo;
	std::vector<double> v_track_kinetic;
	std::vector<double> v_track_kinetic_from_length;

	TBranch *b_num_tracks = friend_tree->Branch("reco_asso_tracks",&v_reco_asso_tracks);


	TBranch *b_num_spacepoints = friend_tree->Branch("reco_track_num_spacepoints",&v_num_spacepoints);
	TBranch *b_straightness_chi2 = friend_tree->Branch("reco_track_straightness_chi2",&v_straightness_chi2);
	TBranch *b_principal_0 = friend_tree->Branch("reco_track_principal_0",&v_principal_0);
	TBranch *b_principal_1 = friend_tree->Branch("reco_track_principal_1",&v_principal_1);
	TBranch *b_principal_2 = friend_tree->Branch("reco_track_principal_2",&v_principal_2);
	TBranch *b_straightness_sumangle = friend_tree->Branch("reco_track_sum_angles",&v_straightness_sumangle);

	TBranch *b_mean_dEdx = friend_tree->Branch("reco_track_mean_dEdx",&v_mean_dEdx);
	TBranch *b_start_mean_dEdx = friend_tree->Branch("reco_track_end_mean_dEdx",&v_end_mean_dEdx);
	TBranch *b_end_mean_dEdx = friend_tree->Branch("reco_track_start_mean_dEdx",&v_start_mean_dEdx);
	TBranch *b_track_length = friend_tree->Branch("reco_track_range",&v_track_range);
	TBranch *b_track_PIDA = friend_tree->Branch("reco_track_PIDA",&v_track_PIDA);
	TBranch *b_bragg_parA = friend_tree->Branch("reco_track_braggA",&v_track_braggA);
	TBranch *b_kinetic = friend_tree->Branch("reco_track_kinetic",&v_track_kinetic);
	TBranch *b_good_calo = friend_tree->Branch("reco_track_good_calo",&v_track_good_calo);

	TBranch *b_kinetic_from_length = friend_tree->Branch("reco_track_kinetic_from_length",&v_track_kinetic_from_length);


	//Truncated Mean Testing
	//TFile* ftest = new TFile("test_trunc.root","recreate");
	bool done = false;
	std::vector<double> bx;
	std::vector<double> by;
	std::vector<double> ay;

	int dra =0;
	int NN = file->tvertex->GetEntries();
	for(int i=0; i< file->tvertex->GetEntries(); i++){

		file->tvertex->GetEntry(i);

		v_bragg_parA = -999;
		v_bragg_parD = -999;		


		v_reco_asso_tracks = reco_asso_tracks;
		//if(i>5000) break;

		if(reco_asso_tracks == 0){
			//welp nothing to do here
		}else if(reco_asso_tracks >0){

			for(int s=0; s<reco_asso_tracks; s++){

				if (i%2000==0)std::cout<<i<<"/"<<NN<<" "<<file->tag<<" "<<vtrk_dEdx->at(s).size()<<" "<<vtrk_resrange->at(s).size()<<" "<<reco_asso_tracks<<" "<<reco_asso_showers<<std::endl;
				//First off lets calculate straightness..etc..
				//

				std::vector<double> startpt = {vtrk_X->at(s).front(), vtrk_Y->at(s).front(), vtrk_Z->at(s).front()};
				std::vector<double> endpt = {vtrk_X->at(s).back(), vtrk_Y->at(s).back(), vtrk_Z->at(s).back()};

				double line_chi2=0;
				double sum_angle=0;				

				if(vtrk_X->at(s).size() != vtrk_Y->at(s).size() && vtrk_X->at(s).size() != vtrk_Z->at(s).size()){
					std::cout<<"ERROR: Track X,Y,Z points of different length X: "<<vtrk_X->at(s).size()<<" Y: "<<vtrk_Y->at(s).size()<<" Z: "<<vtrk_Z->at(s).size()<<std::endl;
					exit(EXIT_FAILURE);
				}

				/*TGraph2D * space_pts_3d = new TGraph2D(vtrk_X->at(s).size(), &vtrk_X->at(s).at(0), &vtrk_Y->at(s).at(0), &vtrk_Z->at(s).at(0));

				  ROOT::Fit::Fitter fitter;
				  SumDistance2 sdist(space_pts_3d);
				  ROOT::Math::Functor fcn(sdist,4);
				  double pStart[4] = {0.1,0.1,0.1,0.1};
				  fitter.SetFCN(fcn,pStart);
				  bool ok = fitter.FitFCN();
				  double fitans = -999;
				  if(ok){
				  const ROOT::Fit::FitResult & result = fitter.Result();
				  fitans = result.MinFcnValue();
				  }
				  delete space_pts_3d;
				  */


				TPrincipal* principal = new TPrincipal(3,"ND");

				for(int x = 0; x< vtrk_X->at(s).size(); x++){
					std::vector<double> tmp = {vtrk_X->at(s).at(x), vtrk_Y->at(s).at(x),vtrk_Z->at(s).at(x)};
					principal->AddRow(&tmp[0]);
				}
				principal->MakePrincipals();
				//principal->Print();
				
				TVectorD * eigen = (TVectorD*) principal->GetEigenValues();
				std::vector<double> veigen = {(*eigen)(0),(*eigen)(1),(* eigen)(2)};		

				delete principal;
				//delete eigen;


				for(int j=0; j< vtrk_X->at(s).size(); j++){


					std::vector<double> pt = {vtrk_X->at(s).at(j), vtrk_Y->at(s).at(j),vtrk_Z->at(s).at(j)};
					double dist = dist_line_point(startpt, endpt, pt);

					line_chi2+= dist*dist;

					if(j==0  || j >= vtrk_X->at(s).size()-2){
						continue;
					}else{
						std::vector<double> last_pt = {vtrk_X->at(s).at(j-1), vtrk_Y->at(s).at(j-1),vtrk_Z->at(s).at(j-1)};
						std::vector<double> next_pt = {vtrk_X->at(s).at(j+1), vtrk_Y->at(s).at(j+1),vtrk_Z->at(s).at(j+1)};

						double costh = cos_angle_3pts(last_pt, pt, next_pt);

						sum_angle += fabs(180.0-fabs(acos(costh)*180.0/3.14159) );


					}


				}

				v_num_spacepoints.push_back(vtrk_X->at(s).size());
				v_straightness_chi2.push_back(line_chi2);
				v_principal_0.push_back(veigen.at(0));
				v_principal_1.push_back(veigen.at(1));
				v_principal_2.push_back(veigen.at(2));
				v_straightness_sumangle.push_back(sum_angle);
				//				std::cout<<line_chi2<<" "<<fitans<<" "<<sum_angle<<" "<<v_num_spacepoints.back()<<std::endl;


				//And now some calorimetry things!

				double sum_dEdx =0;
				double sum_dEdx_start = 0;
				double sum_dEdx_end = 0;
				double num_start =0;
				double num_end =0;

				std::vector<float> start_dEdx;
				std::vector<float> start_resrange;

				std::vector<float> end_dEdx;
				std::vector<float> end_resrange;


				//Calculate the Total Kinetic Energy from length, assuming proton.
				double reco_track_length = sqrt(pow(reco_track_startx[s]-reco_track_endx[s], 2) + pow(reco_track_starty[s]-reco_track_endy[s], 2) + pow(reco_track_startz[s]-reco_track_endz[s], 2));
				v_track_kinetic_from_length.push_back( (gconv->Eval(reco_track_length))/1000.0);






				TruncMean tm;
				std::vector<double> trunc_dEdx;
				std::vector<double> c_resrange;
				std::vector<double> c_dEdx;

				//Hard kill any silly dEdx, Dropping first and last point
				if(vtrk_resrange->at(s).size()>1){
					for(int k=1; k < vtrk_resrange->at(s).size()-1; k++){

						bool is_sensible =vtrk_dEdx->at(s).at(k) < max_realistic_dEdx; 
						bool is_nan =vtrk_dEdx->at(s).at(k) != vtrk_dEdx->at(s).at(k); 
						bool is_inf = std::isinf(vtrk_dEdx->at(s).at(k));
						bool is_nonzero = vtrk_dEdx->at(s).at(k) > min_dEdx;

						if(is_sensible && !is_nan && !is_inf && is_nonzero){
							c_resrange.push_back(vtrk_resrange->at(s).at(k));
							c_dEdx.push_back(vtrk_dEdx->at(s).at(k));
						}
					}
				}
				//If we dont have enough hits, move to next track and tag this track as "BAD"
				if(c_dEdx.size() < min_realistic_hits ){
					v_track_range.push_back(-999);
					v_track_PIDA.push_back(-999);
					v_track_braggA.push_back(-999);
					v_start_mean_dEdx.push_back(-999);	
					v_end_mean_dEdx.push_back(-999);	
					v_mean_dEdx.push_back(-999);	
					v_track_good_calo.push_back(-999);
					v_track_kinetic.push_back(-999);
					continue;
				}

				//otherwise track is OK	
				v_track_good_calo.push_back(1);


				//Calculate appriopiate radius
				double tenth_track = std::max(c_resrange.front(), c_resrange.back())/20.0;
				double min_dx = 999;
				for(int j = c_resrange.size()-1; j>1; j--){
					double dx = fabs(c_resrange.at(j)-c_resrange.at(j-1));
					if(dx < min_dx) min_dx = dx;
				}

				double rad = std::max( min_dx*2, tenth_track); 
				tm.setRadius( rad  );


				tm.CalcTruncMeanProfile(c_resrange,c_dEdx, trunc_dEdx);			


				/*
				   if(false && trunc_dEdx.size()> 150 && trunc_dEdx.back() >  trunc_dEdx.front() && true_track_pdg->at(s) == 13 && dra < 100){
				   bx = c_resrange;
				   by = c_dEdx;
				   ay = trunc_dEdx;
				   std::vector<double> ax = bx;
				   ay.erase(ay.begin());
				   ax.erase(ax.begin());


				   done = true;
				   std::cout<<"Got Done! "<<dra<<std::endl;
				//Some Truncated Mean Testing

				TGraph * gb = new TGraph(bx.size(),&bx[0], &by[0]);
				TGraph * ga = new TGraph(ax.size(),&ax[0], &ay[0]);
				TCanvas *c = new TCanvas();
				gb->Draw("ap");
				gb->SetMarkerStyle(20);
				gb->SetMarkerSize(1.5);


				ga->Draw("l");
				ga->SetLineWidth(2);
				ga->SetLineColor(kRed-7);

				gb->GetXaxis()->SetTitle("Residual Range [cm]");
				gb->GetYaxis()->SetTitle("dE/dx [MeV/cm]");

				dra++;

				c->SaveAs(("truncpics/muon_"+std::to_string(dra)+".png").c_str(),"png");


				}
				*/

				//Check to make sure the vectors are filled in the right order. (Just for bragg variable!)
				std::vector<double> brag_c_resrange = c_resrange;	
				std::vector<double> brag_trunc_dEdx = trunc_dEdx;	

				if(brag_c_resrange.front() < brag_c_resrange.back()){
					std::reverse(brag_c_resrange.begin(),brag_c_resrange.end());
					std::reverse(brag_trunc_dEdx.begin(),brag_trunc_dEdx.end());
				}
				for(int j=0; j< brag_trunc_dEdx.size(); j++){
					//		if(brag_c_resrange.at(j) ==0) brag_c_resrange.at(j)=0.001;
				}
				while(!brag_c_resrange.empty() && brag_c_resrange.back()==0){
					brag_c_resrange.pop_back();
					brag_trunc_dEdx.pop_back();

				}

				if(brag_c_resrange.size()>1){

					std::vector<double> params = {15.0, -0.4};
					TF1 *bragg = new TF1("bragg", bragg_fnc, brag_c_resrange.front(), brag_c_resrange.back(), params.size());
					bragg->SetParameters(&params[0]);
					bragg->SetParNames("bragg_A","bragg_b");
					bragg->SetParLimits(0, 0.01, 30.0);
					//bragg->SetParLimits(1, -1.5, 1);
					bragg->FixParameter(1,-0.42);

					TGraph * pts = new TGraph(brag_trunc_dEdx.size(), &brag_c_resrange[0],&brag_trunc_dEdx[0] );
					//			TFitResultPtr fit_bragg = pts->Fit(bragg,"SQ");	
					TFitResultPtr fit_bragg = pts->Fit(bragg,"SWNQ");	
					v_bragg_parD = fit_bragg->Value(1);		
					v_bragg_parA = fit_bragg->Value(0);		


					if(!fit_bragg->IsValid()){
						std::cout<<"ERROR: fit result is not valid: status code: "<<fit_bragg->Status()<<" "<<std::endl;
						std::cout<<c_resrange.size()<<" "<<trunc_dEdx.size()<<std::endl;
						for(int m=0;m<c_resrange.size(); m++){

							std::cout<<c_resrange.at(m)<<" "<<trunc_dEdx.at(m)<<std::endl;
						}		
						std::cout<<"ERROR: Usually breaks here but for now set as bad_calo: "<<fit_bragg->Status()<<" "<<std::endl;

						//exit(EXIT_FAILURE);
						v_track_good_calo.back() = 0;
					}

					delete bragg;
					delete pts;

				}else{
					v_track_good_calo.back() = 0;
				}

				double pida_sum=0;

				for(int j=0; j<trunc_dEdx.size(); j++){

					//Some error checking
					if( trunc_dEdx.at(j) != trunc_dEdx.at(j) || std::isinf(trunc_dEdx.at(j))){
						std::cout<<"Truncated dedx is either inf or nan @ "<<j<<" "<<trunc_dEdx.at(j)<<std::endl;
						std::cout<<"Vector Length : "<<trunc_dEdx.size()<<std::endl;
						std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
						for(int m=0; m<trunc_dEdx.size(); m++){
							std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
						}
						std::cout<<"Using Radius: "<<rad<<std::endl;
						exit(EXIT_FAILURE);
					}


					if(trunc_dEdx.at(j) > 1000*max_realistic_dEdx){
						//	std::cout<<"LARGE dEdx: "<<c_resrange.at(j)<<" "<<trunc_dEdx.at(j)<<" size: "<<trunc_dEdx.size()<<" pretrunc: "<<vtrk_dEdx->at(s).at(j)<<std::endl;	
					}
					if(trunc_dEdx.at(j) < 0){
						std::cout<<"Negative dEdx: "<<c_resrange.at(j)<<" "<<trunc_dEdx.at(j)<<std::endl;	
					}


					pida_sum += trunc_dEdx.at(j)/(pow(c_resrange.at(j),-0.42));
					sum_dEdx += trunc_dEdx.at(j);

					if( j < trunc_dEdx.size()/2){
						sum_dEdx_start += trunc_dEdx.at(j);
						num_start+=1.0;
					}

					if( j >= trunc_dEdx.size()/2){
						sum_dEdx_end += trunc_dEdx.at(j);
						num_end+=1.0;
					}

				}


				double previous_resrange = c_resrange.front();
				double energy_sum = 0;
				for(int j=0; j< c_resrange.size(); j++){
					double current_resrange = 0;
					if(j+1 < c_resrange.size()) current_resrange = c_resrange.at(j+1);
					double dx = fabs(previous_resrange-current_resrange);
					energy_sum += trunc_dEdx.at(j)*dx/1000.0;	

					if(trunc_dEdx.at(j)<0) std::cout<<"Negative dEdx!! "<<trunc_dEdx.at(j)<<std::endl;
					if(dx<0) std::cout<<"Negative dx!! "<<dx<<" "<<j<<" "<<previous_resrange<<" "<<current_resrange<<std::endl;

					previous_resrange = current_resrange;	
				}	

				v_track_range.push_back(std::max(vtrk_resrange->at(s).front(), vtrk_resrange->at(s).back()));
				v_track_kinetic.push_back(energy_sum);
				v_track_braggA.push_back(v_bragg_parA);
				v_track_PIDA.push_back(pida_sum/(double)c_resrange.size());
				v_mean_dEdx.push_back(sum_dEdx/(double)c_resrange.size());
				v_start_mean_dEdx.push_back(sum_dEdx_start/num_start);
				v_end_mean_dEdx.push_back(sum_dEdx_end/num_end);
			
			
			}

		}

		// OK now some error checking!..
		for(int j=0; j<v_mean_dEdx.size();j++){
			if(v_mean_dEdx.at(j) != v_mean_dEdx.at(j)) 	{
				std::cout<<"Yarp: "<<" "<<v_mean_dEdx.at(j)<<" "<<v_start_mean_dEdx.at(j)<<std::endl;
				v_mean_dEdx.at(j) =-999;
			}
			if(v_start_mean_dEdx.at(j) != v_start_mean_dEdx.at(j)) v_start_mean_dEdx.at(j) =-999;
		}


		friend_tree->Fill();
		v_track_range.clear();
		v_track_kinetic.clear();
		v_track_kinetic_from_length.clear();
		v_track_PIDA.clear();
		v_track_braggA.clear();
		v_start_mean_dEdx.clear();	
		v_end_mean_dEdx.clear();	
		v_mean_dEdx.clear();	
		v_num_spacepoints.clear();	
		v_straightness_sumangle.clear();	
		v_straightness_chi2.clear();	
		v_principal_1.clear();	
		v_principal_0.clear();	
		v_principal_2.clear();	
		v_track_good_calo.clear();

	}
	friend_file_out->cd();
	friend_tree->Write();

	delete friend_file_out;
	//Some Truncated Mean Testing
	/*	ftest->cd();
		TGraph * gb = new TGraph(bx.size(),&bx[0], &by[0]);
		TGraph * ga = new TGraph(bx.size(),&bx[0], &ay[0]);
		TCanvas *c = new TCanvas();
		gb->Draw("ap");
		gb->SetMarkerStyle(20);
		gb->SetMarkerSize(2);

		ga->Draw("l");
		ga->SetLineColor(kRed-7);

		c->Write();


		delete gb;
		delete ga;
		delete c;
		*/
	return 0;
}


int bdt_precalc::genShowerInfo(){
	std::cout<<"bdt_precalc::genShowerInfo() || Starting. "<<std::endl;
	TFile *friend_file_out = new TFile(friend_file_out_name.c_str(),"update");
	TTree * friend_tree = new TTree("shower_info","shower_info");


	file->tvertex->ResetBranchAddresses();

	//Some branches for some basic info
	//	double asso_r
	int reco_asso_tracks = 0;
	int reco_asso_showers = 0;

	double reco_shower_dedx_plane1[100];
	double reco_shower_dedx_plane2[100];
	double reco_shower_dedx_plane0[100];

	double reco_shower_dedxnew_plane1[100];
	double reco_shower_dedxnew_plane2[100];
	double reco_shower_dedxnew_plane0[100];

	double reco_shower_dirx[100];
	double reco_shower_diry[100];
	double reco_shower_dirz[100];

	file->tvertex->SetBranchAddress("reco_asso_tracks", &reco_asso_tracks);
	file->tvertex->SetBranchAddress("reco_asso_tracks", &reco_asso_tracks);
	file->tvertex->SetBranchAddress("reco_asso_showers", &reco_asso_showers);

	file->tvertex->SetBranchAddress("reco_shower_dirx", reco_shower_dirx);
	file->tvertex->SetBranchAddress("reco_shower_diry", reco_shower_diry);
	file->tvertex->SetBranchAddress("reco_shower_dirz", reco_shower_dirz);

	file->tvertex->SetBranchAddress("reco_shower_dedx_plane1", reco_shower_dedx_plane1);
	file->tvertex->SetBranchAddress("reco_shower_dedx_plane2", reco_shower_dedx_plane2);
	file->tvertex->SetBranchAddress("reco_shower_dedx_plane0", reco_shower_dedx_plane0);

	file->tvertex->SetBranchAddress("reco_shower_dedxnew_plane1", reco_shower_dedxnew_plane1);
	file->tvertex->SetBranchAddress("reco_shower_dedxnew_plane0", reco_shower_dedxnew_plane0);
	file->tvertex->SetBranchAddress("reco_shower_dedxnew_plane2", reco_shower_dedxnew_plane2);

	// New branches for FRIEND TREEEE 
	std::vector<double> vec_amal_dEdx;

	TBranch *b_amal_dEdx = friend_tree->Branch("amalgamated_shower_dEdx",&vec_amal_dEdx);

	TH2D* h2 = new TH2D("h2","h2",50,0,6,50,0,90);
	TH2D * h1 = new TH2D("h1","h1",50,0,6,50,-180,180);
	TH2D * h0 = new TH2D("h0","h0",50,0,6,50,-180,180);

	TH2D* hnew2 = new TH2D("hnew2","hnew2",50,0,6,50,0,90);
	TH2D * hnew1 = new TH2D("hnew1","hnew1",50,0,6,50,-180,180);
	TH2D * hnew0 = new TH2D("hnew0","hnew0",50,0,6,50,-180,180);


	vec_amal_dEdx.clear();
	int NN = file->tvertex->GetEntries();
	for(int i=0; i< NN; i++){
		file->tvertex->GetEntry(i);
		if (i%10000==0)std::cout<<i<<"/"<<NN<<" "<<file->tag<<" "<<std::endl;

		for(int s=0; s< reco_asso_showers; s++){		
			double dedx = -999;
			double shower_theta = atan2(reco_shower_diry[s],reco_shower_dirz[s])*180.0/3.14159;
			//First lets check the collection plane:
			//is it with 30 degrees of wire?

			double mina2 = std::min( fabs(shower_theta-90), fabs(shower_theta+90));
			double mina0 = std::min( fabs(shower_theta-30), fabs(shower_theta+150));
			double mina1 =std::min( fabs(shower_theta-150), fabs(shower_theta+30));


			if( mina2 > 35 && ( reco_shower_dedxnew_plane2[s] > 0 || reco_shower_dedx_plane2[s]>0) ){
				dedx = reco_shower_dedxnew_plane2[s];
				if(dedx < 0) dedx = reco_shower_dedx_plane2[s];

			}else{

				if( mina1 > mina2 && ( reco_shower_dedxnew_plane1[s] > 0 || reco_shower_dedx_plane1[s] > 0) ){ // use plane1, or V plane
					dedx =  reco_shower_dedxnew_plane1[s];

					if(dedx < 0) dedx = reco_shower_dedx_plane1[s];



				}else {//use plane0 or U plane
					dedx = reco_shower_dedxnew_plane0[s];

					if(dedx < 0) dedx = reco_shower_dedx_plane0[s];




				}
			}

			if(dedx < 0){
				dedx = std::max(std::max( reco_shower_dedxnew_plane2[s], reco_shower_dedx_plane2[s]), std::max(std::max(reco_shower_dedxnew_plane1[s], reco_shower_dedx_plane1[s]), std::max(reco_shower_dedxnew_plane0[s], reco_shower_dedx_plane0[s]) ) );

			}

			h0->Fill(reco_shower_dedx_plane0[s], mina0);
			h1->Fill(reco_shower_dedx_plane1[s], mina1);
			h2->Fill(reco_shower_dedx_plane2[s], mina2);

			hnew0->Fill(reco_shower_dedxnew_plane0[s], mina0);
			hnew1->Fill(reco_shower_dedxnew_plane1[s], mina1);
			hnew2->Fill(reco_shower_dedxnew_plane2[s], mina2);

			vec_amal_dEdx.push_back(dedx);
			if(dedx < 0) std::cout<<dedx<<" "<<reco_shower_dedx_plane0[s]<<" "<<reco_shower_dedx_plane1[s]<<" "<<reco_shower_dedx_plane2[s]<<std::endl;

		}
		friend_tree->Fill();

		vec_amal_dEdx.clear();
	}


	friend_file_out->cd();
	friend_tree->Write();

	TFile *ftest = new TFile("shower.root","recreate");
	ftest->cd();
	h2->Write();
	h0->Write();
	h1->Write();

	hnew2->Write();
	hnew0->Write();
	hnew1->Write();

	ftest->Close();

	delete friend_file_out;
	return 0;
}



int bdt_precalc::genPi0Info(){
	std::cout<<"bdt_precalc::genPi0Info() || Starting. "<<std::endl;
	TFile *friend_file_out = new TFile(friend_file_out_name.c_str(),"update");

	file->tvertex->ResetBranchAddresses();
	friend_file_out->cd();
	TTree * friend_tree = new TTree("pi0_info","pi0_info");

	std::cout<<"bdt_precalc: genPi0Info"<<std::endl;

	//Some branches for some basic info
	//	double asso_r
	int reco_asso_tracks = 0;
	int reco_asso_showers = 0;
	int pi0_class_number = 0;
	//std::vector<double> longest_asso_track_displacement = 0;	
	std::vector<double> *vall_reco_tracks= 0;
	std::vector<double> *vall_reco_showers= 0;
	std::vector<double> *vall_bp_reco_tracks= 0;
	std::vector<double> *vall_bp_reco_showers= 0;


	TBranch *ballt = 0;
	TBranch *balls = 0;

	TBranch *ballbpt = 0;
	TBranch *ballbps = 0;

	file->tvertex->SetBranchAddress("reco_asso_tracks", &reco_asso_tracks);
	file->tvertex->SetBranchAddress("reco_asso_showers", &reco_asso_showers);
	file->tvertex->SetBranchAddress("all_reco_tracks_dist_from_vertex",&vall_reco_tracks,&ballt);
	file->tvertex->SetBranchAddress("all_reco_showers_dist_from_vertex",&vall_reco_showers,&balls);

	file->tvertex->SetBranchAddress("all_reco_tracks_bp_dist_from_vertex",&vall_bp_reco_tracks,&ballbpt);
	file->tvertex->SetBranchAddress("all_reco_showers_bp_dist_from_vertex",&vall_bp_reco_showers,&ballbps);

	file->tvertex->SetBranchAddress("pi0_class_number",&pi0_class_number);
	// New branches for FRIEND TREEEE 


	double v_pi0_class_number=0;
	std::vector<int> vec_reco_showers_within;
	std::vector<int> vec_reco_tracks_within;


	int reco_showers_bp_within_10 = 0;
	int reco_showers_bp_within_20 = 0;
	int reco_showers_bp_within_30 = 0;

	int reco_tracks_bp_within_10 = 0;
	int reco_tracks_bp_within_20 = 0;
	int reco_tracks_bp_within_30 = 0;

	int collated_showers = 0;
	int collated_tracks = 0;

	int reco_showers_within_10 = 0;
	int reco_showers_within_20 = 0;
	int reco_showers_within_30 = 0;

	int reco_tracks_within_10 = 0;
	int reco_tracks_within_20 = 0;
	int reco_tracks_within_30 = 0;

	TBranch *b_vec_showers = friend_tree->Branch("num_reco_showers_within_Xcm_vertex",&vec_reco_showers_within);
	TBranch *b_vec_tracks = friend_tree->Branch("num_reco_tracks_within_Xcm_vertex",&vec_reco_tracks_within);

	TBranch *b_reco_showers_within_10 = friend_tree->Branch("num_reco_showers_within_10cm_vertex",&reco_showers_within_10);
	TBranch *b_reco_showers_within_20 = friend_tree->Branch("num_reco_showers_within_20cm_vertex",&reco_showers_within_20);
	TBranch *b_reco_showers_within_30 = friend_tree->Branch("num_reco_showers_within_30cm_vertex",&reco_showers_within_30);
	TBranch *b_reco_tracks_within_10 = friend_tree->Branch("num_reco_tracks_within_10cm_vertex",&reco_tracks_within_10);
	TBranch *b_reco_tracks_within_20 = friend_tree->Branch("num_reco_tracks_within_20cm_vertex",&reco_tracks_within_20);
	TBranch *b_reco_tracks_within_30 = friend_tree->Branch("num_reco_tracks_within_30cm_vertex",&reco_tracks_within_30);

	TBranch *b_reco_showers_bp_within_10 = friend_tree->Branch("num_reco_showers_bp_within_10cm_vertex",&reco_showers_bp_within_10);
	TBranch *b_reco_showers_bp_within_20 = friend_tree->Branch("num_reco_showers_bp_within_20cm_vertex",&reco_showers_bp_within_20);
	TBranch *b_reco_showers_bp_within_30 = friend_tree->Branch("num_reco_showers_bp_within_30cm_vertex",&reco_showers_bp_within_30);
	TBranch *b_reco_tracks_bp_within_10 = friend_tree->Branch("num_reco_tracks_bp_within_10cm_vertex",&reco_tracks_bp_within_10);
	TBranch *b_reco_tracks_bp_within_20 = friend_tree->Branch("num_reco_tracks_bp_within_20cm_vertex",&reco_tracks_bp_within_20);
	TBranch *b_reco_tracks_bp_within_30 = friend_tree->Branch("num_reco_tracks_bp_within_30cm_vertex",&reco_tracks_bp_within_30);
	

	TBranch *b_reco_collatedt = friend_tree->Branch("num_other_asso_tracks",&collated_tracks);
	TBranch *b_reco_collateds = friend_tree->Branch("num_other_asso_showers",&collated_showers);


	TBranch *b_pi0_class_number = friend_tree->Branch("pi0_class_number",&v_pi0_class_number);

	int NN = file->tvertex->GetEntries();
	for(int i=0; i< file->tvertex->GetEntries(); i++){

		if (i%10000==0)std::cout<<i<<"/"<<NN<<" "<<file->tag<<" "<<" piClassInfo "<<std::endl;
		reco_showers_within_10 = 0;
		reco_showers_within_20 = 0;
		reco_showers_within_30 = 0;

		reco_tracks_within_10 = 0;
		reco_tracks_within_20 = 0;
		reco_tracks_within_30 = 0;

		collated_showers = 0;
		collated_tracks = 0;

		reco_showers_bp_within_10 = 0;
		reco_showers_bp_within_20 = 0;
		reco_showers_bp_within_30 = 0;

		reco_tracks_bp_within_10 = 0;
		reco_tracks_bp_within_20 = 0;
		reco_tracks_bp_within_30 = 0;


		file->tvertex->GetEntry(i);
		v_pi0_class_number = pi0_class_number;




		for(double c=1; c<32.0; c=c+2.0){
			vec_reco_showers_within.push_back(0);
			vec_reco_tracks_within.push_back(0);


			for(int j=0; j< vall_reco_tracks->size(); j++){

				if(vall_reco_tracks->at(j)<c){
					vec_reco_tracks_within.back()++;

					
				}
			}


			for(int j=0; j< vall_reco_showers->size(); j++){

				if(vall_reco_showers->at(j)<c){
					vec_reco_showers_within.back()++;
				}	
			}
		}


		if(vall_bp_reco_tracks->size() != vall_reco_tracks->size()){
			std::cout<<"Error in track: "<<vall_bp_reco_tracks->size()<<" != "<<vall_reco_tracks->size()<<std::endl;
		}	

		if(vall_bp_reco_showers->size() != vall_reco_showers->size()){
			std::cout<<"Error in shower: "<<vall_bp_reco_showers->size()<<" != "<<vall_reco_showers->size()<<std::endl;
		}	

		for(int j=0; j< vall_bp_reco_showers->size(); j++){

			if(vall_bp_reco_showers->at(j) < 10){
				reco_showers_bp_within_10 ++;
				reco_showers_bp_within_20 ++;
				reco_showers_bp_within_30 ++;
			}else if(vall_bp_reco_showers->at(j) < 20){
				reco_showers_bp_within_20 ++;
				reco_showers_bp_within_30 ++;
			}else if(vall_bp_reco_showers->at(j) < 30){
				reco_showers_bp_within_30 ++;
			}	
		}

		for(int j=0; j< vall_reco_tracks->size(); j++){

			if(vall_reco_tracks->at(j) < 10){
				reco_tracks_within_10 ++;
				reco_tracks_within_20 ++;
				reco_tracks_within_30 ++;
			}else if(vall_reco_tracks->at(j) < 20){
				reco_tracks_within_20 ++;
				reco_tracks_within_30 ++;
			}else if(vall_reco_tracks->at(j) < 30){
				reco_tracks_within_30 ++;
			}	
		}



		for(int j=0; j< vall_bp_reco_showers->size(); j++){
			if(vall_bp_reco_showers->at(j) < 10){
				reco_showers_bp_within_10 ++;
				reco_showers_bp_within_20 ++;
				reco_showers_bp_within_30 ++;


				if(vall_reco_showers->at(j) < 30){
					collated_showers++;
				}

			}else if(vall_bp_reco_showers->at(j) < 20){
				reco_showers_bp_within_20 ++;
				reco_showers_bp_within_30 ++;
			}else if(vall_bp_reco_showers->at(j) < 30){
				reco_showers_bp_within_30 ++;


			}	
		}


		for(int j=0; j< vall_reco_tracks->size(); j++){
			if(vall_reco_tracks->at(j) < 10){
				reco_tracks_within_10 ++;
				reco_tracks_within_20 ++;
				reco_tracks_within_30 ++;

				if(vall_reco_tracks->at(j) < 30){
					collated_tracks++;
				}


			}else if(vall_reco_tracks->at(j) < 20){
				reco_tracks_within_20 ++;
				reco_tracks_within_30 ++;
			}else if(vall_reco_tracks->at(j) < 30){
				reco_tracks_within_30 ++;
			}	
		}

		friend_tree->Fill();

		vec_reco_showers_within.clear();
		vec_reco_tracks_within.clear();

	}
	friend_file_out->cd();
	friend_tree->Write();


	std::cout<<"Done with Pi0ClassInfo"<<std::endl;
	delete friend_file_out;
	return 0;
}



int bdt_precalc::genBNBcorrectionInfo(){
	std::cout<<"bdt_precalc:: Starting BNBcorrectionInfo()"<<std::endl;
	TFile *friend_file_out = new TFile(friend_file_out_name.c_str(),"update");

	file->tvertex->ResetBranchAddresses();

	//first Load up all the flux files
	TFile *f_old = new TFile((bnbcorrection_dir+"/bnbcorrection/bnb_oldflux_volAVTPC.root").c_str(),"read");
	TH1F* h_mu_nue_old = (TH1F*)f_old->Get("h_mu_nue");
	TH1F* h_mu_nuebar_old = (TH1F*)f_old->Get("h_mu_nuebar");
	TH1F* h_mu_numu_old = (TH1F*)f_old->Get("h_mu_numu");
	TH1F* h_mu_numubar_old = (TH1F*)f_old->Get("h_mu_numubar");

	TH1F* h_pi_nue_old = (TH1F*)f_old->Get("h_pi_nue");
	TH1F* h_pi_nuebar_old = (TH1F*)f_old->Get("h_pi_nuebar");
	TH1F* h_pi_numu_old = (TH1F*)f_old->Get("h_pi_numu");
	TH1F* h_pi_numubar_old = (TH1F*)f_old->Get("h_pi_numubar");

	TH1F* h_k0_nue_old = (TH1F*)f_old->Get("h_k0_nue");
	TH1F* h_k0_nuebar_old = (TH1F*)f_old->Get("h_k0_nuebar");
	TH1F* h_k0_numu_old = (TH1F*)f_old->Get("h_k0_numu");
	TH1F* h_k0_numubar_old = (TH1F*)f_old->Get("h_k0_numubar");

	TH1F* h_k_nue_old = (TH1F*)f_old->Get("h_k_nue");
	TH1F* h_k_nuebar_old = (TH1F*)f_old->Get("h_k_nuebar");
	TH1F* h_k_numu_old = (TH1F*)f_old->Get("h_k_numu");
	TH1F* h_k_numubar_old = (TH1F*)f_old->Get("h_k_numubar");

	TH1F* h_nue_old = (TH1F*)h_mu_nue_old->Clone("h_nue_old");
	h_nue_old->Add(h_pi_nue_old);
	h_nue_old->Add(h_k0_nue_old);
	h_nue_old->Add(h_k_nue_old);

	TH1F* h_nuebar_old = (TH1F*)h_mu_nuebar_old->Clone("h_nuebar_old");
	h_nuebar_old->Add(h_pi_nuebar_old);
	h_nuebar_old->Add(h_k0_nuebar_old);
	h_nuebar_old->Add(h_k_nuebar_old);

	TH1F* h_numu_old = (TH1F*)h_mu_numu_old->Clone("h_numu_old");
	h_numu_old->Add(h_pi_numu_old);
	h_numu_old->Add(h_k0_numu_old);
	h_numu_old->Add(h_k_numu_old);

	TH1F* h_numubar_old = (TH1F*)h_mu_numubar_old->Clone("h_numubar_old");
	h_numubar_old->Add(h_pi_numubar_old);
	h_numubar_old->Add(h_k0_numubar_old);
	h_numubar_old->Add(h_k_numubar_old);


	TFile *f_new = new TFile((bnbcorrection_dir+"/bnbcorrection/bnb_newflux_volAVTPC.root").c_str(),"read");
	TH1F* h_mu_nue_new = (TH1F*)f_new->Get("h_mu_nue");
	TH1F* h_mu_nuebar_new = (TH1F*)f_new->Get("h_mu_nuebar");
	TH1F* h_mu_numu_new = (TH1F*)f_new->Get("h_mu_numu");
	TH1F* h_mu_numubar_new = (TH1F*)f_new->Get("h_mu_numubar");

	TH1F* h_pi_nue_new = (TH1F*)f_new->Get("h_pi_nue");
	TH1F* h_pi_nuebar_new = (TH1F*)f_new->Get("h_pi_nuebar");
	TH1F* h_pi_numu_new = (TH1F*)f_new->Get("h_pi_numu");
	TH1F* h_pi_numubar_new = (TH1F*)f_new->Get("h_pi_numubar");

	TH1F* h_k0_nue_new = (TH1F*)f_new->Get("h_k0_nue");
	TH1F* h_k0_nuebar_new = (TH1F*)f_new->Get("h_k0_nuebar");
	TH1F* h_k0_numu_new = (TH1F*)f_new->Get("h_k0_numu");
	TH1F* h_k0_numubar_new = (TH1F*)f_new->Get("h_k0_numubar");

	TH1F* h_k_nue_new = (TH1F*)f_new->Get("h_k_nue");
	TH1F* h_k_nuebar_new = (TH1F*)f_new->Get("h_k_nuebar");
	TH1F* h_k_numu_new = (TH1F*)f_new->Get("h_k_numu");
	TH1F* h_k_numubar_new = (TH1F*)f_new->Get("h_k_numubar");

	TH1F* h_nue_new = (TH1F*)h_mu_nue_new->Clone("h_nue_new");
	h_nue_new->Add(h_pi_nue_new);
	h_nue_new->Add(h_k0_nue_new);
	h_nue_new->Add(h_k_nue_new);

	TH1F* h_nuebar_new = (TH1F*)h_mu_nuebar_new->Clone("h_nuebar_new");
	h_nuebar_new->Add(h_pi_nuebar_new);
	h_nuebar_new->Add(h_k0_nuebar_new);
	h_nuebar_new->Add(h_k_nuebar_new);

	TH1F* h_numu_new = (TH1F*)h_mu_numu_new->Clone("h_numu_new");
	h_numu_new->Add(h_pi_numu_new);
	h_numu_new->Add(h_k0_numu_new);
	h_numu_new->Add(h_k_numu_new);

	TH1F* h_numubar_new = (TH1F*)h_mu_numubar_new->Clone("h_numubar_new");
	h_numubar_new->Add(h_pi_numubar_new);
	h_numubar_new->Add(h_k0_numubar_new);
	h_numubar_new->Add(h_k_numubar_new);



	TH1F * h_nue_ratio = (TH1F*)h_nue_new->Clone("h_nue_ratio");
	h_nue_ratio->Divide(h_nue_old);

	TH1F * h_nuebar_ratio = (TH1F*)h_nuebar_new->Clone("h_nuebar_ratio");
	h_nuebar_ratio->Divide(h_nuebar_old);

	TH1F * h_numu_ratio = (TH1F*)h_numu_new->Clone("h_numu_ratio");
	h_numu_ratio->Divide(h_numu_old);

	TH1F * h_numubar_ratio = (TH1F*)h_numubar_new->Clone("h_numubar_ratio");
	h_numubar_ratio->Divide(h_numubar_old);




	file->f->cd();
	int nu_pdg = 0;
	double true_nu_E = 0;
	std::cout<<"1"<<std::endl;
	file->tvertex->SetBranchAddress("nu_pdg", &nu_pdg);
	std::cout<<"2"<<std::endl;
	file->tvertex->SetBranchAddress("true_nu_E", &true_nu_E);

	// New branches for FRIEND TREEEE 
	std::cout<<"Set up new bnbcorrection_info TTree"<<std::endl;
	friend_file_out->cd();
	TTree * friend_tree = new TTree("bnbcorrection_info","bnbcorrection_info");


	double weight=0;
	TBranch *b_weight = friend_tree->Branch("weight",&weight);



	int NN = file->tvertex->GetEntries();
	for(int i=0; i< NN; i++){
		file->tvertex->GetEntry(i);
		if (i%10000==0)std::cout<<i<<"/"<<NN<<" "<<file->tag<<" "<<std::endl;

		weight = 1.0;



		if(nu_pdg == 12){
			int thisbin = h_nue_ratio->GetXaxis()->FindBin(true_nu_E);	
			weight = h_nue_ratio->GetBinContent(thisbin);
		}else if(nu_pdg == -12){
			int thisbin = h_nuebar_ratio->GetXaxis()->FindBin(true_nu_E);	
			weight = h_nuebar_ratio->GetBinContent(thisbin);
		}else if(nu_pdg == 14){
			int thisbin = h_numu_ratio->GetXaxis()->FindBin(true_nu_E);	
			weight = h_numu_ratio->GetBinContent(thisbin);
		}else if(nu_pdg == -14){
			int thisbin = h_numubar_ratio->GetXaxis()->FindBin(true_nu_E);	
			weight = h_numubar_ratio->GetBinContent(thisbin);
		}else{
			weight = 1.0;
		}



		friend_tree->Fill();


	}

	friend_file_out->cd();
	friend_tree->Write();



	/*	TFile *f_bnbcorr = new TFile("bnbcorr.root","recreate");
		f_bnbcorr->cd();
		h_nue_ratio->Write();
		h_nuebar_ratio->Write();
		h_numu_ratio->Write();
		h_numubar_ratio->Write();
		f_bnbcorr->Close();
		*/

	//delete friend_tree;
	delete friend_file_out;
	std::cout<<"Closed file"<<std::endl;
	//	delete friend_tree;
	//	delete friend_file_out;
	std::cout<<"Done with bnbcorrection_info"<<std::endl;
	return 0;
}



//line between X1 and X2, point X0
double dist_line_point( std::vector<double>X1, std::vector<double> X2, std::vector<double> X0){
	double x1 =X1.at(0);
	double y1 =X1.at(1);
	double z1 =X1.at(2);

	double x2 =X2.at(0);
	double y2 =X2.at(1);
	double z2 =X2.at(2);

	double x0 =X0.at(0);
	double y0 =X0.at(1);
	double z0 =X0.at(2);

	double x10 = x1-x0;
	double y10 = y1-y0;
	double z10 = z1-z0;

	double x21 = x2-x1;
	double y21 = y2-y1;
	double z21 = z2-z1;

	double t = -(x10*x21+y10*y21+z10*z21)/fabs(x21*x21+y21*y21+z21*z21 );

	double d2 = pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0)+(z2-z1)*(z1-z0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2));


	return sqrt(d2);

}

double cos_angle_3pts(std::vector<double> last, std::vector<double> next, std::vector<double> mid){
	double lx = mid.at(0)-last.at(0);		
	double ly = mid.at(1)-last.at(1);		
	double lz = mid.at(2)-last.at(2);		

	double nx = mid.at(0)-next.at(0);		
	double ny = mid.at(1)-next.at(1);		
	double nz = mid.at(2)-next.at(2);		

	double cos = (lx*nx+ly*ny+lz*nz)/(sqrt(lx*lx+ly*ly+lz*lz)*sqrt(nx*nx+ny*ny+nz*nz));

	return cos;
}




