#include "bdt_eff.h"

/*
bdt_efficiency::bdt_efficiency(bdt_file* filein, std::string denomin, double c1, double c2) : file(filein), denominator(denomin){


	//First step, find event entrylist. In future we must actually track this from the event_tree
	int  event_number=0;
	int  run_number=0;
	int  subrun_number=0;
	double wei=0;
	int pass_denom = 0;

	double weighted_num = 0;
	double base_num = 0;

	file->tvertex->ResetBranchAddresses();

	file->tvertex->SetBranchAddress("event_number",&event_number);
	file->tvertex->SetBranchAddress("bnbcorrection_info.weight",&wei);

	//TPCActive
	//Z: 0 to 1036.8
	//X: 0 to 256.35
	//Y: -116.5 to 116.5

	std::set<int> eventIDs;
	event_entry_list = new TEntryList(file->tvertex);

	for(int k=0; k< file->tvertex->GetEntries();k++){
		file->tvertex->GetEntry(k);
		if(k%10000==0) std::cout<<k<<std::endl;	
		if(eventIDs.count(event_number)==0){
			eventIDs.insert(event_number);
			event_entry_list->Enter(k,file->tvertex);
			weighted_num+=wei;
			base_num++;
		}	
	}

	double MOD =file->scale_data*6.6e20/file->pot;
	double volCryo = 199668.427885;
	double volTPC =  101510.0;
	double volTPCActive=  86698.6;

	file->tvertex->SetEntryList(event_entry_list);


	double Ndenominator = file->GetEntries(denominator.c_str());

	std::cout<<"====================Raw Numbers of Events==================="<<std::endl;
	std::cout<<"1: Number of events in FillLightFiles : "<<file->numberofevents_raw<<std::endl;
	std::cout<<"2: Number of events in FillLightFiles scaled to TPCActive: "<<file->numberofevents_raw*volTPCActive/volTPC<<std::endl;
	std::cout<<"3: Number of unique events in vertexed_files : "<<base_num<<std::endl;
	std::cout<<"4: Number of unique events in vertexed_files with BNB_correction : "<<weighted_num<<std::endl;
	std::cout<<"5: Number of approximated events in FillLightFiles with BNB_correction : "<<file->numberofevents_raw*weighted_num/base_num<<std::endl;
	std::cout<<"6: Same scaled to TPCActive: "<<file->numberofevents_raw*weighted_num/base_num*volTPCActive/volTPC<<std::endl;
	double vertex_eff = (weighted_num)/(file->numberofevents_raw*weighted_num/base_num*volTPCActive/volTPC);
	std::cout<<"7: Ratio of (6) and (4): This is a measure of vertexing efficiency : "<<vertex_eff<<std::endl;

	std::cout<<"==================== Denominator events (Here on x3.1) ==================="<<std::endl;
	std::cout<<"Number of events in vertexed_files that pass denominator with BNB_correction : "<<Ndenominator<<std::endl;
	std::cout<<"Number of events in vertexed_files that pass denominator with BNB_correction with vertex eff removed : "<<Ndenominator/vertex_eff<<std::endl;
	double true_denominator = Ndenominator/vertex_eff*MOD;
	std::cout<<"Same scaled to 6.6e20 "<<true_denominator<<std::endl; 



	double nverticies;
	for(int s=0; s<4; s++){
		if(s==2) file->calcCosmicBDTEntryList(c1, c2);
		if(s==3) file->calcBNBBDTEntryList(c1, c2);
		file->setStageEntryList(s);

		double stage_entries = file->GetEntries("1")*MOD;
		if(s==0) nverticies = stage_entries;
		std::cout<<"Stage: "<<s<<" Verticies: "<<stage_entries<<" Efficiency: "<<stage_entries/true_denominator<<std::endl;
	}

	std::cout<<"==================== Precuts - Individual  ==================="<<std::endl;
	file->setStageEntryList(0);

	for(int m=0; m< file->flow.vec_pre_cuts.size(); m++){
		std::string tmpcut = file->flow.vec_pre_cuts.at(m);
		double np = file->GetEntries(tmpcut.c_str())*MOD;
		std::cout<<" + "<<file->flow.vec_pre_cuts.at(m)<<"\t||\t"<<np<<"\t("<<np/nverticies*100<<")\%"<<std::endl;
	}




	std::cout<<"==================== Precuts - One by One  ==================="<<std::endl;
	file->setStageEntryList(0);

	std::string thiscut = "1";
	for(int m=0; m< file->flow.vec_pre_cuts.size(); m++){
		thiscut += "&&"+ file->flow.vec_pre_cuts.at(m);
		double np = file->GetEntries(thiscut.c_str())*MOD;
		std::cout<<" + "<<file->flow.vec_pre_cuts.at(m)<<"\t||\t"<<np<<"\t("<<np/nverticies*100<<")\%"<<std::endl;
	}




}
*/


//THIS IS WITH A BOOL.
bdt_efficiency::bdt_efficiency(bdt_file* filein, std::string denomin, double c1, double c2,bool in) : file(filein), denominator(denomin){
	
	TFile * feff = new TFile("eff.1g1p.root","recreate");
	//First step, find event entrylist. In future we must actually track this from the event_tree
	int  event_number=0;
	int  run_number=0;
	int  subrun_number=0;
	double wei=0;
	int pass_denom = 0;

	double weighted_num = 0;
	double base_num = 0;

	file->tvertex->ResetBranchAddresses();

	file->tvertex->SetBranchAddress("event_number",&event_number);

	if(filein->name.compare(9,6,"bnbext",0,6) != 0){
	//the following is not for BNBext
	file->tvertex->SetBranchAddress("bnbcorrection_info.weight",&wei);
	}

	//TPCActive
	//Z: 0 to 1036.8
	//X: 0 to 256.35
	//Y: -116.5 to 116.5

	std::set<int> eventIDs;
	event_entry_list = new TEntryList(file->tvertex);
	
	int total_number = file->tvertex->GetEntries();

	for(int k=0; k< total_number; k++){
		file->tvertex->GetEntry(k);

		//if(k%10000==0) std::cout<<k<<std::endl;	
		if( k%(total_number/100) == 0){
		    std::cout.precision(3);
		    std::cout<<"\r"<< 100*k/total_number+1 <<"\% of "<< total_number <<" entries complete.";//+1 is for the display 100%
		    std::cout.flush();
		}

		if(eventIDs.count(event_number)==0){
			eventIDs.insert(event_number);
			event_entry_list->Enter(k,file->tvertex);
	    
	    //bnb ext for special treatment: ()? BNBEXT TAKE THIS : non BNBEXT TAKE THIS!
	(filein->name.compare(9,6,"bnbext",0,6) == 0)? weighted_num+= 1 : weighted_num+=wei;

			base_num++;
		}	
	}

	double MOD =file->scale_data*6.6e20/file->pot;
	double volCryo = 199668.427885;
	double volTPC =  101510.0;
	double volTPCActive=  86698.6;

	file->tvertex->SetEntryList(event_entry_list);


	double Ndenominator = file->GetEntries(denominator.c_str());

	std::cout<<"====================Raw Numbers of Events==================="<<std::endl;
	std::cout<<"1: Number of events in FillLightFiles : "<<file->numberofevents_raw<<std::endl;
	std::cout<<"2: Number of events in FillLightFiles scaled to TPCActive: "<<file->numberofevents_raw*volTPCActive/volTPC<<std::endl;
	std::cout<<"3: Number of unique events in vertexed_files : "<<base_num<<std::endl;
	std::cout<<"4: Number of unique events in vertexed_files with BNB_correction : "<<weighted_num<<std::endl;
	std::cout<<"5: Number of approximated events in FillLightFiles with BNB_correction : "<<file->numberofevents_raw*weighted_num/base_num<<std::endl;
	std::cout<<"6: Same scaled to TPCActive: "<<file->numberofevents_raw*weighted_num/base_num*volTPCActive/volTPC<<std::endl;

	double vertex_eff = (weighted_num)/(file->numberofevents_raw*weighted_num/base_num*volTPCActive/volTPC);
	
	std::cout<<"7: Ratio of (6) and (4): This is a measure of vertexing efficiency : "<<vertex_eff<<std::endl;

	std::cout<<"==================== Denominator events (Here on x3.1) ==================="<<std::endl;
	std::cout<<"Number of events in vertexed_files that pass denominator with BNB_correction : "<<Ndenominator<<std::endl;
	std::cout<<"Number of events in vertexed_files that pass denominator with BNB_correction with vertex eff removed : "<<Ndenominator/vertex_eff<<std::endl;
	double true_denominator = Ndenominator/vertex_eff*MOD;
	std::cout<<"Same scaled to 6.6e20 "<<true_denominator<<std::endl; 


//	bdt_variable true_photon("delta_photon_energy","(20, 0 , 1.0)","True Photon Energy [GeV]",false,"d");
//	bdt_variable true_proton("delta_proton_energy-0.938272","(20, 0 , 1.0)","True Proton Kinetic Energy [GeV]",false,"d");
//
//	TH1* h_true_photon_denom = (TH1*)file->getTH1(true_photon, denominator, "photon_true_denom", 6.6e20);
//	TH1* h_true_proton_denom = (TH1*)file->getTH1(true_proton, denominator, "proton_true_denom", 6.6e20);
//
//	TH2* h2_true_photon_proton_denom = (TH2*)file->getTH2(true_proton, true_photon, denominator,"2d proton photon denom",6.6e20);
//
//CHANGE photon turns into electron
//CHANGE proton is unchanged:
//

	

	
	bdt_variable true_photon("true_shower_energy[0]","(20, 0 , 1.0)","True Photon Energy [GeV]",false,"d");
	bdt_variable true_proton("true_track_energy[0]-0.938272","(20, 0 , 1.0)","True Proton Kinetic Energy [GeV]",false,"d");

	if(filein->name.compare(9,6,"bnbext",0,6) == 0){//for bnb external
	true_photon.name = "reco_shower_helpernew_energy[0]";
	true_photon.unit = "Reco Photon Energy [GeV]";
	true_proton.name = "reco_track_energy[0]";
	true_photon.unit = "Reco Proton Kinetic Energy [GeV]";
	}

	TH1* h_true_photon_denom = (TH1*)file->getTH1(true_photon, denominator, "photon_true_denom", 6.6e20);
	TH1* h_true_proton_denom = (TH1*)file->getTH1(true_proton, denominator, "proton_true_denom", 6.6e20);

	TH2* h2_true_photon_proton_denom = (TH2*)file->getTH2(true_proton, true_photon, denominator,"2d proton photon denom",6.6e20);
	


	TH1* h_true_photon_numer;
	TH1* h_true_proton_numer;
	TH2* h2_true_photon_proton_numer;

	double nverticies;
	double finaleff;
	for(int s=0; s<4; s++){
		if(s==2) file->calcCosmicBDTEntryList(c1, c2);
		if(s==3) file->calcBNBBDTEntryList(c1, c2);
		file->setStageEntryList(s);

		double stage_entries = file->GetEntries("1")*MOD;
		if(s==0) nverticies = stage_entries;
		std::cout<<"Stage: "<<s<<" Verticies: "<<stage_entries<<" Efficiency: ";
//		std::cout<<stage_entries/true_denominator<<std::endl;
		std::cout<<stage_entries/nverticies<<std::endl;

		if(s==3){
//CHANGE    		std::cout<<"Getting true photon and proton energies"<<std::endl;
			std::cout<<"Getting true electron and proton energies"<<std::endl;
			h_true_photon_numer = (TH1*)file->getTH1(true_photon, denominator, "photon_true_numer", 6.6e20);
			h_true_proton_numer = (TH1*)file->getTH1(true_proton, denominator, "proton_true_numer", 6.6e20);
			h2_true_photon_proton_numer = (TH2*)file->getTH2(true_proton, true_photon, denominator,"2d proton photon numer",6.6e20);
			finaleff = stage_entries/true_denominator*100.0;

		}

	}

	h_true_photon_denom->Scale(MOD/vertex_eff);
	h_true_proton_denom->Scale(MOD/vertex_eff);
	h_true_photon_numer->Scale(MOD/vertex_eff);
	h_true_proton_numer->Scale(MOD/vertex_eff);

	h2_true_photon_proton_numer->Scale(MOD/vertex_eff);
	h2_true_photon_proton_denom->Scale(MOD/vertex_eff);
	std::cout<<"Writing"<<std::endl;

	feff->cd();
	h_true_photon_numer->Write();
	h_true_proton_numer->Write();
	h_true_photon_denom->Write();
	h_true_proton_denom->Write();
	h2_true_photon_proton_numer->Write();
	h2_true_photon_proton_denom->Write();


	TH1* h_true_photon_ratio = (TH1*)h_true_photon_numer->Clone("h_true_photon_ratio");
	TH1* h_true_proton_ratio = (TH1*)h_true_proton_numer->Clone("h_true_proton_ratio");
	TH2* h2_true_photon_proton_ratio = (TH2*)h2_true_photon_proton_numer->Clone("h2_true_photon_proton_ratio");

	h_true_photon_ratio->Divide(h_true_photon_denom);
	h_true_proton_ratio->Divide(h_true_proton_denom);
	h2_true_photon_proton_ratio->Divide(h2_true_photon_proton_denom);
	feff->cd();

	h_true_photon_ratio->Write();
	h_true_proton_ratio->Write();
	h2_true_photon_proton_ratio->Write();


	//Std::plotting
	//
	TCanvas * c = new TCanvas();
	c->cd();

	h_true_photon_ratio->Scale(100);
	h_true_proton_ratio->Scale(100);

	h_true_photon_ratio->Draw("E1");

	h_true_photon_denom->Scale(80/h_true_photon_denom->Integral());
	h_true_photon_denom->SetFillStyle(3454);
	h_true_photon_denom->SetFillColor(kRed-6);
	h_true_photon_denom->SetLineColor(kRed-6);
	h_true_photon_denom->Draw("hist same");

	h_true_proton_denom->Scale(80/h_true_proton_denom->Integral());
	h_true_proton_denom->SetFillStyle(3444);
	h_true_proton_denom->SetFillColor(kBlue-6);
	h_true_proton_denom->SetLineColor(kBlue-6);
	h_true_proton_denom->Draw("hist same");

	
	//sig_notrack->Scale(200);
	//sig_notrack->Draw("hist same");



	h_true_photon_ratio->DrawCopy("E1 same");
	h_true_photon_ratio->SetTitle("");
	h_true_proton_ratio->Draw("E1 same");

	h_true_photon_ratio->SetLineColor(kRed);
	h_true_photon_ratio->SetLineWidth(2);
	h_true_photon_ratio->SetMarkerStyle(20);
	h_true_photon_ratio->SetMarkerColor(kRed);

	h_true_proton_ratio->SetLineColor(kBlue);
	h_true_proton_ratio->SetMarkerColor(kBlue);
	h_true_proton_ratio->SetLineWidth(2);
	h_true_proton_ratio->SetMarkerStyle(20);

	h_true_photon_ratio->GetYaxis()->SetTitle("Efficiency [%]");
	h_true_photon_ratio->GetXaxis()->SetTitle("True Photon/Proton Energy [GeV]");

	h_true_photon_ratio->SetMaximum(32);
	h_true_photon_ratio->SetMinimum(0);
	h_true_photon_ratio->GetXaxis()->SetRangeUser(0,1);

	TLegend *l = new TLegend(0.13,0.79,0.89,0.89);
	l->SetLineColor(kWhite);
	l->SetLineWidth(0);
	l->SetNColumns(2);
	l->AddEntry(h_true_proton_ratio,"Efficiency: True Proton Energy","lp");
	l->AddEntry(h_true_photon_ratio,"Efficiency: True Photon Energy","lp");
	l->AddEntry(h_true_proton_denom,"#splitline{}{True Proton Spectum (au)}","fl");
	l->AddEntry(h_true_photon_denom,"#splitline{}{True Photon Spectum (au)}","fl");
	l->Draw();


		TLatex txt;
		txt.SetTextSize(0.05);
		txt.SetTextAlign(13);  //align at top
		txt.SetNDC();
		txt.DrawLatex(.60,.70,("#splitline{Total Efficiency:}{"+to_string_prec(finaleff,2)+"%}").c_str());
		



	c->Write();
	//c->SaveAs("eff_a.pdf","pdf");
	std::string namea(filein->tag+"_a.pdf");

	c->SaveAs(namea.c_str() ,"pdf");

	TCanvas * c2d = new TCanvas();
	c2d->cd();

	gStyle->SetPalette(kDarkBodyRadiator);	
	h2_true_photon_proton_ratio->Draw("colz");
	h2_true_photon_proton_ratio->GetXaxis()->SetTitle("True Photon Energy [GeV]");
	h2_true_photon_proton_ratio->GetYaxis()->SetTitle("True Proton Kinetic Energy [GeV]");
	
	std::string nameb(filein->tag+"_b.pdf");
	//c2d->SaveAs("eff_b.pdf","pdf");
	c2d->SaveAs(nameb.c_str(),"pdf");

	


	feff->Close();

	//return;
	
//	nverticies = true_denominator;

	std::cout<<"==================== Precuts - Individual  ==================="<<std::endl;
    	file->setStageEntryList(0);
	
	//introduce the "fiducial volumne" cut: std::string denominator.c_str();
	for(int m=0; m< file->flow.vec_pre_cuts.size(); m++){
		std::string tmpcut = file->flow.vec_pre_cuts.at(m);// +"&&"+ denominator.c_str();

		double np = file->GetEntries(tmpcut.c_str())*MOD;
		std::cout<<" + "<<file->flow.vec_pre_cuts.at(m)<<"\t||\t"<<np<<"\t("<<np/nverticies*100<<")\%"<<std::endl;
	}


	std::cout<<"==================== Precuts - One by One  ==================="<<std::endl;
	file->setStageEntryList(0);

//	std::string thiscut = denominator.c_str();
	std::string thiscut = "1";
	for(int m=0; m< file->flow.vec_pre_cuts.size(); m++){
		thiscut += "&&"+ file->flow.vec_pre_cuts.at(m);
		double np = file->GetEntries(thiscut.c_str())*MOD;
		std::cout<<" + "<<file->flow.vec_pre_cuts.at(m)<<"\t||\t"<<np<<"\t("<<np/nverticies*100<<")\%"<<std::endl;
	}
}



