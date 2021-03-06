


void sig(TFile * file, std::string method) {

  std::string signal_definition = "is_delta_rad == 1 && true_nu_vtx_fid_contained == 1";

  double pot_sp = 2.23455e+22;
  double pot_bnb_cosmic = 1.19874e+19;
  double ngenbnbcosmic = 11750;
  double ngencosmic = 10936;

  double run_pot = 6.6e20;

  TTree * tree_sp = (TTree*)file->Get(("ncdeltarad_"+method).c_str()); 
  TTree * tree_bnb_cosmic = (TTree*)file->Get(("bnbcosmic_"+method).c_str());
  TTree * tree_cosmic = (TTree*)file->Get(("cosmic_"+method).c_str());

  double total_scaled_sp = tree_sp->GetEntries(signal_definition.c_str()) * run_pot / pot_sp;
  double total_scaled_bnb_cosmic = tree_bnb_cosmic->GetEntries() * run_pot / pot_bnb_cosmic;
  double total_scaled_cosmic = tree_cosmic->GetEntries() * run_pot * ngenbnbcosmic / ngencosmic * 10.729 / pot_bnb_cosmic;

  double largest_significance = 0;
  int isp = 0;
  int ibc = 0;
  int ic = 0;
  double lsp = 0;
  double lbc = 0;
  double lc = 0;

  cout << method << "\n";

  for(double d = -1; d < 1; d += 0.05) {
  
    double scaled_sp = tree_sp->GetEntries(("mva > "+std::to_string(d)+"&&"+signal_definition).c_str()) * run_pot / pot_sp;
    double scaled_bnb_cosmic = tree_bnb_cosmic->GetEntries(("mva > "+std::to_string(d)).c_str()) * run_pot / pot_bnb_cosmic;
    double scaled_cosmic = tree_cosmic->GetEntries(("mva > "+std::to_string(d)).c_str()) * run_pot * ngenbnbcosmic / ngencosmic * 10.729 / pot_bnb_cosmic;

    if(scaled_bnb_cosmic+scaled_cosmic) {
      if(scaled_sp / sqrt(scaled_bnb_cosmic+scaled_cosmic) > largest_significance) {
	largest_significance = scaled_sp / sqrt(scaled_bnb_cosmic+scaled_cosmic);
	isp = tree_sp->GetEntries(("mva > "+std::to_string(d)+"&&"+signal_definition).c_str());
	lsp = scaled_sp;
	ibc = tree_bnb_cosmic->GetEntries(("mva > "+std::to_string(d)).c_str());
	lbc = scaled_bnb_cosmic;
	ic = tree_cosmic->GetEntries(("mva > "+std::to_string(d)).c_str());
	lc = scaled_cosmic;
      } 
    }
    else if(scaled_sp) {
      cout << "Background eliminated scaled signal: " << scaled_sp << " efficiency: " << scaled_sp / total_scaled_sp * 100 << " %\n";
    }

  }

  cout << "total - sp: " << tree_sp->GetEntries(signal_definition.c_str()) << " scaled: " << total_scaled_sp << "\n"
       << "        bnb_cosmic: " << tree_bnb_cosmic->GetEntries() << " scaled: " << total_scaled_bnb_cosmic << "\n"
       << "        cosmic: " << tree_cosmic->GetEntries() << " scaled: " << total_scaled_cosmic << "\n"
       << "after - sp: " << isp << " scaled: " << lsp << "\n"
       << "        bnb_cosmic: " << ibc << " scaled: " << lbc << "\n"
       << "        cosmic: " << ic << " scaled: " << lc << "\n"
       << "largest significance: " << largest_significance << "\n"
       << "seff: " << lsp / total_scaled_sp * 100 << " % beff: " << (lbc + lc) / (total_scaled_bnb_cosmic + total_scaled_cosmic) * 100 << " %\n\n";

}



void getsig(std::string const fname) {
  
  TFile * file = TFile::Open(fname.c_str());

  sig(file, "BDTG");
  sig(file, "BDT");
  sig(file, "BDTB");
  sig(file, "BDTD");
  sig(file, "RuleFit");

}
