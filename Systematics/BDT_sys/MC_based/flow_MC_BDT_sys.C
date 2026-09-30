// =============================================================================
//  flow_MC_BDT_sys.C
//
//  BDT systematic study for prompt-D0 v2 / v3 with the scalar-product method
//  on the 2023 PbPb prompt-D0 MC.
//
//  Method
//  ------
//  In data the flow comes from a mass fit (yield method).  In MC the signal is
//  simply the gen-matched, non-swap, prompt reconstructed D0 candidates
//  (matchGEN==1 && isSwap==0 && Dgen_isPrompt==1) -- no mass fit needed.
//
//  The 3-sub-event resolution is computed in a SEPARATE step
//  (Calculate_Resolution_MC.C -> hadd), q-inclusive and centrality-inclusive
//  over the whole MC sample.  This macro reads that file, turns it into
//    res_plus  = sqrt( <HFm.HFp> <HFm.Trk> / <HFp.Trk> )   (D0 with y>0, ref HF-)
//    res_minus = sqrt( <HFm.HFp> <HFp.Trk> / <HFm.Trk> )   (D0 with y<0, ref HF+)
//  and does the whole vn calculation in one pass -- no intermediate ntuple,
//  exactly like charge_vn/for_20Qbin/flow_Analysis_chg.C.
//
//  Centrality
//  ----------
//  There is no meaningful centrality in the MB-embedded MC, so:
//    * events of ALL centralities are used (only a 0-200 sanity bound);
//    * the vn numerator and the resolution are fully centrality-inclusive;
//    * the BDT cut is applied the getDCA_fromMC.C / getEfficiency.C way: the
//      candidate's true event centrality is ignored, and the result is produced
//      once per analysis class, applying that class's BDT cut
//      (bdtHandler.getBDTCut(y, cent_bdt_arg[g], pT)) to the full sample.
//
//  For every selected candidate:  vn = Re( u_n(D0) . Q*_n(HF_ref) ) / res
//     "_woBDT"       : one profile, all gen-matched prompt candidates (no cut)
//     "_wBDT_<cen>"  : one profile per class, candidates passing that class cut
//
//  Build / run:
//     g++ flow_MC_BDT_sys.C $(root-config --cflags --libs) -Wall -O2 -o <exe>
//     ./<exe>  <input_txt>  <output_path>  <istart>  <iend>  [<resolution_file>]
//  output:  <output_path>/ROOT/flow_MC_BDT_sys_<istart>_<iend>.root
//
//  The per-job outputs are hadd-ed and turned into the systematic graphs by
//  make_BDT_systematic.C.
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/BDTHandler.h"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/BDTHandler.cc"
#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/Analysis_bin_BDT.h"

using namespace std;

// safe 3-sub resolution sqrt( AB * AC / BC ); -1 if ill-defined
static double res3sub(double AB, double AC, double BC)
{
  if (BC == 0.0)
    return -1.0;
  double r2 = AB * AC / BC;
  return (r2 > 0.0) ? sqrt(r2) : -1.0;
}

// Check the cut table before handing it to BDTHandler::loadCuts, which does no
// checking itself (unopenable file -> silent return; malformed line -> write at
// index -9; absent bin -> uninitialised cut).  Every |y|<1 bin must be present,
// since that is the only y bin this macro selects.
static void validateBDTCuts(const TString &fname)
{
  ifstream f(fname.Data());
  if (!f.is_open())
  {
    cerr << "FATAL: cannot open BDT cut file " << fname << endl;
    exit(1);
  }
  bool seen[y_bins][cent_bins][pT_bins] = {};
  string line;
  while (getline(f, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    stringstream ss(line);
    int iy, ic, ip;
    double cut;
    if (!(ss >> iy >> ic >> ip >> cut) || iy < 0 || iy >= y_bins ||
        ic < 0 || ic >= cent_bins || ip < 0 || ip >= pT_bins)
    {
      cerr << "FATAL: bad line in " << fname << ": '" << line << "'" << endl;
      exit(1);
    }
    seen[iy][ic][ip] = true;
  }
  for (int ic = 0; ic < cent_bins; ++ic)
    for (int ip = 0; ip < pT_bins; ++ip)
      if (!seen[0][ic][ip])
      {
        cerr << "FATAL: " << fname << " has no cut for |y|<1 bin (cent " << ic
             << ", pT " << ip << ")" << endl;
        exit(1);
      }
}

// BDT pass fraction per v2 pT bin (cent0to10 cut), all samples pooled
// unweighted, from the reference output -> keep probability of the random cut
static void loadRandFrac(const TString &fname, double *frac)
{
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie())
  {
    cerr << "FATAL: cannot open random-cut reference " << fname << endl;
    exit(1);
  }
  for (int ib = 1; ib <= N_PTBINS_V2; ++ib)
  {
    double nwo = 0, nw = 0;
    for (int s = 0; s < N_DPT_SAMPLES; ++s)
    {
      int thr = (int)Dpt_threshold[s];
      TProfile *pwo = (TProfile *)f->Get(Form("pv2_woBDT_Dpt%d", thr));
      TProfile *pw = (TProfile *)f->Get(Form("pv2_wBDT_%s_Dpt%d", cen_name[0], thr));
      if (!pwo || !pw)
      {
        cerr << "FATAL: " << fname << " has no per-sample profiles" << endl;
        exit(1);
      }
      nwo += pwo->GetBinEntries(ib);
      nw += pw->GetBinEntries(ib);
    }
    frac[ib - 1] = (nwo > 0) ? nw / nwo : 0.5;
    printf(">>> random-cut keep fraction pT %g-%g: %.3f\n", pt_edges_v2[ib - 1], pt_edges_v2[ib], frac[ib - 1]);
  }
  f->Close();
}

void flow_MC_BDT_sys(TString input_txt, TString output_path, int istart, int iend,
                     TString res_file = "")
{
  TH1::SetDefaultSumw2();
  TH1::StatOverflows(kTRUE); // match the resolution step when reading its GetMean()
  BDTHandler bdtHandler;
  // reload the validated file so the cuts in use are guaranteed to be the checked ones
  validateBDTCuts(BDT_CUTS_FILE);
  bdtHandler.loadCuts(BDT_CUTS_FILE.Data());

  double rand_frac[N_PTBINS_V2];
  loadRandFrac(RAND_FRAC_REF, rand_frac);
  TRandom3 rng(4357 + istart); // reproducible, independent per job

  if (res_file.Length() == 0)
    res_file = RES_FILE_DEFAULT;

  // -------------------------------------------------------------------------
  //  1. resolution (q-inclusive, centrality-inclusive)  -> 4 scalars
  // -------------------------------------------------------------------------
  TFile *fres = TFile::Open(res_file);
  if (!fres || fres->IsZombie())
  {
    cerr << "FATAL: cannot open resolution file " << res_file << endl;
    return;
  }
  auto meanOf = [&](const char *n) -> double
  {
    TH1D *h = (TH1D *)fres->Get(n);
    if (!h) { cerr << "FATAL: missing " << n << " in " << res_file << endl; exit(1); }
    return h->GetMean();
  };
  double m2mp = meanOf("Q2Q2_HFmHFp_Re");
  double m2mt = meanOf("Q2Q2_HFmTrk_Re");
  double m2pt = meanOf("Q2Q2_HFpTrk_Re");
  double m3mp = meanOf("Q3Q3_HFmHFp_Re");
  double m3mt = meanOf("Q3Q3_HFmTrk_Re");
  double m3pt = meanOf("Q3Q3_HFpTrk_Re");
  fres->Close();

  const double res2_plus  = res3sub(m2mp, m2mt, m2pt); // D0 y>0 -> ref HF-minus
  const double res2_minus = res3sub(m2mp, m2pt, m2mt); // D0 y<0 -> ref HF-plus
  const double res3_plus  = res3sub(m3mp, m3mt, m3pt);
  const double res3_minus = res3sub(m3mp, m3pt, m3mt);
  printf(">>> resolution (incl):  res2_plus=%.4f res2_minus=%.4f  res3_plus=%.4f res3_minus=%.4f\n",
         res2_plus, res2_minus, res3_plus, res3_minus);

  // -------------------------------------------------------------------------
  //  2. output profiles
  // -------------------------------------------------------------------------
  gSystem->Exec(TString::Format("mkdir -p %s/ROOT", output_path.Data()));
  TString outfile = TString::Format("%s/ROOT/flow_MC_BDT_sys_%d_%d.root",
                                    output_path.Data(), istart, iend);
  TFile *fout = new TFile(outfile, "RECREATE");

  // reference: no BDT cut (single, centrality-inclusive)
  TProfile *pv2_wo = new TProfile("pv2_woBDT", "v2 no BDT;p_{T} (GeV);v_{2}", N_PTBINS_V2, pt_edges_v2);
  TProfile *pv3_wo = new TProfile("pv3_woBDT", "v3 no BDT;p_{T} (GeV);v_{3}", N_PTBINS_V3, pt_edges_v3);
  TH1D *hpt_wo = new TH1D("hpt_woBDT", "cand pT no BDT;p_{T} (GeV)", N_PTBINS_V2, pt_edges_v2);

  // one profile per analysis centrality class -- that class's BDT cut applied to
  // the full inclusive candidate sample
  TProfile *pv2_w[N_CENTBINS], *pv3_w[N_CENTBINS];
  TH1D *hpt_w[N_CENTBINS];
  for (int g = 0; g < N_CENTBINS; ++g)
  {
    pv2_w[g] = new TProfile(Form("pv2_wBDT_%s", cen_name[g]), Form("v2 BDT(%s cut);p_{T} (GeV);v_{2}", cen_name[g]), N_PTBINS_V2, pt_edges_v2);
    pv3_w[g] = new TProfile(Form("pv3_wBDT_%s", cen_name[g]), Form("v3 BDT(%s cut);p_{T} (GeV);v_{3}", cen_name[g]), N_PTBINS_V3, pt_edges_v3);
    hpt_w[g] = new TH1D(Form("hpt_wBDT_%s", cen_name[g]), Form("cand pT BDT(%s cut);p_{T} (GeV)", cen_name[g]), N_PTBINS_V2, pt_edges_v2);
  }

  // same profiles per Dpt sample, UNWEIGHTED: each sample is its own closure
  // test (true vn = 0), combined per pT bin by inverse variance downstream
  // (SampleCombination.h).  Avoids the stitch-weight outliers from gen-pT
  // migration across the sample thresholds.
  TProfile *pv2_wo_s[N_DPT_SAMPLES], *pv3_wo_s[N_DPT_SAMPLES];
  TProfile *pv2_w_s[N_CENTBINS][N_DPT_SAMPLES], *pv3_w_s[N_CENTBINS][N_DPT_SAMPLES];
  for (int s = 0; s < N_DPT_SAMPLES; ++s)
  {
    int thr = (int)Dpt_threshold[s];
    pv2_wo_s[s] = new TProfile(Form("pv2_woBDT_Dpt%d", thr), Form("v2 no BDT, Dpt>%d;p_{T} (GeV);v_{2}", thr), N_PTBINS_V2, pt_edges_v2);
    pv3_wo_s[s] = new TProfile(Form("pv3_woBDT_Dpt%d", thr), Form("v3 no BDT, Dpt>%d;p_{T} (GeV);v_{3}", thr), N_PTBINS_V3, pt_edges_v3);
    for (int g = 0; g < N_CENTBINS; ++g)
    {
      pv2_w_s[g][s] = new TProfile(Form("pv2_wBDT_%s_Dpt%d", cen_name[g], thr), Form("v2 BDT(%s cut), Dpt>%d;p_{T} (GeV);v_{2}", cen_name[g], thr), N_PTBINS_V2, pt_edges_v2);
      pv3_w_s[g][s] = new TProfile(Form("pv3_wBDT_%s_Dpt%d", cen_name[g], thr), Form("v3 BDT(%s cut), Dpt>%d;p_{T} (GeV);v_{3}", cen_name[g], thr), N_PTBINS_V3, pt_edges_v3);
    }
  }

  // random-cut control replicas, per sample, unweighted (see N_RAND)
  TProfile *pv2_r_s[N_RAND][N_DPT_SAMPLES], *pv3_r_s[N_RAND][N_DPT_SAMPLES];
  for (int r = 0; r < N_RAND; ++r)
    for (int s = 0; s < N_DPT_SAMPLES; ++s)
    {
      int thr = (int)Dpt_threshold[s];
      pv2_r_s[r][s] = new TProfile(Form("pv2_wRAND%d_Dpt%d", r, thr), Form("v2 random cut %d, Dpt>%d;p_{T} (GeV);v_{2}", r, thr), N_PTBINS_V2, pt_edges_v2);
      pv3_r_s[r][s] = new TProfile(Form("pv3_wRAND%d_Dpt%d", r, thr), Form("v3 random cut %d, Dpt>%d;p_{T} (GeV);v_{3}", r, thr), N_PTBINS_V3, pt_edges_v3);
    }

  // -------------------------------------------------------------------------
  //  3. file loop
  // -------------------------------------------------------------------------
  ifstream file_stream(input_txt.Data());
  string filename;
  int ifile = 0;

  while (file_stream >> filename)
  {
    if (ifile < istart) { ifile++; continue; }
    if (ifile >= iend) break;

    TFile *fin = TFile::Open(filename.c_str());
    if (!fin || fin->IsZombie())
    {
      cout << "Warning: skipping bad file: " << filename << endl;
      if (fin) { fin->Close(); delete fin; }
      ifile++;
      continue;
    }
    cout << ">>> Processing ifile=" << ifile << " : " << filename << endl;

    const int isample = dptSampleIndex(filename.c_str());
    if (isample < 0)
    {
      cerr << "FATAL: cannot identify Dpt sample of " << filename << endl;
      exit(1);
    }

    TTree *tree      = (TTree *)fin->Get("d0Analyzer/VCNtuple_D02kpi");
    TTree *tree_gen  = (TTree *)fin->Get("d0Analyzer/AllGensNtuple");
    TTree *t_evtinfo = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    if (!tree || !tree_gen || !t_evtinfo)
    {
      cout << "Skipping: missing tree(s) in " << filename << endl;
      fin->Close(); delete fin; ifile++;
      continue;
    }
    tree->AddFriend(tree_gen);
    tree->AddFriend(t_evtinfo);

    Int_t candSize = 0, centrality = 0, N_genD0s = 0;
    vector<float> *pT = nullptr, *phi = nullptr, *y = nullptr, *mva = nullptr;
    vector<bool>  *matchGEN = nullptr, *isSwap = nullptr, *Dgen_isPrompt = nullptr;
    vector<float> *genD0_pt = nullptr;
    Float_t genWeight = 1.f;
    Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
    Float_t ephfpAngle[3], ephfmAngle[3];

    tree->SetBranchStatus("*", 0);
    for (const auto &p : {"candSize", "centrality", "pT", "phi", "y", "mva",
                          "matchGEN", "isSwap", "Dgen_isPrompt",
                          "N_genD0s", "genD0_pt", "genWeight",
                          "ephfpAngle", "ephfmAngle", "ephfpQ", "ephfmQ",
                          "ephfpSumW", "ephfmSumW"})
      tree->SetBranchStatus(p, 1);

    tree->SetBranchAddress("candSize", &candSize);
    tree->SetBranchAddress("centrality", &centrality);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("mva", &mva);
    tree->SetBranchAddress("matchGEN", &matchGEN);
    tree->SetBranchAddress("isSwap", &isSwap);
    tree->SetBranchAddress("Dgen_isPrompt", &Dgen_isPrompt);
    tree->SetBranchAddress("N_genD0s", &N_genD0s);
    tree->SetBranchAddress("genD0_pt", &genD0_pt);
    tree->SetBranchAddress("genWeight", &genWeight);
    tree->SetBranchAddress("ephfpAngle", ephfpAngle);
    tree->SetBranchAddress("ephfmAngle", ephfmAngle);
    tree->SetBranchAddress("ephfmQ", ephfmQ);
    tree->SetBranchAddress("ephfpQ", ephfpQ);
    tree->SetBranchAddress("ephfmSumW", ephfmSumW);
    tree->SetBranchAddress("ephfpSumW", ephfpSumW);

    Long64_t nevt = tree->GetEntries();
    for (Long64_t ii = 0; ii < nevt; ++ii)
    {
      tree->GetEntry(ii);
      if (ii % 20000 == 0)
        printf("  entry %lld / %lld (%.1f%%)\n", ii, nevt, 100. * ii / nevt);

      if (centrality < 0 || centrality >= 200) // sanity only -- centrality inclusive
        continue;

      double w = 1.0;
      if (USE_MC_WEIGHT)
      {
        double p_gen_max = 0.0;
        for (int ig = 0; ig < N_genD0s; ++ig)
          p_gen_max = max(p_gen_max, (double)genD0_pt->at(ig));
        w = genWeight * getStitchWeight(p_gen_max);
        if (w <= 0.0)
          continue;
      }

      float sumW2 = ephfpSumW[1] + ephfmSumW[1];
      float sumW3 = ephfpSumW[2] + ephfmSumW[2];
      if (sumW2 <= 0 || sumW3 <= 0)
        continue;

      TComplex Q2_HFm(ephfmQ[1] * TMath::Cos(2. * ephfmAngle[1]), ephfmQ[1] * TMath::Sin(2. * ephfmAngle[1]));
      TComplex Q2_HFp(ephfpQ[1] * TMath::Cos(2. * ephfpAngle[1]), ephfpQ[1] * TMath::Sin(2. * ephfpAngle[1]));
      TComplex Q3_HFm(ephfmQ[2] * TMath::Cos(3. * ephfmAngle[2]), ephfmQ[2] * TMath::Sin(3. * ephfmAngle[2]));
      TComplex Q3_HFp(ephfpQ[2] * TMath::Cos(3. * ephfpAngle[2]), ephfpQ[2] * TMath::Sin(3. * ephfpAngle[2]));

      for (int ic = 0; ic < candSize; ++ic)
      {
        if (!matchGEN->at(ic) || isSwap->at(ic))
          continue;
        if (!Dgen_isPrompt->at(ic))
          continue;

        float Pt  = pT->at(ic);
        float Y   = y->at(ic);
        float Phi = phi->at(ic);
        float Mva = mva->at(ic);

        if (Pt < MIN_PT_ANA || Pt >= MAX_PT_ANA || fabs(Y) >= MAX_Y_ANA)
          continue;

        TComplex u2(TMath::Cos(2. * Phi), TMath::Sin(2. * Phi));
        TComplex u3(TMath::Cos(3. * Phi), TMath::Sin(3. * Phi));

        double v2_obs, v3_obs, r2, r3;
        if (Y > 0)
        {
          v2_obs = (u2 * TComplex::Conjugate(Q2_HFm)).Re();
          v3_obs = (u3 * TComplex::Conjugate(Q3_HFm)).Re();
          r2 = res2_plus;  r3 = res3_plus;
        }
        else
        {
          v2_obs = (u2 * TComplex::Conjugate(Q2_HFp)).Re();
          v3_obs = (u3 * TComplex::Conjugate(Q3_HFp)).Re();
          r2 = res2_minus; r3 = res3_minus;
        }

        double v2 = (r2 > 0.0) ? v2_obs / r2 : 0.0;
        double v3 = (r3 > 0.0) ? v3_obs / r3 : 0.0;

        // reference: no BDT
        if (r2 > 0.0) { pv2_wo->Fill(Pt, v2, w); hpt_wo->Fill(Pt, w); }
        if (r3 > 0.0) {  pv3_wo->Fill(Pt, v3, w); }
        if (r2 > 0.0) pv2_wo_s[isample]->Fill(Pt, v2);
        if (r3 > 0.0) pv3_wo_s[isample]->Fill(Pt, v3);

        // per-class BDT cut (candidate's true centrality is not used)
        for (int g = 0; g < N_CENTBINS; ++g)
        {
          if (Mva <= bdtHandler.getBDTCut(Y, cent_bdt_arg[g], Pt))
            continue;
          if (r2 > 0.0) { pv2_w[g]->Fill(Pt, v2, w); hpt_w[g]->Fill(Pt, w); }
          if (r3 > 0.0) {  pv3_w[g]->Fill(Pt, v3, w); }
          if (r2 > 0.0) pv2_w_s[g][isample]->Fill(Pt, v2);
          if (r3 > 0.0) pv3_w_s[g][isample]->Fill(Pt, v3);
        }

        // random-cut control: one keep decision per candidate and replica,
        // shared by v2 and v3 like the BDT decision
        const double keep = rand_frac[TMath::BinarySearch(N_PTBINS_V2 + 1, pt_edges_v2, (double)Pt)];
        for (int ir = 0; ir < N_RAND; ++ir)
        {
          if (rng.Rndm() >= keep)
            continue;
          if (r2 > 0.0) pv2_r_s[ir][isample]->Fill(Pt, v2);
          if (r3 > 0.0) pv3_r_s[ir][isample]->Fill(Pt, v3);
        }
      } // candidate loop
    } // event loop

    fin->Close();
    delete fin;
    ifile++;
  } // file loop

  fout->cd();
  fout->Write(0, TObject::kOverwrite);
  fout->Close();
  cout << ">>> wrote " << outfile << endl;
}

int main(int argc, char *argv[])
{
  if (argc != 5 && argc != 6)
  {
    cout << "Usage: " << argv[0] << " input_txt output_path istart iend [resolution_file]" << endl;
    return 1;
  }
  TString res = (argc == 6) ? argv[5] : "";
  flow_MC_BDT_sys(argv[1], argv[2], std::stoi(argv[3]), std::stoi(argv[4]), res);
  return 0;
}
