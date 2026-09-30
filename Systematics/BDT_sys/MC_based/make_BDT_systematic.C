// =============================================================================
//  make_BDT_systematic.C
//
//  Final step of the BDT systematic study.  Reads the hadd-ed output of
//  flow_MC_BDT_sys.C and produces, per pT bin and per analysis centrality class:
//        dvn      = vn(class-g BDT cut) - vn(no BDT)
//        dvn/vn   = dvn / vn(no BDT)             <-- the relative BDT systematic
//
//  The reference vn(no BDT) is a single centrality-inclusive profile; the BDT
//  cut is the one for class g applied to the full inclusive sample (this is the
//  getDCA_fromMC.C / getEfficiency.C convention -- the candidate's true event
//  centrality is not used).
//
//  Writes:
//     vN_woBDT                 (one graph, the common reference)
//     vN_wBDT_<cen>            (per class)
//     vN_absdiff_<cen>         (per class)  = vn(BDT) - vn(noBDT)
//     vN_reldiff_<cen>         (per class)  = dvn / vn(noBDT)
//  and prints a summary table.
//
//  Run:
//     root -l -b -q 'make_BDT_systematic.C("<flow_combined.root>","BDT_sys_output.root")'
//  or compiled:
//     g++ make_BDT_systematic.C $(root-config --cflags --libs) -Wall -O2 -o <exe>
//     ./<exe> <flow_combined.root> BDT_sys_output.root
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TString.h"

#include "/home/saha115/D0_ESE/CMSSW_13_2_11/src/Systematics/BDT_sys/MC_based/Analysis_bin_BDT.h"

using namespace std;

// build the reference (no-BDT) graph from a profile
static TGraphErrors *ref_graph(TProfile *pwo, int npt, const double *edges, const char *vtag)
{
  TGraphErrors *g = new TGraphErrors();
  g->SetName(Form("%s_woBDT", vtag));
  g->SetTitle(Form("%s no BDT (cent-incl);p_{T} (GeV);%s", vtag, vtag));
  for (int ip = 0; ip < npt; ++ip)
  {
    double xc = 0.5 * (edges[ip] + edges[ip + 1]);
    double xe = 0.5 * (edges[ip + 1] - edges[ip]);
    g->SetPoint(ip, xc, pwo->GetBinContent(ip + 1));
    g->SetPointError(ip, xe, pwo->GetBinError(ip + 1));
  }
  return g;
}

// per-class comparison: wBDT graph + absolute & relative systematic graphs
static void one_class(TFile *fout, TProfile *pwo, TProfile *pw, int npt,
                      const double *edges, const char *vtag, const char *ctag)
{
  if (!pwo || !pw)
  {
    cerr << "[skip] missing profile for " << vtag << " " << ctag << endl;
    return;
  }

  TGraphErrors *g_w = new TGraphErrors();
  TGraphErrors *g_abs = new TGraphErrors();
  TGraphErrors *g_rel = new TGraphErrors();
  g_w->SetName(Form("%s_wBDT_%s", vtag, ctag));
  g_abs->SetName(Form("%s_absdiff_%s", vtag, ctag));
  g_rel->SetName(Form("%s_reldiff_%s", vtag, ctag));
  g_w->SetTitle(Form("%s BDT(%s cut);p_{T} (GeV);%s", vtag, ctag, vtag));
  g_abs->SetTitle(Form("%s(BDT)-%s(noBDT) %s;p_{T} (GeV);#Delta%s", vtag, vtag, ctag, vtag));
  g_rel->SetTitle(Form("BDT systematic %s %s;p_{T} (GeV);#Delta%s/%s", vtag, ctag, vtag, vtag));

  cout << "\n==== " << vtag << "   BDT cut = " << ctag << "  (ref = no BDT, cent-incl) ====\n";
  cout << left << setw(12) << "pT"
       << setw(20) << "vn(noBDT)" << setw(20) << "vn(BDT)"
       << setw(16) << "diff" << setw(12) << "rel.diff" << "\n";

  for (int ip = 0; ip < npt; ++ip)
  {
    double xc = 0.5 * (edges[ip] + edges[ip + 1]);
    double xe = 0.5 * (edges[ip + 1] - edges[ip]);
    double vwo = pwo->GetBinContent(ip + 1), ewo = pwo->GetBinError(ip + 1);
    double vw = pw->GetBinContent(ip + 1), ew = pw->GetBinError(ip + 1);

    g_w->SetPoint(ip, xc, vw);
    g_w->SetPointError(ip, xe, ew);

    double d = vw - vwo;
    // wBDT sample is a subset of woBDT -> strongly correlated; Barlow error
    // on the difference between a sample and its subsample.
    double de = sqrt(fabs(ew * ew - ewo * ewo));
    g_abs->SetPoint(ip, xc, d);
    g_abs->SetPointError(ip, xe, de);

    double rel = (vwo != 0.0) ? d / vwo : 0.0;
    double rele = (vwo != 0.0) ? de / fabs(vwo) : 0.0;
    g_rel->SetPoint(ip, xc, rel);
    g_rel->SetPointError(ip, xe, rele);

    cout << left << setw(12) << Form("%.0f-%.0f", edges[ip], edges[ip + 1])
         << setw(20) << Form("%.5f +- %.5f", vwo, ewo)
         << setw(20) << Form("%.5f +- %.5f", vw, ew)
         << setw(16) << Form("%.5f", d)
         << setw(12) << Form("%.1f%%", 100. * rel) << "\n";
  }

  fout->cd();
  g_w->Write();
  g_abs->Write();
  g_rel->Write();
}

int make_BDT_systematic(TString fin_name =
                            "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
                            "BDT_sys_Flow_MC_prompt_Sept10_v2/ROOT/flow_MC_BDT_out_combined.root",
                        TString fout_name = "BDT_sys_output.root")
{
  TH1::StatOverflows(kTRUE);

  TFile *fin = TFile::Open(fin_name);
  if (!fin || fin->IsZombie())
  {
    cerr << "cannot open " << fin_name << endl;
    return 1;
  }
  TFile *fout = new TFile(fout_name, "RECREATE");

  TProfile *pv2_wo = (TProfile *)fin->Get("pv2_woBDT");
  TProfile *pv3_wo = (TProfile *)fin->Get("pv3_woBDT");
  if (!pv2_wo || !pv3_wo)
  {
    cerr << "FATAL: missing pv2_woBDT / pv3_woBDT in " << fin_name << endl;
    return 1;
  }

  fout->cd();
  ref_graph(pv2_wo, N_PTBINS_V2, pt_edges_v2, "v2")->Write();
  ref_graph(pv3_wo, N_PTBINS_V3, pt_edges_v3, "v3")->Write();

  for (int g = 0; g < N_CENTBINS; ++g)
  {
    one_class(fout, pv2_wo, (TProfile *)fin->Get(Form("pv2_wBDT_%s", cen_name[g])),
              N_PTBINS_V2, pt_edges_v2, "v2", cen_name[g]);
    one_class(fout, pv3_wo, (TProfile *)fin->Get(Form("pv3_wBDT_%s", cen_name[g])),
              N_PTBINS_V3, pt_edges_v3, "v3", cen_name[g]);
  }

  // copy the profiles + pT spectra through for convenience
  fout->cd();
  pv2_wo->Write();
  pv3_wo->Write();
  for (int g = 0; g < N_CENTBINS; ++g)
    for (const char *pfx : {"pv2_wBDT_", "pv3_wBDT_", "hpt_wBDT_"})
    {
      TObject *o = fin->Get(Form("%s%s", pfx, cen_name[g]));
      if (o) o->Write();
    }
  { TObject *o = fin->Get("hpt_woBDT"); if (o) o->Write(); }

  fout->Write(0, TObject::kOverwrite);
  fout->Close();
  fin->Close();
  cout << "\n>>> wrote " << fout_name << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  TString fin_name = (argc > 1) ? argv[1]
                                : "/scratch/negishi/saha115/D0_ESE_out/CMSSW_13_2_11/src/"
                                  "BDT_sys_Flow_MC_prompt_Sept10_v2/ROOT/flow_MC_BDT_out_combined.root";
  TString fout_name = (argc > 2) ? argv[2] : "BDT_sys_output.root";
  return make_BDT_systematic(fin_name, fout_name);
}
