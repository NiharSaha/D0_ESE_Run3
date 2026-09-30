#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "TNtuple.h"

using namespace std;

// Helper function to get 1% centrality bin from hiHF.

int get_cent1_from_hiHF(Float_t hiHF, const std::vector<double> &hf_boundaries)
{
  int n_bins = hf_boundaries.size() - 1;

  if (hiHF >= hf_boundaries.back())
    return 0;

  for (int i = 0; i < n_bins; ++i)
  {
    if (hiHF >= hf_boundaries[i] && hiHF < hf_boundaries[i + 1])
    {
      // Invert the index: the highest hiHF bracket (i=199) maps to centrality 0 (most central)
      int cent = (n_bins - 1) - i;
      return cent / 2;
    }
  }
  return -1;
}

void make_q2_slices(TString input_txt, TString output_path, int istart, int iend)
{

  TH1::SetDefaultSumw2();

  // 1 (UP), -1 (DOWN), 0 (Nominal)
  int sys_mode = -1;

  std::vector<double> hf_table_nominal = {0, 10.3417, 11.0768, 11.7884, 12.5014, 13.223, 13.979, 14.7243, 15.4805, 16.2881, 17.0953, 17.9161, 18.7604, 19.6304, 20.551, 21.5208, 22.5138, 23.5369, 24.5918, 25.6714, 26.8101, 28.006, 29.2436, 30.5117, 31.8466, 33.228, 34.6998, 36.2222, 37.8447, 39.5591, 41.3376, 43.1468, 45.0235, 47.0435, 49.0904, 51.2926, 53.6121, 56.0542, 58.5133, 61.1585, 63.8113, 66.639, 69.5748, 72.568, 75.676, 78.9026, 82.3018, 85.7734, 89.444, 93.2175, 97.2133, 101.253, 105.319, 109.546, 113.867, 118.346, 123.019, 127.884, 132.907, 138.133, 143.553, 149.142, 154.94, 160.848, 166.99, 173.359, 179.998, 186.842, 193.842, 201.117, 208.532, 216.211, 224.158, 232.3, 240.736, 249.435, 258.302, 267.463, 277.014, 286.845, 296.83, 307.189, 317.783, 328.756, 339.992, 351.395, 363.16, 375.455, 387.794, 400.501, 413.561, 426.97, 440.651, 454.496, 468.886, 483.371, 498.4, 513.754, 529.457, 545.441, 561.988, 578.779, 595.858, 613.597, 631.446, 649.849, 668.467, 687.351, 706.767, 726.71, 747.018, 767.66, 788.701, 810.28, 832.266, 854.725, 877.463, 900.464, 924.078, 948.465, 973.01, 998.014, 1023.62, 1049.62, 1075.96, 1102.82, 1129.95, 1157.78, 1186.03, 1215.15, 1244.29, 1274.23, 1304.67, 1335.36, 1366.46, 1398.45, 1431.19, 1464.02, 1497.65, 1531.96, 1566.77, 1602.1, 1637.92, 1674.27, 1711.1, 1749.07, 1787.21, 1826.02, 1865.25, 1905.66, 1946.56, 1987.95, 2030.81, 2073.58, 2117.19, 2161.79, 2206.96, 2252.79, 2300.01, 2347.25, 2395.35, 2444.79, 2494.51, 2544.84, 2596.27, 2649.26, 2703.27, 2758.12, 2813.87, 2870.23, 2927.2, 2985.36, 3045.09, 3105.31, 3166.48, 3229.14, 3293.26, 3359.1, 3425.63, 3493.44, 3562.41, 3633.68, 3706.37, 3780.38, 3856.7, 3934.4, 4013.53, 4095.4, 4178.3, 4263.29, 4350.7, 4440.56, 4532.84, 4628.66, 4727.32, 4827.83, 4933.06, 5042.94, 5161.8, 5305.56, 8171.19};
  std::vector<double> hf_table_up = {0, 10.3379, 11.0725, 11.7802, 12.4919, 13.21, 13.9601, 14.7025, 15.4608, 16.2614, 17.0642, 17.8839, 18.7243, 19.5915, 20.501, 21.4676, 22.4592, 23.4754, 24.5231, 25.597, 26.7278, 27.914, 29.1421, 30.4015, 31.7342, 33.0996, 34.549, 36.0668, 37.6735, 39.3801, 41.1438, 42.9432, 44.8014, 46.7956, 48.853, 51.0064, 53.3073, 55.7264, 58.1748, 60.78, 63.436, 66.2411, 69.112, 72.1054, 75.1783, 78.3803, 81.7626, 85.1756, 88.7616, 92.5521, 96.492, 100.55, 104.72, 109.084, 113.514, 118.197, 123.071, 128.061, 133.341, 138.724, 144.257, 149.905, 155.906, 162.144, 168.39, 174.967, 181.708, 188.672, 196.01, 203.584, 211.181, 219.091, 227.066, 235.409, 244.073, 253.075, 262.124, 271.539, 281.203, 290.991, 301.253, 311.645, 322.338, 333.437, 344.657, 356.169, 368.062, 380.364, 392.808, 405.651, 418.805, 432.246, 445.973, 459.933, 474.24, 488.99, 504.095, 519.521, 535.239, 551.417, 568.002, 584.782, 601.983, 619.827, 637.684, 656.132, 674.811, 693.728, 713.198, 733.254, 753.661, 774.269, 795.462, 817.018, 839.074, 861.556, 884.288, 907.465, 931.318, 955.542, 980.139, 1005.18, 1030.77, 1056.79, 1083.18, 1110.08, 1137.3, 1165.1, 1193.47, 1222.35, 1251.81, 1281.66, 1311.96, 1342.79, 1373.97, 1406.01, 1438.58, 1471.65, 1505.24, 1539.46, 1574.28, 1609.71, 1645.39, 1681.8, 1718.61, 1756.42, 1794.69, 1833.55, 1872.86, 1913.08, 1953.8, 1995.29, 2038.11, 2080.97, 2124.46, 2168.93, 2214.04, 2260.09, 2307.06, 2354.41, 2402.35, 2451.68, 2501.29, 2551.71, 2603.17, 2655.91, 2709.96, 2764.55, 2820.31, 2876.54, 2933.35, 2991.78, 3051.22, 3111.18, 3172.27, 3234.87, 3298.85, 3364.68, 3430.89, 3498.44, 3567.44, 3638.38, 3711.01, 3784.97, 3861.04, 3938.56, 4017.48, 4099.05, 4181.81, 4266.65, 4353.98, 4443.42, 4535.56, 4630.97, 4729.41, 4829.54, 4934.49, 5044.11, 5162.64, 5306.21, 8171.19};
  std::vector<double> hf_table_down = {0, 10.3463, 11.0861, 11.7998, 12.5172, 13.2458, 14.0069, 14.7568, 15.519, 16.3334, 17.1454, 17.9715, 18.8204, 19.7036, 20.6264, 21.605, 22.6104, 23.6387, 24.7063, 25.796, 26.9534, 28.1605, 29.4073, 30.6918, 32.0377, 33.4416, 34.9364, 36.4788, 38.1179, 39.8435, 41.6583, 43.4757, 45.4068, 47.4261, 49.5297, 51.6006, 53.7062, 55.9019, 58.199, 60.5762, 63.0623, 65.6426, 68.3145, 71.1183, 74.0371, 77.0526, 80.1847, 83.4274, 86.8108, 90.3357, 93.9759, 97.7734, 101.745, 105.854, 110.126, 114.482, 119.024, 123.747, 128.668, 133.745, 139.036, 144.527, 150.19, 156.04, 162.051, 168.268, 174.704, 181.437, 188.373, 195.443, 202.785, 210.302, 218.137, 226.175, 234.398, 242.965, 251.767, 260.763, 270.085, 279.822, 289.708, 299.86, 310.338, 321.122, 332.245, 343.563, 355.174, 367.154, 379.524, 392.109, 405.031, 418.299, 431.845, 445.699, 459.798, 474.234, 489.117, 504.386, 519.952, 535.84, 552.219, 569.001, 585.895, 603.3, 621.38, 639.438, 658.099, 676.89, 696.136, 715.903, 736.115, 756.791, 777.69, 799.191, 821.022, 843.312, 866.052, 889.133, 912.586, 936.843, 961.302, 986.171, 1011.74, 1037.65, 1064.03, 1090.8, 1118.01, 1145.77, 1173.82, 1202.86, 1232.18, 1262.05, 1292.32, 1323.16, 1354.33, 1386.08, 1418.62, 1451.56, 1485.4, 1519.47, 1554.19, 1589.76, 1625.49, 1661.89, 1698.83, 1736.63, 1775, 1813.88, 1853.12, 1893.44, 1934.42, 1975.76, 2018.46, 2061.74, 2105.37, 2149.7, 2194.95, 2240.94, 2288.16, 2335.73, 2383.92, 2433.41, 2483.29, 2533.65, 2585.27, 2638.35, 2692.45, 2747.33, 2803.27, 2859.75, 2916.98, 2975.45, 3035.05, 3095.57, 3156.99, 3219.67, 3284.06, 3349.98, 3417.04, 3485.29, 3554.01, 3625.58, 3698.49, 3772.71, 3849.45, 3927.34, 4006.82, 4089.11, 4172.24, 4257.81, 4345.44, 4435.8, 4528.43, 4624.6, 4723.91, 4824.73, 4930.41, 5040.92, 5160.37, 5304.54, 8171.19};

  // Select the active table based on the flag
  std::vector<double> active_hf_table;
  TString sys_suffix = "";
  if (sys_mode == 1)
  {
    active_hf_table = hf_table_up;
    sys_suffix = "_sysUp";
  }
  else if (sys_mode == -1)
  {
    active_hf_table = hf_table_down;
    sys_suffix = "_sysDown";
  }
  else
  {
    active_hf_table = hf_table_nominal;
    sys_suffix = "_nominal";
  }

  ifstream file_stream(input_txt.Data());
  TString outfile_ntup = TString::Format("%s/ROOT/Quantiles_ESE_ntuples_%d_%d%s.root", output_path.Data(), istart, iend, sys_suffix.Data());
  TString outfile_hist = TString::Format("%s/ROOT/Quantiles_ESE_hists_%d_%d%s.root", output_path.Data(), istart, iend, sys_suffix.Data());
  TFile *fout_ntup = new TFile(outfile_ntup, "RECREATE");
  TFile *fout_hist = new TFile(outfile_hist, "RECREATE");
  string filename;
  int ifile = 0;

  // const int N_CENTBIN = 4;
  // Int_t min_centbin[N_CENTBIN] = {0, 10, 30, 50};
  // Int_t max_centbin[N_CENTBIN] = {10, 30, 50, 90};
  const int N_CENTBINS_1PC = 50;

  // Define your 1% histograms exactly as before
  TH1D *hist_q2_tot[N_CENTBINS_1PC];
  TH1D *hist_q3_tot[N_CENTBINS_1PC];

  for (int j = 0; j < N_CENTBINS_1PC; ++j)
  {
    hist_q2_tot[j] = new TH1D(Form("hist_q2_tot_cent%i_%i", j, j + 1), Form("hist_q2_tot_cent%i_%i", j, j + 1), 35000, 0.0, 0.70);
    hist_q2_tot[j]->SetDirectory(0);
    hist_q3_tot[j] = new TH1D(Form("hist_q3_tot_cent%i_%i", j, j + 1), Form("hist_q3_tot_cent%i_%i", j, j + 1), 35000, 0.0, 0.70);
    hist_q3_tot[j]->SetDirectory(0);
  }

  // TNtuple *nt_ese_global = new TNtuple("nt_global", "nt_global", "centrality:q2_hfp:q2_hfm:q3_hfp:q3_hfm:q2_total:q3_total");

  TNtuple *nt[N_CENTBINS_1PC];
  for (int i_cen = 0; i_cen < N_CENTBINS_1PC; i_cen++)
  {
    nt[i_cen] = new TNtuple(Form("nt_qn_cent%i_%i", i_cen, i_cen + 1), Form("nt_qn_cent%i_%i", i_cen, i_cen + 1), "q2_hfp:q2_hfm:q3_hfp:q3_hfm:q2_hf_total:q3_hf_total:q2_hfp_w:q2_hfm_w:q3_hfp_w:q3_hfm_w:q2_hf_total_w:q3_hf_total_w");
    nt[i_cen]->SetDirectory(0);
  }

  while (file_stream >> filename)
  {
    if (ifile < istart)
    {
      ifile++;
      continue;
    }
    if (ifile >= iend)
      break;

    TFile *fin = TFile::Open(filename.c_str());
    if (!fin || fin->IsZombie())
    {
      ifile++;
      continue;
    }

    std::cout << ">>> Processing ifile=" << ifile << " : " << filename << std::endl;

    TTree *t_eventinfoana = (TTree *)fin->Get("eventinfoana/EventInfoNtuple");
    Short_t centrality;
    Float_t ephfmQ[3], ephfpQ[3], ephfmSumW[3], ephfpSumW[3];
    Float_t hiHF;

    t_eventinfoana->SetBranchAddress("HFsumET", &hiHF);
    t_eventinfoana->SetBranchAddress("centrality", &centrality);
    t_eventinfoana->SetBranchAddress("ephfmQ", ephfmQ);
    t_eventinfoana->SetBranchAddress("ephfpQ", ephfpQ);
    t_eventinfoana->SetBranchAddress("ephfmSumW", ephfmSumW);
    t_eventinfoana->SetBranchAddress("ephfpSumW", ephfpSumW);

    Int_t n_entries = t_eventinfoana->GetEntries();
    std::cout << "Nevents = " << n_entries << std::endl;

    for (Long64_t ii = 0; ii < n_entries; ii++)
    {
      t_eventinfoana->GetEntry(ii);

      if (ii % 10000 == 0)
        printf("Current entry of the loop = %lld out of %d : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

      Int_t cent = get_cent1_from_hiHF(hiHF, active_hf_table);
      if (cent < 0 || cent >= N_CENTBINS_1PC)
        continue;

      Double_t q2_hfp = (ephfpSumW[1] > 0) ? ephfpQ[1] / ephfpSumW[1] : -999;
      Double_t q2_hfm = (ephfmSumW[1] > 0) ? ephfmQ[1] / ephfmSumW[1] : -999;
      Double_t q3_hfp = (ephfpSumW[2] > 0) ? ephfpQ[2] / ephfpSumW[2] : -999;
      Double_t q3_hfm = (ephfmSumW[2] > 0) ? ephfmQ[2] / ephfmSumW[2] : -999;
      Double_t q2_total_norm = (ephfmSumW[1] + ephfpSumW[1] > 0) ? (ephfmQ[1] + ephfpQ[1]) / (ephfmSumW[1] + ephfpSumW[1]) : -999;
      Double_t q3_total_norm = (ephfmSumW[2] + ephfpSumW[2] > 0) ? (ephfmQ[2] + ephfpQ[2]) / (ephfmSumW[2] + ephfpSumW[2]) : -999;

      // nt_ese_global->Fill(cent, q2_hfp, q2_hfm, q3_hfp, q3_hfm, q2_total_norm, q3_total_norm);

      if (cent >= 0 && cent < N_CENTBINS_1PC)
      {
        hist_q2_tot[cent]->Fill(q2_total_norm);
        hist_q3_tot[cent]->Fill(q3_total_norm);
        nt[cent]->Fill(q2_hfp, q2_hfm, q3_hfp, q3_hfm, q2_total_norm, q3_total_norm, ephfpSumW[1], ephfmSumW[1], ephfpSumW[2], ephfmSumW[2], (ephfmSumW[1] + ephfpSumW[1]), (ephfmSumW[2] + ephfpSumW[2]));
      }

    } // -- evt loop --
    fin->Close();
    delete fin;
    ifile++;
  }

  fout_hist->cd();
  fout_hist->mkdir("q2_raw_1pc");
  fout_hist->cd("q2_raw_1pc");

  for (int j = 0; j < N_CENTBINS_1PC; ++j)
  {
    if (hist_q2_tot[j]->GetEntries() > 0)
      hist_q2_tot[j]->Write();
  }

  fout_hist->cd();
  fout_hist->mkdir("q3_raw_1pc");
  fout_hist->cd("q3_raw_1pc");

  for (int j = 0; j < N_CENTBINS_1PC; ++j)
  {
    if (hist_q3_tot[j]->GetEntries() > 0)
      hist_q3_tot[j]->Write();
  }

  fout_ntup->cd();
  fout_ntup->mkdir("Q_ntuple");
  fout_ntup->cd("Q_ntuple");
  // nt_ese_global->Write();
  for (int i_cen = 0; i_cen < N_CENTBINS_1PC; i_cen++)
  {
    nt[i_cen]->Write();
  }

  fout_hist->Close();
  fout_ntup->Close();
}

int main(int argc, char *argv[])
{
  if (argc == 5) // Expecting 4 arguments: input_txt, output_path, istart, iend

  {
    TString input_txt = argv[1];
    TString output_path = argv[2];
    int istart = std::stoi(argv[3]);
    int iend = std::stoi(argv[4]);

    make_q2_slices(input_txt, output_path, istart, iend);
  }
  else
  {
    std::cout << "Usage: ./your_program input_txt output_path istart iend" << std::endl;
    return 1;
  }
  return 0;
}
