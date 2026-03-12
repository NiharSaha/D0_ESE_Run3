#pragma once

#include "vnbinning_generated_Mar9.h"


const int N_CENTBINS = 3;
const int cen_edges[N_CENTBINS+1] = {0, 10, 30, 50};
const char *cen_name[N_CENTBINS] = {"cent0to10", "cent10to30", "cent30to50"};


const int    N_PTBINS  = 12;
const double pt_edges[N_PTBINS+1] = {1,2,3,4,5,6,8,10,15,20,40,60,100};
const char  *pt_name[N_PTBINS]    = {"pT1to2","pT2to3","pT3to4","pT4to5","pT5to6",
                                      "pT6to8","pT8to10","pT10to15","pT15to20",
                                      "pT20to40","pT40to60","pT60to100"};


const int    N_MASSBINS    = 48;
const double fit_range_low  = 1.75;
const double fit_range_high = 1.99;
const double width = (fit_range_high - fit_range_low) / N_MASSBINS;

//--- vn binning array ------
const int N_VBINS_V2 = 62;   // 63 edges
const int N_VBINS_V3 = 66;   // 67 edges

static const double* vnb_v2_table[3][12] = {
    { vnbinning_v2_cent0to10_pT1to2,   vnbinning_v2_cent0to10_pT2to3,
      vnbinning_v2_cent0to10_pT3to4,   vnbinning_v2_cent0to10_pT4to5,
      vnbinning_v2_cent0to10_pT5to6,   vnbinning_v2_cent0to10_pT6to8,
      vnbinning_v2_cent0to10_pT8to10,  vnbinning_v2_cent0to10_pT10to15,
      vnbinning_v2_cent0to10_pT15to20, vnbinning_v2_cent0to10_pT20to40,
      vnbinning_v2_cent0to10_pT40to60, vnbinning_v2_cent0to10_pT60to100 },
    { vnbinning_v2_cent10to30_pT1to2,   vnbinning_v2_cent10to30_pT2to3,
      vnbinning_v2_cent10to30_pT3to4,   vnbinning_v2_cent10to30_pT4to5,
      vnbinning_v2_cent10to30_pT5to6,   vnbinning_v2_cent10to30_pT6to8,
      vnbinning_v2_cent10to30_pT8to10,  vnbinning_v2_cent10to30_pT10to15,
      vnbinning_v2_cent10to30_pT15to20, vnbinning_v2_cent10to30_pT20to40,
      vnbinning_v2_cent10to30_pT40to60, vnbinning_v2_cent10to30_pT60to100 },
    { vnbinning_v2_cent30to50_pT1to2,   vnbinning_v2_cent30to50_pT2to3,
      vnbinning_v2_cent30to50_pT3to4,   vnbinning_v2_cent30to50_pT4to5,
      vnbinning_v2_cent30to50_pT5to6,   vnbinning_v2_cent30to50_pT6to8,
      vnbinning_v2_cent30to50_pT8to10,  vnbinning_v2_cent30to50_pT10to15,
      vnbinning_v2_cent30to50_pT15to20, vnbinning_v2_cent30to50_pT20to40,
      vnbinning_v2_cent30to50_pT40to60, vnbinning_v2_cent30to50_pT60to100 }
};

static const double* vnb_v3_table[3][12] = {
    { vnbinning_v3_cent0to10_pT1to2,   vnbinning_v3_cent0to10_pT2to3,
      vnbinning_v3_cent0to10_pT3to4,   vnbinning_v3_cent0to10_pT4to5,
      vnbinning_v3_cent0to10_pT5to6,   vnbinning_v3_cent0to10_pT6to8,
      vnbinning_v3_cent0to10_pT8to10,  vnbinning_v3_cent0to10_pT10to15,
      vnbinning_v3_cent0to10_pT15to20, vnbinning_v3_cent0to10_pT20to40,
      vnbinning_v3_cent0to10_pT40to60, vnbinning_v3_cent0to10_pT60to100 },
    { vnbinning_v3_cent10to30_pT1to2,   vnbinning_v3_cent10to30_pT2to3,
      vnbinning_v3_cent10to30_pT3to4,   vnbinning_v3_cent10to30_pT4to5,
      vnbinning_v3_cent10to30_pT5to6,   vnbinning_v3_cent10to30_pT6to8,
      vnbinning_v3_cent10to30_pT8to10,  vnbinning_v3_cent10to30_pT10to15,
      vnbinning_v3_cent10to30_pT15to20, vnbinning_v3_cent10to30_pT20to40,
      vnbinning_v3_cent10to30_pT40to60, vnbinning_v3_cent10to30_pT60to100 },
    { vnbinning_v3_cent30to50_pT1to2,   vnbinning_v3_cent30to50_pT2to3,
      vnbinning_v3_cent30to50_pT3to4,   vnbinning_v3_cent30to50_pT4to5,
      vnbinning_v3_cent30to50_pT5to6,   vnbinning_v3_cent30to50_pT6to8,
      vnbinning_v3_cent30to50_pT8to10,  vnbinning_v3_cent30to50_pT10to15,
      vnbinning_v3_cent30to50_pT15to20, vnbinning_v3_cent30to50_pT20to40,
      vnbinning_v3_cent30to50_pT40to60, vnbinning_v3_cent30to50_pT60to100 }
};
