#pragma once

#include "vnbinning_generated.h" 
#include "ESE_Cuts_MB11to21_Jan22.h"

// --- Mass distribution fit settings ---
const int N_MASSBINS = 48; //48 
Float_t fit_range_low = 1.75;
Float_t fit_range_high = 1.99;
Double_t D0_mass = 1.8648;
Float_t hist_range_low = 1.75;
Float_t hist_range_high = 1.99;
Float_t width = (hist_range_high-hist_range_low)/N_MASSBINS; //-->0.005

//--- Bin count constants ---
static const int N_QBINS      = 10;
static const int N_CENTBINS_1 = 80;
static const int N_CENTBINS   = 6;
static const int N_PTBINS     = 12;
//static const int N_VBINS      = 62;
static const int N_VBINS_V2   = 62;   // <-- v2 specific
static const int N_VBINS_V3   = 66;   // <--- Equal to VNBINNING array size -1 
static const int N_BINS       = 2300;

// ---- Centrality bin definitions ---

static const int min_centbin[N_CENTBINS] = {0, 10, 20, 30, 40, 50};
static const int max_centbin[N_CENTBINS] = {10, 20, 30, 40, 50, 80};
static const int cen_edges[N_CENTBINS+1] = {0, 10, 20, 30, 40, 50, 80};

static const char* cen_name[N_CENTBINS] = {
  "cent0to10", "cent10to20", "cent20to30", "cent30to40",
  "cent40to50", "cent50to80"
};
TString cen_label[N_CENTBINS]={"cent0to10","cent10to20","cent20to30","cent30to40","cent40to50","cent50to80"};

static const char* centbin[N_CENTBINS] = {
    "0 <= cent < 10 %",  "10 <= cent < 20 %", "20 <= cent < 30 %",
    "30 <= cent < 40 %", "40 <= cent < 50 %",  "50 <= cent < 80 %"
};



// --- pT bin definitions ----
static const double pt_edges[N_PTBINS+1] = {
    1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 40, 60, 100
};
static const char* pt_name[N_PTBINS] = {
    "pT1to2",   "pT2to3",   "pT3to4",   "pT4to5",
    "pT5to6",   "pT6to8",   "pT8to10",  "pT10to15",
    "pT15to20", "pT20to40", "pT40to60", "pT60to100"
};
static const char* pTbin[N_PTBINS] = {
    "1.0 < pT < 2.0 GeV",   "2.0 < pT < 3.0 GeV",
    "3.0 < pT < 4.0 GeV",   "4.0 < pT < 5.0 GeV",
    "5.0 < pT < 6.0 GeV",   "6.0 < pT < 8.0 GeV",
    "8.0 < pT < 10.0 GeV",  "10.0 < pT < 15.0 GeV",
    "15.0 < pT < 20.0 GeV", "20.0 < pT < 40.0 GeV",
    "40.0 < pT < 60.0 GeV", "60.0 < pT < 100.0 GeV"
};

// --- vnbinning v2 lookup table [N_CENTBINS][N_PTBINS] ----

static const double* vnb_v2_table[N_CENTBINS][N_PTBINS] = {
    // cent0to10
    { vnbinning_v2_cent0to10_pT1to2,   vnbinning_v2_cent0to10_pT2to3,
      vnbinning_v2_cent0to10_pT3to4,   vnbinning_v2_cent0to10_pT4to5,
      vnbinning_v2_cent0to10_pT5to6,   vnbinning_v2_cent0to10_pT6to8,
      vnbinning_v2_cent0to10_pT8to10,  vnbinning_v2_cent0to10_pT10to15,
      vnbinning_v2_cent0to10_pT15to20, vnbinning_v2_cent0to10_pT20to40,
      vnbinning_v2_cent0to10_pT40to60, vnbinning_v2_cent0to10_pT60to100 },
    // cent10to20
    { vnbinning_v2_cent10to20_pT1to2,   vnbinning_v2_cent10to20_pT2to3,
      vnbinning_v2_cent10to20_pT3to4,   vnbinning_v2_cent10to20_pT4to5,
      vnbinning_v2_cent10to20_pT5to6,   vnbinning_v2_cent10to20_pT6to8,
      vnbinning_v2_cent10to20_pT8to10,  vnbinning_v2_cent10to20_pT10to15,
      vnbinning_v2_cent10to20_pT15to20, vnbinning_v2_cent10to20_pT20to40,
      vnbinning_v2_cent10to20_pT40to60, vnbinning_v2_cent10to20_pT60to100 },
    // cent20to30
    { vnbinning_v2_cent20to30_pT1to2,   vnbinning_v2_cent20to30_pT2to3,
      vnbinning_v2_cent20to30_pT3to4,   vnbinning_v2_cent20to30_pT4to5,
      vnbinning_v2_cent20to30_pT5to6,   vnbinning_v2_cent20to30_pT6to8,
      vnbinning_v2_cent20to30_pT8to10,  vnbinning_v2_cent20to30_pT10to15,
      vnbinning_v2_cent20to30_pT15to20, vnbinning_v2_cent20to30_pT20to40,
      vnbinning_v2_cent20to30_pT40to60, vnbinning_v2_cent20to30_pT60to100 },
    // cent30to40
    { vnbinning_v2_cent30to40_pT1to2,   vnbinning_v2_cent30to40_pT2to3,
      vnbinning_v2_cent30to40_pT3to4,   vnbinning_v2_cent30to40_pT4to5,
      vnbinning_v2_cent30to40_pT5to6,   vnbinning_v2_cent30to40_pT6to8,
      vnbinning_v2_cent30to40_pT8to10,  vnbinning_v2_cent30to40_pT10to15,
      vnbinning_v2_cent30to40_pT15to20, vnbinning_v2_cent30to40_pT20to40,
      vnbinning_v2_cent30to40_pT40to60, vnbinning_v2_cent30to40_pT60to100 },
    // cent40to50
    { vnbinning_v2_cent40to50_pT1to2,   vnbinning_v2_cent40to50_pT2to3,
      vnbinning_v2_cent40to50_pT3to4,   vnbinning_v2_cent40to50_pT4to5,
      vnbinning_v2_cent40to50_pT5to6,   vnbinning_v2_cent40to50_pT6to8,
      vnbinning_v2_cent40to50_pT8to10,  vnbinning_v2_cent40to50_pT10to15,
      vnbinning_v2_cent40to50_pT15to20, vnbinning_v2_cent40to50_pT20to40,
      vnbinning_v2_cent40to50_pT40to60, vnbinning_v2_cent40to50_pT60to100 },
    // cent50to80
    { vnbinning_v2_cent50to80_pT1to2,   vnbinning_v2_cent50to80_pT2to3,
      vnbinning_v2_cent50to80_pT3to4,   vnbinning_v2_cent50to80_pT4to5,
      vnbinning_v2_cent50to80_pT5to6,   vnbinning_v2_cent50to80_pT6to8,
      vnbinning_v2_cent50to80_pT8to10,  vnbinning_v2_cent50to80_pT10to15,
      vnbinning_v2_cent50to80_pT15to20, vnbinning_v2_cent50to80_pT20to40,
      vnbinning_v2_cent50to80_pT40to60, vnbinning_v2_cent50to80_pT60to100 }

};

//------ vnbinning v3 lookup table [N_CENTBINS][N_PTBINS] -------

static const double* vnb_v3_table[N_CENTBINS][N_PTBINS] = {
    // cent0to10
    { vnbinning_v3_cent0to10_pT1to2,   vnbinning_v3_cent0to10_pT2to3,
      vnbinning_v3_cent0to10_pT3to4,   vnbinning_v3_cent0to10_pT4to5,
      vnbinning_v3_cent0to10_pT5to6,   vnbinning_v3_cent0to10_pT6to8,
      vnbinning_v3_cent0to10_pT8to10,  vnbinning_v3_cent0to10_pT10to15,
      vnbinning_v3_cent0to10_pT15to20, vnbinning_v3_cent0to10_pT20to40,
      vnbinning_v3_cent0to10_pT40to60, vnbinning_v3_cent0to10_pT60to100 },
    // cent10to20
    { vnbinning_v3_cent10to20_pT1to2,   vnbinning_v3_cent10to20_pT2to3,
      vnbinning_v3_cent10to20_pT3to4,   vnbinning_v3_cent10to20_pT4to5,
      vnbinning_v3_cent10to20_pT5to6,   vnbinning_v3_cent10to20_pT6to8,
      vnbinning_v3_cent10to20_pT8to10,  vnbinning_v3_cent10to20_pT10to15,
      vnbinning_v3_cent10to20_pT15to20, vnbinning_v3_cent10to20_pT20to40,
      vnbinning_v3_cent10to20_pT40to60, vnbinning_v3_cent10to20_pT60to100 },
    // cent20to30
    { vnbinning_v3_cent20to30_pT1to2,   vnbinning_v3_cent20to30_pT2to3,
      vnbinning_v3_cent20to30_pT3to4,   vnbinning_v3_cent20to30_pT4to5,
      vnbinning_v3_cent20to30_pT5to6,   vnbinning_v3_cent20to30_pT6to8,
      vnbinning_v3_cent20to30_pT8to10,  vnbinning_v3_cent20to30_pT10to15,
      vnbinning_v3_cent20to30_pT15to20, vnbinning_v3_cent20to30_pT20to40,
      vnbinning_v3_cent20to30_pT40to60, vnbinning_v3_cent20to30_pT60to100 },
    // cent30to40
    { vnbinning_v3_cent30to40_pT1to2,   vnbinning_v3_cent30to40_pT2to3,
      vnbinning_v3_cent30to40_pT3to4,   vnbinning_v3_cent30to40_pT4to5,
      vnbinning_v3_cent30to40_pT5to6,   vnbinning_v3_cent30to40_pT6to8,
      vnbinning_v3_cent30to40_pT8to10,  vnbinning_v3_cent30to40_pT10to15,
      vnbinning_v3_cent30to40_pT15to20, vnbinning_v3_cent30to40_pT20to40,
      vnbinning_v3_cent30to40_pT40to60, vnbinning_v3_cent30to40_pT60to100 },
    // cent40to50
    { vnbinning_v3_cent40to50_pT1to2,   vnbinning_v3_cent40to50_pT2to3,
      vnbinning_v3_cent40to50_pT3to4,   vnbinning_v3_cent40to50_pT4to5,
      vnbinning_v3_cent40to50_pT5to6,   vnbinning_v3_cent40to50_pT6to8,
      vnbinning_v3_cent40to50_pT8to10,  vnbinning_v3_cent40to50_pT10to15,
      vnbinning_v3_cent40to50_pT15to20, vnbinning_v3_cent40to50_pT20to40,
      vnbinning_v3_cent40to50_pT40to60, vnbinning_v3_cent40to50_pT60to100 },
    // cent50to80
    { vnbinning_v3_cent50to80_pT1to2,   vnbinning_v3_cent50to80_pT2to3,
      vnbinning_v3_cent50to80_pT3to4,   vnbinning_v3_cent50to80_pT4to5,
      vnbinning_v3_cent50to80_pT5to6,   vnbinning_v3_cent50to80_pT6to8,
      vnbinning_v3_cent50to80_pT8to10,  vnbinning_v3_cent50to80_pT10to15,
      vnbinning_v3_cent50to80_pT15to20, vnbinning_v3_cent50to80_pT20to40,
      vnbinning_v3_cent50to80_pT40to60, vnbinning_v3_cent50to80_pT60to100 }

};
