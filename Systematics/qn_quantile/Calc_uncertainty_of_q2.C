/bin/bash: q: command not found
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"

#include "TMath.h"
#include "TComplex.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include <TNtuple.h>

#include "/home/chand140/ESE_analysis_May5_2024/CMSSW_10_3_3_patch1/src/include.h"
using namespace std;

void Calc_uncertainty_of_q2()
{
    Float_t cen_min[100], cen_max[100], evt_plane[6][100], evt_plane_err;
    Float_t maximum_uncer[100], max_diff, diff, min_evt_plane;
    
    ifstream inf1 ("RescorSave/Rescor_1_326544/trackmid2_1.dat");
    ifstream inf2 ("RescorSave/Rescor_326545_326619/trackmid2_1.dat");
    ifstream inf3 ("RescorSave/Rescor_326620_326886/trackmid2_1.dat");
    ifstream inf4 ("RescorSave/Rescor_326887_327146/trackmid2_1.dat");
    ifstream inf5 ("RescorSave/Rescor_327147_327229/trackmid2_1.dat");
    ifstream inf6 ("RescorSave/Rescor_327230_327999/trackmid2_1.dat");

    for (int i=0; i<100; i++)
    {
        inf1 >> cen_min[i] >> cen_max[i] >> evt_plane[0][i] >> evt_plane_err;
        inf2 >> cen_min[i] >> cen_max[i] >> evt_plane[1][i] >> evt_plane_err;
        inf3 >> cen_min[i] >> cen_max[i] >> evt_plane[2][i] >> evt_plane_err;
        inf4 >> cen_min[i] >> cen_max[i] >> evt_plane[3][i] >> evt_plane_err;
        inf5 >> cen_min[i] >> cen_max[i] >> evt_plane[4][i] >> evt_plane_err;
        inf6 >> cen_min[i] >> cen_max[i] >> evt_plane[5][i] >> evt_plane_err;
    }

    for (int i=0; i<100; i++)
    {
        //Find the maximum difference between the mean values
        max_diff = 0.0;    
        for (int j=0; j<5; j++) {
            diff = fabs(evt_plane[0][i]-evt_plane[j+1][i]);
            if (diff>max_diff) max_diff=diff;
        }
            
        for (int k=1; k<5; k++) {
            diff = fabs(evt_plane[1][i]-evt_plane[k+1][i]);
            if (diff>max_diff) max_diff=diff;
        }

        for (int l=2; l<5; l++) {
            diff = fabs(evt_plane[2][i]-evt_plane[l+1][i]);
            if (diff>max_diff) max_diff=diff;
        }

        for (int m=3; m<5; m++) {
            diff = fabs(evt_plane[3][i]-evt_plane[m+1][i]);
            if (diff>max_diff) max_diff=diff;
        }

        diff = fabs(evt_plane[4][i]-evt_plane[5][i]);
        if (diff>max_diff) max_diff=diff;

        //Find the least mean value
        min_evt_plane = 10000.0;
        for (int j=0; j<6; j++) {
            if (evt_plane[j][i]<min_evt_plane) min_evt_plane=evt_plane[j][i];
        }

        //Calculate the maximum uncertainty
        maximum_uncer[i]=1.0*max_diff/min_evt_plane;
        cout << "Relative uncertainty for cent bin " << i << " is: " << 100.0*maximum_uncer[i] << "\%" << endl;
    }

    //Calculate new q2 boundaries based on the variations:
    //Change the q2 edges by +/- relative amount 

    //Define new edges:
    Float_t q2_binning_up[N_CENTBINS_1][N_Q2BINS+1], q2_binning_down[N_CENTBINS_1][N_Q2BINS+1];

    for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
        for (int i_q2=0; i_q2<N_Q2BINS+1; i_q2++) {
            q2_binning_up[i_cen][i_q2] = q2_binning[i_cen][i_q2]*(1.0+maximum_uncer[i_cen]); //increase q2 cuts
            q2_binning_down[i_cen][i_q2] = q2_binning[i_cen][i_q2]*(1.0-maximum_uncer[i_cen]); //decrease q2 cuts
        }
    }

    cout << endl << endl << "Float_t q2_binning_up[N_CENTBINS_1][N_Q2BINS+1]={";
    for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
        cout << "{";
        for (int i_q2=0; i_q2<N_Q2BINS+1; i_q2++) {
            cout << q2_binning_up[i_cen][i_q2];
            if (i_q2==N_Q2BINS) continue;
            cout << ",";
        }
        cout << "}";
        if (i_cen==N_CENTBINS_1-1) continue;
        cout << ",";
    }
    cout << "};" << endl << endl;

    cout << "Float_t q2_binning_down[N_CENTBINS_1][N_Q2BINS+1]={";
    for (int i_cen=0; i_cen<N_CENTBINS_1; i_cen++) {
        cout << "{";
        for (int i_q2=0; i_q2<N_Q2BINS+1; i_q2++) {
            cout << q2_binning_down[i_cen][i_q2];
            if (i_q2==N_Q2BINS) continue;
            cout << ",";
        }
        cout << "}";
        if (i_cen==N_CENTBINS_1-1) continue;
        cout << ",";
    }
    cout << "};" << endl << endl;
    
}
