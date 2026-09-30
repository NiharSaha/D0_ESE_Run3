#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TList.h>
#include <TParameter.h>
#include <iostream>
#include <fstream>

using namespace std;

// Read-only logger: pulls the already-fitted Yield/YieldError/Significance out of
// HIST_DATA_DCA_Aug26.root (the FitInfo TList written by get_Data_DCA_fromMass.C for
// every cent/pT/DCA bin) and prints/saves a plain table. Does not touch or re-run
// get_Data_DCA_fromMass.C.
void SignalSignificanceLog()
{
	const char *inFileName = "HIST_DATA_DCA_Aug26.root";
	TFile *f = TFile::Open(inFileName);
	if (!f || f->IsZombie())
	{
		cout << "ERROR: could not open " << inFileName << endl;
		return;
	}

	const int N_CENTBIN = 5;
	const char *cent[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};

	const int PTBIN = 9;
	const char *pt[PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};

	const int dca_bins = 10;
	Double_t DCA[dca_bins + 1] = {0, 0.0011, 0.0024, 0.0039, 0.0059, 0.0085, 0.01179, 0.016, 0.0214, 0.0366, 0.135};
	const char *label_dca[dca_bins] = {"dca_0to0p0011", "dca_0p0011to0p0024", "dca_0p0024to0p0039", "dca_0p0039to0p0059", "dca_0p0059to0p0085", "dca_0p0085to0p01179", "dca_0p01179to0p016", "dca_0p016to0p0214", "dca_0p0214to0p0366", "dca_0p0366to0p135"};

	ofstream logFile("SignalSignificance_Log.txt");

	for (int ipt = 0; ipt < PTBIN; ipt++)
	{
		TString hdr = Form("\n===== pT: %s =====", pt[ipt]);
		cout << hdr << endl;
		logFile << hdr << "\n";

		for (int icent = 0; icent < N_CENTBIN; icent++)
		{
			TDirectory *dir = (TDirectory *)f->Get(Form("MassFits_%s_%s", cent[icent], pt[ipt]));
			if (!dir)
				continue;

			TString sub = Form("-- %s --", cent[icent]);
			cout << sub << endl;
			logFile << sub << "\n";

			for (int idca = 0; idca < dca_bins; idca++)
			{
				TH1F *h = (TH1F *)dir->Get(Form("dmass_fit_%s_%s_%s", cent[icent], pt[ipt], label_dca[idca]));
				if (!h)
					continue;
				TList *info = (TList *)h->GetListOfFunctions()->FindObject("FitInfo");
				if (!info)
					continue;

				auto getp = [&](const char *name) -> Double_t
				{
					TParameter<Double_t> *p = (TParameter<Double_t> *)info->FindObject(name);
					return p ? p->GetVal() : 0.0;
				};

				Double_t yield = getp("Yield");
				Double_t yerr = getp("YieldError");
				Double_t sbSig = getp("Significance"); // sideband S/sqrt(S+B), from the tiering logic
				Double_t fitSig = (yerr > 0) ? yield / yerr : 0.0; // Yield/YieldError, as requested

				char buf[256];
				snprintf(buf, sizeof(buf),
				         "  [%d] DCA %.4f-%.4f : Yield = %10.1f +/- %8.1f   Sig(Yield/Err) = %8.2f   Sig(sideband) = %8.2f",
				         idca, DCA[idca], DCA[idca + 1], yield, yerr, fitSig, sbSig);
				cout << buf << endl;
				logFile << buf << "\n";
			}
		}
	}

	logFile.close();
	f->Close();
	cout << "\n>>> Wrote SignalSignificance_Log.txt" << endl;
}
