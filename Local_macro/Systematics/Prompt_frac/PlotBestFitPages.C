#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

void PlotBestFitPages()
{
	gROOT->SetBatch(kTRUE);

	const char *inFileName = "Output_Chi2_vs_SF_AllBins_Sept10.root";
	const char *tmpDir = "/tmp/PlotBestFitPages_tmp";
	const char *mergeScript = "merge_bestfit_grid.py";

	const Int_t N_CENTBIN = 5;
	const Int_t N_PTBIN = 9;
	const char *cent[N_CENTBIN] = {"cent0to10", "cent10to20", "cent20to30", "cent30to40", "cent40to50"};
	const char *pt[N_PTBIN] = {"pT2to3", "pT3to4", "pT4to5", "pT5to6", "pT6to8", "pT8to10", "pT10to15", "pT15to30", "pT30to100"};

	// One output PDF per canvas type, each pulled from the same 45 per-bin directories.
	const Int_t N_PDF = 4;
	const char *canvasName[N_PDF] = {"Canvas_BestFit", "Canvas_Chi2_vs_SF", "Canvas_Overlay_Shape", "Canvas_Overlay_Shape_Gen"};
	const char *outPdfName[N_PDF] = {"BestFit_AllBins.pdf", "Chi2vsSF_AllBins.pdf", "OverlayShape_AllBins.pdf", "OverlayShapeGen_AllBins.pdf"};
	Bool_t pdfOk[N_PDF];

	TFile *f = TFile::Open(inFileName);
	if (!f || f->IsZombie())
	{
		cout << "ERROR: could not open " << inFileName << endl;
		return;
	}

	// Each cell is saved straight from the source canvas as its own single-page vector PDF
	// (PDF, unlike PNG, has no rasterization loss). merge_bestfit_grid.py then arranges the
	// 45 cells into the final 3x3-grid-per-centrality, 5-page PDF (no re-rastering at any
	// step). This runs once per canvas type, producing one independent PDF each time.
	for (int ipdf = 0; ipdf < N_PDF; ipdf++)
	{
		gSystem->mkdir(tmpDir, kTRUE);

		for (int icent = 0; icent < N_CENTBIN; icent++)
		{
			for (int ipt = 0; ipt < N_PTBIN; ipt++)
			{
				TDirectory *dir = (TDirectory *)f->Get(Form("%s_%s", cent[icent], pt[ipt]));
				if (!dir)
				{
					cout << "Missing directory for " << cent[icent] << " " << pt[ipt] << endl;
					continue;
				}
				TCanvas *src = (TCanvas *)dir->Get(canvasName[ipdf]);
				if (!src)
				{
					cout << "Missing " << canvasName[ipdf] << " for " << cent[icent] << " " << pt[ipt] << endl;
					continue;
				}
				src->SaveAs(Form("%s/cell_%d_%d.pdf", tmpDir, icent, ipt));
			}
			cout << "Saved cell PDFs for " << cent[icent] << " (" << canvasName[ipdf] << ")" << endl;
		}

		Int_t rc = gSystem->Exec(Form("/usr/local/bin/python3.11 %s %s %s", mergeScript, tmpDir, outPdfName[ipdf]));
		gSystem->Exec(Form("rm -rf %s", tmpDir));

		pdfOk[ipdf] = (rc == 0);
		if (rc != 0)
			cout << "ERROR: merge script failed for " << outPdfName[ipdf] << " (exit code " << rc << ")" << endl;
		else
			cout << ">>> Done! Saved " << outPdfName[ipdf] << endl;
	}

	f->Close();

	// Final summary, printed with absolute paths so they can be copied straight into
	// "open <path>" or a file browser.
	TString outDir = gSystem->WorkingDirectory();
	cout << "\n>>> Output PDFs:" << endl;
	for (int ipdf = 0; ipdf < N_PDF; ipdf++)
		if (pdfOk[ipdf])
			cout << outDir << "/" << outPdfName[ipdf] << endl;
}
