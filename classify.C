void classify(){
  TFile *A=TFile::Open("noEff/Flow_output_cen0_noEff.root"), *B=TFile::Open("withEff/Flow_output_cen0_withEff.root");
  const char* v2pt[9]={"pT2to3","pT3to4","pT4to5","pT5to6","pT6to8","pT8to10","pT10to15","pT15to30","pT30to100"};
  // per-position counters: 0=outer third, 1=middle, 2=central
  int nOK[3]={0},nFail[3]={0},nMid[3]={0};
  printf("%-10s %2s %3s %8s %8s %6s %7s %7s  %s\n","pT","q","iv","Y_noE","relE_noE","1/k","ratio","class","");
  for(int ip=0;ip<9;ip++) for(int iq=0;iq<12;iq++){
    TH1D *ma=(TH1D*)B->Get(Form("mass_fitted/h_mass_cent0to10_q2bin%d_%s",iq,v2pt[ip]));
    TH1D *ha=(TH1D*)A->Get(Form("histograms_SP/h_v2_hist_cent0to10_q2bin%d_%s",iq,v2pt[ip]));
    TH1D *hb=(TH1D*)B->Get(Form("histograms_SP/h_v2_hist_cent0to10_q2bin%d_%s",iq,v2pt[ip]));
    if(!ma||!ha||!hb) continue;
    double sw=0,sw2=0; for(int b=ma->FindBin(1.73);b<=ma->FindBin(2.0);b++){sw+=ma->GetBinContent(b); sw2+=pow(ma->GetBinError(b),2);}
    if(sw<=0) continue; double invk=sqrt(sw/sw2);   // expected ratio if WL failed
    int n=ha->GetNbinsX();
    for(int i=1;i<=n;i++){
      double ya=ha->GetBinContent(i),ea=ha->GetBinError(i),yb=hb->GetBinContent(i),eb=hb->GetBinError(i);
      if(!(ya>0&&yb>0&&ea>0&&eb>0)) continue;
      double r=(eb/yb)/(ea/ya);
      // log-distance to the two hypotheses
      double dF=fabs(log(r/invk)), dO=fabs(log(r));
      const char* cls = (dF<0.3*fabs(log(invk))) ? "FAILED" : (dO<0.3*fabs(log(invk)) ? "OK" : "between");
      double pos=fabs((i-0.5)-n/2.0)/(n/2.0); int zone = pos>0.66?0:(pos>0.33?1:2);
      if(!strcmp(cls,"FAILED")) nFail[zone]++; else if(!strcmp(cls,"OK")) nOK[zone]++; else nMid[zone]++;
      if(ip<3 && iq==0) printf("%-10s %2d %3d %8.1f %8.3f %6.3f %7.3f %7s\n",v2pt[ip],iq,i-1,ya,ea/ya,invk,r,cls);
    }
  }
  const char* zn[3]={"outer SP bins  ","middle SP bins ","central SP bins"};
  printf("\nSummary over all v2 (q, pT) bins, cent0to10, 1M entries:\n");
  for(int z=0;z<3;z++) printf("  %s  WL OK: %3d   WL FAILED (L-like error): %3d   in between: %3d\n",zn[z],nOK[z],nFail[z],nMid[z]);
}
