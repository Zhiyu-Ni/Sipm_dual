
{

    
    
    const int entry=6;
    double energy[]={1 , 5, 10, 30, 60, 100};
    double energy1[]={1, 5, 10, 30, 60, 100};
    double u330_u330[]={};
      double sipm_sipm[]={};
    
    TCanvas * cc=new TCanvas();
    
    TGraph *g2 = new TGraph(6,energy,u330_u330);
    g2->SetName("g2");
    g2->SetMarkerStyle(20);
    //g2->Draw("AP+");
    g2->SetMarkerColor(kBlue);
    g2->SetTitle("U330/U330");
    g2->SetDrawOption("AP+");
    
    TGraph *g = new TGraph(5,energy1, sipm_sipm);
    //g->Draw("AP");
    //g->GetXaxis()->SetTitle("#theta (degrees)");
    //g->GetYaxis()->SetTitle("C/S (a.u.)");
    g->SetName("g");
    g->SetDrawOption("AP+");
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(20);
    g->SetTitle("Sipm/Sipm");
    

    TMultiGraph *gg=new TMultiGraph();
    gg->Add(g);
    gg->Add(g2);
    gg->Draw("AP");
    gg->GetXaxis()->SetTitle("Energy(Gev)");
    gg->GetYaxis()->SetTitle("Mean number of C photons per Gev");
    gg->GetXaxis()->SetLimits(0.,110.);
    gg->SetMinimum(15.);
    gg->SetMaximum(50.);
    //gg->SetTitle(Form("S_%g_C%g",s_in,c_in));
    cc->BuildLegend();


    //cc->SaveAs(Form("S_%g_C%g.png",s_in,c_in));
    //delete cc;
    //}
    //}
    
    
}
