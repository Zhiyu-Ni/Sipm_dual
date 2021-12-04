
{

    
    
    const int entry=6;
    double energy[]={1 , 5, 10, 30, 60, 100};
    double energy1[]={ 5, 10, 30, 60, 100};
    double u330_u330[]={23.98 , 22.74, 22.5, 22.3, 22.3, 22.36};
      double u330_ug5[]={28.172,27.749,27.467,27.496,27.530};
    
    TCanvas * cc=new TCanvas();
    
    TGraph *g2 = new TGraph(6,energy,u330_u330);
    g2->SetName("g2");
    g2->SetMarkerStyle(20);
    //g2->Draw("AP+");
    g2->SetMarkerColor(kBlue);
    g2->SetTitle("U330/U330");
    g2->SetDrawOption("AP+");
    
    TGraph *g = new TGraph(5,energy1,u330_ug5 );
    //g->Draw("AP");
    //g->GetXaxis()->SetTitle("#theta (degrees)");
    //g->GetYaxis()->SetTitle("C/S (a.u.)");
    g->SetName("g");
    g->SetDrawOption("AP+");
    g->SetMarkerColor(kRed);
    g->SetMarkerStyle(20);
    g->SetTitle("U330/UG5");
    

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
