
{
//=========Macro generated from canvas: c1_n3/c1_n3
//=========  (Wed Jul 15 09:45:47 2020) by ROOT version 6.16/00
    /*
   TCanvas *c1_n3 = new TCanvas("c1_n3", "c1_n3",101,26,700,500);
   c1_n3->Range(-75.625,-14.08847,75.625,126.7963);
   c1_n3->SetFillColor(0);
   c1_n3->SetBorderMode(0);
   c1_n3->SetBorderSize(2);
   c1_n3->SetFrameBorderMode(0);
   c1_n3->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("C photons detected at both ends");
   
   Double_t Graph_fx1[23] = {
   -55,
   -50,
   -45,
   -40,
   -35,
   -30,
   -25,
   -20,
   -15,
   -10,
   -5,
   0,
   5,
   10,
   15,
   20,
   25,
   30,
   35,
   40,
   45,
   50,
   55};
   Double_t Graph_fy1[23] = {
   107.5474,
   88.9,
   74.7,
   59.512,
   49.644,
   42.196,
   35.828,
   32.596,
   26.16,
   23.032,
   19.76,
   17.416,
   14.116,
   10.472,
   7.692,
   8.412,
   8.288,
   10.344,
   12.012,
   13.592,
   15.612,
   23.916,
   32.888};
   TGraph *graph = new TGraph(23,Graph_fx1,Graph_fy1);
   graph->SetName("Graph");
   graph->SetTitle("Left side C (noise)");
   graph->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Left side C (noise)",100,-66,66);
   Graph_Graph1->SetMinimum(6.9228);
   Graph_Graph1->SetMaximum(117.533);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx2[23] = {
   -55,
   -50,
   -45,
   -40,
   -35,
   -30,
   -25,
   -20,
   -15,
   -10,
   -5,
   0,
   5,
   10,
   15,
   20,
   25,
   30,
   35,
   40,
   45,
   50,
   55};
   Double_t Graph_fy2[23] = {
   17.064,
   12.436,
   9.128,
   7.392,
   5.772,
   5.388,
   4.924,
   4.976,
   4.34,
   5.176,
   8.128,
   9.788,
   12.828,
   14.54,
   16.676,
   20.064,
   22.288,
   28.14,
   30.616,
   36.348,
   39.588,
   51.996,
   59.8};
   graph = new TGraph(23,Graph_fx2,Graph_fy2);
   graph->SetName("Graph");
   graph->SetTitle("Right side C (signal)");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Right side C (signal)",100,-66,66);
   Graph_Graph2->SetMinimum(3.906);
   Graph_Graph2->SetMaximum(65.346);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleOffset(1);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   multigraph->Add(graph,"");
   multigraph->Draw("AP");
   multigraph->GetXaxis()->SetTitle("#theta (degrees)");
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelSize(0.035);
   multigraph->GetXaxis()->SetTitleSize(0.035);
   multigraph->GetXaxis()->SetTitleOffset(1);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("C photons");
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelSize(0.035);
   multigraph->GetYaxis()->SetTitleSize(0.035);
   multigraph->GetYaxis()->SetTitleFont(42);
   
   TLegend *leg = new TLegend(0.1014493,0.106383,0.4014493,0.316383,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","Left side C (noise)","lpf");
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","Right side C (signal)","lpf");
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   c1_n3->Modified();
   c1_n3->cd();
   c1_n3->SetSelected(c1_n3);
    */
    /*
    //old emission spec
    Double_t L_C[23] = {
    107.5474, 88.9, 74.7, 59.512, 49.644, 42.196, 35.828, 32.596, 26.16, 23.032, 19.76, 17.416, 14.116, 10.472, 7.692, 8.412, 8.288, 10.344, 12.012, 13.592, 15.612, 23.916, 32.888};
    Double_t R_C[23] = {
    17.064, 12.436, 9.128, 7.392, 5.772, 5.388, 4.924, 4.976, 4.34, 5.176, 8.128, 9.788, 12.828, 14.54, 16.676, 20.064, 22.288, 28.14, 30.616, 36.348, 39.588, 51.996, 59.8};
    Double_t L_S[23] = {
    251.324, 171.18, 125.952, 94.196, 76.272, 62.152, 53.076, 51.004, 44.236, 39.532, 39.468, 38.516, 41.364, 41.244, 42.692, 50.06, 53.8, 66.616, 76.12, 92.596, 114.412, 167.796, 241.636};
    Double_t R_S[23] = {
    828.44, 576.332, 428.356, 325.508, 262.936, 216.408, 188.32, 178.568, 159.588, 141.372, 144.212, 139.468, 148.588, 151.304, 157.04, 185.636, 198.624, 248.868, 288.268, 355.596, 441.196, 660.252, 939.0121};
    */
    Double_t L_C[23] = {   117.8614,
        96.89604,
        77.07426,
        59.81188,
        52.34653,
        43.4604,
        38.31188,
        32.37624,
        26.9505,
        24.65842,
        21.64851,
        18.55446,
        12.77228,
        11.25248,
        10.0198,
        8.980198,
        10.60396,
        12.16337,
        13.09406,
        16.41089,
        20.79208,
        28.72277,
        38.89109};
    Double_t L_C_e[23] = {};
    Double_t R_C[23] = {   19.4802,
        13.56931,
        10.22277,
        7.881188,
        6.668317,
        5.925743,
        5.643564,
        5.019802,
        5.009901,
        5.89604,
        8.351485,
        10.36634,
        11.87129,
        14.07426,
        17.69307,
        18.67327,
        23.13861,
        27.21782,
        30.69307,
        35.97525,
        43.73762,
        55.5495,
        62.88119};
    Double_t R_C_e[23] = { };
    Double_t L_S[23] = {   1113.495,
        806.8911,
        569.4455,
        409.8218,
        345.0941,
        278.9703,
        245.3168,
        218.7327,
        196.4406,
        182.4257,
        184.3713,
        170.6089,
        162.8267,
        183.7871,
        205.9505,
        195.4307,
        244.3267,
        286.4703,
        327.8614,
        413.7574,
        552.0644,
        785.3762,
        1117.609};
    Double_t R_S[23] = {   343.5743,
        252.3614,
        182.1733,
        133.1188,
        112.3614,
        91.33168,
        81.02475,
        73.77723,
        64.54455,
        61.9901,
        61.92079,
        57.81188,
        55.29703,
        61.54455,
        70.65842,
        68.17327,
        85.99505,
        99.66832,
        114.6436,
        146.104,
        198.2624,
        285.0198,
        409.6782};
    
    
    
    
    const int entry=23;
    double angle_the[]={-54.69 , -50.139 , -44.981 , -40.126 , -35.12 , -29.81 , -25.259 , -19.798 , -14.791 , -10.088 , -5.082 , -0.228 , 4.779 , 10.088 , 14.943 , 19.949 , 25.259 , 29.962 , 34.665 , 39.671 , 44.981 , 49.987 , 54.69};
    double C_S_ratio[]={0.90844 , 0.90413 , 0.90642 , 0.957 , 0.98125 , 0.99451 , 1.01436 , 1.02104 , 1.03211 , 1.02122 , 1.0762 , 1.22777 , 1.56595 , 1.79437 , 2.07547 , 2.31266 , 2.39838 , 2.32821 , 2.09778 , 1.90907 , 1.60181 , 1.34504 , 1.23975};

    double angle_y[23];
    double s_in;
    double c_in;
    s_in = 0.1;
    c_in = 0.2;
    double sum_1=0,sum_2=0;;
    for(int i=0;i<23;i++){
        angle_y[i] = (R_C[i] + s_in*R_S[i]) / (L_S[i]*2 + c_in*L_C[i]);
        //angle_y[i] = (R_S[i]) / (L_S[i]);
        if(i<9){
          sum_1+=angle_y[i];
          sum_2+=C_S_ratio[i];
        }
    }
    
    for(int i=0;i<23;i++){
        angle_y[i] = angle_y[i] *(sum_2/sum_1);
    }
    
    
    TCanvas * cc=new TCanvas();
    
    TGraph *g2 = new TGraph(23,angle_the,C_S_ratio);
    g2->SetMarkerStyle(20);
    g2->Draw("AP+");
    g2->SetMarkerColor(kRed);
    g2->SetTitle("paper");

    TGraph *g = new TGraph(23,angle_the, angle_y);
    g->Draw("AP");
    g->GetXaxis()->SetTitle("#theta (degrees)");
    g->GetYaxis()->SetTitle("C/S (a.u.)");
    g->SetMarkerColor(kBlue);
    g->SetMarkerStyle(20);
    g->SetTitle("Simu");
    

    TMultiGraph *gg=new TMultiGraph();
    gg->Add(g);
    gg->Add(g2);
    gg->Draw("AP");
    gg->GetXaxis()->SetTitle("#theta (degrees)");
    gg->GetYaxis()->SetTitle("C/S (a.u.)");
    gg->SetTitle(Form("S_%g_C%g",s_in,c_in));
    cc->BuildLegend();
    //cc->SaveAs(Form("S_%g_C%g.png",s_in,c_in));
    //delete cc;
    //}
    //}
    
    
}
