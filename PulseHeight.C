#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include "TNtuple.h"
#include "TStyle.h"
#include <iostream>
#include <tuple>
#include "FairLogger.h"
#include "TRDBase/Digit.h"
#include "TRDBase/TRDSimParam.h"
#include "TRDBase/TRDCommonParam.h"
#include "DataFormatsTRD/Constants.h"
#endif
using namespace o2::trd;
using namespace std;

void PulseHeight(string pathE="./Cosmics/cosmicsimages")
{
  
  vector<string> digifile = {"./Cosmics/data/trddigits_502238_0_00.root","./Cosmics/data/trddigits_502238_0_01.root",
  "./Cosmics/data/trddigits_502238_0_02.root","./Cosmics/data/trddigits_502238_0_03.root",
  "./Cosmics/data/trddigits_502238_0_04.root","./Cosmics/data/trddigits_502238_0_05.root",
  "./Cosmics/data/trddigits_502238_0_06.root","./Cosmics/data/trddigits_502238_0_07.root",
  "./Cosmics/data/trddigits_502238_0_08.root","./Cosmics/data/trddigits_502238_0_09.root",
  "./Cosmics/data/trddigits_502238_0_10.root","./Cosmics/data/trddigits_502238_0_11.root",
  "./Cosmics/data/trddigits_502238_0_12.root","./Cosmics/data/trddigits_502238_0_13.root",
  "./Cosmics/data/trddigits_502238_0_14.root","./Cosmics/data/trddigits_502238_0_15.root",
  "./Cosmics/data/trddigits_502238_0_16.root","./Cosmics/data/trddigits_502238_0_17.root",
  "./Cosmics/data/trddigits_502238_0_18.root","./Cosmics/data/trddigits_502238_0_19.root",
  "./Cosmics/data/trddigits_502238_0_20.root","./Cosmics/data/trddigits_502238_0_21.root","./Cosmics/data/trddigits_502238_1_00.root","./Cosmics/data/trddigits_502238_1_01.root",
  "./Cosmics/data/trddigits_502238_1_02.root","./Cosmics/data/trddigits_502238_1_03.root",
  "./Cosmics/data/trddigits_502238_1_04.root","./Cosmics/data/trddigits_502238_1_05.root",
  "./Cosmics/data/trddigits_502238_1_06.root","./Cosmics/data/trddigits_502238_1_07.root",
  "./Cosmics/data/trddigits_502238_1_08.root","./Cosmics/data/trddigits_502238_1_09.root",
  "./Cosmics/data/trddigits_502238_1_10.root","./Cosmics/data/trddigits_502238_1_11.root",
  "./Cosmics/data/trddigits_502238_1_12.root","./Cosmics/data/trddigits_502238_1_13.root",
  "./Cosmics/data/trddigits_502238_1_14.root","./Cosmics/data/trddigits_502238_1_15.root",
  "./Cosmics/data/trddigits_502238_1_16.root","./Cosmics/data/trddigits_502238_1_17.root",
  "./Cosmics/data/trddigits_502238_1_18.root","./Cosmics/data/trddigits_502238_1_19.root",
  "./Cosmics/data/trddigits_502238_1_20.root","./Cosmics/data/trddigits_502238_1_21.root"};

//vector<string> digifile = {"./muon_gun_10GeV/generated_files/trddigits.root"};

//vector<string> digifile = {"./Cosmics/data/trddigits_502238_0_00.root"};
  TFile* fin;
  TTree* digitTree;
  std::vector<Digit>* digitCont;

  int nev = 0;

  TProfile* ph = new TProfile("average pulse height", "Average Pulse Height", 30, -0.5, 29.5);

  TH2F* ph2D = new TH2F("Total Pulse Height", "Pulse Height Spectrum",30,-0.5,29.5,400,0.,400.);
    
  LOG(INFO) << nev << " entries found";

 
  int tbmax = 0;
  int tbhi = 0;
  int tblo = 0;

  int det = 0;
  int row = 0;
  int pad = 0;
  int channel = 0;
  map<tuple<int,int,int>,ArrayADC> dataMap;
  int tbsum[540][16][144]={{{0}}};

auto f = TFile::Open("ntuples.root","RECREATE"); //root file where tntuples are being stored
TNtuple *t = new TNtuple("nt","nt","sm:tb:ph");
  
for (const auto &item : digifile){
  fin = TFile::Open(item.data());
  digitTree = (TTree*)fin->Get("o2sim");
  digitCont = nullptr;
  digitTree->SetBranchAddress("TRDDigit", &digitCont);
  nev = digitTree->GetEntries();
  LOG(INFO) << nev << " entries found";

  
  for (int iev = 0; iev < nev; iev++) {

  for (int i =0; i<540;i++){
    for (int j = 0; j<16;j++){
      for (int k = 0 ; k<144;k++){
        tbsum[i][j][k] = 0;
      }
    }
  }
  digitTree->GetEvent(iev);
    for (const auto& digit : *digitCont) {
        // loop over det, pad, row?
        channel = digit.getChannel();
        if (channel == 0 || channel == 1 || channel ==20){continue;}

        auto adcs = digit.getADC();
        det = digit.getDetector();
        row = digit.getPadRow();
        pad = digit.getPadCol();
        int supermod = det/(o2::trd::constants::NSTACK*o2::trd::constants::NLAYER);
        tbsum[det][row][pad] = 0;

        for (int tb = 0; tb < o2::trd::constants::TIMEBINS; tb++) {
          ADC_t adc = adcs[tb];
          
          if (adc > (ADC_t)SimParam::instance()->getADCoutRange()) {
            LOG(INFO) << "Out of range ADC " << adc;
            continue;
          }
          tbsum[det][row][pad] += adc;

        }
          
        if (tbsum[det][row][pad]>=400){
          dataMap.insert(make_pair(make_tuple(det,row,pad), adcs));

        }
        else{
          tbsum[det][row][pad] = 0;
        }
      }// end digitcont
    
      for (int d=0;d<540;d++) {
        for (int r=0;r<16;r++) {
          for (int c=2;c<142;c++) {
            int s = d/(o2::trd::constants::NSTACK*o2::trd::constants::NLAYER);
            if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1]) {
              if (tbsum[d][r][c-1] > tbsum[d][r][c+1]) {
                tbmax = tbsum[d][r][c];
                tbhi = tbsum[d][r][c-1];
                tblo = tbsum[d][r][c+1];
                auto adcMax = dataMap.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMap.find(make_tuple(d,r,c-1));
                auto adcLo = dataMap.find(make_tuple(d,r,c+1));
                
                if (dataMap.find(make_tuple(d,r,c-2)) == dataMap.end()){
                  if (tblo > 400 ){
                    int phVal = 0;
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      t->Fill(s,tb,phVal);
                      ph->Fill(tb,phVal,1);
                      ph2D->Fill(tb,phVal);
                    }
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c-2));
                  if (tblo > 400){
                    int phVal = 0;
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      t->Fill(s,tb,phVal);
                      ph->Fill(tb,phVal,1);
                      ph2D->Fill(tb,phVal);
                    }
                  }
                }
                
              } 
              else {
                tbmax = tbsum[d][r][c];
                tbhi = tbsum[d][r][c+1];
                tblo = tbsum[d][r][c-1];
                auto adcMax = dataMap.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMap.find(make_tuple(d,r,c+1));
                auto adcLo = dataMap.find(make_tuple(d,r,c-1));
                if (dataMap.find(make_tuple(d,r,c+2))== dataMap.end()){
                  
                  if (tblo > 400){
                    int phVal = 0;
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      ph->Fill(tb,phVal,1);
                      t->Fill(s,tb,phVal);
                      ph2D->Fill(tb,phVal);
                    }
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c+2));
                  if (tblo > 400){
                    int phVal = 0;
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      t->Fill(s,tb,phVal);
                      ph->Fill(tb,phVal,1);
                      ph2D->Fill(tb,phVal);
                    }                    
                  }
                }
              }//end else
            }// end if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1])
          }  // end for c
        }//end for r
      }// end for d 
      dataMap.clear();
    } //end event 
    fin->Close();
  }
  
  t->Write();
  
  //myFile.close();


  TCanvas* ph2d_cnv = new TCanvas("ph2d", "", 800,600);
  
  ph2D->Draw("colz");
  ph->Draw("SAME");
 gStyle->SetOptStat(11);
  //ph->SetLineColor(kRed);
 // double scale = (ph->GetEntries()) / 30;
  //ph->Scale(1/scale);
  ph->SetMarkerColor(2);
   ph->SetMarkerSize(0.7);
   ph->SetMarkerStyle(21);
   ph2D->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2D->GetYaxis()->SetTitle("ADC value");
  TLegend* l = new TLegend(0.1,0.7,0.48,0.8);
  l->AddEntry(ph,"Average Pulse Height","lep");
  l->Draw();

  TLegend* border2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  border2->SetBorderSize(0); // no border
  border2->SetFillStyle(0);
  border2->SetFillColor(0); // Legend background should be white
  border2->SetTextFont(42);
  border2->SetTextSize(0.03); // Increase entry font size!
  //border2->AddEntry(ph, "pulse height", "C");
  border2->Draw();
  ph2d_cnv->SaveAs("./Cosmics/ph2d.pdf");
/**
 TCanvas* ph_cnv = new TCanvas("ph", "ph", 800,600);
  //ph2D->Draw("colz");
  ph->Draw("SAME");
  ph->SetLineColor(kBlack);
  double scale2 = (ph->GetEntries()) / 30;
  ph->Scale(1/scale2);
*/

TH2F* ph2DSM0 = new TH2F("ph2DSM0", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM1 = new TH2F("ph2DSM1", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM2 = new TH2F("ph2DSM2", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM3 = new TH2F("ph2DSM3", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM4 = new TH2F("ph2DSM4", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM5 = new TH2F("ph2DSM5", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM6 = new TH2F("ph2DSM6", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM7 = new TH2F("ph2DSM7", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM8 = new TH2F("ph2DSM8", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM9 = new TH2F("ph2DSM9", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM10 = new TH2F("ph2DSM10", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM11 = new TH2F("ph2DSM11", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM12 = new TH2F("ph2DSM12", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM13 = new TH2F("ph2DSM13", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM14 = new TH2F("ph2DSM14", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM15 = new TH2F("ph2DSM15", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM16 = new TH2F("ph2DSM16", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);
TH2F* ph2DSM17 = new TH2F("ph2DSM17", "Pulse Height Spectrum",30, -0.5, 29.5,400,0.,400.);

TProfile *hprof0 = new TProfile("hprof0","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof1 = new TProfile("hprof1","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof2 = new TProfile("hprof2","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof3 = new TProfile("hprof3","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof4 = new TProfile("hprof4","Profile of ph versus tb",30,-0.5,29.5,0,400);   
TProfile *hprof5 = new TProfile("hprof5","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof6 = new TProfile("hprof6","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof7 = new TProfile("hprof7","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof8 = new TProfile("hprof8","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof9 = new TProfile("hprof9","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof10 = new TProfile("hprof10","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof11 = new TProfile("hprof11","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof12 = new TProfile("hprof12","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof13 = new TProfile("hprof13","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof14 = new TProfile("hprof14","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof15 = new TProfile("hprof15","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof16 = new TProfile("hprof16","Profile of ph versus tb",30,-0.5,29.5,0,400);
TProfile *hprof17 = new TProfile("hprof17","Profile of ph versus tb",30,-0.5,29.5,0,400);

  //ph2dsm_cnv->Divide(3,3);
  //ph2dsm_cnv->cd(1);
  TCanvas* ph2dsm_cnv0 = new TCanvas("Pulse Height Spectrum: Supermodule 0", "Pulse Height Spectrum: Supermodule 0", 800,600);
  t->Draw("ph:tb>>ph2DSM0","sm==0","colz");  
  t->Draw("ph:tb>>hprof0","sm==0","prof hist SAME");
  //t->Draw("ph:tb>>ph2DSM0_ave","sm==0","SAME");
   //double scale = (hprof0->GetEntries()) / 30;
   //hprof0->Scale(1/scale);
   hprof0->SetMarkerColor(2);
   hprof0->SetMarkerSize(0.7);
   hprof0->SetMarkerStyle(21);

   //gStyle->SetOptStat(11);
   //gStyle->SetStats(0);
    ph2DSM0->SetTitle("Pulse Height Spectrum: Supermodule 0");
   ph2DSM0->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM0->GetYaxis()->SetTitle("ADC value");
  TLegend* l0 = new TLegend(0.1,0.7,0.48,0.8);
  l0->AddEntry(ph,"Average Pulse Height","lep");
  l0->Draw();
   //ph2dsm_cnv->Modified();
   //ph2dsm_cnv->Update();
 

  //ph2dsm_cnv->cd(2);
  TCanvas* ph2dsm_cnv1 = new TCanvas("Pulse Height Spectrum: Supermodule 1", "Pulse Height Spectrum: Supermodule 1", 800,600);
  t->Draw("ph:tb>>ph2DSM1","sm==1","colz");
  t->Draw("ph:tb>>hprof1","sm==1","prof hist SAME");
  gStyle->SetOptStat(0);
  hprof1->SetMarkerColor(2);
   hprof1->SetMarkerSize(0.7);
   hprof1->SetMarkerStyle(21);
   ph2DSM1->SetTitle("Pulse Height Spectrum: Supermodule 1");
   ph2DSM1->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM1->GetYaxis()->SetTitle("ADC value");
  TLegend* l1 = new TLegend(0.1,0.7,0.48,0.8);
  l1->AddEntry(ph,"Average Pulse Height","lep");
  l1->Draw();
  //ph2dsm_cnv->Modified();
  //ph2dsm_cnv->Update();

  //ph2dsm_cnv->cd(3);
TCanvas* ph2dsm_cnv2 = new TCanvas("Pulse Height Spectrum: Supermodule 2", "Pulse Height Spectrum: Supermodule 2", 800,600);
  t->Draw("ph:tb>>ph2DSM2","sm==2","colz");
  t->Draw("ph:tb>>hprof2","sm==2","prof hist SAME");
  hprof2->SetMarkerColor(2);
   hprof2->SetMarkerSize(0.7);
   hprof2->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM2->SetTitle("Pulse Height Spectrum: Supermodule 2");
  ph2DSM2->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM2->GetYaxis()->SetTitle("ADC value");
  TLegend* l2 = new TLegend(0.1,0.7,0.48,0.8);
  l2->AddEntry(ph,"Average Pulse Height","lep");
  l2->Draw();

  //ph2dsm_cnv->cd(4);
  TCanvas* ph2dsm_cnv3 = new TCanvas("Pulse Height Spectrum: Supermodule 3", "Pulse Height Spectrum: Supermodule 3", 800,600);
  t->Draw("ph:tb>>ph2DSM3","sm==3","colz");
  t->Draw("ph:tb>>hprof3","sm==3","prof hist SAME");
  hprof3->SetMarkerColor(2);
   hprof3->SetMarkerSize(0.7);
   hprof3->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM3->SetTitle("Pulse Height Spectrum: Supermodule 3");
  ph2DSM3->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM3->GetYaxis()->SetTitle("ADC value");
  TLegend* l3 = new TLegend(0.1,0.7,0.48,0.8);
  l3->AddEntry(ph,"Average Pulse Height","lep");
  l3->Draw();
  

  //ph2dsm_cnv->cd(5);
TCanvas* ph2dsm_cnv4 = new TCanvas("Pulse Height Spectrum: Supermodule 4", "Pulse Height Spectrum: Supermodule 4", 800,600);
  t->Draw("ph:tb>>ph2DSM4","sm==4","colz");
  t->Draw("ph:tb>>hprof4","sm==4","prof hist SAME");
  hprof4->SetMarkerColor(2);
   hprof4->SetMarkerSize(0.7);
   hprof4->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM4->SetTitle("Pulse Height Spectrum: Supermodule 4");
  ph2DSM4->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM4->GetYaxis()->SetTitle("ADC value");
  TLegend* l4 = new TLegend(0.1,0.7,0.48,0.8);
  l4->AddEntry(ph,"Average Pulse Height","lep");
  l4->Draw();

  //ph2dsm_cnv->cd(6);
TCanvas* ph2dsm_cnv5 = new TCanvas("Pulse Height Spectrum: Supermodule 5", "Pulse Height Spectrum: Supermodule 5", 800,600);
  t->Draw("ph:tb>>ph2DSM5","sm==5","colz");
  t->Draw("ph:tb>>hprof5","sm==5","prof hist SAME");
  hprof5->SetMarkerColor(2);
   hprof5->SetMarkerSize(0.7);
   hprof5->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM5->SetTitle("Pulse Height Spectrum: Supermodule 5");
  ph2DSM5->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM5->GetYaxis()->SetTitle("ADC value");
  TLegend* l5 = new TLegend(0.1,0.7,0.48,0.8);
  l5->AddEntry(ph,"Average Pulse Height","lep");
  l5->Draw();

  //ph2dsm_cnv->cd(7);
TCanvas* ph2dsm_cnv6 = new TCanvas("Pulse Height Spectrum: Supermodule 6", "Pulse Height Spectrum: Supermodule 6", 800,600);
  t->Draw("ph:tb>>ph2DSM6","sm==6","colz");
  t->Draw("ph:tb>>hprof6","sm==6","prof hist SAME");
  hprof6->SetMarkerColor(2);
   hprof6->SetMarkerSize(0.7);
   hprof6->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM6->SetTitle("Pulse Height Spectrum: Supermodule 6");
  ph2DSM6->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM6->GetYaxis()->SetTitle("ADC value");
  TLegend* l6 = new TLegend(0.1,0.7,0.48,0.8);
  l6->AddEntry(ph,"Average Pulse Height","lep");
  l6->Draw();

  //ph2dsm_cnv->cd(8);
TCanvas* ph2dsm_cnv7 = new TCanvas( "Pulse Height Spectrum: Supermodule 7", "Pulse Height Spectrum: Supermodule 7", 800,600);
  t->Draw("ph:tb>>ph2DSM7","sm==7","colz");
  t->Draw("ph:tb>>hprof7","sm==7","prof hist SAME");
  hprof7->SetMarkerColor(2);
   hprof7->SetMarkerSize(0.7);
   hprof7->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM7->SetTitle("Pulse Height Spectrum: Supermodule 7");
  ph2DSM7->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM7->GetYaxis()->SetTitle("ADC value");
  TLegend* l7 = new TLegend(0.1,0.7,0.48,0.8);
  l7->AddEntry(ph,"Average Pulse Height","lep");
  l7->Draw();

  //ph2dsm_cnv->cd(9);
TCanvas* ph2dsm_cnv8 = new TCanvas("Pulse Height Spectrum: Supermodule 8", "Pulse Height Spectrum: Supermodule 8", 800,600);
  t->Draw("ph:tb>>ph2DSM8","sm==8","colz");
  t->Draw("ph:tb>>hprof8","sm==8","prof hist SAME");
  hprof8->SetMarkerColor(2);
   hprof8->SetMarkerSize(0.7);
   hprof8->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM8->SetTitle("Pulse Height Spectrum: Supermodule 8");
  ph2DSM8->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM8->GetYaxis()->SetTitle("ADC value");
  TLegend* l8 = new TLegend(0.1,0.7,0.48,0.8);
  l8->AddEntry(ph,"Average Pulse Height","lep");
  l8->Draw();

 //TCanvas* ph2dsm_cnv2 = new TCanvas("ph2dSM:SM9-SM17", "ph2dSM:SM9-SM17", 800,600);
  //ph2dsm_cnv2->Divide(3,3);
  //ph2dsm_cnv2->cd(1);

TCanvas* ph2dsm_cnv9 = new TCanvas("Pulse Height Spectrum: Supermodule 9", "Pulse Height Spectrum: Supermodule 9", 800,600);
  t->Draw("ph:tb>>ph2DSM9","sm==9","colz");
  t->Draw("ph:tb>>hprof9","sm==9","prof hist SAME");
  hprof9->SetMarkerColor(2);
   hprof9->SetMarkerSize(0.7);
   hprof9->SetMarkerStyle(21);
gStyle->SetOptStat(0);
ph2DSM9->SetTitle("Pulse Height Spectrum: Supermodule 9");
  ph2DSM9->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM9->GetYaxis()->SetTitle("ADC value");
  TLegend* l9 = new TLegend(0.1,0.7,0.48,0.8);
  l9->AddEntry(ph,"Average Pulse Height","lep");
  l9->Draw();

  //ph2dsm_cnv2->cd(2);
TCanvas* ph2dsm_cnv10 = new TCanvas("Pulse Height Spectrum: Supermodule 10", "Pulse Height Spectrum: Supermodule 10", 800,600);
  t->Draw("ph:tb>>ph2DSM10","sm==10","colz");
  t->Draw("ph:tb>>hprof10","sm==10","prof hist SAME");
  hprof10->SetMarkerColor(2);
   hprof10->SetMarkerSize(0.7);
   hprof10->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM10->SetTitle("Pulse Height Spectrum: Supermodule 10");
  ph2DSM10->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM10->GetYaxis()->SetTitle("ADC value");
  TLegend* l10 = new TLegend(0.1,0.7,0.48,0.8);
  l10->AddEntry(ph,"Average Pulse Height","lep");
  l10->Draw();

  //ph2dsm_cnv2->cd(3);
TCanvas* ph2dsm_cnv11 = new TCanvas("Pulse Height Spectrum: Supermodule 11", "Pulse Height Spectrum: Supermodule 11", 800,600);
  t->Draw("ph:tb>>ph2DSM11","sm==11","colz");
  t->Draw("ph:tb>>hprof11","sm==11","prof hist SAME");
  hprof11->SetMarkerColor(2);
   hprof11->SetMarkerSize(0.7);
   hprof11->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM11->SetTitle("Pulse Height Spectrum: Supermodule 11");
  ph2DSM11->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM11->GetYaxis()->SetTitle("ADC value");
  TLegend* l11 = new TLegend(0.1,0.7,0.48,0.8);
  l11->AddEntry(ph,"Average Pulse Height","lep");
  l11->Draw();

  //ph2dsm_cnv2->cd(4);
TCanvas* ph2dsm_cnv12 = new TCanvas("Pulse Height Spectrum: Supermodule 12", "Pulse Height Spectrum: Supermodule 12", 800,600);
  t->Draw("ph:tb>>ph2DSM12","sm==12","colz");
  t->Draw("ph:tb>>hprof12","sm==12","prof hist SAME");
  hprof12->SetMarkerColor(2);
   hprof12->SetMarkerSize(0.7);
   hprof12->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM12->SetTitle("Pulse Height Spectrum: Supermodule 12");
  ph2DSM12->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM12->GetYaxis()->SetTitle("ADC value");
  TLegend* l12 = new TLegend(0.1,0.7,0.48,0.8);
  l12->AddEntry(ph,"Average Pulse Height","lep");
  l12->Draw();

  //ph2dsm_cnv2->cd(5);
TCanvas* ph2dsm_cnv13 = new TCanvas("Pulse Height Spectrum: Supermodule 13", "Pulse Height Spectrum: Supermodule 13", 800,600);
  t->Draw("ph:tb>>ph2DSM13","sm==13","colz");
  t->Draw("ph:tb>>hprof13","sm==13","prof hist SAME");
  hprof13->SetMarkerColor(2);
   hprof13->SetMarkerSize(0.7);
   hprof13->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM13->SetTitle("Pulse Height Spectrum: Supermodule 13");
  ph2DSM13->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM13->GetYaxis()->SetTitle("ADC value");
  TLegend* l13 = new TLegend(0.1,0.7,0.48,0.8);
  l13->AddEntry(ph,"Average Pulse Height","lep");
  l13->Draw();

  //ph2dsm_cnv2->cd(6);
TCanvas* ph2dsm_cnv14 = new TCanvas("Pulse Height Spectrum: Supermodule 14", "Pulse Height Spectrum: Supermodule 14", 800,600);
  t->Draw("ph:tb>>ph2DSM14","sm==14","colz");
  t->Draw("ph:tb>>hprof14","sm==14","prof hist SAME");
  hprof14->SetMarkerColor(2);
   hprof14->SetMarkerSize(0.7);
   hprof14->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM14->SetTitle("Pulse Height Spectrum: Supermodule 14");
  ph2DSM14->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM14->GetYaxis()->SetTitle("ADC value");
  TLegend* l14 = new TLegend(0.1,0.7,0.48,0.8);
  l14->AddEntry(ph,"Average Pulse Height","lep");
  l14->Draw();

  //ph2dsm_cnv2->cd(7);
TCanvas* ph2dsm_cnv15 = new TCanvas("Pulse Height Spectrum: Supermodule 15", "Pulse Height Spectrum: Supermodule 15", 800,600);
  t->Draw("ph:tb>>ph2DSM15","sm==15","colz");
  t->Draw("ph:tb>>hprof15","sm==15","prof hist SAME");
  hprof15->SetMarkerColor(2);
   hprof15->SetMarkerSize(0.7);
   hprof15->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM15->SetTitle("Pulse Height Spectrum: Supermodule 15");
  ph2DSM15->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM15->GetYaxis()->SetTitle("ADC value");
  TLegend* l15 = new TLegend(0.1,0.7,0.48,0.8);
  l15->AddEntry(ph,"Average Pulse Height","lep");
  l15->Draw();

  //ph2dsm_cnv2->cd(8);
TCanvas* ph2dsm_cnv16 = new TCanvas("Pulse Height Spectrum: Supermodule 16", "Pulse Height Spectrum: Supermodule 16", 800,600);
  t->Draw("ph:tb>>ph2DSM16","sm==16","colz");
  t->Draw("ph:tb>>hprof16","sm==16","prof hist SAME");
  hprof16->SetMarkerColor(2);
   hprof16->SetMarkerSize(0.7);
   hprof16->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM16->SetTitle("Pulse Height Spectrum: Supermodule 16");
  ph2DSM16->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM16->GetYaxis()->SetTitle("ADC value");
  TLegend* l16 = new TLegend(0.1,0.7,0.48,0.8);
  l16->AddEntry(ph,"Average Pulse Height","lep");
  l16->Draw();

  //ph2dsm_cnv2->cd(9);
  TCanvas* ph2dsm_cnv17 = new TCanvas("Pulse Height Spectrum: Supermodule 17", "Pulse Height Spectrum: Supermodule 17", 800,600);
  t->Draw("ph:tb>>ph2DSM17","sm==17","colz");
  t->Draw("ph:tb>>hprof17","sm==17","prof hist SAME");
  hprof17->SetMarkerColor(2);
   hprof17->SetMarkerSize(0.7);
   hprof17->SetMarkerStyle(21);
  gStyle->SetOptStat(0);
  ph2DSM17->SetTitle("Pulse Height Spectrum: Supermodule 17");
  ph2DSM17->GetXaxis()->SetTitle("time (time bin = 100ns)");
   ph2DSM17->GetYaxis()->SetTitle("ADC value");
  TLegend* l17 = new TLegend(0.1,0.7,0.48,0.8);
  l17->AddEntry(ph,"Average Pulse Height","lep");
  l17->Draw();
  
  f->Close();

}// end of macro