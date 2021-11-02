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

//constexpr int kMINENTRIES = 100;

void tbsumDigits(string pathE="./Cosmics/cosmicsimages")
{
/**
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
  "./Cosmics/data/trddigits_502238_1_20.root","./Cosmics/data/trddigits_502238_1_21.root"}; */
  //TFile* fin = TFile::Open(digifile.data());
  //vector<string> digifile = {"./Cosmics/data/trddigits_502238_0_00.root"};
 vector<string> digifile = {"./muon_gun_10GeV/generated_files/trddigits.root"};
  auto f = TFile::Open("ntuples.root","RECREATE"); //root file where tntuples are being stored

  TFile* fin;
  //TTree* digitTree = (TTree*)fin->Get("o2sim"); //= "./e_gun_2GeV/generated_files/trddigits.root" string digifileP="Cosmics/data/trddigits_502238_0_01.root"
  TTree* digitTree;
  ofstream myFile (pathE+".csv"); //= "./p_gun_2GeV/generated_files/trddigits.root" , string pathP="./Cosmics/cosmics2"
  myFile<<"adcMax,adcHi,adcLo,adcHiNeighbour"<<"\n";



  //auto f2 = TFile::Open("./Cosmics/padrows.root","RECREATE");

  //std::vector<Digit>* digitCont = nullptr;
  std::vector<Digit>* digitCont;
  //digitTree->SetBranchAddress("TRDDigit", &digitCont);
  //int nev = digitTree->GetEntries();
  int nev = 0;

  TH1F* htbsum = new TH1F("tbsum", "Histogram of Sums Over Time Bins", 100, 0, 10000);
  TH1F* htbhi = new TH1F("htbhi", "Histogram of Sums Over Time Bins", 100, 0, 10000);
  TH1F* htblo = new TH1F("htblo", "Histogram of Sums Over Time Bins", 100, 0, 10000);
  TH1F* htbmax = new TH1F("htbmax", "Histogram of Sums Over Time Bins", 100, 0, 10000);
  TH1F* ADCSpectrum = new TH1F("ADC Spectrum", "ADC Spectrum (7400<tbsum<7600)", 40, 230, 270);

  TH2F* padrow = new TH2F("padrow", "padrow",144,0.,144.,30,0.,30.);
  
  TCanvas* cnv = new TCanvas("cnv_padrow", "cnv_padrow", 800,600);
  TCanvas* spectrum = new TCanvas("ADC Spectrum", "ADC Spectrum", 800,600);
  
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
        if(tbsum[det][row][pad]>=7400 && tbsum[det][row][pad]<=7600){
          for (int tb = 0; tb<o2::trd::constants::TIMEBINS; tb++){
            ADCSpectrum->Fill(adcs[tb]);
          }
          
        }  
        if (tbsum[det][row][pad]>=400){
          dataMap.insert(make_pair(make_tuple(det,row,pad), adcs));
          htbsum->Fill(tbsum[det][row][pad]);
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
                    htbmax->Fill(tbmax);
                    htbhi->Fill(tbhi);
                    htblo->Fill(tblo);
                    int phVal = 0;
                    padrow->Reset();
                    //myFile<<d<<" "<<r<<" "<<c<<" "<<tbmax<<" "<<tbhi<<" "<<tblo<<"\n";
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c-1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c+1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c-2, tb, 0); 
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                    }
                    myFile<<"\n";
                    
                    string title = "det: "+to_string(d)+", row: "+to_string(r);
                    padrow->SetTitle(title.c_str());
                    padrow->SetStats(0);
                    //padrow->Write();
                    padrow->Draw("colz");
                    padrow->GetXaxis()->SetRange(c-3,c+3);
                    
                    padrow->Draw("text,same");
                    cnv->Modified();
                    cnv->Update();
                    cnv->SaveAs(Form("./muon_gun_10GeV/padrows/padrow_%04d_%03d_%02d.pdf",iev,d,r));
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c-2));
                  if (tblo > 400){
                    htbmax->Fill(tbmax);
                    htbhi->Fill(tbhi);
                    htblo->Fill(tblo);
                    int phVal = 0;
                    padrow->Reset();
                    //myFile<<d<<" "<<r<<" "<<c<<" "<<tbmax<<" "<<tbhi<<" "<<tblo<<"\n";
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c-1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c+1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c-2, tb, (adcHiNeighbour->second)[tb]); 
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                    }
                    myFile<<"\n";
                    string title = "det: "+to_string(d)+", row: "+to_string(r);
                    padrow->SetTitle(title.c_str());
                    padrow->SetStats(0);
                    //padrow->Write();
                    padrow->Draw("colz");
                    //padrow->GetXaxis()->SetLimits(c-2,c+20);
                    padrow->GetXaxis()->SetRange(c-3,c+3);

                    padrow->Draw("text,same");
                    cnv->Modified();
                    cnv->Update();
                    cnv->SaveAs(Form("./muon_gun_10GeV/padrows/padrow_%04d_%03d_%02d.pdf",iev,d,r));
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
                    htbmax->Fill(tbmax);
                    htbhi->Fill(tbhi);
                    htblo->Fill(tblo);
                    
                    int phVal = 0;
                    padrow->Reset();
                   // myFile<<d<<" "<<r<<" "<<c<<" "<<tbmax<<" "<<tbhi<<" "<<tblo<<"\n";
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c+1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c-1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c+2, tb,0); 
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                    }
                    //cout<<endl;
                    myFile << "\n";
                    string title = "det: "+to_string(d)+", row: "+to_string(r);
                    padrow->SetTitle(title.c_str());
                    padrow->SetStats(0);
                    //padrow->Write();
                    padrow->Draw("colz");
                    //padrow->GetXaxis()->SetLimits(c-2,c+20);
                    padrow->GetXaxis()->SetRange(c-3,c+3);
                    padrow->Draw("text,same");
                    cnv->Modified();
                    cnv->Update();
                    cnv->SaveAs(Form("./muon_gun_10GeV/padrows/padrow_%04d_%03d_%02d.pdf",iev,d,r));
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c+2));
                  if (tblo > 400){
                    //cout<<tbmax<<endl;
                    htbmax->Fill(tbmax);
                    htbhi->Fill(tbhi);
                    htblo->Fill(tblo);
                    int phVal = 0;
                    padrow->Reset();
                    //myFile<<d<<" "<<r<<" "<<c<<" "<<tbmax<<" "<<tbhi<<" "<<tblo<<"\n";
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c+1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c-1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c+2, tb, (adcHiNeighbour->second)[tb]); 
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                    }
                    myFile << "\n";
                    
                    
                    //padrow->GetXaxis()->SetLimits(c-2,c+20);
                    string title = "det: "+to_string(d)+", row: "+to_string(r);
                    padrow->SetTitle(title.c_str());
                    padrow->SetStats(0);
                    //padrow->Write();
                    padrow->Draw("colz");
                    padrow->GetXaxis()->SetRange(c-3,c+3);
                    padrow->Draw("text,same");
                    cnv->Modified();
                    cnv->Update();
                    cnv->SaveAs(Form("./muon_gun_10GeV/padrows/padrow_%04d_%03d_%02d.pdf",iev,d,r));
                  }
                }

              }//end else
            }// end if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1])
          }  // end for c
        }//end for r
      }// end for d 
      dataMap.clear();
    } //end event 
  }
  ADCSpectrum->Draw();
  ADCSpectrum->GetXaxis()->SetTitle("ADC Value");
  ADCSpectrum->GetYaxis()->SetTitle("Number of counts");
  ADCSpectrum->SaveAs("./Cosmics/ADCSpectrum.pdf");

  //spectrum->Close();
  myFile.close();

  TCanvas* c3 = new TCanvas("c3", "TB Sum", 600, 600);
  gPad->SetLogy();
  //c3->Divide(1,2);

  htbmax->SetLineColor(kRed);
  htblo->SetLineColor(kBlue);
  htbhi->SetLineColor(kGreen); 
  htbsum->SetLineColor(kBlack);

  //c3->cd(1);
  
  htbsum->Draw();
  htbmax->Draw("SAME");
  htbhi->Draw("SAME");
  htblo->Draw("SAME"); 
  //c3->cd(2); 


  TLegend* border = new TLegend(0.1,0.7,0.48,0.8);
  //border->SetBorderSize(0); // no border
  border->SetFillStyle(0);
  border->SetFillColor(0); // Legend background should be white
  border->SetTextFont(42);
  border->SetTextSize(0.03); // Increase entry font size!
  border->AddEntry(htbsum, "Total Time Bin Sum", "l");
  border->AddEntry(htbmax, "Maximum Time Bin Sum", "l");
  border->AddEntry(htbhi, "Time Bin Sum: High Neighbour", "l");
  border->AddEntry(htblo, "Time Bin Sum: Low Neigbour", "l");
  border->Draw();
  htbsum->GetXaxis()->SetTitle("sum of ADC Value over 30 time bins");
  htbsum->GetYaxis()->SetTitle("Number of counts");
  gStyle->SetOptStat(0);
  c3->SaveAs("./Cosmics/tbsum.pdf");

}// end of macro
