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

void tbsumDigits(std::string digifile, string digifileP, string pathE, string pathP)
{

  TFile* fin = TFile::Open(digifile.data());
  TTree* digitTree = (TTree*)fin->Get("o2sim"); //= "./e_gun_2GeV/generated_files/trddigits.root"

  ofstream myFile (pathE+".csv");//= "./p_gun_2GeV/generated_files/trddigits.root"
  myFile<<"adcMax,adcHi,adcLo,adcHiNeighbour"<<"\n";

  auto f = TFile::Open("./e_gun_2GeV/ntuples.root","RECREATE"); //root file where tntuples are being stored

  std::vector<Digit>* digitCont = nullptr;
  digitTree->SetBranchAddress("TRDDigit", &digitCont);
  int nev = digitTree->GetEntries();


  TH1F* htbsum = new TH1F("tbsum", "Tbsum", 100, 0, 3000);
  TH1F* htbhi = new TH1F("htbhi", "Tbsum", 100, 0, 3000);
  TH1F* htblo = new TH1F("htblo", "Tbsum", 100, 0, 3000);
  TH1F* htbmax = new TH1F("htbmax", "Tbsum", 100, 0, 3000);
  TH1F* ph = new TH1F("pulse height", "Pulse height", 30, -0.5, 29.5);

  TH2F* padrow = new TH2F("padrow", "padrow",144,0.,144.,30,0.,30.);
  TCanvas* cnv = new TCanvas("cnv_padrow", "cnv_padrow", 800,600);
  LOG(INFO) << nev << " entries found";

  //TNtuple *t = new TNtuple("nt","nt","d:r:c:m:h:l");
  int tbmax = 0;
  int tbhi = 0;
  int tblo = 0;

  int det = 0;
  int row = 0;
  int pad = 0;
  int channel = 0;
   map<tuple<int,int,int>,ArrayADC> dataMap;
  int tbsum[540][16][144]={{{0}}};
  for (int iev = 0; iev < nev; iev++) {
    //tbsum = {{{0}}};
    digitTree->GetEvent(iev);
      for (const auto& digit : *digitCont) {
          // loop over det, pad, row?
          auto adcs = digit.getADC();
          det = digit.getDetector();
          row = digit.getPadRow();
          pad = digit.getPadCol();
          
          tbsum[det][row][pad] = 0;
          channel = digit.getChannel();
          dataMap.insert(make_pair(make_tuple(det,row,pad), adcs));
          if (channel == 0 || channel == 1 || channel ==20){continue;}
          else{
            for (int tb = 0; tb < o2::trd::constants::TIMEBINS; tb++) {
              ADC_t adc = adcs[tb];

              //cout<<adc<<endl;
              if (adc > (ADC_t)SimParam::instance()->getADCoutRange()) {
                LOG(INFO) << "Out of range ADC " << adc;
                continue;
              }
              tbsum[det][row][pad] += adc;
            }
            //cout<<endl;
            htbsum->Fill(tbsum[det][row][pad]);
          }
        }// end digitcont
      for (int d=0;d<540;d++) {
        for (int r=0;r<16;r++) {
          for (int c=2;c<142;c++) {
            
            if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1]) {
              if (tbsum[d][r][c-1] > tbsum[d][r][c+1]) {
                tbmax = tbsum[d][r][c];
                tbhi = tbsum[d][r][c-1];
                tblo = tbsum[d][r][c+1];
                //cout<<c-2<<" "<<c-1<<" "<<c<<" "<<" "<<c+1<<endl;
                auto adcMax = dataMap.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMap.find(make_tuple(d,r,c-1));
                auto adcLo = dataMap.find(make_tuple(d,r,c+1));
                
                if (dataMap.find(make_tuple(d,r,c-2))== dataMap.end()){
                  if (tblo != 0){
                    htbmax->Fill(tbsum[d][r][c]);
                    htbhi->Fill(tbsum[d][r][c-1]);
                    htblo->Fill(tbsum[d][r][c+1]);
                    //t->Fill(d,r,c,tbmax,tbhi,tblo);
                    int phVal = 0;
                    //padrow->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      //cout<<tbsum[d][r][c-2]<<endl;
                      //cout<<d<<" "<<r<<" "<<c-2<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      //cout<<(adcHiNeighbour->second)[tb]<<endl;
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      
                      ph->Fill(tb,phVal);
                      /**
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c-1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c+1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c-2, tb, (adcHiNeighbour->second)[tb]); */
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                    }
                    //cout<<endl;
                    myFile<<"\n";
                    
                    //padrow->Draw("colz");
                    //cnv->SaveAs(Form("./pion_gun_sim/padrow_%03d_%02d.pdf",d,r));
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c-2));
                  if (tblo != 0){
                    htbmax->Fill(tbsum[d][r][c]);
                    htbhi->Fill(tbsum[d][r][c-1]);
                    htblo->Fill(tbsum[d][r][c+1]);
                    //t->Fill(d,r,c,tbmax,tbhi,tblo);
                    int phVal = 0;
                    //padrow->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      //cout<<tbsum[d][r][c-2]<<endl;
                      //cout<<d<<" "<<r<<" "<<c-2<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      //cout<<(adcHiNeighbour->second)[tb]<<endl;
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      
                      ph->Fill(tb,phVal);
                      /**
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c-1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c+1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c-2, tb, (adcHiNeighbour->second)[tb]); */
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                    }
                    //cout<<endl;
                    myFile<<"\n";
                    
                    //padrow->Draw("colz");
                    //cnv->SaveAs(Form("./pion_gun_sim/padrow_%03d_%02d.pdf",d,r));
                  }
                }
                
              } 
              else {
                tbmax = tbsum[d][r][c];
                tbhi = tbsum[d][r][c+1];
                tblo = tbsum[d][r][c-1];
                //cout<<c-1<<" "<<c<<" "<<c+1<<" "<<" "<<c+2<<endl;
                auto adcMax = dataMap.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMap.find(make_tuple(d,r,c+1));
                auto adcLo = dataMap.find(make_tuple(d,r,c-1));
                //auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c+2));
                if (dataMap.find(make_tuple(d,r,c+2))== dataMap.end()){
                  
                  if (tblo != 0){
                    
                    htbmax->Fill(tbsum[d][r][c]);
                    htbhi->Fill(tbsum[d][r][c+1]);
                    htblo->Fill(tbsum[d][r][c-1]);
                    
                    //t->Fill(d,r,c,tbmax,tbhi,tblo);
                    int phVal = 0;
                    //padrow->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      //cout<<(adcHiNeighbour->second)[tb]<<endl;;
                      //cout<<tbsum[d][r][c-2]<<endl;
                      //cout<<d<<" "<<r<<" "<<c-2<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      ph->Fill(tb,phVal);
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<" "<<(adcHiNeighbour->second)[tb]<<endl;

                      /**
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c+1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c-1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c+2, tb, (adcHiNeighbour->second)[tb]); */
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                    }
                    //cout<<endl;
                    myFile << "\n";
                    //padrow->Draw("colz");
                    //cnv->SaveAs(Form("./pion_gun_sim/padrow_%03d_%02d.pdf",d,r));
                  }
                }
                else{
                  auto adcHiNeighbour = dataMap.find(make_tuple(d,r,c+2));
                  if (tblo != 0){
                    
                    htbmax->Fill(tbsum[d][r][c]);
                    htbhi->Fill(tbsum[d][r][c+1]);
                    htblo->Fill(tbsum[d][r][c-1]);
                    
                    //t->Fill(d,r,c,tbmax,tbhi,tblo);
                    int phVal = 0;
                    //padrow->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      //cout<<(adcHiNeighbour->second)[tb]<<endl;;
                      //cout<<tbsum[d][r][c-2]<<endl;
                      //cout<<d<<" "<<r<<" "<<c-2<<" "<<(adcHiNeighbour->second)[tb]<<endl;
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      ph->Fill(tb,phVal);
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<" "<<(adcHiNeighbour->second)[tb]<<endl;

                      /**
                      padrow->Fill(c, tb, (adcMax->second)[tb]);
                      padrow->Fill(c+1, tb, (adcHi->second)[tb]);
                      padrow->Fill(c-1, tb, (adcLo->second)[tb]);
                      padrow->Fill(c+2, tb, (adcHiNeighbour->second)[tb]); */
                      myFile<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                    }
                    //cout<<endl;
                    myFile << "\n";
                    //padrow->Draw("colz");
                    //cnv->SaveAs(Form("./pion_gun_sim/padrow_%03d_%02d.pdf",d,r));
                  }
                }

              }//end else
            }// end if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1])
          }  // end for c
          
        }//end for r
      }// end for d 
    } //end event 

  myFile.close();

  //t->Write();
/**
  TCanvas* c3 = new TCanvas("c3", "e-gun TB Sum", 600, 600);
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


  TLegend* border = new TLegend(0.7, 0.7, 0.9, 0.9);
  border->SetBorderSize(0); // no border
  border->SetFillStyle(0);
  border->SetFillColor(0); // Legend background should be white
  border->SetTextFont(42);
  border->SetTextSize(0.03); // Increase entry font size!
  border->AddEntry(htbsum, "htbsum", "l");
  border->AddEntry(htbmax, "htbmax", "l");
  border->AddEntry(htbhi, "htbhi", "l");
  border->AddEntry(htblo, "htblo", "l");
  border->Draw();

  c3->SaveAs("./e_gun_2GeV/tbsum.pdf");
*/
  //TCanvas* c4 = new TCanvas("c4", "Pulse Height", 600, 600);
  //ph->Draw();
  
  //c4->SaveAs("pulseheight.pdf");*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* finP = TFile::Open(digifileP.data());
  TTree* digitTreeP = (TTree*)finP->Get("o2sim");

  ofstream myFileP (pathP+".csv");
  myFileP<<"adcMax,adcHi,adcLo,adcHiNeighbour"<<"\n";

  //auto fP = TFile::Open("./p_gun_2GeV/ntuples.root","RECREATE"); //root file where tntuples are being stored

  std::vector<Digit>* digitContP = nullptr;
  digitTreeP->SetBranchAddress("TRDDigit", &digitContP);
  int nevP = digitTreeP->GetEntries();

  TH1F* htbsumP = new TH1F("tbsum", "Tbsum", 100, 0, 3000);
  TH1F* htbhiP = new TH1F("htbhi", "Tbsum", 100, 0, 3000);
  TH1F* htbloP = new TH1F("htblo", "Tbsum", 100, 0, 3000);
  TH1F* htbmaxP = new TH1F("htbmax", "Tbsum", 100, 0, 3000);
  TH1F* phP = new TH1F("pulse height", "Pulse height", 30, -0.5, 29.5);

  //TH2F* padrowP = new TH2F("pion padrow", "pion padrow",144,0.,144.,30,0.,30.);
  //TCanvas* cnvP = new TCanvas("cnv_padrowP", "cnv_padrowP", 800,600);
  
  LOG(INFO) << nevP << " entries found";

  //TNtuple *tP = new TNtuple("nt","nt","d:r:c:m:h:l");
  int tbmaxP = 0;
  int tbhiP = 0;
  int tbloP = 0;

  int detP = 0;
  int rowP = 0;
  int padP = 0;
  int channelP = 0;

  map<tuple<int,int,int>,ArrayADC> dataMapP;

  for (int i =0; i<540;i++){
    for (int j = 0; j<16;j++){
      for (int k = 0 ; k<144;k++){
        tbsum[i][j][k] = 0;
      }
    }
  }

  for (int iev = 0; iev < nevP; iev++) {
    digitTreeP->GetEvent(iev);
      for (const auto& digit : *digitContP) {

          // loop over det, pad, row?
          auto adcs = digit.getADC();
          detP = digit.getDetector();
          rowP = digit.getPadRow();
          padP = digit.getPadCol();
          tbsum[detP][rowP][padP] = 0;
          channelP = digit.getChannel();
          dataMapP.insert(make_pair(make_tuple(detP,rowP,padP), adcs));
          if (channelP == 0 || channelP == 1 || channelP ==20){continue;}
          else{
            for (int tb = 0; tb < o2::trd::constants::TIMEBINS; tb++) {
              ADC_t adc = adcs[tb];
              if (adc > (ADC_t)SimParam::instance()->getADCoutRange()) {
                LOG(INFO) << "Out of range ADC " << adc;
                continue;
              }
              tbsum[detP][rowP][padP] += adc;
            }
            //htbsumP->Fill(tbsum[detP][rowP][padP]);
          }
        }// end digitcont
        
      for (int d=0;d<540;d++) {
        for (int r=0;r<16;r++) {
          for (int c=2;c<142;c++) {
            
            if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1]) {
              if (tbsum[d][r][c-1] > tbsum[d][r][c+1]) {
                tbmaxP = tbsum[d][r][c];
                tbhiP = tbsum[d][r][c-1];
                tbloP = tbsum[d][r][c+1];
                //cout<<c-2<<" "<<c-1<<" "<<c<<" "<<" "<<c+1<<endl;
                auto adcMax = dataMapP.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMapP.find(make_tuple(d,r,c-1));
                auto adcLo = dataMapP.find(make_tuple(d,r,c+1));
                if (dataMapP.find(make_tuple(d,r,c-2))== dataMapP.end()){
                  if (tbloP != 0){
                  
                    htbmaxP->Fill(tbsum[d][r][c]);
                    htbhiP->Fill(tbsum[d][r][c-1]);
                    htbloP->Fill(tbsum[d][r][c+1]); 
                    //tP->Fill(d,r,c,tbmaxP,tbhiP,tbloP);
                    int phVal = 0;
                    //padrowP->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<endl;
                      
                      phP->Fill(tb,phVal);
                      /**
                      padrowP->Fill(c, tb, (adcMax->second)[tb]);
                      padrowP->Fill(c-1, tb, (adcHi->second)[tb]);
                      padrowP->Fill(c+1, tb, (adcLo->second)[tb]);
                      padrowP->Fill(c-2, tb, (adcHiNeighbour->second)[tb]); */
                      
                      myFileP<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                    }
                    myFileP<<"\n";
                    
                    //padrowP->Draw("colz");
                  }
                }
                else{
                  auto adcHiNeighbour = dataMapP.find(make_tuple(d,r,c-2));
                  if (tbloP != 0){
                      htbmaxP->Fill(tbsum[d][r][c]);
                      htbhiP->Fill(tbsum[d][r][c-1]);
                      htbloP->Fill(tbsum[d][r][c+1]); 
                      //tP->Fill(d,r,c,tbmaxP,tbhiP,tbloP);
                      int phVal = 0;
                      //padrowP->Reset();
                      for (int tb = 0 ; tb<30;tb++){
                        phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                        //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<endl;
                        
                        phP->Fill(tb,phVal);
                        /**
                        padrowP->Fill(c, tb, (adcMax->second)[tb]);
                        padrowP->Fill(c-1, tb, (adcHi->second)[tb]);
                        padrowP->Fill(c+1, tb, (adcLo->second)[tb]);
                        padrowP->Fill(c-2, tb, (adcHiNeighbour->second)[tb]); */
                        
                        myFileP<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                      }
                      myFileP<<"\n";
                      
                      //padrowP->Draw("colz");
                    }
                  }
                  
              } 
              else {
                //padrow = new TH2F("padrow", "padrow",(c+2)-(c-1),c-1,c+2,30,0.,30.);
                tbmaxP = tbsum[d][r][c];
                tbhiP = tbsum[d][r][c+1];
                tbloP = tbsum[d][r][c-1];
                //cout<<c-1<<" "<<c<<" "<<c+1<<" "<<" "<<c+2<<endl;
                auto adcMax = dataMapP.find(make_tuple(d,r,c)) ;
                auto adcHi = dataMapP.find(make_tuple(d,r,c+1));
                auto adcLo = dataMapP.find(make_tuple(d,r,c-1));
                
                if (dataMapP.find(make_tuple(d,r,c+2))== dataMapP.end()){
                  if (tbloP != 0){
                  
                  htbmaxP->Fill(tbsum[d][r][c]);
                  htbhiP->Fill(tbsum[d][r][c+1]);
                  htbloP->Fill(tbsum[d][r][c-1]);
                  
                  //tP->Fill(d,r,c,tbmaxP,tbhiP,tbloP);
                  int phVal = 0;
                  //padrowP->Reset();
                  for (int tb = 0 ; tb<30;tb++){
                    phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                    //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<endl;
                    
                    phP->Fill(tb,phVal);
                    /**
                    padrowP->Fill(c, tb, (adcMax->second)[tb]);
                    padrowP->Fill(c+1, tb, (adcHi->second)[tb]);
                    padrowP->Fill(c-1, tb, (adcLo->second)[tb]);
                    padrowP->Fill(c+2, tb, (adcHiNeighbour->second)[tb]); */
                  
                    myFileP<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<0<<"\n";
                  }
                  myFileP << "\n";
                  //padrowP->Draw("colz");
                }
              }
                else{
                  auto adcHiNeighbour = dataMapP.find(make_tuple(d,r,c+2));
                  if (tbloP != 0){
                    htbmaxP->Fill(tbsum[d][r][c]);
                    htbhiP->Fill(tbsum[d][r][c+1]);
                    htbloP->Fill(tbsum[d][r][c-1]);
                    
                    //tP->Fill(d,r,c,tbmaxP,tbhiP,tbloP);
                    int phVal = 0;
                    //padrowP->Reset();
                    for (int tb = 0 ; tb<30;tb++){
                      phVal = ((adcMax->second)[tb] + (adcHi->second)[tb] + (adcLo->second)[tb]);
                      //cout<<(adcMax->second)[tb]<<" "<<(adcHi->second)[tb]<<" "<<(adcLo->second)[tb]<<endl;
                      
                      phP->Fill(tb,phVal);
                      /**
                      padrowP->Fill(c, tb, (adcMax->second)[tb]);
                      padrowP->Fill(c+1, tb, (adcHi->second)[tb]);
                      padrowP->Fill(c-1, tb, (adcLo->second)[tb]);
                      padrowP->Fill(c+2, tb, (adcHiNeighbour->second)[tb]); */
                    
                      myFileP<<(adcMax->second)[tb]<<","<<(adcHi->second)[tb]<<","<<(adcLo->second)[tb]<<","<<(adcHiNeighbour->second)[tb]<<"\n";
                    }
                    myFileP << "\n";
                    //padrowP->Draw("colz");
                  }
                }
              }//end else
            }// end if (tbsum[d][r][c]>tbsum[d][r][c-1] && tbsum[d][r][c]>tbsum[d][r][c+1])
          }  // end for c
        }//end for r
      }// end for d 
    } //end event 

  myFileP.close();

  //tP->Write();
  /**
  TCanvas* c3P = new TCanvas("c3P", "p-gun TB Sum", 600, 600);
  gPad->SetLogy();
  //c3->Divide(1,2);

  htbmaxP->SetLineColor(kRed);
  htbloP->SetLineColor(kBlue);
  htbhiP->SetLineColor(kGreen); 
  htbsumP->SetLineColor(kBlack);

  //c3->cd(1);
  htbsumP->Draw();
  htbmaxP->Draw("SAME");
  htbhiP->Draw("SAME");
  htbloP->Draw("SAME"); 
  //c3->cd(2);


  TLegend* borderP = new TLegend(0.7, 0.7, 0.9, 0.9);
  borderP->SetBorderSize(0); // no border
  borderP->SetFillStyle(0);
  borderP->SetFillColor(0); // Legend background should be white
  borderP->SetTextFont(42);
  borderP->SetTextSize(0.03); // Increase entry font size!
  borderP->AddEntry(htbsumP, "htbsum", "l");
  borderP->AddEntry(htbmaxP, "htbmax", "l");
  borderP->AddEntry(htbhiP, "htbhi", "l");
  borderP->AddEntry(htbloP, "htblo", "l");
  borderP->Draw();

  c3->SaveAs("./p_gun_2GeV/tbsum.pdf");

  TCanvas* c4 = new TCanvas("c4", "Pulse Height", 600, 600);
  ph->Draw();
  ph->SetLineColor(kRed);
  double scale = (ph->GetEntries()) / 30;
  ph->Scale(1/scale);
  
  TLegend* border2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  border2->SetBorderSize(0); // no border
  border2->SetFillStyle(0);
  border2->SetFillColor(0); // Legend background should be white
  border2->SetTextFont(42);
  border2->SetTextSize(0.03); // Increase entry font size!
  border2->AddEntry(ph, "electron gun pulse height", "C");
  border2->Draw();


  phP->Draw("SAME");
  phP->SetLineColor(kBlue);
  double scaleP = (phP->GetEntries()) / 30;
  phP->Scale(1/scaleP);
  TLegend* border2P = new TLegend(0.7, 0.7, 0.9, 0.9);
  border2P->SetBorderSize(0); // no border
  border2P->SetFillStyle(0);
  border2P->SetFillColor(0); // Legend background should be white
  border2P->SetTextFont(42);
  border2P->SetTextSize(0.03); // Increase entry font size!
  //border2P->SetHeader("Pulse Height Spectra","C");
  border2P->AddEntry(phP, "Pion gun pulse height", "lep");
   border2P->AddEntry(ph, "Electron gun pulse height", "lep");
  border2P->Draw();
  c4->SaveAs("pulseheight.pdf");
*/
}// end of macro
