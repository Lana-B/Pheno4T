#include "SampleAnalyzer/User/Analyzer/CMS_SUS_13_013.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool CMS_SUS_13_013::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    INFO << "        <>  Analysis: CMS-SUS-13-016, JHEP 01 (2014) 163                              <>" << endmsg;
    INFO << "        <>   (Search for new physics in events with same-sign                         <>" << endmsg;
    INFO << "        <>     dileptons and jets in pp collisions at √s = 8 TeV)                     <>"<< endmsg;
    INFO << "        <>  Recasted by: L.Beck, D.Dobur, B.Fuks, J.Keaveney, K.Mawatari, F.Blekman   <>" << endmsg;
    INFO << "        <>  Contact: lana.beck@cern.ch                                                <>" << endmsg;
    INFO << "        <>           fuks@cern.ch                                                     <>" << endmsg;
    INFO << "        <>           ddidar@mail.cern.ch                                              <>" << endmsg;
    INFO << "        <>  Based on MadAnalysis 5 v1.1.11                                            <>" << endmsg;
    INFO << "        <>  For more information, see                                                 <>" << endmsg;
    INFO << "        <>  http://madanalysis.irmp.ucl.ac.be/wiki/PhysicsAnalysisDatabase            <>" << endmsg;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    
    cout << "BEGIN Initialization" << endl;
    
    Manager()->AddRegionSelection("SR28");
    
    Manager()->AddCut("2 leptons");
    Manager()->AddCut("same sign leptons");
    Manager()->AddCut("veto low mass");
    Manager()->AddCut("veto Z");

    Manager()->AddCut("Njets>4", "SR28");
    Manager()->AddCut("NBjets>2", "SR28");
    Manager()->AddCut("MET>120", "SR28");
    Manager()->AddCut("HT>400", "SR28");

    Manager()->AddHisto("njet",15,0,15);
    
    dCounterSelectionEff = 0;

    
    cout << "END   Initialization" << endl;
    return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void CMS_SUS_13_013::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
    cout << "BEGIN Finalization" << endl;
    
    //////// Plotting efficiency curves /////////
    
    cout<<"Selection efficiency = "<<dCounterSelectionEff<<endl;
    
    TFile* myOutput = new TFile("MyOutput.root", "RECREATE");
    
    double fauxHT, fauxHTEff;
    for (int i=0; i<=800; i+=20)
    {
        fauxHT = (double)i;
        fauxHTEff = dFunctionHT(fauxHT, "SR28");
        vecHT.push_back(fauxHT);
        vecHTEff.push_back(fauxHTEff);
    }
    
    double fauxMET, fauxMETEff;
    for (int i=0; i<=200; i+=5)
    {
        fauxMET = (double)i;
        fauxMETEff = dFunctionMET(fauxMET, "SR28");
        vecMET.push_back(fauxMET);
        vecMETEff.push_back(fauxMETEff);
    }
    
    double fauxJetReco, fauxJetRecoEff;
    for (int i=0; i<=100; i+=2)
    {
        fauxJetReco = (double)i;
        fauxJetRecoEff = dFunctionJetReco(fauxJetReco, "SR28");
        vecJetReco.push_back(fauxJetReco);
        vecJetRecoEff.push_back(fauxJetRecoEff);
    }
    
    double fauxBTag, fauxBTagEff;
    for (int i=40; i<601; i+=1)
    {
        fauxBTag = (double)i;
        fauxBTagEff = dFunctionBTag(fauxBTag, "SR28");
        vecBTag.push_back(fauxBTag);
        vecBTagEff.push_back(fauxBTagEff);
    }
    
    double fauxMuon, fauxMuonEff;
    for (int i=10; i<120; i+=5)
    {
        fauxMuon = (double)i;
        fauxMuonEff = dFunctionLepton(fauxMuon, 13, "SR28");
        vecMuon.push_back(fauxMuon);
        vecMuonEff.push_back(fauxMuonEff);
    }
    double fauxElectron, fauxElectronEff;
    for (int i=10; i<120; i+=5)
    {
        fauxElectron = (double)i;
        fauxElectronEff = dFunctionLepton(fauxElectron, 11, "SR28");
        vecElectron.push_back(fauxElectron);
        vecElectronEff.push_back(fauxElectronEff);
    }

    //myOutput->cd();
    
    TCanvas *c1 = new TCanvas();
    c1->cd();
    EffBTag = new TGraph(vecBTag.size(), &(vecBTag[0]), &(vecBTagEff[0]));
    EffBTag->SetLineColor(kBlue);
    EffBTag->SetLineWidth(2);
    EffBTag->Draw("APL");
    EffBTag->GetYaxis()->SetRangeUser(0, 1);
    EffBTag->SetTitle("BTag");
    EffBTag->Draw("APL");
    EffBTag->GetXaxis()->SetTitle("gen Jet PT");
    EffBTag->GetYaxis()->SetTitle("Efficiency");
    EffBTag->Draw("APL");
    EffBTag->Write("BTag");
    c1->SaveAs("BTagEff.png");
    
    TCanvas *c2 = new TCanvas();
    c2->cd();
    EffHT = new TGraph(vecHT.size(), &(vecHT[0]), &(vecHTEff[0]));
    EffHT->SetLineColor(kRed);
    EffHT->SetLineWidth(2);
    EffHT->SetTitle("HT");
    EffHT->Draw("APL");
    EffHT->GetXaxis()->SetTitle("gen HT");
    EffHT->GetYaxis()->SetTitle("Efficiency");
    EffHT->Draw("APL");
    EffHT->Write("HT");
    c2->SaveAs("HTEff.png");
    
    TCanvas *c3 = new TCanvas();
    c3->cd();
    EffMET = new TGraph(vecMET.size(), &(vecMET[0]), &(vecMETEff[0]));
    EffMET->SetLineColor(kGreen);
    EffMET->SetLineWidth(2);
    EffMET->SetTitle("MET");
    EffMET->Draw("APL");
    EffMET->GetXaxis()->SetTitle("gen ETmiss");
    EffMET->GetYaxis()->SetTitle("Efficiency");
    EffMET->Draw("APL");
    
    EffMET->Write("MET");
    c3->SaveAs("METEff.png");
    
    TCanvas *c4 = new TCanvas();
    c4->cd();
    EffJetReco = new TGraph(vecJetReco.size(), &(vecJetReco[0]), &(vecJetRecoEff[0]));
    
    EffJetReco->SetLineColor(kBlue);
    EffJetReco->SetLineWidth(2);
    EffJetReco->SetTitle("Jet reconstruction");
    EffJetReco->Draw("APL");
    
    EffJetReco->GetXaxis()->SetTitle("gen Jet PT");
    EffJetReco->GetYaxis()->SetTitle("Efficiency");
    EffJetReco->Draw("APL");
    
    EffJetReco->Write("JetReco");
    c4->SaveAs("JetRecoEff.png");
    
    
    TCanvas *c5 = new TCanvas();
    c5->cd();
    EffMuon = new TGraph(vecMuon.size(), &(vecMuon[0]), &(vecMuonEff[0]));
    EffElectron = new TGraph(vecElectron.size(), &(vecElectron[0]), &(vecElectronEff[0]));
    EffElectron->SetLineColor(kBlue);
    EffElectron->SetLineWidth(2);
    EffMuon->GetYaxis()->SetRangeUser(0, 0.7);
    
    EffMuon->SetLineColor(kRed);
    EffMuon->SetLineWidth(2);
    EffMuon->SetTitle("Jet reconstruction");
    EffMuon->Draw("APL");
    EffMuon->GetXaxis()->SetTitle("lepton PT");
    EffMuon->GetYaxis()->SetTitle("Efficiency");
    EffMuon->Draw("APL");
    EffElectron->Draw("same");
    
    //  EffMuon->Write("Muon");
    
    c5->SaveAs("LeptonEff.png");
    
    
    delete EffMET;
    delete EffHT;
    delete EffBTag;
    delete EffJetReco;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    
    cout << "END   Finalization" << endl;
}


// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool CMS_SUS_13_013::Execute(SampleFormat& sample, const EventFormat& event)
{
  // ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ***************************************************************************

    if (event.mc() !=0)
     {
         dBTagEff = 1;
         dLeptonEff = 1;
         
         double myEventWeight;
         if(Configuration().IsNoEventWeight()) myEventWeight=1.;
         else if(event.mc()->weight()!=0.) myEventWeight = event.mc()->weight();
         else
         {
             WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
             return false;
         }
         Manager()->InitializeForNewEvent(myEventWeight);
         
         
         //cout << "---------------NEW EVENT-------------------" << endl;

         std::vector<const MCParticleFormat*> electrons, muons, positrons, antimuons, jets, btags, MCMET;
         std::vector<const MCParticleFormat*> leptons; //electrons and muons
         std::vector<const MCParticleFormat*> posileptons; //positrons and antimuons
         std::vector<const MCParticleFormat*> negaleptons; //electrons and muons
         
         PHYSICS->mcConfig().AddHadronicId(5);  //identifying bjets as hadronic
         PHYSICS->mcConfig().AddHadronicId(21);  //identifying jets as hadronic
         PHYSICS->mcConfig().AddInvisibleId(12); //identifying met as invisible

         for (unsigned int i=0;i<event.mc()->particles().size();i++)
         {
             const MCParticleFormat* part = &event.mc()->particles()[i];
             
             /*cout<<"ID  :  "<<part->pdgid();
             if (PHYSICS->Id->IsInitialState(part)) cout << " (Initial state) ";
             else if (PHYSICS->Id->IsFinalState(part)) cout << " (Final state) ";
             else cout << " (Intermediate state) ";
             cout<<""<<endl;*/
             
             //---------------------------------------------------------------------------------------------//
             //-------------------------------Add particle to vector collections----------------------------//
             //---------------------------------------------------------------------------------------------//

             
             if(part->statuscode() != 1) continue;
             if(part->pdgid() == 11) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) electrons.push_back(part);
             }
             else if(part->pdgid() == 13) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) muons.push_back(part);
             }
             else if(part->pdgid() == -11) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) positrons.push_back(part);
             }
             else if(part->pdgid() == -13) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) antimuons.push_back(part);
             }
             else if(std::abs(part->pdgid()) == 5) {
                 if(part->momentum().Pt()>40 && std::abs(part->momentum().Eta())<2.5) btags.push_back(part);
             }
             else if(std::abs(part->pdgid()) == 21) {
                 if(part->momentum().Pt()>40 && std::abs(part->momentum().Eta())<2.5) jets.push_back(part);
             }
             else if(std::abs(part->pdgid()) == 12) {
                 MCMET.push_back(part);
             }
             
             if(std::abs(part->pdgid()) == 11 || std::abs(part->pdgid()) == 13) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) leptons.push_back(part);
             }
             if(part->pdgid() == 11 || part->pdgid() == 13) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) negaleptons.push_back(part);
             }
             if(part->pdgid() == -11 || part->pdgid() == -13) {
                 if(part->momentum().Pt()>20 && std::abs(part->momentum().Eta())<2.5) posileptons.push_back(part);
             }

         }
         
         //---------------------------------------------------------------------------------------------//
         //-----------------------Apply baseline cuts 2 same sign leptons-------------------------------//
         //---------------------------------------------------------------------------------------------//


         if ( !Manager()->ApplyCut( (leptons.size() > 1),"2 leptons")) return true;

         bool SSlep = false;
         if(posileptons.size() > 1 || negaleptons.size() > 1){
             SSlep = true;
         }
         if ( !Manager()->ApplyCut( SSlep,"same sign leptons")) return true;
         
         //---------------------------------------------------------------------------------------------//
         //---------------------------------Making selection region cuts--------------------------------//
         //---------------------------------------------------------------------------------------------//
         
         if( !Manager()->ApplyCut((jets.size() > 4), "Njets>4"))  return true;
         if( !Manager()->ApplyCut((btags.size() > 2), "NBjets>2"))  return true;
         if( !Manager()->ApplyCut((event.mc()->MET().pt() >120    ), "MET>120"))  return true;
         if( !Manager()->ApplyCut((event.mc()->THT() >400   ), "HT>400"))  return true;
         
         if( leptons.size() > 3 ){
             if (leptons.size() > 4){
                 cout<<"4 leptons!!"<<endl;
             }
         }
         
         //---------------------------------------------------------------------------------------------//
         //----------------------------------------------veto-------------------------------------------//
         //--------------------low-mass bound-state or γ∗ → l+l− in the final state---------------------//
         //-------------------as well as multiboson (WZ, ZZ and tribosons) production-------------------//
         //---------------------------------------------------------------------------------------------//
         
         if( muons.size()>0 && antimuons.size()>0 ){  //for 2 opposite sign same flavour muons
             if ( electrons.size() > 0 ){
                 if (antimuons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and antimuon <12GeV |eta|<2.4) return true;
                     if (antimuons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
             if ( muons.size() > 1 ){
                 if (antimuons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of either muon and antimuon <12GeV |eta|<2.4) return true;
                     if (antimuons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
                 
             }
             if (positrons.size() > 0 ){
                 if (muons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and antimuon <12GeV |eta|<2.4) return true;
                     if (muons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
             if ( antimuons.size()>1 ){
                 if (muons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and either antimuon <12GeV |eta|<2.4) return true;
                     if (muons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
         }
         
         //---------------------------------------------------------------------------------------------//
         
         if( electrons.size()>0 && positrons.size()>0 ){  //for 2 opposite sign same flavour electrons
             if ( muons.size() > 0 ){
                 if (positrons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and antimuon <12GeV |eta|<2.4) return true;
                     if (positrons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
             if ( electrons.size() > 1 ){
                 if (positrons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of either muon and antimuon <12GeV |eta|<2.4) return true;
                     if (positrons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
             if (antimuons.size() > 0 ){
                 if (electrons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and antimuon <12GeV |eta|<2.4) return true;
                     if (electrons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
             if ( positrons.size()>1 ){
                 if (electrons[0]->pt()>10){
                     cout<<"VETO"<<endl;
                     //if (invariant mass of muon and either antimuon <12GeV |eta|<2.4) return true;
                     if (electrons[0]->pt()>5){
                         cout<<"VETO"<<endl;
                         //if (invariant mass of muon and antimuon 76<Mll<106GeV |eta|<2.4) return true;
                     }
                 }
             }
         }
         
         
         dMETEff = dFunctionMET(event.mc()->MET().pt(), "SR28");
         dHTEff = dFunctionHT(event.mc()->THT(), "SR28");
         
         
         for( int i=0; i<btags.size(); i++ )
         {
             dBTagEff *= dFunctionJetReco(btags[i]->momentum().Pt(), "SR28") * dFunctionBTag(btags[i]->momentum().Pt(), "SR28");
         }
         
         for ( int i=0; i<leptons.size(); i++){
             dLeptonEff *= dFunctionLepton(leptons[i]->momentum().Pt(), leptons[i]->pdgid(), "SR28");
         }

         dSelectionEff = dHTEff * dMETEff * dBTagEff * dLeptonEff;
         dCounterSelectionEff += dSelectionEff;
         cout<<"selection eff:  "<<dSelectionEff<<endl;
         return true;
     }//end of event.mc()
     
     return true;
}//end of event loop

double CMS_SUS_13_013::dFunctionMET(double dGenMET, std::string SearchReg){
    double E_inf=0;
    double numxhalf=0;
    double sigma=1;
    if (SearchReg == "SR28"){
        E_inf = 0.999;      // +\- 0.001
        numxhalf = 117.85;  // +\-0.09
        sigma = 36.90;      // +\-0.14
    }
    double efficiency  = 0.5 * E_inf * (erf((dGenMET - numxhalf)/sigma) + 1.);
    return efficiency;
}

double CMS_SUS_13_013::dFunctionHT(double dGenHT, std::string SearchReg){
    double E_inf=0;
    double numxhalf=0;
    double sigma=1;
    if (SearchReg == "SR28"){
        E_inf = 0.999;      // +\- 0.001
        numxhalf = 378.69;  // +\-0.17
        sigma = 59.41;      // +\-0.26
    }
    double efficiency  = 0.5 * E_inf * (erf((dGenHT - numxhalf)/sigma) + 1.);
    return efficiency;
}

double CMS_SUS_13_013::dFunctionLepton(double dGenLepPt, int leptFlav, std::string SearchReg){
    double E_inf=0;
    double E_10=0;
    double sigma=1;
    if (SearchReg == "SR28"){
        if (leptFlav == 11 || leptFlav == -11){ //ELECTRON
            E_inf = 0.640;      // +\- 0.001
            E_10 = 0.170;  // +\-0.002
            sigma = 36.94;      // +\-0.320
        }
        if (leptFlav == 13 || leptFlav == -13){ //MUON
            E_inf = 0.673;      // +\- 0.001
            E_10 = 0.332;  // +\-0.003
            sigma = 29.65;      // +\-0.382
        }
    }
    double partofeff = (dGenLepPt - 10)/sigma ;
    double efficiency  = ( E_inf * ( erf(partofeff) ) ) + (E_10*(1 - erf(partofeff)));
    return efficiency;
}

double CMS_SUS_13_013::dFunctionBTag(double dGenJetPt, std::string SearchReg){
    ///////Efficiency for a reconstructed jet that matchers to a b-quark to be btagged by the CMS algorithm/////
    double dA = 1.55e-06; // +/-0.05e-07
    double dB = -4.26e-04; // +/- 0.12e-04
    double dC = 0.0391; // +/- 0.0008
    double dD = -0.496; // +/- 0.020
    double dE = -3.26e-04; // +/- 0.01e-04
    double dF = 0.7681; // +/- 0.0016
    double effbTagged = 0;
    
    if (dGenJetPt > 120){
        effbTagged = (dE*dGenJetPt) + dF;
    }
    else{
        effbTagged = (dA * pow(dGenJetPt,3) )  + (dB * pow(dGenJetPt, 2)) + (dC * dGenJetPt) + dD;
    }
    
    return effbTagged;
}

double CMS_SUS_13_013::dFunctionJetReco(double dGenJetPt, std::string SearchReg){
    ///////Efficiency for a given GenJet to be reconstructed by CMS detector with PT>40///////
    double E_inf = 1.0;
    double numxhalf = 29.8;
    double sigma = 18.8;
    
    double effJetReco  = 0.5 * E_inf * ((erf((dGenJetPt - numxhalf)/sigma)) + 1.);
    
    return effJetReco;
}

bool CMS_SUS_13_013::dFunctionVETO(std::vector leptons1, std::vector leptons2){
    if (leptons1.size()>1 && leptons2.size() && leptons2[1]->pt()>5){//if there are two SS leptons and
                                                                        //a third OS lepton which has pt>5
        for (i=0; i<2, i++){                                            //for each of the SS lepton check the following
            if(leptons1[i]->pdgid() + leptons2[1]->pdgid() == 0)        //if the leptons have the same flavour
                TLorentzVector combinedM = leptons1[i].momentum + leptons2[1].momentum;
                if (combinedM.M()<12){return true;}                     //mass<12 => gamma veto
                if (leptons2[1].pt()>10 && 76<combinedM.M()<106){return true}              //76<mass<106 => Z veto
                else if {return false;}
        }
    }
    
}

