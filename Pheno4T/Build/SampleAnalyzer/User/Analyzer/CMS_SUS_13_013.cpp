#include "SampleAnalyzer/User/Analyzer/CMS_SUS_13_013.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool CMS_SUS_13_013::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
    cout<<""<<endl;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    INFO << "        <>  Analysis: CMS-SUS-13-016, JHEP 01 (2014) 163                              <>" << endmsg;
    INFO << "        <>   (Search for new physics in events with same-sign                         <>" << endmsg;
    INFO << "        <>     dileptons and jets in pp collisions at √s = 8 TeV)                     <>" << endmsg;
    INFO << "        <>  Recasted by: L.Beck, D.Dobur, B.Fuks, J.Keaveney, K.Mawatari, F.Blekman   <>" << endmsg;
    INFO << "        <>  Contact: lana.beck@cern.ch                                                <>" << endmsg;
    INFO << "        <>           fuks@cern.ch                                                     <>" << endmsg;
    INFO << "        <>           ddidar@mail.cern.ch                                              <>" << endmsg;
    INFO << "        <>  Based on MadAnalysis 5 v1.1.11                                            <>" << endmsg;
    INFO << "        <>  For more information, see                                                 <>" << endmsg;
    INFO << "        <>  http://madanalysis.irmp.ucl.ac.be/wiki/PhysicsAnalysisDatabase            <>" << endmsg;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    cout<<""<<endl;
    cout << "BEGIN Initialization" << endl;
    
    Manager()->AddRegionSelection("SR28");
    
    Manager()->AddCut("2 leptons");
    Manager()->AddCut("same sign leptons");
    Manager()->AddCut("Njets>=4", "SR28");
    Manager()->AddCut("3rd lepton veto");
    
    dCounterSelectionEff = 0;
    dCounterPassedEvents = 0;

    METcounter = 0;
    HTcounter = 0;
    leptoncounter = 0;
    btagcounter = 0;
    
    cout << "END Initialization" << endl;
    cout<<""<<endl;

    return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void CMS_SUS_13_013::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
    cout<<""<<endl;
    cout << "BEGIN Finalization" << endl;
    cout<<""<<endl;
    cout<<"        <><><><>><><><><><><><><><><><><><><><><><><><>"<<endl; 
    cout<<"        <> Events passed = "<<dCounterPassedEvents<<"                       <>"<<endl;
    
    cout<<"        <> Weighted events = "<<dCounterSelectionEff<<"                 <>"<<endl;
    cout<<"        <> Selection efficiency for SR28= "<<dCounterSelectionEff*100/10000<<"%  <>"<<endl;
    cout<<"        <><><><>><><><><><><><><><><><><><><><><><><><>"<<endl; 

         //cout<<"average METeff = "<<METcounter/dCounterPassedEvents<"  HTeff = "<<HTcounter/dCounterPassedEvents<<"  leptoneff = "<<leptoncounter/dCounterPassedEvents<<"  btagcounter = "<<btagcounter/dCounterPassedEvents<<endl;

    cout<<""<<endl;
    
    //////// Plotting efficiency curves /////////
    
    for (int i=0; i<=800; i+=20){   //Filling HT curve vector
        vecHT.push_back((double)i);
        vecHTEff.push_back(dFunctionHT((double)i, "SR28"));
    }
    
    for (int i=0; i<=200; i+=5){   //Filling MET curve vector
        vecMET.push_back((double)i);
        vecMETEff.push_back(dFunctionMET((double)i, "SR28"));
    }
    
    for (int i=0; i<=100; i+=2){   //filling Jet Reco efficiency curve vector
        vecJetReco.push_back((double)i);
        vecJetRecoEff.push_back(dFunctionJetReco((double)i, "SR28"));
    }
    
    for (int i=0; i<601; i+=2){   //Filling bTag efficiency curve vector
        vecBTag.push_back((double)i);
        vecBTagEff.push_back(dFunctionBTag((double)i, "SR28"));
    }
    
    for (int i=10; i<120; i+=5){    //Filling muon efficiency curve vector (pdgid = 13)
        vecMuon.push_back((double)i);
        vecMuonEff.push_back(dFunctionLepton((double)i, 13, "SR28"));
    }

    for (int i=10; i<120; i+=5){    //Filling electron efficiency curve vector (pdgid = 11)
        vecElectron.push_back((double)i);
        vecElectronEff.push_back(dFunctionLepton((double)i, 11, "SR28"));
    }

    //Graph for the efficiency of tagging a b jet as a b
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
    c1->SaveAs("BTagEff.png");
    
    //Graph for HT efficiency
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
    c2->SaveAs("HTEff.png");
    
    //Graph for MET efficiency
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
    c3->SaveAs("METEff.png");
    
    //Graph for the efficiency of reconstructing a jet
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
    c4->SaveAs("JetRecoEff.png");
    
    //Graph for the efficiency of reconstructing leptons - muons and electrons
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

    cout<<""<<endl;
    cout << "END   Finalization" << endl;
    cout<<""<<endl;

}


// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool CMS_SUS_13_013::Execute(SampleFormat& sample, const EventFormat& event)
{
    if (event.mc() !=0)
     {
         double myEventWeight;
         if(Configuration().IsNoEventWeight()) myEventWeight=1.;
         else if(event.mc()->weight()!=0.) myEventWeight = event.mc()->weight();
         else
         {
             WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
             return false;
         }
         Manager()->InitializeForNewEvent(myEventWeight);
         if (myEventWeight != 1)
         {
           cout<<myEventWeight<<endl;
         }
        
         
         //cout << "---------------NEW EVENT-------------------" << endl;

         std::vector<const MCParticleFormat*> electrons, muons, positrons, antimuons, jets, bjets, MCMET;
         std::vector<const MCParticleFormat*> leptons; //electrons and muons of either charge
         std::vector<const MCParticleFormat*> posleptons; //positrons and antimuons
         std::vector<const MCParticleFormat*> negleptons; //electrons and muons
         std::vector<const MCParticleFormat*> lightsnbs; //lights and bs
         
         PHYSICS->mcConfig().AddHadronicId(5);   //identifying bjets as hadronic
         PHYSICS->mcConfig().AddHadronicId(21);  //identifying jets as hadronic
         PHYSICS->mcConfig().AddHadronicId(15);  //hadronically decaying taus
         PHYSICS->mcConfig().AddHadronicId(-15); //hadronically decaying anti-taus
         PHYSICS->mcConfig().AddInvisibleId(12); //identifying met as invisible

         for (unsigned int i=0;i<event.mc()->particles().size();i++)
         {
             const MCParticleFormat* part = &event.mc()->particles()[i];
             
             //---------------------------------------------------------------------------------------------//
             //-------------------------------Add particles to vector collections---------------------------//
             //---------------------------------------------------------------------------------------------//

             if(part->statuscode() != 1) continue; //ie. skip if not a final state particle
             //if(part->pdgid() == 15 || part->pdgid() == -15 ){ /*cout<<"!!!!!!!!  TAU !!!!!!!"<<endl;*/}
             if(part->pdgid() == 11) {
                if(std::abs(part->momentum().Eta())<2.5 && !( 1.4442<std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566  )){
                    electrons.push_back(part);
                    leptons.push_back(part);
                    negleptons.push_back(part);
                }
             }
             else if(part->pdgid() == 13) {
                if(std::abs(part->momentum().Eta())<2.5){
                    muons.push_back(part);
                    leptons.push_back(part);
                    negleptons.push_back(part);
                }
             }
             else if(part->pdgid() == -11) {
                if(std::abs(part->momentum().Eta())<2.5 && !( 1.4442<std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566  )){
                    positrons.push_back(part);
                    leptons.push_back(part);
                    posleptons.push_back(part);
                }
             }
             else if(part->pdgid() == -13) {
                if(std::abs(part->momentum().Eta())<2.5 ){
                    antimuons.push_back(part);
                    leptons.push_back(part);
                    posleptons.push_back(part);
                }
             }
             else if(std::abs(part->pdgid()) == 5) {
                if(std::abs(part->momentum().Eta())<2.5) bjets.push_back(part);
             }
             else if(std::abs(part->pdgid()) == 12) {
                MCMET.push_back(part);
             }
             
             if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 5) { //light quarks and b quarks for btagging weight
                 if(std::abs(part->momentum().Eta())<2.5) lightsnbs.push_back(part);
             }
             if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 5 || std::abs(part->pdgid()) == 15) { 
                //lights, bs, taus...ie. all jets
                if(std::abs(part->momentum().Eta())<2.5) jets.push_back(part);
             }

         }


         //---------------------------------------------------------------------------------------------//
         //-----------------------Apply baseline cuts 2 same sign leptons-------------------------------//
         //---------------------------------------------------------------------------------------------//

         if ( !Manager()->ApplyCut( (leptons.size() > 1),"2 leptons")) return true; //there are at least 2 leptons with |eta|>2.5       

         SSlep = false;  //is there a same sign pair?
         if(posleptons.size() > 1 || negleptons.size() > 1){
            SSlep = true;  //if there are either 2 positive or 2 negative leptons, then there exists a same sign pair
         }
         if ( !Manager()->ApplyCut( SSlep,"same sign leptons")) return true;
         
         //---------------------------------------------------------------------------------------------//
         //---------------------------------Making selection region cuts--------------------------------//
         //---------------------------------------------------------------------------------------------//

         if( !Manager()->ApplyCut((jets.size() > 3), "Njets>=4"))  return true;

         //---------------------------------------------------------------------------------------------//
         //----------------------------------------------veto-------------------------------------------//
         //--------------------low-mass bound-state or γ∗ → l+l− in the final state---------------------//
         //-------------------as well as multiboson (WZ, ZZ and tribosons) production-------------------//
         //---------------------------------------------------------------------------------------------//
 
         VETObool=true;  //default is that it is not vetoed

         bool negbool = dFunctionVETO(negleptons, posleptons); //if two SS negative leptons with pT>20, check positive leptons with no pT cut
         bool posbool = dFunctionVETO(posleptons, negleptons); //if two SS positive leptons with pT>20, check negative leptons with no pT cut

         if( negbool == true || posbool == true){
             VETObool = false; //if the veto criteria are satisified the event does not pass the selection
         }

         if ( !Manager()->ApplyCut( VETObool,"3rd lepton veto")) return true;
         //---------------------------------------------------------------------------------------------//
         //------------------------------Calculate and combine efficiences------------------------------//
         //---------------------------------------------------------------------------------------------//
         
         dMETEff = dFunctionMET(event.mc()->MET().pt(), "SR28");
         dHTEff = dFunctionHT(event.mc()->THT(), "SR28");
         dLeptonEff = dFunctionTotalLepton(posleptons, negleptons);
         dBTagEff = dFunctionTotalBTag(lightsnbs); //includes lights being mistagged as bs

         METcounter += dMETEff;
         HTcounter += dHTEff;
         leptoncounter += dLeptonEff;
         btagcounter += dBTagEff;


         dSelectionEff = dHTEff * dMETEff * dBTagEff * dLeptonEff;

         dCounterPassedEvents += 1;             //counts number of Events
         dCounterSelectionEff += dSelectionEff; //counts weighted number of Events

         return true;
     }//end of event.mc()
     
     return true;
}//end of event loop


//////////////////////////////////////////////////////////////////////////////
//////*********Custom functions for calculating the efficiences********///////
//////////////////////////////////////////////////////////////////////////////

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

double CMS_SUS_13_013::dFunctionLepton(double dGenLepPt, int leptFlav, std::string SearchReg){ //Function for individual lepton ID efficiency
    double E_inf=0;
    double E_10=0;
    double sigma=1;
    if (SearchReg == "SR28"){
        if (leptFlav == 11 || leptFlav == -11){ //ELECTRON
            E_inf = 0.597;      // +\- 0.001
            E_10 = 0.133;       // +\-0.002
            sigma = 37.75;      // +\-0.300
        }
        if (leptFlav == 13 || leptFlav == -13){ //MUON
            E_inf = 0.617;      // +\- 0.001
            E_10 = 0.291;       // +\-0.002
            sigma = 29.949;     // +\-0.377
        }
    }
    double partofeff = (dGenLepPt - 10)/sigma ;
    double efficiency  = ( E_inf * ( erf(partofeff) ) ) + (E_10*(1 - erf(partofeff)));

    return efficiency;
}

double CMS_SUS_13_013::dFunctionBTag(double dGenJetPt, std::string SearchReg){
    ///////Efficiency for a reconstructed jet that matches to a b-quark to be btagged by the CMS algorithm/////
    double dA = 1.55e-06;   // +/-0.05e-07
    double dB = -4.26e-04;  // +/- 0.12e-04
    double dC = 0.0391;     // +/- 0.0008
    double dD = -0.496;     // +/- 0.020
    double dE = -3.26e-04;  // +/- 0.01e-04
    double dF = 0.7681;     // +/- 0.0016
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
    double numxhalf = 34.9; //29.8;
    double sigma = 14.0;//18.8;
    
    double effJetReco  = 0.5 * E_inf * ((erf((dGenJetPt - numxhalf)/sigma)) + 1.);

    return effJetReco;
}

double CMS_SUS_13_013::dFunctionBTagCombined(double dGenJetPt, std::string SearchReg){
    ///////efficiency for a single jet to be reconstructed AND tagged as a b
    double Combined  = dFunctionJetReco(dGenJetPt, SearchReg) * dFunctionBTag(dGenJetPt, SearchReg);
    return Combined;
}

double CMS_SUS_13_013::dFunctionTotalBTag(std::vector<const MCParticleFormat*> lightsnbs){
    ///////efficiency for >= 2 reco level btags to be found given there are n gen level bjets
    double total = dFunctionGE2("b", lightsnbs);
    return total;
}

bool CMS_SUS_13_013::dFunctionVETO(std::vector<const MCParticleFormat*> leptons1, std::vector<const MCParticleFormat*> leptons2){
    bool bResult = false;
    if (leptons1.size()>1 && leptons2.size()>0){//if there are two SS leptons and
        for (int j=0; j<leptons2.size();j++){
            if(leptons2[j]->momentum().Pt()>5){//a third OS lepton which has pt>5
                for (int i=0; i<2; i++){       //for each of the SS lepton check the following
                    if(leptons1[i]->pdgid() + leptons2[j]->pdgid() == 0){        //if the leptons have the same flavour
                        TLorentzVector combinedM;
                        combinedM = leptons1[i]->momentum() + leptons2[j]->momentum();
                        if (combinedM.M()<12){//mass<12 => gamma veto
                            bResult = true;
                        }
                        else if (leptons2[j]->momentum().Pt()>10 && 76<combinedM.M() && combinedM.M()<106){//76<mass<106 => Z veto
                            bResult = true;
                        }
                    }

                }
            }
        }
    }
    return bResult;
}

double CMS_SUS_13_013::dFunctionTotalLepton(std::vector<const MCParticleFormat*> posleptons, std::vector<const MCParticleFormat*> negleptons){

    double w_0pairs;
    double w_GE1pair; //efficiency for >= one pair
    
    if ( posleptons.size() > 1 && negleptons.size() > 1){
        w_0pairs = (1 - dFunctionGE2("lepton", posleptons)) * (1 - dFunctionGE2("lepton", negleptons));
        w_GE1pair = 1 -  w_0pairs;
    }
    else if ( posleptons.size() > 1 && negleptons.size() <2){
        w_GE1pair = dFunctionGE2("lepton", posleptons);
    }
    else if ( negleptons.size() > 1 && posleptons.size() <2){
        w_GE1pair = dFunctionGE2("lepton", negleptons);
    }
    else{
        throw std::invalid_argument("there do not seem to be 2 same sign leptons");
    }
    return w_GE1pair;
}


double CMS_SUS_13_013::dFunctionGE2(std::string particleType, std::vector<const MCParticleFormat*> vec_particles){ //function for calculating getting >=2 objects at reco level given there are n at gen level
    double w_GE2;
    
    if (vec_particles.size() < 2){
        w_GE2 = 0;
    }
    
    else{  // if there are at least 2 particle in the vector we can calculate eff >=2
        double w_0bTags = 1;
        double w_1bTag = 0;
        double product, eff_i, eff_j;
        
        for (int i =0; i<vec_particles.size(); i++){
            if (particleType == "b"){
                if(vec_particles[i]->pdgid() == 5){
                    eff_i = dFunctionBTagCombined(vec_particles[i]->momentum().Pt(), "SR28");
                }
                else if(vec_particles[i]->pdgid() ==21){
                    eff_i = 0.01;
                }
            }
            else if ( particleType == "lepton" ){
                eff_i = dFunctionLepton( vec_particles[i]->momentum().Pt(), vec_particles[i]->pdgid(), "SR28");
            }
            else{
                throw std::invalid_argument("particle type not available");
            }
            if (eff_i < 0){
                eff_i = 0;
            }
            //cout<<" i: "<<i<<"  pdgid: "<<vec_particles[i]->pdgid()<<"   eff: "<<eff_i<<endl;
            w_0bTags *= (1 - eff_i);
            
            product = eff_i;
            for (int j=0; j< vec_particles.size(); j++){
                if (j == i){continue;}
                else{
                    if ( particleType == "b"){
                        if(vec_particles[j]->pdgid() == 5){
                            eff_j = dFunctionBTagCombined(vec_particles[j]->momentum().Pt(), "SR28");
                        }
                        else if (vec_particles[j]->pdgid() ==21){
                            eff_j = 0.01;
                        }
                    }
                    else if ( particleType == "lepton" ){
                        eff_j = dFunctionLepton( vec_particles[j]->momentum().Pt(), vec_particles[j]->pdgid(), "SR28");
                    }
                    else{
                        throw std::invalid_argument("particle type not available");
                    }
                    if (eff_j < 0){
                        eff_j = 0;
                    }
                    product *= ( 1 - eff_j );
                        //cout<<" j: "<<j<<"  pdgid: "<<vec_particles[j]->pdgid()<<"   eff: "<<eff_j<<endl;
                }
            }//end of j for loop
            w_1bTag += product;
        }//end of i for loop
        w_GE2 = 1 - w_0bTags - w_1bTag;
        //cout<<"w0: "<<w_0bTags<<"  w1: "<<w_1bTag<<"  w_GE2: "<<w_GE2<<endl;
    }
        
    return w_GE2;
}

