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
    INFO << "        <>  Analysis: CMS-SUS-13-016, JHEP 01 (2014) 163	 			 <>" << endmsg;
    INFO << "        <>   (Search for new physics in events with same-sign dileptons and jets in pp collisions at s√ = 8 TeV ) "<< endmsg;
    INFO << "        <>  Recasted by: L.Beck, D.Dobur, B.Fuks, J.Keaveney, K.Mawatari, F.Blekman   <>" << endmsg;
    INFO << "        <>  Contact: lana.beck@cern.ch             		 			 <>" << endmsg;
    INFO << "        <>           fuks@cern.ch                     	 			 <>" << endmsg;
    INFO << "        <>           ddidar@mail.cern.ch             				 <>" << endmsg;
    INFO << "        <>  Based on MadAnalysis 5 v1.1.11            	 			 <>" << endmsg;
    INFO << "        <>  For more information, see                 				 <>" << endmsg;
    INFO << "        <>  http://madanalysis.irmp.ucl.ac.be/wiki/PhysicsAnalysisDatabase		 <>" << endmsg;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    
  cout << "BEGIN Initialization" << endl;

    Manager()->AddRegionSelection("SR28");
    
    Manager()->AddCut("2 leptons");
    Manager()->AddCut("same sign leptons");
    Manager()->AddCut("lepton PT>20");
    Manager()->AddCut("lepton |eta|<2.5");
    
    Manager()->AddCut("MET>120", "SR28");
    Manager()->AddCut("HT>400", "SR28");
    Manager()->AddCut("Njets>4", "SR28");
    Manager()->AddCut("NBjets>2", "SR28");
    
    
    // how to add histo : Manager()->AddHisto("",10,0,300);
    //EffMET = plots.Add_TGraph();
    
    eventCounter = 0;
    
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

    // saving histos
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
        cout<<"dGenJetPt: "<<fauxJetReco<<"  Eff: "<<fauxJetRecoEff<<endl;
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
    //myOutput->cd();
    double fauxElectron, fauxElectronEff;
    for (int i=10; i<120; i+=5)
    {
        fauxElectron = (double)i;
        fauxElectronEff = dFunctionLepton(fauxElectron, 11, "SR28");
        vecElectron.push_back(fauxElectron);
        vecElectronEff.push_back(fauxElectronEff);
    }
    
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
    
    //plots.Write_TextFormat(out());
    //plots.Finalize();

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
        // Initial state
       /* for (unsigned int i=0;i<event.mc()->particles().size();i++)
         {
         const MCParticleFormat& part = event.mc()->particles()[i];
         TLorentzVector sumInvisible; //vectorial sum of invisible momentum
         if (PHYSICS->Id->IsInvisible(part))
         {
         sumInvisible += part->momentum();
         }
         }
         cout<<"sumInvisible: "<< sumInvisible <<endl;
        std::string sSR = "SR28";
        
        double dMET = event.mc()->MET().momentum().Pt();
        double dMETEff = dFunctionMET(dMET, sSR);
        //cout<<"MET: "<<dMET<<"  efficiency: "<<dMETEff<<endl;
        */
        double dHT = (double)event.mc()->TET();//.momentum().Pt();
        cout<<dHT<<endl;
        //double dHTEff = dFunctionHT(dHT, sSR);
        //cout<<"HT: "<<dHT<<"  efficiency: "<<dHTEff<<endl;
        //vecHT.push_back(dHT);
        //vecHTEff.push_back(dHTEff);
        //eventCounter += 1;
/*
        
        std::vector<const MCParticleFormat*> electrons;
        for(unsigned int i=0; i<event.mc()->particles().size(); i++)
        {
            const MCParticleFormat* prt =
            &event.mc()->particles()[i]; if(prt->statuscode() != 1) continue;
            if(std::abs(prt->pdgid()) == 11) {
                if(prt->momentum().Pt()>50) electrons.push_back(prt);
            }
        }
        
        vector<const RecLeptonFormat*> CurrentElectron;
*/
        
        //std::vector<const MCParticleFormat*> muons;

        
        
        return true;
        
    }
    

    
  /*
  if (event.mc()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Initial state
    for (unsigned int i=0;i<event.mc()->particles().size();i++)
    {
      const MCParticleFormat& part = event.mc()->particles()[i];

      cout << "----------------------------------" << endl;
      cout << "MC particle" << endl;
      cout << "----------------------------------" << endl;

      // display index particle
      cout << "index=" << i+1;

      // display the status code
      cout << "Status Code=" << part.statuscode() << endl;
      if (PHYSICS->Id->IsInitialState(part)) cout << " (Initial state) ";
      else if (PHYSICS->Id->IsFinalState(part)) cout << " (Intermediate state) ";
      else cout << " (Final state) ";
      cout << endl;

      // pdgid
      cout << "pdg id=" << part.pdgid() << endl;
      if (PHYSICS->Id->IsInvisible(part)) cout << " (invisible particle) ";
      else cout << " (visible particle) ";
      cout << endl;

      // display kinematics information
      cout << "px=" << part.px()
                << " py=" << part.py()
                << " pz=" << part.pz()
                << " e="  << part.e()
                << " m="  << part.m() << endl;
      cout << "pt=" << part.pt() 
                << " eta=" << part.eta() 
                << " phi=" << part.phi() << endl;

      // display particle mother id
      if (part.mother1()==0) 
      {
        cout << "particle with no mother." << endl;
      }
      else if (part.mother2()==0 || part.mother1()==part.mother2())
      {
        const MCParticleFormat* mother = part.mother1();
        cout << "particle coming from the decay of " 
             << mother->pdgid() << "." << endl;
      }
      else
      {
        const MCParticleFormat* mother1 = part.mother1();
        const MCParticleFormat* mother2 = part.mother2();
        cout << "particle coming from interaction between "
             << mother1->pdgid() << " and " << mother2->pdgid()
             << "." << endl;
      }
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.mc()->MET().pt()
         << " phi=" << event.mc()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.mc()->MHT().pt()
              << " phi=" << event.mc()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.mc()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.mc()->THT() << endl; 
   cout << endl;

  return true;
  }

*/

  // ***************************************************************************
  // Example of analysis with reconstructed objects
  // Concerned samples : 
  //   - LHCO samples
  //   - LHE/STDHEP/HEPMC samples after applying jet-clustering algorithm
  // ***************************************************************************

  if (event.rec()!=0)
  {
    cout << "---------------NEW EVENT-------------------" << endl;

    // Looking through the reconstructed electron collection
    for (unsigned int i=0;i<event.rec()->electrons().size();i++)
    {
      const RecLeptonFormat& elec = event.rec()->electrons()[i];
      cout << "----------------------------------" << endl;
      cout << "Electron" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << elec.charge() << endl;
      cout << "px=" << elec.px()
                << " py=" << elec.py()
                << " pz=" << elec.pz()
                << " e="  << elec.e()
                << " m="  << elec.m() << endl;
      cout << "pt=" << elec.pt() 
                << " eta=" << elec.eta() 
                << " phi=" << elec.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
                << elec.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed muon collection
    for (unsigned int i=0;i<event.rec()->muons().size();i++)
    {
      const RecLeptonFormat& mu = event.rec()->muons()[i];
      cout << "----------------------------------" << endl;
      cout << "Muon" << endl;
      cout << "----------------------------------" << endl;
      cout << "index=" << i+1 
                << " charge=" << mu.charge() << endl;
      cout << "px=" << mu.px()
                << " py=" << mu.py()
                << " pz=" << mu.pz()
                << " e="  << mu.e()
                << " m="  << mu.m() << endl;
      cout << "pt=" << mu.pt() 
                << " eta=" << mu.eta() 
                << " phi=" << mu.phi() << endl;
      cout << "ET/PT isolation criterion =" << mu.ET_PT_isol() << endl;
      cout << "pointer address to the matching MC particle: " 
           << mu.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed hadronic tau collection
    for (unsigned int i=0;i<event.rec()->taus().size();i++)
    {
      const RecTauFormat& tau = event.rec()->taus()[i];
      cout << "----------------------------------" << endl;
      cout << "Tau" << endl;
      cout << "----------------------------------" << endl;
      cout << "tau: index=" << i+1 
                << " charge=" << tau.charge() << endl;
      cout << "px=" << tau.px()
                << " py=" << tau.py()
                << " pz=" << tau.pz()
                << " e="  << tau.e()
                << " m="  << tau.m() << endl;
      cout << "pt=" << tau.pt() 
                << " eta=" << tau.eta() 
                << " phi=" << tau.phi() << endl;
      cout << "pointer address to the matching MC particle: " 
           << tau.mc() << endl;
      cout << endl;
    }

    // Looking through the reconstructed jet collection
    for (unsigned int i=0;i<event.rec()->jets().size();i++)
    {
      const RecJetFormat& jet = event.rec()->jets()[i];
      cout << "----------------------------------" << endl;
      cout << "Jet" << endl;
      cout << "----------------------------------" << endl;
      cout << "jet: index=" << i+1 
           << " charge=" << jet.charge() << endl;
      cout << "px=" << jet.px()
           << " py=" << jet.py()
           << " pz=" << jet.pz()
           << " e="  << jet.e()
           << " m="  << jet.m() << endl;
      cout << "pt=" << jet.pt() 
           << " eta=" << jet.eta() 
           << " phi=" << jet.phi() << endl;
      cout << "b-tag=" << jet.btag()
           << " true b-tag (before eventual efficiency)=" 
           << jet.true_btag() << endl;
      cout << "EE/HE=" << jet.EEoverHE()
           << " ntracks=" << jet.ntracks() << endl;
      cout << endl;
    }

    // Transverse missing energy (MET)
    cout << "MET pt=" << event.rec()->MET().pt()
         << " phi=" << event.rec()->MET().phi() << endl;
    cout << endl;

    // Transverse missing hadronic energy (MHT)
    cout << "MHT pt=" << event.rec()->MHT().pt()
              << " phi=" << event.rec()->MHT().phi() << endl;
    cout << endl;

    // Total transverse energy (TET)
    cout << "TET=" << event.rec()->TET() << endl;
    cout << endl;

    // Total transverse hadronic energy (THT)
    cout << "THT=" << event.rec()->THT() << endl;
    cout << endl;
  }
  
  return true;
}

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

