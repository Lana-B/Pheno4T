#include "SampleAnalyzer/User/Analyzer/B2G_12_023.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool B2G_12_023::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{

    cout<<""<<endl;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    INFO << "        <>  Analysis: CMS-B2G-12-023, Physics Letters B 731 (2014) 173      <>" << endmsg;
    INFO << "        <>   (Search for baryon number violating top-quark decays           <>" << endmsg;
    INFO << "        <>          in pp collisions at sqrt(s) = 8 TeV)                    <>" << endmsg;
    INFO << "        <>  Recasted by: L.Beck                                             <>" << endmsg;
    INFO << "        <>  Contact: lana.beck@cern.ch                                      <>" << endmsg;
    INFO << "        <>  Based on MadAnalysis 5 v1.1.11                                  <>" << endmsg;
    INFO << "        <>  For more information, see                                       <>" << endmsg;
    INFO << "        <>  http://madanalysis.irmp.ucl.ac.be/wiki/PhysicsAnalysisDatabase  <>" << endmsg;
    INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;
    cout<<""<<endl;

    cout << "BEGIN Initialization" << endl;
    // initialize variables, histos
    cout << "END   Initialization" << endl;
    Manager()->AddRegionSelection("BNVdecays");

    Manager()->AddCut("exactly 1 lepton");
    Manager()->AddCut("Njets>=5");
    Manager()->AddCut("bjets>=1");

    dCounterPassedEvents = 0;

    return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void B2G_12_023::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
    cout << "BEGIN Finalization" << endl;

    cout<<""<<endl;
    cout<<"        <><><><><><><><><><><><><><><><><><><><><><><><>"<<endl; 
    cout<<"        <> Events passed = "<<dCounterPassedEvents<<"                        <>"<<endl;
        //cout<<"        <> Weighted events = "<<dCounterSelectionEff<<"                  <>"<<endl;
    //cout<<"        <> Selection efficiency for SR28= "<<dCounterSelectionEff*100/10000<<"%  <>"<<endl;
    cout<<"        <><><><><><><><><><><><><><><><><><><><><><><><>"<<endl; 

    cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool B2G_12_023::Execute(SampleFormat& sample, const EventFormat& event)
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
                if(std::abs(part->momentum().Eta())<2.4 && !( 1.4442<std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566) && part->momentum().Pt()>20){
                    electrons.push_back(part);
                    leptons.push_back(part);
                    negleptons.push_back(part);
                }
            }
            else if(part->pdgid() == 13) {
                if(std::abs(part->momentum().Eta())<2.4 && part->momentum().Pt()>20){
                    muons.push_back(part);
                    leptons.push_back(part);
                    negleptons.push_back(part);
                }
             }
            else if(part->pdgid() == -11) {
                if(std::abs(part->momentum().Eta())<2.4 && !( 1.4442<std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566) && part->momentum().Pt()>20){
                    positrons.push_back(part);
                    leptons.push_back(part);
                    posleptons.push_back(part);
                }
            }
            else if(part->pdgid() == -13) {
                if(std::abs(part->momentum().Eta())<2.4 && part->momentum().Pt()>20){
                    antimuons.push_back(part);
                    leptons.push_back(part);
                    posleptons.push_back(part);
                }
            }
            else if(std::abs(part->pdgid()) == 5) {
                if(std::abs(part->momentum().Eta())<2.4) bjets.push_back(part);
            }
            else if(std::abs(part->pdgid()) == 12) {
                MCMET.push_back(part);
            }
         
            if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 5) { //light quarks and b quarks for btagging weight
                if(std::abs(part->momentum().Eta())<2.4) lightsnbs.push_back(part);
            }
            if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 5 || std::abs(part->pdgid()) == 15) { 
                //lights, bs, taus...ie. all jets
                if(std::abs(part->momentum().Eta())<2.4) jets.push_back(part);
            }
        }


        //---------------------------------------------------------------------------------------------//
        //-------------------------------- Apply baseline cuts ----------------------------------------//
        //---------------------------------------------------------------------------------------------//

        if ( !Manager()->ApplyCut( (leptons.size() == 1), "exactly 1 lepton")) return true; //there are at least 2 leptons with |eta|>2.5       
        if ( !Manager()->ApplyCut((jets.size() > 4), "Njets>=5"))  return true;
        if ( !Manager()->ApplyCut((bjets.size() > 0), "bjets>=1"))  return true;

        //---------------------------------------------------------------------------------------------//
        //------------------------------Calculate and combine efficiences------------------------------//
        //---------------------------------------------------------------------------------------------//

        // dMETEff = dFunctionMET(event.mc()->MET().pt(), "SR28");
        // dHTEff = dFunctionHT(event.mc()->THT(), "SR28");
        // dLeptonEff = dFunctionTotalLepton(posleptons, negleptons);
        // dBTagEff = dFunctionTotalBTag(lightsnbs); //includes lights being mistagged as bs

        // METcounter += dMETEff;
        // HTcounter += dHTEff;
        // leptoncounter += dLeptonEff;
        // btagcounter += dBTagEff;


        //dSelectionEff = dHTEff * dMETEff * dBTagEff * dLeptonEff;
        //Manager()->SetCurrentEventWeight((float)dSelectionEff);



        dCounterPassedEvents += 1;             //counts number of Events
        //dCounterSelectionEff += dSelectionEff; //counts weighted number of Events

        return true;
    }//end of event.mc()
     









  // ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ***************************************************************************
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
      if (PHYSICS->IsInitialState(part)) cout << " (Initial state) ";
      else if (PHYSICS->IsFinalState(part)) cout << " (Intermediate state) ";
      else cout << " (Final state) ";
      cout << endl;

      // pdgid
      cout << "pdg id=" << part.pdgid() << endl;
      if (PHYSICS->IsInvisible(part)) cout << " (invisible particle) ";
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
  /*
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
  */
  return true;
}

void B2G_12_023::dFunctionBeta(double dBeta, string lepton){


    double dFactorA = 1 + (( dSigma_tW * dE_B_tW )/( dSigma_tt * dE_B_tt ));
    double dFactorB = 1 + (( dSigma_tt * dE_B_tt )/( dSigma_tW * dE_B_tW ));

    double dFactor1 = ( 1/dFactorA )  *  ( dE_T_tt/dE_B_tt );
    double dFactor2 = ( 1/dFactorB )  *  ( dE_T_tW/dE_B_tW );

    double dN_T_exp = (  ( dN_B_obs - dN_B_bck ) * (Factor1 + Factor2)  ) + dN_T_bck;

    return dN_T_exp;

}

void B2G_12_023::dFunctionEfficiencies(double dBeta, string lepton, string channel, string TorB){

    if(lepton=)



}

