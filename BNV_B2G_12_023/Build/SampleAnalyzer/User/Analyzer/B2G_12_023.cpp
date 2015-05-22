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
    Manager()->AddRegionSelection("Basic");
    Manager()->AddRegionSelection("Tight");

    Manager()->AddCut("exactly 1 lepton");
    Manager()->AddCut("Njets>=5");
    Manager()->AddCut("bjets>=1");
    //Manager()->AddCut("MET<20","Tight");
    //Manager()->AddCut("X^2<20","Tight");    

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
    //double output = dFunctionBeta(0, "Muon");
    //cout<<output<<endl;
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
        if ( !Manager()->ApplyCut( ( jets.size() > 4), "Njets>=5"))  return true;
        if ( !Manager()->ApplyCut( (bjets.size() > 0), "bjets>=1"))  return true;

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
     
  return true;
}

double B2G_12_023::dFunctionBeta(double dBeta, std::string lepton){

    dE_B_tt = dFunctionEfficiencies(dBeta, lepton, "Basic", "tt");
    dE_T_tt = dFunctionEfficiencies(dBeta, lepton, "Tight", "tt");
    dE_B_tW = dFunctionEfficiencies(dBeta, lepton, "Basic", "tW");
    dE_T_tW = dFunctionEfficiencies(dBeta, lepton, "Tight", "tW");

    double dFactorA = 1 + (( dSigma_tW * dE_B_tW )/( dSigma_tt * dE_B_tt ));
    double dFactorB = 1 + (( dSigma_tt * dE_B_tt )/( dSigma_tW * dE_B_tW ));

    double dFactor1 = ( 1/dFactorA )  *  ( dE_T_tt/dE_B_tt );
    double dFactor2 = ( 1/dFactorB )  *  ( dE_T_tW/dE_B_tW );

    double dN_T_exp = (  ( dN_B_obs - dN_B_bck ) * (dFactor1 + dFactor2)  ) + dN_T_bck;

    return dN_T_exp;

}

double B2G_12_023::dFunctionEfficiencies(double dBeta, std::string lepton, std::string TorB, std::string channel){

    if(lepton=="Muon"){
        if (TorB == "Basic"){
            dEpsilon_X_SM_SM   = 8.1E-3;   // ±1.5E3
            dEpsilon_X_BNV_SM  = 7.37E-02; // ±0.89E-2            
            dEpsilon_X_BNV_BNV = 1E-02;    // ±0.16E-2
            dEpsilon_X_SM      = 2.68E-03; // ±0.32E-3
            dEpsilon_X_BNV     = 2.72E-02; // ±0.42E-2
        }
        else if (TorB == "Tight"){
            dEpsilon_X_SM_SM   = 4.62E-04; // ±0.93E-2
            dEpsilon_X_BNV_SM  = 1.86E-02; // ±0.32E-2           
            dEpsilon_X_BNV_BNV = 1.74E-03; // ±0.32E-3
            dEpsilon_X_SM      = 1.13E-04; // ±0.14E-4
            dEpsilon_X_BNV     = 5.38E-03; // ± 0.84E-03    
        }
    }

    else if(lepton=="Electron"){
        if (TorB == "Basic"){
            dEpsilon_X_SM_SM   = 8.0E-03;  // ±1.3E-03
            dEpsilon_X_BNV_SM  = 7.33E-02; // ±0.88E-02 
            dEpsilon_X_BNV_BNV = 1.55E-02; // ±0.25E-02   
            dEpsilon_X_SM      = 2.64E-03; // ±0.31E-03
            dEpsilon_X_BNV     = 2.80E-02; // ±0.42E-02           
        }
        else if (TorB == "Tight"){
            dEpsilon_X_SM_SM   = 4.24E-04; // ±0.85E-04
            dEpsilon_X_BNV_SM  = 1.62E-02; // ±0.27E-02           
            dEpsilon_X_BNV_BNV = 2.64E-03; // ±0.55E-03
            dEpsilon_X_SM      = 8.21E-05; // ±0.99E-05 
            dEpsilon_X_BNV     = 5.84E-03; // ±0.82E-03            
        }
    }
    else throw std::invalid_argument(" Wrong lepton name selected ");

    if(channel == "tt"){
        dEpsilon_X_chan = 2*dBeta*(1-dBeta)*dEpsilon_X_BNV_SM + pow((1-dBeta),2)*dEpsilon_X_SM_SM + pow(dBeta,2)*dEpsilon_X_BNV_BNV;
    }
    else if (channel == "tW"){
        dEpsilon_X_chan = (1-dBeta)*dEpsilon_X_SM + dBeta*dEpsilon_X_BNV;
    }
    else throw std::invalid_argument("*** Wrong channel selected ***");

    return dEpsilon_X_chan;
}

