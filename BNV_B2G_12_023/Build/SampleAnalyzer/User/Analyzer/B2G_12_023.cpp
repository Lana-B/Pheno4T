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
    Manager()->AddRegionSelection("Basic");
    Manager()->AddRegionSelection("Tight");

    Manager()->AddCut("exactly 1 lepton");
    Manager()->AddCut("Njets>=5");
    Manager()->AddCut("bjets>=1");
    Manager()->AddCut("MET<20","Tight");
    Manager()->AddCut("Chi^2<20","Tight"); 

    Manager()->AddHisto("ET_miss", 40, 0, 200, "Basic");
    Manager()->AddHisto("Chi2", 20, 0, 200, "Basic");
    
    c1 = new TCanvas();
    c2 = new TCanvas();
    hET_miss = new TH1D("ET_miss", "ET_miss", 40, 0, 200);
    hChi2 = new TH1D("Chi2", "Chi2", 20, 0, 200);

    dCounterPassedEvents = 0;

    cout << "END   Initialization" << endl;


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
    cout<<"        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"<<endl; 
    cout<<"        <> Events passed in surviving selection region = "<<dCounterPassedEvents<<" <>"<<endl;
    cout<<"        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"<<endl; 
    double output = dFunctionBeta(0, "Muon");
    cout<<output<<endl;
    c1->cd();
    hET_miss->Draw();
    c1->SaveAs("../Output/ET_miss.png");
    c2->cd();
    hChi2->Draw();
    c2->SaveAs("../Output/Chi2.png");
    cout<<""<<endl;
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

        std::vector<const MCParticleFormat*> electrons, muons, positrons, antimuons, jets, bjets, taus, MCMET;
        std::vector<const MCParticleFormat*> leptons; //electrons and muons of either charge
        std::vector<const MCParticleFormat*> lights_taus; //lights and bs
         
        PHYSICS->mcConfig().AddHadronicId(5);   //identifying bjets as hadronic
        PHYSICS->mcConfig().AddHadronicId(21);  //identifying jets as hadronic
        PHYSICS->mcConfig().AddHadronicId(15);  //hadronically decaying taus
        PHYSICS->mcConfig().AddHadronicId(-15); //hadronically decaying anti-taus
        PHYSICS->mcConfig().AddInvisibleId(12); //identifying met as invisible

        int jet_40_counter = 0;
        int jet_55_counter = 0;
        int jet_70_counter = 0;

        for (unsigned int i=0;i<event.mc()->particles().size();i++)
        {
            const MCParticleFormat* part = &event.mc()->particles()[i];
         
            //---------------------------------------------------------------------------------------------//
            //-------------------------------Add particles to vector collections---------------------------//
            //---------------------------------------------------------------------------------------------//

            if(part->statuscode() != 1) continue; //ie. skip if not a final state particle

            if(part->pdgid() == 15 || part->pdgid() == -15 ){ 
                if(std::abs(part->momentum().Eta())<2.4) taus.push_back(part);
            }

            if(part->pdgid() == 11) {
                if(std::abs(part->momentum().Eta())<2.5 && !( 1.444 <std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566) && part->momentum().Pt()>30){
                    electrons.push_back(part);
                    leptons.push_back(part);
                }
            }
            else if(part->pdgid() == 13) {
                if(std::abs(part->momentum().Eta())<2.1 && part->momentum().Pt()>24){
                    muons.push_back(part);
                    leptons.push_back(part);
                }
             }
            else if(part->pdgid() == -11) {
                if(std::abs(part->momentum().Eta())<2.5 && !( 1.444 <std::abs(part->momentum().Eta()) && std::abs(part->momentum().Eta())< 1.566) && part->momentum().Pt()>30){
                    positrons.push_back(part);
                    leptons.push_back(part);
                }
            }
            else if(part->pdgid() == -13) {
                if(std::abs(part->momentum().Eta())<2.1 && part->momentum().Pt()>24){
                    antimuons.push_back(part);
                    leptons.push_back(part);
                }
            }
            else if(std::abs(part->pdgid()) == 5) {
                if(std::abs(part->momentum().Eta())<2.4 && part->momentum().Pt()>30 ) bjets.push_back(part);
            }
            else if(std::abs(part->pdgid()) == 12) {
                MCMET.push_back(part);
            }
         
            if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 15) { //light quarks and taus for W invariant mass
                if(std::abs(part->momentum().Eta())<2.4 && part->momentum().Pt()>30 ) {
                    jets.push_back(part);
                    lights_taus.push_back(part);
                }
            }
            if(std::abs(part->pdgid()) == 21 || std::abs(part->pdgid()) == 15 || std::abs(part->pdgid()) == 5) { //light quarks, taus and b quarks jet momentum requirements
                if(std::abs(part->momentum().Eta())<2.4 && part->momentum().Pt()>40) {
                    jet_40_counter++;
                    if(part->momentum().Pt()>55){
                        jet_55_counter++;
                        if( part->momentum().Pt()>70){
                            jet_70_counter++;
                        }
                    }
                }
            }
        }

        if( jet_70_counter < 1 || jet_55_counter < 2 || jet_40_counter < 3) return true; 
        //This ensures jets of momenta >70, >55, >40. The later cut of Njets>=5 with pt>30 takes care of the 4th and 5th jets momenta requirement.

        jets.insert( jets.end(), bjets.begin(), bjets.end());

        //---------------------------------------------------------------------------------------------//
        //-------------------------------- Apply Basic Selection cuts ---------------------------------//
        //---------------------------------------------------------------------------------------------//

        if ( !Manager()->ApplyCut( (leptons.size() == 1), "exactly 1 lepton")) return true; //there are at least 2 leptons with |eta|>2.5       
        if ( !Manager()->ApplyCut( ( jets.size() > 4), "Njets>=5"))  return true;
        if ( !Manager()->ApplyCut( (bjets.size() > 0), "bjets>=1"))  return true;


        //---------------------------------------------------------------------------------------------//
        //-------------------------------- Apply Tight Selection cuts ---------------------------------//
        //---------------------------------------------------------------------------------------------//
        double dChi2 = dJetCombiner(jets, lights_taus, leptons);
        // cout<<dChi2<<endl;

        Manager()->FillHisto("ET_miss", event.mc()->MET().pt());
        Manager()->FillHisto("Chi2", dChi2);
        hET_miss->Fill(event.mc()->MET().pt());
        hChi2->Fill(dChi2);

        if ( !Manager()->ApplyCut( ( event.mc()->MET().pt() < 20 ), "MET<20"))  return true;
        if ( !Manager()->ApplyCut( ( dChi2 < 20), "Chi^2<20"))  return true;

        dCounterPassedEvents += 1;             //counts number of Events

        return true;
    }//end of event.mc()
     
  return true;
}

double B2G_12_023::dJetCombiner(std::vector<const MCParticleFormat*> jets, std::vector<const MCParticleFormat*> lights_taus, std::vector<const MCParticleFormat*> leptons){
    double dMass_W = 82.4; //GeV
    double dSigma_W = 9.1; //GeV
    double dMass_HadTop = 171.94; //GeV
    double dSigma_HadTop = 14.8; //GeV
    double dMass_BNVTop = 174.8; //GeV
    double dSigma_BNVTop = 17.2; //GeV
    double dChi2_Total_Smallest = 1000; //Smallest chi^2 must be less than 20. This is a starting value for later logic.

    if ( lights_taus.size() < 2 ) return dChi2_Total_Smallest;

    //Loop for calculating invariant mass of W originating from hadronically decaying top. b-tagged jets not included.
    for(int j1 = 0; j1<lights_taus.size()-1; j1++){
        for(int j2 = j1+1; j2<lights_taus.size(); j2++){

            TLorentzVector LVector_W;
            LVector_W = lights_taus[j1]->momentum() + lights_taus[j2]->momentum();
            double dInvMass_W = LVector_W.M();
            double dChi2_W = pow((dInvMass_W - dMass_W),2) / pow(dSigma_W, 2);

            // Loop for hadronic top
            for(int j3 = 0; j3<jets.size(); j3++){
                if (j3 == j2 || j3 == j1) continue;

                TLorentzVector LVector_HadTop;
                LVector_HadTop = jets[j1]->momentum() + jets[j2]->momentum() + jets[j3]->momentum();
                double dInvMass_HadTop = LVector_HadTop.M();                
                double dChi2_HadTop =  pow((dInvMass_HadTop - dMass_HadTop),2) / pow(dSigma_HadTop, 2);

                // Loop for BNV decaying top
                for(int j4=0; j4<jets.size()-1; j4++){
                    if(j4 == j3 || j4 == j2 || j4 == j1) continue;

                    for(int j5=j4+1; j5<jets.size(); j5++){
                        if (j5 == j3 || j5 == j2 || j5 == j1) continue;

                        TLorentzVector LVector_BNVTop;
                        LVector_BNVTop = jets[j4]->momentum() + jets[j5]->momentum() + leptons[0]->momentum();
                        double dInvMass_BNVTop = LVector_BNVTop.M();                
                        double dChi2_BNVTop =  pow((dInvMass_BNVTop - dMass_BNVTop),2) / pow(dSigma_BNVTop, 2); 

                        double dChi2_Total = dChi2_W + dChi2_HadTop + dChi2_BNVTop;
                        if (dChi2_Total < dChi2_Total_Smallest){
                            dChi2_Total_Smallest = dChi2_Total;
                        }
                        // cout<<j1<<" : "<<j2<<" : "<<j3<<" : "<<j4<<" : "<<j5<<endl;   
                        // cout<<"Wmass: "<<dInvMass_W<<"  WChi2: "<<dChi2_W<<"  hadtop: "<<dInvMass_HadTop<<"  hadtopChi2: "<<dChi2_HadTop<<"  dInvMass_BNVTop: "<<dInvMass_BNVTop<<"  BNVChi2: "<<dChi2_BNVTop<<"  TotalChi2: "<<dChi2_Total<<endl;

                    }
                }
            }
        }
    }
    return dChi2_Total_Smallest;
}

double B2G_12_023::dFunctionBeta(double dBeta, std::string lepton){

    dSigma_tt = 22.2; // ±1.5pb
    dSigma_tW = 246; // ± 12pb
    dN_T_bck =8100;
    dN_B_bck = 290;
    dN_B_obs = 47950;
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
