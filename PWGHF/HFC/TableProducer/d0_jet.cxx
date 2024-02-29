// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \author Xu wang <wangxuwx@mails.ccnu.edu.cn>.
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "d0minijet.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::constants::physics;
using namespace o2::framework;
//  using namespace o2::do_mini_jet;

// histogram binning definition
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = o2::constants::math::TwoPI;
const int massAxisBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;

double getDeltaPhi(double phiHadron, double phiD)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

struct A_D0minijet
{
    SliceCache cache;
    Produces<aod::DHadronPair> entryD0HadronPair;
    //Produces<aod::DHadronRecoInfo> entryD0HadronRecoInfo;
   // Produces<aod::dominijet> rowdominijet;
    std::vector<int> vtriggerFlag;
    Configurable<double> cutforwardjet{"cutforwardjet", 0.9, "the range of forwardjet"};
    Configurable<double> cutbackwardjet{"cutbackwardjet", 1.4, "the range of backwardjet"};
    Configurable<double> cuttriggerpt{"cuttriggerpt", 1., "the range of trigger particle pt"};
    Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
    Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};
    Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
    Configurable<double> ptCandMin{"ptCandMin", -1., "min. cand. pT"};
    Configurable<int> triggerFlag{"triggerFlag", 1, "trigger Flag"};
    Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
    Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
    Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};

    HfHelper hfHelper;

    HistogramRegistry registry
    {
        "registry",
        {
          {"hPtCand", "D0,D0bar candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
          {"hPtProng0", "D0,D0bar candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
          {"hPtProng1", "D0,D0bar candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
          {"hEta", "D0,D0bar candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
          {"hPhi", "D0,D0bar candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
          {"hY", "D0,D0bar candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
        }
    };
    Partition<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
    void init(InitContext const&)
    {
       auto vbins = (std::vector<double>)bins;
       //registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
       registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
       registry.add("hMassD0", "D0,D0bar candidates;inv. mass D0 only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
       registry.add("hMassD0bar", "D0,D0bar candidates;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    }
    void processData(aod::Collision const& collision,aod::TracksWDca const& tracks,soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&)
    {
      int nTracks = 0;
      int ntriggers = 0;
      double setparticlestatus = 0;
      if (collision.numContrib() > 1) {
        for (const auto& track : tracks) {
          if (track.eta() < -4.0 || track.eta() > 4.0) {
            vtriggerFlag.push_back(0);
            continue;
          }
          if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) {
            vtriggerFlag.push_back(0);
            continue;
          }
          if( track.pt() > cuttriggerpt)
          {
            ntriggers++;
          }
          nTracks++;
        }
      }
      registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
      if (nTracks < multMin || nTracks > multMax) {
        return;
      }
      registry.fill(HIST("hMultiplicity"), nTracks);

    auto selectedD0CandidatesGrouped = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache); 

    for (const auto& candidate1 : selectedD0CandidatesGrouped)
    {
      if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate1)) > yCandMax) {
        continue;
      }
      if (ptCandMin >= 0. && candidate1.pt() < ptCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }    
       // ========================== Fill mass histo  ================================
      if (candidate1.isSelD0() >= selectionFlagD0) 
      {
         registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt());
         registry.fill(HIST("hMassD0"), hfHelper.invMassD0ToPiK(candidate1), candidate1.pt());
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) 
      {
          registry.fill(HIST("hMass"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt());
          registry.fill(HIST("hMassD0bar"), hfHelper.invMassD0barToKPi(candidate1), candidate1.pt());
      }
      // ========================== Fill general histos ================================
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), hfHelper.yD0(candidate1));
      for (const auto& track : tracks) 
      {
        registry.fill(HIST("hTrackCounter"), 1); // fill total no. of tracks
        // Remove D0 daughters by checking track indices
        if ((candidate1.prong0Id() == track.globalIndex()) || (candidate1.prong1Id() == track.globalIndex())) {
          continue;
        }
        if (std::abs(track.dcaXY()) >= 1. || std::abs(track.dcaZ()) >= 1.)
          {continue;} // Remove secondary tracks

        registry.fill(HIST("hTrackCounter"), 2); // fill no. of tracks before soft pion removal

        if(track.pt() < cuttriggerpt) {continue;}
        double deltaPhi = getDeltaPhi(track.phi(), candidate1.phi());
        //double deltaEta = track.eta() - candidate1.eta();
        setparticlestatus = 0;
        if( std::abs(deltaPhi) <= cutforwardjet)
        {
          setparticlestatus = 1;
        }
        if( std::abs(deltaPhi - o2::constants::math::PI) <= cutbackwardjet)
        {
          setparticlestatus = 2;
        }

        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),track.eta() - candidate1.eta(),candidate1.pt(),track.pt(),setparticlestatus);
      }
    }
  }
  PROCESS_SWITCH(A_D0minijet, processData, "d0 mini jet", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<A_D0minijet>(cfgc)};
  return workflow;
}

