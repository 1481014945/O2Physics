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

#ifndef PWGUD_CORE_DGSELECTOR_H_
#define PWGUD_CORE_DGSELECTOR_H_

#include <vector>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "Framework/Logger.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/DGCutparHolder.h"

// -----------------------------------------------------------------------------
// add here Selectors for different types of diffractive events
// Selector for Double Gap events
class DGSelector
{
 public:
  // constructor/destructor
  DGSelector() { fPDG = TDatabasePDG::Instance(); }
  ~DGSelector() { delete fPDG; }

  template <typename CC, typename BCs, typename TCs, typename FWs>
  int Print(DGCutparHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(info, "Size of array %i", collision.size());
    return 1;
  }

  // Function to check if collision passes DG filter
  template <typename CC, typename BCs, typename TCs, typename FWs>
  int IsSelected(DGCutparHolder diffCuts, CC& collision, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    LOGF(debug, "Collision %f", collision.collisionTime());
    LOGF(debug, "Number of close BCs: %i", bcRange.size());

    // return if FIT veto is found in any of the compatible BCs
    // Double Gap (DG) condition
    // 4 types of vetoes:
    //  0 TVX
    //  1 TSC
    //  2 TCE
    //  3 TOR
    for (auto const& bc : bcRange) {
      /* for debuging
      auto isVetoed = udhelpers::FITveto(bc, diffCuts);
      auto isClean = udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());
      LOGF(info, "<IsSelected> isVetoed: %d isClean: %d", isVetoed, isClean);
      if (isVetoed) {
        return 1;
      }
      */

      if (udhelpers::FITveto(bc, diffCuts)) {
        return 1;
      }
    }

    // forward tracks
    LOGF(debug, "FwdTracks %i", fwdtracks.size());
    if (!diffCuts.withFwdTracks()) {
      for (auto& fwdtrack : fwdtracks) {
        LOGF(debug, "  %i / %f / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p(), fwdtrack.trackTimeRes());
        // only consider tracks with MID (good timing)
        if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) {
          return 2;
        }
      }
    }

    // no global tracks which are not vtx tracks
    // no vtx tracks which are not global tracks
    // no PV tracks with ITS only
    auto rgtrwTOF = 0.; // fraction of PV tracks with TOF hit
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        return 3;
      }
      if (diffCuts.globalTracksOnly() && track.isPVContributor() && !track.isGlobalTrack()) {
        return 4;
      }
      if (!diffCuts.ITSOnlyTracks() && track.isPVContributor() && !track.hasTPC()) {
        return 5;
      }

      // update fraction of PV tracks with TOF hit
      if (track.isPVContributor() && track.hasTOF()) {
        rgtrwTOF += 1.;
      }
    }
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    if (rgtrwTOF < diffCuts.minRgtrwTOF()) {
      return 6;
    }

    // number of vertex tracks
    if (collision.numContrib() < diffCuts.minNTracks() || collision.numContrib() > diffCuts.maxNTracks()) {
      return 7;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // consider only vertex tracks

    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      if (track.isPVContributor()) {

        // PID
        // if (!udhelpers::hasGoodPID(diffCuts, track)) {
        //   return 8;
        // }

        // pt
        lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
        if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
          return 9;
        }

        // eta
        if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
          return 10;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }

    // net charge
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      return 11;
    }
    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 12;
    }

    // good collision selection according to event-selection task
    /* if (!goodCollision(collision, diffCuts)){
      return 13;
    }*/
    // good collision selection per-partes
    std::vector<int> sels = diffCuts.collisionSel();
    // timeframe border
    if (sels[0] && !udhelpers::cutNoTimeFrameBorder(collision)) {
      return 13;
    }
    // same bunch pileup
    if (sels[1] && !udhelpers::cutNoSameBunchPileup(collision)) {
      return 14;
    }
    // ITS ROF border
    if (sels[2] && !udhelpers::cutNoITSROFrameBorder(collision)) {
      return 15;
    }
    // z-vtx from PV and FT0 agrees
    if (sels[3] && !udhelpers::cutIsGoodZvtxFT0vsPV(collision)) {
      return 16;
    }
    // is ITS-TPC track
    if (sels[4] && !udhelpers::cutIsVertexITSTPC(collision)) {
      return 17;
    }

    // if we arrive here then the event is good!
    return 0;
  };

  // Function to check if BC passes DG filter (without associated collision)
  template <typename BCs, typename TCs, typename FWs>
  int IsSelected(DGCutparHolder diffCuts, BCs& bcRange, TCs& tracks, FWs& fwdtracks)
  {
    // return if FIT veto is found in any of the compatible BCs
    // Double Gap (DG) condition
    // 4 types of vetoes:
    //  0 TVX
    //  1 TSC
    //  2 TCE
    //  3 TOR
    for (auto const& bc : bcRange) {
      if (udhelpers::FITveto(bc, diffCuts)) {
        return 1;
      }
    }

    // no activity in muon arm
    if (!diffCuts.withFwdTracks()) {
      for (auto& fwdtrack : fwdtracks) {
        LOGF(info, "  %i / %f / %f / %f / %f", fwdtrack.trackType(), fwdtrack.eta(), fwdtrack.pt(), fwdtrack.p(), fwdtrack.trackTimeRes());
        // only consider tracks with MID (good timing)
        if (fwdtrack.trackType() == 0 || fwdtrack.trackType() == 3) {
          return 2;
        }
      }
    }

    // number of tracks
    if (static_cast<int>(tracks.size()) < diffCuts.minNTracks() || static_cast<int>(tracks.size()) > diffCuts.maxNTracks()) {
      return 7;
    }

    // PID, pt, and eta of tracks, invariant mass, and net charge
    // which particle hypothesis?
    auto mass2Use = 0.;
    TParticlePDG* pdgparticle = fPDG->GetParticle(diffCuts.pidHypothesis());
    if (pdgparticle != nullptr) {
      mass2Use = pdgparticle->Mass();
    }

    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    for (auto& track : tracks) {
      // PID
      if (!udhelpers::hasGoodPID(diffCuts, track)) {
        return 8;
      }

      // pt
      lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
      if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
        return 9;
      }

      // eta
      if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
        return 10;
      }
      netCharge += track.sign();
      ivm += lvtmp;
    }

    // net charge
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      return 11;
    }

    // invariant mass
    if (ivm.M() < diffCuts.minIVM() || ivm.M() > diffCuts.maxIVM()) {
      return 12;
    }

    // if we arrive here then the event is good!
    return 0;
  };

 private:
  TDatabasePDG* fPDG;

  ClassDefNV(DGSelector, 1);
};

#endif // PWGUD_CORE_DGSELECTOR_H_
