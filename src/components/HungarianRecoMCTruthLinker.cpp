/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "k4FWCore/Transformer.h"

#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "IoUTruthMatcher.h"

#include <string>
#include <vector>

using RecoTruthLinkCollection = edm4hep::RecoMCParticleLinkCollection;
using CaloTruthLinks = std::vector<const edm4hep::CaloHitMCParticleLinkCollection*>;
using CaloSimRelations = std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>;
using TrackTruthLinks = std::vector<const edm4hep::TrackMCParticleLinkCollection*>;
using TrackerHitSimRelations = std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>;

/** Links each ReconstructedParticle (PandoraPFOs, HitPF, ...) to at most one MCParticle.
 *
 *  Gaudi counterpart of MarlinReco's RecoMCTruthLinker, using the same parameter names for the collections it
 *  shares, but with a different matching: the IoU between the particle's constituents (calorimeter hits of all its
 *  clusters + its tracks) and each MCParticle's, resolved one-to-one by maximising the summed IoU. Weight = IoU.
 *
 *  Hit/track truth can come from either of two interchangeable sources (configure exactly one of each pair):
 *   - calorimeter: CalohitMCTruthLinkName (CaloHit -> MCParticle links)
 *                  or SimCalorimeterHitRelationNames (CaloHit -> SimCaloHit relations, MCParticle contributions)
 *   - tracks:      TrackMCTruthLinkName (Track -> MCParticle links)
 *                  or TrackerHitsRelInputCollections (TrackerHit -> SimTrackerHit relations, majority of hits)
 *  Calorimeter secondaries are attributed to their ancestor produced outside the calorimeter. Only MCParticles
 *  contained in MCParticleCollection are considered as truth.
 */
struct HungarianRecoMCTruthLinker final
    : k4FWCore::Transformer<RecoTruthLinkCollection(
          const edm4hep::MCParticleCollection&, const edm4hep::ReconstructedParticleCollection&, const CaloTruthLinks&,
          const CaloSimRelations&, const TrackTruthLinks&, const TrackerHitSimRelations&)> {

  HungarianRecoMCTruthLinker(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,
                    {KeyValues("MCParticleCollection", {"MCParticles"}),
                     KeyValues("RecoParticleCollection", {"PandoraPFOs"}),
                     KeyValues("CalohitMCTruthLinkName", {"CalohitMCTruthLink"}),
                     KeyValues("SimCalorimeterHitRelationNames", {}),
                     KeyValues("TrackMCTruthLinkName", {"SiTracksMCTruthLink"}),
                     KeyValues("TrackerHitsRelInputCollections", {})},
                    KeyValues("RecoMCTruthLinkName", {"RecoMCTruthLink"})) {}

  StatusCode initialize() override {
    if (const auto sc = Transformer::initialize(); sc.isFailure()) {
      return sc;
    }
    auto exactlyOne = [this](const std::string& a, const std::string& b) {
      if (inputLocations(a).empty() == inputLocations(b).empty()) {
        error() << "Exactly one of " << a << " and " << b << " must be set" << endmsg;
        return false;
      }
      return true;
    };
    if (!exactlyOne("CalohitMCTruthLinkName", "SimCalorimeterHitRelationNames") ||
        !exactlyOne("TrackMCTruthLinkName", "TrackerHitsRelInputCollections")) {
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }

  RecoTruthLinkCollection operator()(const edm4hep::MCParticleCollection& mcParticles,
                                     const edm4hep::ReconstructedParticleCollection& recoParticles,
                                     const CaloTruthLinks& caloTruthLinks, const CaloSimRelations& caloSimRelations,
                                     const TrackTruthLinks& trackTruthLinks,
                                     const TrackerHitSimRelations& trackerHitSimRelations) const override {
    iou_truth::MatchConfig cfg;
    cfg.iouThreshold = m_iouThreshold.value();
    cfg.barrelRadius = m_caloBarrelRadius.value();
    cfg.nBarrelSides = m_caloBarrelSides.value();
    cfg.endCapZ = m_caloEndcapZ.value();

    iou_truth::KeySet allowedMC;
    allowedMC.reserve(mcParticles.size());
    for (const auto& mc : mcParticles) {
      allowedMC.insert(iou_truth::makeKey(mc.getObjectID()));
    }

    std::vector<iou_truth::Constituents> constituents;
    constituents.reserve(recoParticles.size());
    for (const auto& reco : recoParticles) {
      constituents.push_back(iou_truth::collectConstituents(reco));
    }

    const auto caloHitTruth = caloSimRelations.empty() ? iou_truth::bestTruthFromLinks(caloTruthLinks)
                                                       : iou_truth::bestCaloTruthFromSimRelations(caloSimRelations);
    const auto trackTruth = trackerHitSimRelations.empty()
                                ? iou_truth::bestTruthFromLinks(trackTruthLinks)
                                : iou_truth::bestTrackTruthFromHitRelations(constituents, trackerHitSimRelations);

    const auto matches = iou_truth::matchByIoU(constituents, caloHitTruth, trackTruth, allowedMC, cfg);

    auto links = RecoTruthLinkCollection{};
    for (std::size_t i = 0; i < recoParticles.size(); ++i) {
      if (!matches[i].matched) {
        continue;
      }
      auto link = links.create();
      link.setFrom(recoParticles[i]);
      link.setTo(matches[i].mc);
      link.setWeight(matches[i].iou);
    }

    debug() << "RecoParticles: " << recoParticles.size() << ", MCParticles: " << mcParticles.size()
            << ", links: " << links.size() << endmsg;

    return links;
  }

private:
  Gaudi::Property<float> m_iouThreshold{this, "IoUThreshold", 0.25f,
                                        "Minimum IoU for a ReconstructedParticle-to-MCParticle link"};

  Gaudi::Property<float> m_caloBarrelRadius{
      this, "CaloBarrelRadius", 2150.f,
      "Calorimeter inner barrel radius [mm], used to decide if an MCParticle vertex is inside the calorimeter"};

  Gaudi::Property<int> m_caloBarrelSides{this, "CaloBarrelSides", 12,
                                         "Number of sides of the calorimeter barrel polygon"};

  Gaudi::Property<float> m_caloEndcapZ{
      this, "CaloEndcapZ", 2307.f,
      "Calorimeter endcap inner |z| [mm], used to decide if an MCParticle vertex is inside the calorimeter"};
};

DECLARE_COMPONENT(HungarianRecoMCTruthLinker)
