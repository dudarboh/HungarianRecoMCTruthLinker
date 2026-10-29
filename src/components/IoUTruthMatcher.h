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
#pragma once

// IoU-based truth matching of reconstructed particles to MCParticles, ported from k4PFHitML's TruthMatcher.h
// (https://github.com/key4hep/k4PFHitML). The one-to-one assignment uses Boost.Graph's maximum weighted matching.

#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "podio/ObjectID.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numbers>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace iou_truth {

struct MatchConfig {
  float iouThreshold{0.25f};
  float barrelRadius{2150.f};
  int nBarrelSides{12};
  float endCapZ{2307.f};
};

struct ObjectIDKey {
  int index{-1};
  unsigned collectionID{0};

  bool operator==(const ObjectIDKey& other) const { return index == other.index && collectionID == other.collectionID; }
  bool operator<(const ObjectIDKey& other) const {
    return index != other.index ? index < other.index : collectionID < other.collectionID;
  }

  bool valid() const { return index >= 0; }
};

struct ObjectIDKeyHash {
  std::size_t operator()(const ObjectIDKey& key) const {
    const auto a = static_cast<std::uint64_t>(key.collectionID);
    const auto b = static_cast<std::uint32_t>(key.index);
    return std::hash<std::uint64_t>{}((a << 32) ^ b);
  }
};

template <typename T>
using KeyMap = std::unordered_map<ObjectIDKey, T, ObjectIDKeyHash>;
using KeySet = std::unordered_set<ObjectIDKey, ObjectIDKeyHash>;

inline ObjectIDKey makeKey(const podio::ObjectID& id) { return ObjectIDKey{id.index, id.collectionID}; }

// Hit/track -> the MCParticle that contributed most to it.
using BestTruthMap = KeyMap<edm4hep::MCParticle>;

struct Match {
  bool matched{false};
  edm4hep::MCParticle mc{};
  float iou{0.f};
};

// Calorimeter hits (from all clusters, including sub-clusters) and tracks of a reconstructed particle, each once.
struct Constituents {
  std::vector<edm4hep::CalorimeterHit> caloHits;
  std::vector<edm4hep::Track> tracks;
};

inline Constituents collectConstituents(const edm4hep::ReconstructedParticle& reco) {
  Constituents out;
  KeySet seenHits, seenTracks, seenClusters;

  std::function<void(const edm4hep::Cluster&)> addCluster = [&](const edm4hep::Cluster& cluster) {
    if (!seenClusters.insert(makeKey(cluster.getObjectID())).second) {
      return;
    }
    for (const auto& hit : cluster.getHits()) {
      if (seenHits.insert(makeKey(hit.getObjectID())).second) {
        out.caloHits.push_back(hit);
      }
    }
    for (const auto& sub : cluster.getClusters()) {
      addCluster(sub);
    }
  };

  for (const auto& cluster : reco.getClusters()) {
    addCluster(cluster);
  }
  for (const auto& trk : reco.getTracks()) {
    if (seenTracks.insert(makeKey(trk.getObjectID())).second) {
      out.tracks.push_back(trk);
    }
  }
  return out;
}

// ---- Truth sources: each builds a BestTruthMap, so they are interchangeable ----

// From ready-made object -> MCParticle links (e.g. CalohitMCTruthLink, SiTracksMCTruthLink): highest weight wins.
template <typename LinkCollectionT>
BestTruthMap bestTruthFromLinks(const std::vector<const LinkCollectionT*>& linkCollections) {
  KeyMap<float> bestWeight;
  BestTruthMap result;
  for (const auto* links : linkCollections) {
    for (const auto& link : *links) {
      const auto key = makeKey(link.getFrom().getObjectID());
      const auto [it, inserted] = bestWeight.try_emplace(key, link.getWeight());
      if (inserted || link.getWeight() > it->second) {
        it->second = link.getWeight();
        result[key] = link.getTo();
      }
    }
  }
  return result;
}

// From CalorimeterHit -> SimCalorimeterHit relations: the MCParticle with the largest summed contribution energy.
inline BestTruthMap
bestCaloTruthFromSimRelations(const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& relations) {
  BestTruthMap result;
  for (const auto* coll : relations) {
    for (const auto& rel : *coll) {
      KeyMap<float> energyPerMC;
      KeyMap<edm4hep::MCParticle> mcByKey;
      for (const auto& contrib : rel.getTo().getContributions()) {
        const auto mc = contrib.getParticle();
        if (!mc.isAvailable()) {
          continue;
        }
        const auto mcKey = makeKey(mc.getObjectID());
        energyPerMC[mcKey] += contrib.getEnergy();
        mcByKey.emplace(mcKey, mc);
      }
      const auto best = std::ranges::max_element(energyPerMC, {}, [](const auto& kv) { return kv.second; });
      if (best != energyPerMC.end()) {
        result[makeKey(rel.getFrom().getObjectID())] = mcByKey.at(best->first);
      }
    }
  }
  return result;
}

// From TrackerHit -> SimTrackerHit relations: per track, the MCParticle that produced most of its hits.
inline BestTruthMap
bestTrackTruthFromHitRelations(const std::vector<Constituents>& particles,
                               const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& relations) {
  BestTruthMap hitTruth;
  for (const auto* coll : relations) {
    for (const auto& rel : *coll) {
      const auto mc = rel.getTo().getParticle();
      if (mc.isAvailable()) {
        hitTruth[makeKey(rel.getFrom().getObjectID())] = mc;
      }
    }
  }

  BestTruthMap result;
  for (const auto& particle : particles) {
    for (const auto& trk : particle.tracks) {
      KeyMap<int> hitsPerMC;
      KeyMap<edm4hep::MCParticle> mcByKey;
      for (const auto& hit : trk.getTrackerHits()) {
        const auto it = hitTruth.find(makeKey(hit.getObjectID()));
        if (it == hitTruth.end()) {
          continue;
        }
        const auto mcKey = makeKey(it->second.getObjectID());
        hitsPerMC[mcKey] += 1;
        mcByKey.emplace(mcKey, it->second);
      }
      const auto best = std::ranges::max_element(hitsPerMC, {}, [](const auto& kv) { return kv.second; });
      if (best != hitsPerMC.end()) {
        result[makeKey(trk.getObjectID())] = mcByKey.at(best->first);
      }
    }
  }
  return result;
}

namespace detail {

  constexpr int BIT_BACKSCATTER = 29;
  constexpr int BIT_DECAYED_IN_TRACKER = 27;

  inline bool checkBit(int value, int bit) { return ((value >> bit) & 1) == 1; }

  inline bool backscatteredAndTracker(int simStatus) {
    return checkBit(simStatus, BIT_BACKSCATTER) && !checkBit(simStatus, BIT_DECAYED_IN_TRACKER);
  }

  inline bool isProducedInCalo(const edm4hep::Vector3d& vtx, const MatchConfig& cfg) {
    for (int i = 0; i < cfg.nBarrelSides; ++i) {
      const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(cfg.nBarrelSides);
      if (vtx.x * std::cos(angle) + vtx.y * std::sin(angle) > cfg.barrelRadius) {
        return true;
      }
    }
    return std::abs(vtx.z) > cfg.endCapZ;
  }

  // Secondaries created inside the calorimeter are attributed to their first ancestor produced outside of it.
  inline edm4hep::MCParticle remapToCaloMother(edm4hep::MCParticle mc, const MatchConfig& cfg) {
    while (mc.isAvailable()) {
      if (!isProducedInCalo(mc.getVertex(), cfg) || backscatteredAndTracker(mc.getSimulatorStatus())) {
        break;
      }
      const auto parents = mc.getParents();
      if (parents.empty()) {
        break;
      }
      mc = parents[0];
    }
    return mc;
  }

  // One-to-one assignment maximising the summed weight; weights[r][c] <= 0 means "cannot be matched".
  // Returns, per column, the matched row or -1.
  inline std::vector<int> maximumWeightAssignment(const std::vector<std::vector<float>>& weights, std::size_t nCols) {
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
                                        boost::property<boost::edge_weight_t, double>>;
    const std::size_t nRows = weights.size();
    Graph g(nRows + nCols);
    for (std::size_t r = 0; r < nRows; ++r) {
      for (std::size_t c = 0; c < nCols; ++c) {
        if (weights[r][c] > 0.f) {
          boost::add_edge(r, nRows + c, static_cast<double>(weights[r][c]), g);
        }
      }
    }

    std::vector<boost::graph_traits<Graph>::vertex_descriptor> mate(nRows + nCols);
    boost::maximum_weighted_matching(g, &mate[0]);

    std::vector<int> rowOfCol(nCols, -1);
    for (std::size_t c = 0; c < nCols; ++c) {
      const auto m = mate[nRows + c];
      if (m != boost::graph_traits<Graph>::null_vertex()) {
        rowOfCol[c] = static_cast<int>(m);
      }
    }
    return rowOfCol;
  }

} // namespace detail

// Matches every reconstructed particle (given by its constituents) to at most one MCParticle and vice versa,
// maximising the summed IoU. IoU = |shared hits+tracks| / |union of hits+tracks|; pairs below cfg.iouThreshold are
// never matched. Calorimeter-hit truth is remapped to the ancestor produced outside the calorimeter; truth particles
// not in `allowedMC` count as background. Result is indexed like `particles`.
inline std::vector<Match> matchByIoU(const std::vector<Constituents>& particles, const BestTruthMap& caloHitTruth,
                                     const BestTruthMap& trackTruth, const KeySet& allowedMC,
                                     const MatchConfig& cfg = {}) {
  const ObjectIDKey backgroundKey{-1, 0};

  std::vector<KeyMap<int>> intersections(particles.size());
  KeyMap<int> truthSizes;
  KeyMap<edm4hep::MCParticle> truthParticles;
  std::vector<int> predSizes(particles.size(), 0);

  auto count = [&](std::size_t i, const ObjectIDKey& objKey, const BestTruthMap& truthMap, bool remapCaloMother) {
    ++predSizes[i];
    ObjectIDKey truthKey = backgroundKey;
    if (const auto it = truthMap.find(objKey); it != truthMap.end()) {
      const auto mc = remapCaloMother ? detail::remapToCaloMother(it->second, cfg) : it->second;
      const auto mcKey = makeKey(mc.getObjectID());
      if (mc.isAvailable() && allowedMC.contains(mcKey)) {
        truthKey = mcKey;
        truthParticles.emplace(mcKey, mc);
      }
    }
    intersections[i][truthKey] += 1;
    truthSizes[truthKey] += 1;
  };

  for (std::size_t i = 0; i < particles.size(); ++i) {
    for (const auto& hit : particles[i].caloHits) {
      count(i, makeKey(hit.getObjectID()), caloHitTruth, true);
    }
    for (const auto& trk : particles[i].tracks) {
      count(i, makeKey(trk.getObjectID()), trackTruth, false);
    }
  }

  std::vector<ObjectIDKey> truthKeys;
  truthKeys.reserve(truthSizes.size());
  for (const auto& [k, _] : truthSizes) {
    if (k.valid()) {
      truthKeys.push_back(k);
    }
  }
  std::sort(truthKeys.begin(), truthKeys.end());

  std::vector<std::vector<float>> iou(truthKeys.size(), std::vector<float>(particles.size(), 0.f));
  for (std::size_t t = 0; t < truthKeys.size(); ++t) {
    const int truthSize = truthSizes[truthKeys[t]];
    for (std::size_t s = 0; s < particles.size(); ++s) {
      const auto it = intersections[s].find(truthKeys[t]);
      if (it == intersections[s].end()) {
        continue;
      }
      const float value = static_cast<float>(it->second) / static_cast<float>(predSizes[s] + truthSize - it->second);
      iou[t][s] = value < cfg.iouThreshold ? 0.f : value;
    }
  }

  const auto truthRowOf = detail::maximumWeightAssignment(iou, particles.size());

  std::vector<Match> result(particles.size());
  for (std::size_t s = 0; s < particles.size(); ++s) {
    if (const int t = truthRowOf[s]; t >= 0) {
      result[s] = Match{true, truthParticles.at(truthKeys[t]), iou[t][s]};
    }
  }
  return result;
}

} // namespace iou_truth
