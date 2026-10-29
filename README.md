# HungarianRecoMCTruthLinker

Gaudi (k4FWCore) algorithm that links each `ReconstructedParticle` (e.g.
`PandoraPFOs`, `HitPF`) to at most one `MCParticle`, using IoU and a
one-to-one (Hungarian-style) assignment.

## How to install

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
cd HungarianRecoMCTruthLinker
mkdir build && cd build
cmake .. -G Ninja -DCMAKE_INSTALL_PREFIX=../install
ninja install
cd .. && k4_local_repo        # puts ./install on the environment (run from the package directory)
```

In a new shell, source the Key4hep setup and run `k4_local_repo` from the
package directory again.

## How to run

```bash
k4run options/runHungarianRecoMCTruthLinker.py \
      --inputFiles CLD_REC.edm4hep.root --outputFile out.edm4hep.root \
      --RecoParticleCollection PandoraPFOs HitPF
```

- One linker runs per reco collection. Each writes
  `<collection>MCTruthLinkHungarian`, an
  `edm4hep::RecoMCParticleLinkCollection` whose weight is the IoU.
- By default, hit and track truth comes from `CalohitMCTruthLink` and
  `SiTracksMCTruthLink`. To use sim-hit relations instead:

  ```bash
  ... --SimCalorimeterHitRelationNames RelationCaloHit RelationMuonHit \
      --TrackerHitsRelInputCollections VXDTrackerHitRelations VXDEndcapTrackerHitRelations \
          InnerTrackerBarrelHitsRelations OuterTrackerBarrelHitsRelations \
          InnerTrackerEndcapHitsRelations OuterTrackerEndcapHitsRelations
  ```

## How to test

The tests need an EDM4hep input file with `PandoraPFOs`, `MCParticles`,
`CalohitMCTruthLink`, `SiTracksMCTruthLink` and the tracker-hit relations.
They run on the installed package, so install first and run `k4_local_repo`.

```bash
cd build
cmake .. -DHUNGARIANRECOMCTRUTHLINKER_TEST_INPUT=/path/to/CLD_REC.edm4hep.root
ninja install && ctest
```

---

## Algorithm

The matching is ported from
[k4PFHitML](https://github.com/key4hep/k4PFHitML)'s `TruthMatcher.h`:

1. **Constituents**: a particle's constituents are the calorimeter hits of
   all its clusters (including sub-clusters) and its tracks.
2. **Truth labels**: every hit and track gets one truth label, the
   MCParticle that contributed most to it. For a calorimeter hit, a
   secondary created inside the calorimeter is replaced by its ancestor that
   was produced outside it.
3. **IoU**: for each pair, IoU = |shared hits+tracks| / |union|. Pairs below
   `IoUThreshold` are dropped.
4. **Assignment**: the one-to-one assignment maximises the summed IoU. It
   uses Boost.Graph's header-only `maximum_weighted_matching`, which needs
   Boost >= 1.89 (the version with floating-point weights). Key4hep ships
   1.89 or newer from release 2026-02-01 on.

## Parameters

Parameter names follow MarlinReco's
[RecoMCTruthLinker](https://github.com/iLCSoft/MarlinReco/tree/master/Analysis/RecoMCTruthLink)
for the collections they share. So far only a subset of its inputs is
supported.

| Parameter | Type | Default |
|---|---|---|
| `MCParticleCollection` | `MCParticleCollection`; only these particles are used as truth | `MCParticles` |
| `RecoParticleCollection` | `ReconstructedParticleCollection` | `PandoraPFOs` |
| `CalohitMCTruthLinkName` | `CaloHitMCParticleLinkCollection` (list) | `CalohitMCTruthLink` |
| `SimCalorimeterHitRelationNames` | `CaloHitSimCaloHitLinkCollection` (list) | (empty) |
| `TrackMCTruthLinkName` | `TrackMCParticleLinkCollection` (list) | `SiTracksMCTruthLink` |
| `TrackerHitsRelInputCollections` | `TrackerHitSimTrackerHitLinkCollection` (list) | (empty) |
| `RecoMCTruthLinkName` | output `RecoMCParticleLinkCollection` | `RecoMCTruthLink` |
| `IoUThreshold` | float | 0.25 |
| `CaloBarrelRadius`, `CaloBarrelSides`, `CaloEndcapZ` | calorimeter geometry for remapping secondaries to their ancestors (CLD defaults) | 2150 mm, 12, 2307 mm |

Each detector has two interchangeable truth sources. Exactly one of each
pair must be non-empty, otherwise `initialize()` fails:

- **Calorimeter hits**:
  - `CalohitMCTruthLinkName`: the MCParticle with the highest link weight
    wins.
  - `SimCalorimeterHitRelationNames`: the MCParticle with the largest summed
    contribution energy wins.
- **Tracks**:
  - `TrackMCTruthLinkName`: the MCParticle with the highest link weight
    wins.
  - `TrackerHitsRelInputCollections`: the MCParticle that produced most of
    the track's hits wins.

A few behaviours differ from MarlinReco:

- **Inputs not needed**: MarlinReco's `SimCaloHitCollections`,
  `SimTrackerHitCollections`, `TrackCollection` and `ClusterCollection`
  aren't needed here. The sim hits are reached through the relations, and
  the tracks and clusters through the particles.
- **No MC → reco link collection**: EDM4hep links can be followed in both
  directions, so there is no separate `MCTruthRecoLinkName` output.
- **One link per particle**: MarlinReco (with `FullRecoRelation=true`) keeps
  every reco↔MC pair. This linker keeps one link per particle and one per
  MCParticle.
