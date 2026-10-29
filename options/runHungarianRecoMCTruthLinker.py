#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Configurables import HungarianRecoMCTruthLinker
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser


parser_group = parser.add_argument_group("runHungarianRecoMCTruthLinker.py custom options")
parser_group.add_argument("--inputFiles", nargs="+", metavar=("file1", "file2"), required=True,
                          help="One or multiple input files")
parser_group.add_argument("--outputFile", default="output_HungarianRecoMCTruthLink.edm4hep.root")
parser_group.add_argument("--num_ev", type=int, help="Number of events to process (-1 means all)", default=-1)
parser_group.add_argument("--RecoParticleCollection", nargs="+", default=["PandoraPFOs"],
                          help="ReconstructedParticle collection(s) to link (e.g. PandoraPFOs HitPF); "
                               "one linker is run per collection")
parser_group.add_argument("--RecoMCTruthLinkSuffix", default="MCTruthLinkHungarian",
                          help="Output link collection is named <RecoParticleCollection><suffix>")
parser_group.add_argument("--MCParticleCollection", default="MCParticles")
parser_group.add_argument("--IoUThreshold", type=float, default=0.25,
                          help="Minimum IoU for a ReconstructedParticle-to-MCParticle link")

# Calorimeter hit truth: either CaloHit->MCParticle links, or CaloHit->SimCaloHit relations
calo = parser_group.add_mutually_exclusive_group()
calo.add_argument("--CalohitMCTruthLinkName", nargs="+", help="default: CalohitMCTruthLink")
calo.add_argument("--SimCalorimeterHitRelationNames", nargs="+",
                  help="e.g. RelationCaloHit RelationMuonHit (replaces --CalohitMCTruthLinkName)")
# Track truth: either Track->MCParticle links, or TrackerHit->SimTrackerHit relations
trk = parser_group.add_mutually_exclusive_group()
trk.add_argument("--TrackMCTruthLinkName", nargs="+", help="default: SiTracksMCTruthLink")
trk.add_argument("--TrackerHitsRelInputCollections", nargs="+",
                 help="e.g. VXDTrackerHitRelations VXDEndcapTrackerHitRelations InnerTrackerBarrelHitsRelations "
                      "OuterTrackerBarrelHitsRelations InnerTrackerEndcapHitsRelations "
                      "OuterTrackerEndcapHitsRelations (replaces --TrackMCTruthLinkName)")
args = parser.parse_known_args()[0]

if args.SimCalorimeterHitRelationNames:
    calo_inputs = dict(CalohitMCTruthLinkName=[], SimCalorimeterHitRelationNames=args.SimCalorimeterHitRelationNames)
else:
    calo_inputs = dict(CalohitMCTruthLinkName=args.CalohitMCTruthLinkName or ["CalohitMCTruthLink"],
                       SimCalorimeterHitRelationNames=[])
if args.TrackerHitsRelInputCollections:
    track_inputs = dict(TrackMCTruthLinkName=[], TrackerHitsRelInputCollections=args.TrackerHitsRelInputCollections)
else:
    track_inputs = dict(TrackMCTruthLinkName=args.TrackMCTruthLinkName or ["SiTracksMCTruthLink"],
                        TrackerHitsRelInputCollections=[])

svc = IOSvc("IOSvc")
svc.Input = args.inputFiles
svc.Output = args.outputFile
svc.outputCommands = ["keep *"]

linkers = [
    HungarianRecoMCTruthLinker(f"HungarianRecoMCTruthLinker_{reco}",
                               MCParticleCollection=[args.MCParticleCollection],
                               RecoParticleCollection=[reco],
                               RecoMCTruthLinkName=[reco + args.RecoMCTruthLinkSuffix],
                               IoUThreshold=args.IoUThreshold,
                               **calo_inputs,
                               **track_inputs,
                               )
    for reco in args.RecoParticleCollection
]

ApplicationMgr(TopAlg=linkers,
               EvtSel="NONE",
               EvtMax=args.num_ev,
               ExtSvc=[svc],
               OutputLevel=INFO,
               )
