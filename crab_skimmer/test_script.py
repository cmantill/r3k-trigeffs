#!/usr/bin/env python3
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from dataclasses import dataclass, field
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

# files = [
#     "file:/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023C/c482dcb8-40e3-4b89-bc20-ff60bd9b64b1.root",
# ] 

# # works but slow as a snail
# files = [
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/100000/a1e25d4f-8d3e-4fcf-911a-b704fe0cbb64.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/100000/b455e58a-d507-4106-969f-5b5f1eb50eb2.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/100000/f0705288-d0b4-46d0-aa0e-1249977b602a.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/100000/1b1426d9-58c1-46b9-b3e9-016d3b99cdaf.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/100000/a89ec040-7de7-4e03-a39b-dcf4f8aff51a.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/140000/03583a6e-6f8d-4077-a2b8-b54415b1669a.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/140000/75193ae8-8835-43dc-88c1-41bb7d66bd03.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2023C/ParkingDoubleMuonLowMass0/NANOAOD/NanoAODv15_v4-v1/140000/567eca87-e9e2-42e4-8a6a-6120924c6b4f.root",
# ]

# # 2023, CRISTINA'S NANOs
# files = [
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/06f415b4-a1c0-4533-853c-10ea56c6ae93.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/0e8b2090-2b43-4551-a5e3-96d0654fcb07.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/11208a37-dfdf-4cbc-8550-843a4d9af532.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/2165fe51-7fff-48bf-907d-05434f8d6688.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/2d3aa714-8407-4487-a4b2-55fa371696a8.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/443dd2ed-e80c-4f78-91ac-8d83503e288f.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/44e4ed4f-06c8-4a06-8d28-8b9cea46eb26.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/598c7f6b-fd90-4d8a-a990-e597561b6b43.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/65f23a3b-848b-4762-87e0-5bf13ae86095.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/7b6ee4a2-fe1d-4800-ab97-71826195949f.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/84e77735-8c5a-45d3-9ec3-9180fa85f36e.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/88d2919d-f256-4c6d-b70c-9a91f3c5160f.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/8955439d-1be9-4345-9888-5704cf263929.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/9b185efd-7591-452f-bd93-d783afcfa074.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/9c3c0827-6e5d-438a-88b7-1ae586c7fe73.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/b01a38ee-67e9-4fe6-9f89-63e2c3409ab8.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/20dc130f-f123-465e-8b05-7096e69d5f88.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/336e37ee-a452-4010-9fda-c11097ee6fc2.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/5a5713e1-b93e-4b18-a2da-bb14bcfdad4a.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/6bf805c8-0819-4505-8cfa-c4c454946030.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/8a3d149b-9014-44f4-b6c7-f58ae069f817.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/8bf1e3da-4365-4106-a097-2a2e60be8c25.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/94bcf94a-b3bf-46a5-a50c-c70c000c15c8.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/a0683440-f831-4d56-afd5-5ba9cb12420a.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/a51b9a18-576e-4026-a3d6-300de66dd05f.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/a53a2ed2-a26b-4e07-8e68-1f66e2f8a16a.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/affe815f-b34f-4e44-86a5-fa0dea1ad598.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/b04d651f-a99d-425c-8af3-8f9bc84ff4e4.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/b3f8dfe0-2f90-40e6-af69-929b42475399.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/b6e4d854-faeb-4a76-9b74-b27f40b9af58.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/c6713f03-3690-4287-b399-989ddf43bf0a.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/d9536504-57e0-41a7-988f-c9d2fd08cc73.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/e472526c-c3f7-4ed9-9025-419372306f94.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/e70a5739-9553-48b5-9df0-4a1e0bac0b08.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/f1ae30c4-a79e-4989-8b52-6bf290710eba.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/fe158785-6474-4778-b6ac-db75e012f14e.root",
#     # "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleMuonLowMass0/Run2023D/ff2b4b46-6c09-490d-930c-e7c12bf5eef0.root",
# ]
    
# # 2022, NOAH's NANOs
# files = [
#     f"/eos/cms/store/group/phys_bphys/DiElectronX/nzipper/ParkingDoubleMuonLowMass_NanoAOD/ParkingDoubleMuonLowMass0/Run3_ParkingDoubleMuonLowMass0_2022C_v1_Nano/221215_115329/0000/output_{i}.root" for i in range(500, 999) if os.path.exists(f"/eos/cms/store/group/phys_bphys/DiElectronX/nzipper/ParkingDoubleMuonLowMass_NanoAOD/ParkingDoubleMuonLowMass0/Run3_ParkingDoubleMuonLowMass0_2022C_v1_Nano/221215_115329/0000/output_{i}.root")
# ]

# # BuToKJpsi, 2022preEE
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/BuToKJPsi_JPsiToEE_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen/crab_BuToKJpsi_Toee_2022preEE/250417_151533/0000/DoubleElectronNANO_Run3_2022_mc_noskim_allNano_2025Apr17_{i}.root" for i in range(1, 29)
# ]

# # BuToKJpsi, 2022postEE
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/BuToKJPsi_JPsiToEE_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen/crab_BuToKJpsi_Toee_2022postEE/250417_151546/0000/DoubleElectronNANO_Run3_2022_mc_noskim_allNano_2025Apr17_{i}.root" for i in range(1, 105)
# ]

# # BuToKJpsi, 2023preBpix
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/BuToKJPsi_JPsiToEE_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen/crab_BuToKJpsi_Toee_2023preBPix/250417_151206/0000/DoubleElectronNANO_Run3_2023_mc_noskim_allNano_2025Apr17_{i}.root" for i in range(1, 138)
# ]

# # BuToKJpsi, 2023BPix
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/BuToKJPsi_JPsiToEE_SoftQCD_TuneCP5_13p6TeV_pythia8-evtgen/crab_BuToKJpsi_Toee_2023BPix/250417_151220/0000/DoubleElectronNANO_Run3_2023_mc_noskim_allNano_2025Apr17_{i}.root" for i in range(1, 85)
# # ]

# # prompt jpsi, 2022preEE
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8/crab_JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8_2022preEE/250411_142530/0000/DoubleElectronNANO_Run3_2022_mc_noskim_allNano_2025Apr11_{i}.root" for i in range(1, 57)
# ]

# # prompt jpsi, 2022postEE
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8/crab_JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8_2022postEE/250411_142540/0000/DoubleElectronNANO_Run3_2022_mc_noskim_allNano_2025Apr11_{i}.root" for i in range(1, 145)
# ]

# # prompt jpsi, 2023preBpix
# files = [
#     f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8/crab_JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8_2023preBPix/250411_104418/0000/DoubleElectronNANO_Run3_2023_mc_noskim_allNano_2025Apr11_{i}.root" for i in range(1, 68)
# ]

# prompt jpsi, 2023BPix
files = [
    f"/eos/cms/store/cmst3/group/xee/backgroundSamples/noskim_allnanoColl/JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8/crab_JPsiToEE_pth10toInf_TuneCP5_13p6TeV_pythia8_2023BPix/250411_085148/0000/DoubleElectronNANO_Run3_2023_mc_noskim_allNano_2025Apr11_{i}.root" for i in range(1, 52)
]

# # JETMET TEST (2022D)
# files = [
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/45e25b86-3fe3-487a-aa07-312549f537b5.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/48771e54-654a-47f3-8d80-f5513c1d5a49.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/b5a1076c-1a83-4d5d-8c82-aaad9dda3882.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/b7175047-0b29-4f6c-9f5a-4d6f7b8addb9.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/25f6f54d-5850-4300-a64c-eea4cb7691e9.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/89c115ee-d3f2-4614-89cf-40aacb55d0b7.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/c345cf6a-4228-4b2e-9d97-8f8750a46bb0.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/b4e6dd86-eb33-4248-8562-c3dbc159af3b.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/732070eb-d28f-49e5-886d-903c77039f9a.root",
#     "root://cms-xrd-global.cern.ch///store/data/Run2022D/JetMET/NANOAOD/22Sep2023-v1/2550000/0ecf0d80-d9c6-42b8-802a-73f65c21df95.root",
# # ]

# # DoubleEle 2023 (for reference)
# files = [
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_990.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_991.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_992.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_993.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_994.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_995.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_996.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_997.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_998.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass/crab_Run2023Dv1/251007_142120/0000/DoubleElectronNANO_Run3_2023_data_allNano_2025Oct07_999.root",
# ]

# # DoubleEle 2022 (for reference)
# files = [
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_990.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_991.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_992.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_993.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_994.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_995.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_996.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_997.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_998.root",
#     "/eos/cms/store/cmst3/group/xee/tree_v1/data/allnanoColl/ParkingDoubleElectronLowMass0/crab_Run2022Fv0/251023_172321/0000/DoubleElectronNANO_Run3_2022_data_allNano_2025Oct23_999.root",
# ]


# # check if files are actually there
# for f in files:
#     if not os.path.exists(f) and not f.startswith("root://"):
#         print(f"File {f} does not exist!")

# files = [f for f in files if os.path.exists(f)]

# when running on MC:
# - uncomment cut w/o muon trigger requirements
# - comment out json file
# (opposite when running on data)

@dataclass
class ParamSet:
    # outDir: str = 'skim_BuToKJpsi_2023preBPix'
    # outDir: str = 'skim_BuToKJpsi_2023BPix'
    # outDir: str = 'skim_BuToKJpsi_2022preEE'
    # outDir: str = 'skim_JetMET_2022D'
    # outDir: str = 'skim_DoubleEle_2022'
    # outDir: str = 'skim_DoubleEle_2023'
    # outDir: str = 'skim_promptJpsi_2022preEE'
    # outDir: str = 'skim_promptJpsi_2022postEE'
    # outDir: str = 'skim_promptJpsi_2023preBpix'
    outDir: str = 'skim_promptJpsi_2023BPix'
    inputFiles = files
    # # data
    # cut: str = 'nElectron > 1 && Electron_pt[0] > 4. && Electron_pt[1] > 4. && abs(Electron_eta[0]) < 1.22 && abs(Electron_eta[1]) < 1.22 && abs(Electron_dz[0] - Electron_dz[1]) <= 1. && sqrt(2*Electron_pt[0]*Electron_pt[1]*(cosh(Electron_eta[0]-Electron_eta[1]) - cos(Electron_phi[0]-Electron_phi[1]))) < 5. && sqrt(pow(Electron_eta[0]-Electron_eta[1],2) + pow(acos(cos(Electron_phi[0]-Electron_phi[1])),2)) > 0.03 && Electron_charge[0] + Electron_charge[1] == 0 && (HLT_DoubleMu4_3_Bs || HLT_DoubleMu4_3_Jpsi || HLT_DoubleMu4_3_LowMass || HLT_DoubleMu4_LowMass_Displaced || HLT_Mu0_L1DoubleMu || HLT_Mu4_L1DoubleMu || HLT_DoubleMu3_Trk_Tau3mu || HLT_DoubleMu3_TkMu_DsTau3Mu || HLT_DoubleMu4_MuMuTrk_Displaced || HLT_DoubleMu4_Jpsi_Displaced || HLT_DoubleMu4_Jpsi_NoVertexing || HLT_DoubleMu4_JpsiTrkTrk_Displaced || HLT_DoubleMu4_JpsiTrk_Bc || HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass || HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05)'
    # MC
    cut: str = 'nElectron > 1 && Electron_pt[0] > 4. && Electron_pt[1] > 4. && abs(Electron_eta[0]) < 1.22 && abs(Electron_eta[1]) < 1.22 && abs(Electron_dz[0] - Electron_dz[1]) <= 1. && sqrt(2*Electron_pt[0]*Electron_pt[1]*(cosh(Electron_eta[0]-Electron_eta[1]) - cos(Electron_phi[0]-Electron_phi[1]))) < 5. && sqrt(pow(Electron_eta[0]-Electron_eta[1],2) + pow(acos(cos(Electron_phi[0]-Electron_phi[1])),2)) > 0.03 && Electron_charge[0] + Electron_charge[1] == 0'
    # json: str = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'
    # json: str = '/eos/home-n/npalmeri/DiEleAnalyzer/trigger_sfs/CMSSW_14_0_6/src/r3k-trigeffs/json_files/jay/trigger_OR_2022_jay.json'
    # json: str = '/eos/home-n/npalmeri/DiEleAnalyzer/trigger_sfs/CMSSW_14_0_6/src/r3k-trigeffs/json_files/Eras_CDEFG/trigger_OR.json'
    json : str = None,
    branchesIn: list[str] = field(default_factory=lambda: ["drop *", "keep DiElectron_*", "keep Electron_*", "keep L1_DoubleEG*_er1p2_dR_Max*", "keep HLT_DoubleEle*_mMax6*", "keep PV_*", "keep GenPart_*"])
    branchesOut: str = None
    label: str = '_skim'
    friend: bool = False
    maxEntries: int = None
    haddFileName: str = 'skim_promptJpsi_2023BPix/output_skim.root'

class SkimEvents(Module):
    def __init__(self):
        self.writeHistFile = True
    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
    def analyze(self, event):
        return True


params = ParamSet()
p = PostProcessor(
        params.outDir,
        params.inputFiles,
        cut=params.cut,
        postfix=params.label,
        branchsel=params.branchesIn,
        outputbranchsel=params.branchesOut,
        modules=[SkimEvents()],
        haddFileName=params.haddFileName,
        compression=('LZMA:9'),
        friend=params.friend,
        maxEntries=params.maxEntries,
        # jsonInput=params.json, #comment if running on MC
        provenance=True,
        fwkJobReport=True,
)


p.run()

print("DONE")
