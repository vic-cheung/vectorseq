from dataclasses import dataclass
from itertools import chain
from typing import Tuple, List


@dataclass
class BrainGenes:
    NON_NEURON_MARKERS: Tuple[str] = (
        "Mag",
        "Mal",
        "Mog",
        "Olig1",
        "Olig2",
        "Trf",
        "P2ry12",
        "Cx3cr1",
        "C1qa",
        "Tmem119",
        "Aldh1l1",
        "Aqp4",
        "Gfap",
        "Slc1a2",
        "Flt1",
        "Mbp",
        "Pdgfra",
        "Aldh1a1",
        "Cx3cr1",
        "Cldn5",
        "Foxc1",
    )
    NEURON_MARKERS: Tuple[str] = ("Rbfox3", "Snap25")
    TG_MARKERS: Tuple[str] = (
        "Ef1a-FLPo",
        "Ef1a-mCherry-IRES-Cre",
        "hSyn-Cre",
        "tdTomato",
        "turboRFP",
        "CAV-GFP",
        "Cav2-GFP",
        "NLS-mRuby2",
        "ChR2-EYFP",
        "AAVrg-CAG-tdTomato",
        "AAVrg-CAG-GFP",
        "HSV-Cre",
        "CAV-Cre",
        "DreO",
        "fdio-EYFP",
    )
    EXCITATORY_MARKERS: Tuple[str] = ("Slc17a7", "Slc17a6", "Slc6a3", "Chat")
    INHIBITORY_MARKERS: Tuple[str] = ("Gad2", "Gad1")

    def all_genes(self) -> List[str]:
        return list(set(chain(*[gene_tuple for gene_tuple in self.__dict__.values()])))


@dataclass
class CortexGenes(BrainGenes):
    V1_MARKERS: Tuple[str] = (
        "Mog",
        "Mag",
        "Aqp4",
        "Slc1a2",
        "Slc17a7",
        "Camk2a",
        "Gad1",
        "Cntnap5a",
        "Selplg",
        "Ctss",
        "Flt1",
        "Cldn5",
        "Acta2",
        "Des",
    )
    V1_TG: Tuple[str] = (
        "DreO",
        "Ef1a-FLPo",
        "Ef1a-mCherry-IRES-Cre",
        "turboRFP",
        "fdio-EYFP",
        "tdTomato",
    )
    V1_MARKERS_WITH_TG: Tuple[str] = tuple(
        [*V1_TG, *V1_MARKERS,]
    )

    MICROGLIA_MARKERS: Tuple[str] = ("Selplg", "Ctss")
    ENDOTHELIA_MARKERS: Tuple[str] = ("Flt1", "Cldn5")
    MURAL_MARKERS: Tuple[str] = ("Acta2", "Des")
    OLIGODENDROCYTE_MARKERS: Tuple[str] = ("Mog", "Mag")
    ASTROCYTE_MARKERS: Tuple[str] = ("Aqp4", "Slc1a2")
    EXCITATORY_CLUSTER_MARKERS: Tuple[str] = ("Slc17a7", "Camk2a")
    INHIBITORY_CLUSTER_MARKERS: Tuple[str] = ("Gad1", "Cntnap5a")

    DICT_V1_MARKERS = {
        "microglia": MICROGLIA_MARKERS,
        "endothelia": ENDOTHELIA_MARKERS,
        "mural cells": MURAL_MARKERS,
        "oligodendrocytes": OLIGODENDROCYTE_MARKERS,
        "astrocytes": ASTROCYTE_MARKERS,
        "excitatory_markers": EXCITATORY_CLUSTER_MARKERS,
        "inhibitory markers": INHIBITORY_CLUSTER_MARKERS,
    }


@dataclass
class SuperiorColliculusGenes(BrainGenes):
    SCRG_EXCITATORY_MARKERS: Tuple[str] = (
        "Cmss1",
        "Cdk8",
        "Vwc2l",
        "Adcyap1",
        "Prlr",
        "Cpa6",
        "Slc5a7",
        "Pax5",
        "Slc10a4",
        "Slc18a2",
        "Mybpc1",
        "Hpse2",
        "Etv1",
        "Lrig1",
        "Plce1",
        "Tshz3",
        "Otx2os1",
        "Ndnf",
        "B230323A14Rik",
        "Rxfp2",
        "B130024G19Rik",
        "Ccdc192",
        "Tmem132c",
        "Cbln2",
        "Pitx2",
        "Pmfbp1",
        "Cartpt",
        "D130079A08Rik",
        "Piezo2",
        "Sulf1",
        "Sntb1",
        "Nek11",
        "Rorb",
        "C1ql1",
        "Trpc3",
        "Arhgap6",
        "Megf11",
        "Hmcn1",
        "Nek7",
        "Ttn",
        "Neb",
        "Plekhg1",
        "Synpo2",
        "Smoc2",
        "Glis3",
        "Gm10754",
        "Acss1",
    )

    SCRG_TG: Tuple[str] = (
        "AAVrg-CAG-tdTomato",
        "AAVrg-CAG-GFP",
        "HSV-Cre",
    )

    SCRG_EXCITATORY_MARKERS_WITH_TG: Tuple[str] = tuple(
        [*SCRG_TG, *SCRG_EXCITATORY_MARKERS,]
    )


@dataclass
class VentralMidbrainGenes(BrainGenes):
    VENTRAL_MIDBRAIN_MARKERS: Tuple[str] = (
        "Rln3",
        "Tac1",
        "Tcf7l2",
        "Fras1",
        "Cpa6",
        "Ndnf",
        "Cmss1",
        "Camk1d",
        "Sim1",
        "Gulp1",
        "Ebf2",
        "Trhde",
        "Ebf3",
        "Cdh6",
        "Pdzd2",
        "Nos1",
        "Ecel1",
        "Piezo2",
        "Dlx6os1",
        "Tcf4",
        "Pax6",
        "Cdh23",
        "Cep112",
        "Col25a1",
        "Otx2os1",
        "Pcbp3",
        "Gm12239",
        "Chd7",
        "Tmem132c",
        "Gm2164",
        "Sema3e",
        "Ror1",
        "Gm42951",
        "Gm41333",
    )

    VENTRAL_MIDBRAIN_TG: Tuple[str] = (
        "Ef1a-FLPo",
        "Ef1a-mCherry-IRES-Cre",
        "hSyn-Cre",
        "tdTomato",
        "turboRFP",
        "CAV-GFP",
        "NLS-mRuby2",
        "ChR2-EYFP",
    )

    VENTRAL_MIDBRAIN_MARKERS_WITH_TG: Tuple[str] = tuple(
        [*VENTRAL_MIDBRAIN_TG, *VENTRAL_MIDBRAIN_MARKERS,]
    )

    VENTRAL_MIDBRAIN_TG_SUBSET: Tuple[str] = (
        "Ef1a-FLPo",
        "Ef1a-mCherry-IRES-Cre",
    )

    VENTRAL_MIDBRAIN_MARKERS_WITH_TG_SUBSET: Tuple[str] = tuple(
        [*VENTRAL_MIDBRAIN_TG_SUBSET, *VENTRAL_MIDBRAIN_MARKERS,]
    )
