from typing import Optional
from array import array
import sys


class Gene:
    __slots__ = ("feature_type", "biotype", "length", "has_cds", "has_exon", "has_multiple_exons", "has_children", "category")

    def __init__(self, feature_type: str, biotype: Optional[str], length: int):
        self.feature_type = sys.intern(feature_type)
        self.biotype = sys.intern(biotype) if biotype else None
        self.length = length
        self.has_cds = False
        self.has_exon = False
        self.has_multiple_exons = False
        self.has_children = False
        self.category = None

    def set_category(self):
        if self.feature_type == "pseudogene":
            self.category = sys.intern("pseudogene")
        elif self.has_cds or self.biotype == "protein_coding":
            self.category = sys.intern("coding")
        elif self.has_exon:
            if self.has_multiple_exons or self.length > 200:
                self.category = sys.intern("long_non_coding")
            else:
                self.category = sys.intern("short_non_coding")
        
class Transcript:
    __slots__ = (
        "gene_id",
        "type",
        "exons_lengths",
        "exon_len_sum",
        "exon_count",
        "cds_len_sum",
        "cds_count",
        "cds_lengths",
        "length",
    )

    def __init__(self, gene_id: str, ttype: Optional[str]):
        self.gene_id = sys.intern(gene_id)
        self.type = sys.intern(ttype) if ttype else None
        self.length = 0
        self.exons_lengths = array("i")
        self.exon_len_sum = 0
        self.exon_count = 0
        self.cds_len_sum = 0
        self.cds_count = 0
        self.cds_lengths = array("i")


class OrphanFeature:
    __slots__ = ("feature_id", "feature_type", "length", "parent_id", "biotype", "resolved")

    def __init__(
        self,
        feature_id: str,
        feature_type: str,
        length: int,
        parent_id: Optional[str],
        biotype: Optional[str],
    ):
        self.feature_id = feature_id
        self.feature_type = sys.intern(feature_type)
        self.length = length
        self.parent_id = parent_id
        self.biotype = sys.intern(biotype) if biotype else None
        self.resolved = False
