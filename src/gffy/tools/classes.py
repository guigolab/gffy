from typing import Optional
from array import array


class Gene:
    __slots__ = (
        "feature_type",
        "biotype",
        "length",
        "has_cds",
        "has_exon",
        "has_multiple_exons",
        "category")

    def __init__(self, feature_type: str, biotype: Optional[str], length: int):
        self.feature_type = feature_type
        self.biotype = biotype
        self.length = length
        self.has_cds = False
        self.has_exon = False
        self.has_multiple_exons = False
        self.category = None

    def set_category(self):
        if self.feature_type == "pseudogene":
            self.category = "pseudogene"
            return
        elif self.has_cds or self.biotype == "protein_coding":
            self.category = "coding"
            return
        elif self.has_exon:
            if self.has_multiple_exons or self.length > 200:
                self.category = "long_non_coding"
            else:
                self.category = "short_non_coding"
        
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
        self.gene_id = gene_id
        self.type = ttype
        self.length = 0
        self.exons_lengths = array("i")
        self.exon_len_sum = 0
        self.exon_count = 0
        self.cds_len_sum = 0
        self.cds_count = 0
        self.cds_lengths = array("i")


class OrphanFeature:
    __slots__ = ("feature_id", "feature_type", "length", "parent_ids", "biotype", "resolved")

    def __init__(
        self,
        feature_id: str,
        feature_type: str,
        length: int,
        parent_ids: list,
        biotype: Optional[str],
    ):
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.length = length
        self.parent_ids = parent_ids
        self.biotype = biotype
        self.resolved = False
