from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def concatenator(fs):
    ali_objects = [AlignIO.read(i, "fasta") for i in fs]
    l = [[rec.id for rec in x] for x in ali_objects]
    all_ids = set(([item for sublist in l for item in sublist]))
    d = dict(zip(all_ids,["" for x in range(len(all_ids))]))
    for ali in ali_objects:
        for rec in ali:
            seq = rec.seq
            spp = rec.id
            d[spp] = "".join([d[spp], str(seq)])
    return(MultipleSeqAlignment([SeqRecord(Seq(d[n]), id=n, description="") for n in d.keys()]))
