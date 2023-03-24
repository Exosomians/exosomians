import os.path
from collections import defaultdict
from typing import Optional

import pandas as pd
from Bio import SeqIO


def load_fasta(path: str, seq_key: Optional[str] = 'seq', target_key: Optional[str] = None) -> pd.DataFrame:
    records = defaultdict(list)
    with open(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = getattr(record, seq_key)
            records[seq_key].append(seq)

            if target_key is not None:
                label = getattr(record, target_key)
                records[target_key].append(label)

    df = pd.DataFrame(records)
    return df
