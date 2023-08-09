import sys
import numpy as np
import pandas as pd
from .FeaturesData import *

def write_output(inputData: FeaturesData, results:dict,cutoff=0.2) :

    countpos = {motif: 0 for motif in results}

    protein, res_id, motif_id, probability, negative_5_pos, motifseq, positive_5_pos = [], [], [], [], [], [], []

    for prot in inputData.seqData:
        seq = inputData.seqData[prot]
        seqLength = len(seq)

        for i, aa in enumerate(seq):
            for motif in results:
                #print(len(results[motif]))

                if i >= allMotifs[motif]["windLeft"] and i < seqLength - (allMotifs[motif]["motifSpan"] + allMotifs[motif]["windRight"]):
                    if round(results[motif][countpos[motif]][1], 4) >= cutoff:
                        protein.append(prot)
                        res_id.append(i+1)
                        motif_id.append(motif)
                        probability.append(100 * results[motif][countpos[motif]][1])
                        negative_5_pos.append(seq[i-5:i])
                        motifseq.append(seq[i:i + allMotifs[motif]["motifSpan"]])
                        positive_5_pos.append(seq[i + allMotifs[motif]["motifSpan"]:i + allMotifs[motif]["motifSpan"] + 5])
                
                    countpos[motif] += 1

    df = pd.DataFrame({'protein': protein, 'res_id': res_id, 'motif_id': motif_id, 'probability': probability, 'negative_5_pos': negative_5_pos, 'motifseq': motifseq, 'positive_5_pos': positive_5_pos})
    df.to_csv("test.csv", index=False, sep='\t')



