from src.ModuleData import *
from src.FeaturesData import *
import sys
import os
import argparse
import logging
import pandas as pd

sys.path.insert(0, os.path.abspath("."))

def predict(input: Path, output_directory: Path, threads: int):
    """Predict NLR-related motifs"""

    scriptDir = Path(__file__).resolve().parent

    motifs = allMotifs

    inputData = generateFeatures(inputFasta=Path(input), output_directory = Path(output_directory), threads = threads)
    results = {}

    CCexpress = ModuleData.loadModels(
        modelsPath = { 'extEDVID': str(scriptDir) + '/models/MLP_CC_extEDVID.pkl' })
    for p in CCexpress.predictors:
        X = generateXmat( inputData, p)
        results[p] = CCexpress.predictors[p].model.predict_proba( X )

    TIRexpress = ModuleData.loadModels(
        modelsPath={
                     'bA': str(scriptDir) + '/models/MLP_TIR_bA.pkl',
                     'aA': str(scriptDir) + '/models/MLP_TIR_aA.pkl',
                     'bC': str(scriptDir) + '/models/MLP_TIR_bC.pkl',
                     'aC': str(scriptDir) + '/models/MLP_TIR_aC.pkl',
                     'bDaD1': str(scriptDir) + '/models/MLP_TIR_bD-aD1.pkl',
                     'aD3': str(scriptDir) + '/models/MLP_TIR_aD3.pkl'
                     })
    
    for p in TIRexpress.predictors:
        X = generateXmat(inputData, p)
        results[p] = TIRexpress.predictors[p].model.predict_proba( X )

    NBSexpress = ModuleData.loadModels(
        modelsPath={
                      'VG': str(scriptDir) + '/models/MLP_NBS_VG.pkl',
                     'P-loop': str(scriptDir) + '/models/MLP_NBS_P-loop.pkl',
                     'RNSB-A': str(scriptDir) + '/models/MLP_NBS_RNSB-A.pkl',
                     'RNSB-B': str(scriptDir) + '/models/MLP_NBS_RNSB-B.pkl',
                     'RNSB-C': str(scriptDir) + '/models/MLP_NBS_RNSB-C.pkl',
                     'RNSB-D': str(scriptDir) + '/models/MLP_NBS_RNSB-D.pkl',
                     'Walker-B': str(scriptDir) + '/models/MLP_NBS_Walker-B.pkl',
                     'GLPL': str(scriptDir) + '/models/MLP_NBS_GLPL.pkl',
                     'MHD': str(scriptDir) + '/models/MLP_NBS_MHD.pkl'
                     })
    
    for p in NBSexpress.predictors:
        X = generateXmat(inputData, p)
        results[p] = NBSexpress.predictors[p].model.predict_proba( X )

    LRRexpress = ModuleData.loadModels(
        modelsPath={'LxxLxL': str(scriptDir) + '/models/MLP_LRR_LxxLxL.pkl'})
    for p in LRRexpress.predictors:
        X = generateXmat(inputData, p)
        results[p] = LRRexpress.predictors[p].model.predict_proba( X )

    write_output(inputData, results, output_directory)

def write_output(inputData: FeaturesData, results: dict, output_dir: Path, cutoff=0.2) :

    countpos = {motif: 0 for motif in results}

    protein, res_id, motif_id, probability, negative_5_pos, motifseq, positive_5_pos = [], [], [], [], [], [], []

    for prot in inputData.seqData:
        seq = inputData.seqData[prot]
        seqLength = len(seq)

        for i, aa in enumerate(seq):
            for motif in results:

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
    df.to_csv(str(output_dir) + '/nlrexpress.csv', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='NLRexpress: a tool for NLR protein subfamily assignment')
    parser.add_argument('-i', '--input', help='Input file in FASTA format', required=True)
    parser.add_argument('-o', '--output_directory', help='Output directory', required=True)
    parser.add_argument('-t', '--threads', help='Number of CPUs to use', required=False, default=4, type=int)
    parser.add_argument('--verbose', help='Verbose mode', required=False, action='store_true', default=False)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level = logging.INFO)

    predict(args.input, args.output_directory, args.threads)