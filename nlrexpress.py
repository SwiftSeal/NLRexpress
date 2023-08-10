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

def annotate(input: Path, output_directory: Path):

    CC_motifs = ["extEDVID"]
    NBS_motifs = ["VG", "P-loop", "RNSB-A", "RNSB-B", "RNSB-C", "RNSB-D", "Walker-B", "GLPL", "MHD"]
    TIR_motifs = ["bA", "aA", "bC", "aC", "bDaD1", "aD3"]
    LRR_motifs = ["LxxLxL"]

    threshold = 0.8
    allowed_gaps = 1

    # read input into pandas df
    df = pd.read_csv(input)

    # get unique proteins
    proteins = df['protein'].unique()

    for protein in proteins:
        annotation = []

        # subset df for protein
        df_subset = df[df["protein"] == protein]

        # ensure sorted by residue id
        df_subset = df_subset.sort_values(by=["res_id"])

        # get lists of data
        residues = df_subset["res_id"].tolist()
        motifs = df_subset["motif_id"].tolist()
        probabilities = df_subset["probability"].tolist()

        for i in range(len(residues)):
            # check if residue is a CC motif
            if motifs[i] in CC_motifs and probabilities[i] >= threshold:
                annotation.append("CC")
            elif motifs[i] in NBS_motifs and probabilities[i] >= threshold:
                annotation.append("NBS")
            elif motifs[i] in TIR_motifs and probabilities[i] >= threshold:
                annotation.append("TIR")
            elif motifs[i] in LRR_motifs and probabilities[i] >= threshold:
                annotation.append("LRR")
            else:
                raise Exception(f"Unknown motif at residue {residues[i]} of protein {protein}")
        
        print(annotation)

def main():
    parser = argparse.ArgumentParser("NLRexpress")
    parser.add_argument("--debug", action="store_true", help="Print debug messages")

    subparsers = parser.add_subparsers(dest="command")
    predict_parser = subparsers.add_parser("nlrexpress", help="Predict NLR-related motifs")
    predict_parser.add_argument("--input", type=str, required=True, help="Input FASTA file")
    predict_parser.add_argument("--output_directory", type=str, required=True, help="Output directory")
    predict_parser.add_argument("--threads", type=int, default=4, help="Number of threads")

    annotate_parser = subparsers.add_parser("annotate", help="Annotate proteins with NLR-related motifs")
    annotate_parser.add_argument("--input", type=str, required=True, help="Input FASTA file")
    annotate_parser.add_argument("--output_directory", type=str, required=True, help="Output directory")

    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    if args.command == "predict":
        predict(args.input, args.output_directory, args.threads)
    elif args.command == "annotate":
        annotate(args.input, args.output_directory)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()