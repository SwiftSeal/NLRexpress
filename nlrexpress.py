import click
from src.ModuleData import *
from src.util import *
from datetime import datetime
import os
import pandas as pd
import argparse


sys.path.insert(0, os.path.abspath("."))


def predict( input: Path, outdir: Path, module: str, outformat:str, cpunum:int, skipjhmmer:bool, writeinputfile:bool ):
    """Predict NLR-related motifs"""

    scriptDir = Path(__file__).resolve().parent

    motifs = allMotifs

    inputData = generateFeatures( inputFasta=Path(input), outdir=Path(outdir), motifs=motifs, skipJhmmer=skipjhmmer, cpuNum=cpunum, annotations={}, writeInputFile=writeinputfile)
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

    write_output(inputData, results)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='NLRexpress: a tool for NLR protein subfamily assignment')
    parser.add_argument('-i', '--input', help='Input file in FASTA format', required=True)
    parser.add_argument('-o', '--outdir', help='Output directory', required=True)
    parser.add_argument('-m', '--module', help='Module to run (CC, TIR, NBS, LRR, all)', required=True)
    parser.add_argument('-f', '--outformat', help='Output format (csv, json)', required=False)
    parser.add_argument('-c', '--cpunum', help='Number of CPUs to use', required=False, default=4, type=int)
    parser.add_argument('--writeinputfile', help='Write input file', required=False, default=True, type=bool)
    parser.add_argument('--skipjhmmer', help='Skip jhmmert', required=False, default=False, type=bool)
    args = parser.parse_args()

    predict(args.input, args.outdir, args.module, args.outformat, args.cpunum, args.writeinputfile, args.skipjhmmer)










