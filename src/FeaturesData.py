from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import logging
from datetime import datetime
import subprocess
from Bio import SeqIO

allMotifs = {
    "extEDVID": {"windLeft": 5, "windRight": 5, "motifSpan": 12},
    "bA": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "aA": {"windLeft": 5, "windRight": 5, "motifSpan": 7},
    "bC": {"windLeft": 5, "windRight": 5, "motifSpan": 8},
    "aC": {"windLeft": 5, "windRight": 5, "motifSpan": 6},
    "bDaD1": {"windLeft": 5, "windRight": 5, "motifSpan": 16},
    "aD3": {"windLeft": 5, "windRight": 5, "motifSpan": 13},
    "VG": {"windLeft": 5, "windRight": 5, "motifSpan": 5},
    "P-loop": {"windLeft": 5, "windRight": 5, "motifSpan": 9},
    "RNSB-A": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "Walker-B": {"windLeft": 5, "windRight": 5, "motifSpan": 8},
    "RNSB-B": {"windLeft": 5, "windRight": 5, "motifSpan": 7},
    "RNSB-C": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "RNSB-D": {"windLeft": 5, "windRight": 5, "motifSpan": 9},
    "GLPL": {"windLeft": 5, "windRight": 5, "motifSpan": 5},
    "MHD": {"windLeft": 5, "windRight": 5, "motifSpan": 3},
    "LxxLxL": {"windLeft": 5, "windRight": 5, "motifSpan": 6}
}

@dataclass
class FeaturesData:
    """

    """
    seqData: dict
    hmmData: dict

def generateFeatures(inputFasta: Path, output_directory: Path, threads: int) -> FeaturesData :

    processesInputFasta = Path(str(output_directory) + "/" + str(inputFasta.stem) + '.fasta_proc')
    seqData = processFastaFile( inputFasta, processesInputFasta )

    run_jackhmmer(processesInputFasta, output_directory, threads, 'hmmer_db/targetDB.fasta')

    logging.info('Preparing features: Parsing HMM profile - started')

    try:
        hmmFile1 = str(output_directory) + "/" + str(processesInputFasta.stem) + '-1.hmm'
        hmm_it1 = parse_hmm_multiprot(hmmFile1)
    except FileNotFoundError:
        raise FileNotFoundError('Preparing features: HMM profile iteration 1 was not found at. Execution stopped ')

    try:
        hmmFile2 = str(output_directory) + "/" + str(processesInputFasta.stem) + '-2.hmm'
        hmm_it2 = parse_hmm_multiprot(hmmFile2)

    except FileNotFoundError:
        logging.warning('Preparing features: HMM profile iteration 2 was not found at: ' + hmmFile2 + '. The first iteration HMM profile will be used alone.')
        hmm_it2 = parse_hmm_multiprot(hmmFile1)

    hmmData = generateInputFile(seqData, hmm_it1, hmm_it2)

    return FeaturesData(seqData = seqData, hmmData = hmmData)

def generateXmat (FeaturesData:dict, motif:str):
    X = []

    logging.info('Preparing features: NN input for motif ' + motif + ' started')

    for prot in FeaturesData.hmmData:
        seqLength = len( FeaturesData.seqData[prot] )

        for i in range(seqLength):

                windLeft = allMotifs[motif]["windLeft"]
                windRight = allMotifs[motif]["windRight"]
                motifSpan = allMotifs[motif]["motifSpan"]

                if i >= windLeft and i < seqLength - ( motifSpan + windRight) :
                    features = FeaturesData.hmmData[prot]

                    X.append([])
                    for w in range(windLeft * (-1), motifSpan + windRight + 1):
                        X[-1] += features[i+w]

    logging.info('Preparing features: NN input for motif ' + motif + ' done')

    return X

def generateInputFile( seqData:dict, hmm_it1:dict, hmm_it2:dict) -> dict:
    data = {}

    for name in seqData:
        seq = seqData[name]
        data[name] = []

        for i, aa in enumerate(seq):
            data[name].append([])
            for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                data[name][-1].append(val)

            if name in hmm_it2:
                for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                    data[name][-1].append(val)
            else:
                for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                    data[name][-1].append(val)

    return data

def processFastaFile(input:Path, output:Path) -> dict :
    seqData = {}

    with open(input, 'r') as inputFile, open(output, 'w') as outputFile:
        for record in SeqIO.parse(inputFile, "fasta"):
            seqData[record.id] = str(record.seq)
            SeqIO.write(record, outputFile, "fasta")

    return seqData

def run_jackhmmer(inputFasta: Path, output_directory: Path, threads: int, target_db: str):

    logging.info('jackhmmer - started')
    scriptDir = Path(__file__).resolve().parents[1]
    subprocess.run(["jackhmmer",
                    "--cpu", str(threads),
                    "-o", "/dev/null",
                    "-N", "2",
                    "-E", "1e-5",
                    "--domE", "1e-5",
                    "--noali",
                    "--chkhmm", str(output_directory) + "/" + str(inputFasta.stem),
                    inputFasta,
                    str(scriptDir) + '/' + target_db,
                    ],
                   stdout=subprocess.PIPE)
    logging.info('jackhmmer - done')

def parse_hmm_multiprot( hmmFile:Path ) -> dict:
    """

    :param hmmFile:
    :return:
    """

    # TODO Better try catch system to check for inconsistancies

    header1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    header2 = ["m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d"]

    hmm = {}

    try:
        with open(hmmFile, 'r') as f:
            lines = f.readlines()

            start = 0
            hascomp = 0
            length = 0

            for i, line in enumerate(lines):
                if line[0:4] == "NAME":
                    name = line.split()[1]
                    if name[-3:]=="-i1": name = name[:-3]

                    hmm[ name ] = []
                    start = 0
                    hascomp = 0

                elif line[0:4] == "LENG":
                    length = int( line.split()[1] )

                elif line[0:3] == 'HMM':
                    start = i + 5

                elif "COMPO" in line:
                    hascomp = 1

                elif start > 0 and i >= start and i <= 3*(start + length) :
                    if not hascomp:
                        logging.error("Problem parsing HMM file: no COMP line was found for prot " + name + line)
                        raise Exception("Problem parsing HMM file: no COMP line was found for prot ", name, line)

                    elif length == 0:
                        logging.error("Problem parsing HMM file: no LENGTH line was found for prot " + name)
                        raise Exception("Problem parsing HMM file: no LENGTH line was found for prot ", name)

                    else:
                        l = line.split()
                        if len(l) > 3:

                            if (i - start) % 3 == 0:
                                hmm[name].append({})
                                for i, val in enumerate(l[1:21]):
                                    if val == '*': val = 'inf'
                                    hmm[name][-1][header1[i]] = float(val)

        return hmm

    except FileNotFoundError:
        logging.error('Preparing features: HMM profile not found at: ' + str(hmmFile))
        raise FileNotFoundError('Preparing features: HMM profile not found at: ' + str(hmmFile))

    except:
       logging.error("Something went wrong when parsing the HMM file " + str(hmmFile) )
       raise Exception( "Something went wrong when parsing the HMM file ", str(hmmFile) )
