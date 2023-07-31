import sys
import re
import pandas as pd
from io import StringIO

CC_MOTIFS = ["extEDVID"]
TIR_MOTIFS = ["bA",
              "aA",
              "bC",
              "aC",
              "bDaD1",
              "aD3"]
NBARC_MOTIFS = ["VG",
                "P-loop",
                "RNSB-A",
                "Walker-B",
                "RNSB-B",
                "RNSB-C",
                "GLPL",
                "RNSB-D",
                "MHD"]
LRR_MOTIFS = ["LxxLxL"]

MOTIFS = {"CC": CC_MOTIFS, "TIR": TIR_MOTIFS, "NBARC": NBARC_MOTIFS, "LRR": LRR_MOTIFS}

def clean_input(input_data):
    # The nlrexpress long output is weirdly formatted.
    # Some columns have | to delimit them.
    # Some (not all) columns are separated by tabs.
    # Whitespace is used to pad the columns.
    # To fix this, Deleta all |, replace all whitespaces with a tab, then replace strings of tabs with a single tab.

    input_data = input_data.replace("|", "")

    # Replace all whitespace with a single tab
    input_data = input_data.replace(" ", "\t")

    # Replace all strings of tabs with a single tab
    input_data = re.sub(r'\t+', '\t', input_data)

    # trim tailing tab from each line
    input_data = re.sub(r'\t$', '', input_data, flags=re.MULTILINE)

    # create a pandas dataframe from the tabular data
    df = pd.read_csv(StringIO(input_data), sep="\t")

    # remove any row where 'Proba' is < 80
    df = df[df['Proba'] >= 80]

    return df

def classify_motif(motif):
    # classify motif as CC, TIR, NBARC, or LRR
    for motif_type, motifs in MOTIFS.items():
        if motif in motifs:
            return motif_type

    return "Unknown"

def classify_NLR(df):

    prot_names = df['#ProtName'].unique()

    for prot_name in prot_names:
        print(prot_name)

        df_prot = df[df['#ProtName'] == prot_name]
        
        # sort df_prot by 'ResID' and get list of motifs from 'Motif' column
        motifs = df_prot.sort_values(by=['ResId'])['Motif'].tolist()

        # Get the order of CC, TIR, NBARC, and LRR domains
        # To do this, iterate through the motif list.
        # Identify the MOTIFS class the current motif belongs to.
        # If the current motif classification is not the same as the current classification, then the current motif is the first motif of the current classification.
        #

        domains = []

        current_classification = ""
        for motif in motifs:
            motif_type = classify_motif(motif)

            if motif_type != current_classification:
                domains.append(motif_type)
                current_classification = motif_type
        
        print(domains)



        

def main():
    try:
        # Get the raw tabular output of nlrexpress from file path
        with open("sample\output_ref\zar1_rpp1.short.output.txt", "r") as f:
            input_data = f.read()
        
        # Clean the input data
        df = clean_input(input_data)

        # Classify the NLR
        classify_NLR(df)

    except KeyboardInterrupt:
        # Handling KeyboardInterrupt (Ctrl+C) gracefully
        print("Script terminated by user.")
        sys.exit(0)

if __name__ == "__main__":
    main()
