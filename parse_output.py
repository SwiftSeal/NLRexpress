import sys
import re
import pandas as pd

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

    # create a pandas dataframe from the tabular data
    df = pd.read_csv(input_data, sep="\t")

    return df

def main():
    try:
        # Get the raw tabular output of nlrexpress from stdin
        input_data = sys.stdin.read().strip()

        df = clean_input(input_data)

        print(df)

    except KeyboardInterrupt:
        # Handling KeyboardInterrupt (Ctrl+C) gracefully
        print("Script terminated by user.")
        sys.exit(0)

if __name__ == "__main__":
    main()
