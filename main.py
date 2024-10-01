import pandas as pd


def read_file(path, delimiter):
    try:
        # Read the TSV file
        df = pd.read_csv(path, sep='\t')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def main():
    pathway_path = "string_interactions_short.tsv"
    stringdb_path = "9606.protein.links.v12.0.txt"

    # read in files
    df_pathway = read_file(pathway_path, "\t")
    df_stringdb = read_file(stringdb_path, " ")

    print(df_pathway.head(3))
    print("")
    print(df_stringdb.head(3))

if __name__ == "__main__":
    main()