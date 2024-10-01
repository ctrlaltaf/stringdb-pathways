import pandas as pd


def read_file(path, delimiter):
    try:
        # Read the TSV file
        df = pd.read_csv(path, sep=delimiter)
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

    interaction = []
    for index, row in df_pathway.iterrows():
        node1 = row["node1_string_id"]
        node2 = row["node2_string_id"]
        combined = '_'.join(sorted([node1, node2]))
        interaction.append(combined)
    df_pathway["interaction"] = interaction

    print(df_pathway.head(3))
    print(df_stringdb.head(3))


    interaction = []
    for index, row in df_stringdb.iterrows():
        node1 = row["protein1"]
        node2 = row["protein2"]
        combined = '_'.join(sorted([node1, node2]))
        interaction.append(combined)
    df_stringdb["interaction"] = interaction

if __name__ == "__main__":
    main()