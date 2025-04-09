from Bio import Entrez, SeqIO
import pandas as pd
import time
import random
import os

Entrez.email = "lorcan.brenders@gmail.com"

def fetch_ids(organ, search_term, retmax=10000):
    """Fetch protein IDs from Entrez in batches to avoid server blocking"""
    ids = []
    batch_size = 500
    for start in range(0, retmax, batch_size):
        handle = Entrez.esearch(
            db="protein", 
            term=f"{search_term}[Organism] AND {organ}[All Fields]", 
            retstart=start,
            retmax=batch_size
        )
        record = Entrez.read(handle)
        handle.close()
        ids.extend(record["IdList"])
        time.sleep(0.4)  # Be polite to NCBI
    return ids[:retmax]

def fetch_sequences(ids, organ):
    """Fetch protein sequences given a list of IDs"""
    protein_data = []
    batch_size = 100
    for i in range(0, len(ids), batch_size):
        print(f"Fetching batch {i // batch_size + 1} of {len(ids) // batch_size + 1}")
        batch_ids = ids[i:i + batch_size]
        try:
            handle = Entrez.efetch(db="protein", id=batch_ids, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "genbank")
            for record in records:
                protein_data.append({
                    "Organ": organ,
                    "Protein_Name": record.description.split(" ")[0],
                    "Amino_Acid_Sequence": str(record.seq),
                    "Sequence_Length": len(record.seq)
                })
            handle.close()
            time.sleep(0.4)
            print(f"Fetched {len(batch_ids)} sequences for {organ}")
        except Exception as e:
            print(f"Error processing batch {batch_ids}: {e}")
    return protein_data

organs = ["heart", "brain", "liver", "kidney", "muscle"]
search_term = "Homo sapiens"

all_data = []

for organ in organs:
    print(f"Fetching for {organ}...")
    ids = fetch_ids(organ, search_term, retmax=10000)
    seq_data = fetch_sequences(ids, organ)
    df = pd.DataFrame(seq_data)
    
    if df.shape[0] < 10000:
        print(f"Warning: only {df.shape[0]} sequences retrieved for {organ}")
    
    # Shuffle and split
    df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    train = df_shuffled.iloc[:7000]
    test_pool = df_shuffled.iloc[7000:]
    test_sample = test_pool.sample(n=150, random_state=42)

    # Save intermediate results
    train.to_csv(f"train_{organ}.csv", index=False)
    test_sample.to_csv(f"test_{organ}.csv", index=False)

    all_data.append((organ, test_sample))

# Combine all training data
train_all = pd.concat([pd.read_csv(f"train_{organ}.csv") for organ in organs])
train_all.to_csv("train_sequences.csv", index=False)

# Combine all test samples
test_all = pd.concat([sample for _, sample in all_data])
test_all.to_csv("test_sequences.csv", index=False)

# Generate majorite_{organe}.csv files
for target_organ in organs:
    main_part = test_all[test_all["Organ"] == target_organ].sample(n=random.randint(60, 105), random_state=42)
    other_part = test_all[test_all["Organ"] != target_organ].sample(n=150 - len(main_part), random_state=42)
    mixed_df = pd.concat([main_part, other_part]).sample(frac=1, random_state=42).reset_index(drop=True)
    filename = f"majorite_{target_organ}.csv"
    mixed_df.to_csv(filename, index=False)
    print(f"Saved {filename} with {len(main_part)} sequences from {target_organ}")

print("All done!")
