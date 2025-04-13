import csv

# Step 1: Read the chromosome sizes file into a dictionary
def load_chr_sizes(file_path):
    chr_size_dict = {}
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            chr_name = row["chr"]
            chr_size = int(row["sizes"])
            chr_size_dict[chr_name] = chr_size
    return chr_size_dict
