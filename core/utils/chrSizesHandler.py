import csv

def load_chr_sizes(file_path):
    chr_size_dict = {}
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # skip header row
        for row in reader:
            if len(row) < 2:
                continue  # skip incomplete lines
            chr_name = row[0].strip()
            try:
                chr_size = int(row[1])
                chr_size_dict[chr_name] = chr_size
            except ValueError:
                print(f"Skipping invalid size for {chr_name}: {row[1]}")
    return chr_size_dict
