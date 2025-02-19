import csv

def recmapHandler(rec_map, chr_start, chr_end, chunk_size):
    """
    Processes the recombination map CSV file using the csv module and returns
    the average recombination rate per chunk.
    
    Parameters:
        rec_map (str): Path to the recombination map CSV file.
        chr_start (int): Starting position of the chromosome.
        chr_end (int): Ending position of the chromosome.
        chunk_size (int): Size of each chunk in base pairs.
    
    Returns:
        list: Average recombination rate for each chunk.
        
    Raises:
        ValueError: If the CSV file does not have 'start' and 'rate' as headers,
                    or if any row has an invalid value for these columns.
    """
    rec_map_data = []
    with open(rec_map, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        # Check if header contains 'start' and 'rate'
        if reader.fieldnames is None or 'start' not in reader.fieldnames or 'rate' not in reader.fieldnames:
            raise ValueError("CSV file must have 'start' and 'rate' as header columns.")
        
        for row in reader:
            try:
                start_val = int(row['start'])
            except ValueError:
                raise ValueError(f"Column 'start' value '{row['start']}' is not a valid integer.")
            try:
                rate_val = float(row['rate'])
            except ValueError:
                raise ValueError(f"Column 'rate' value '{row['rate']}' is not a valid float.")
            rec_map_data.append({
                'start': start_val,
                'rate': rate_val
            })

    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size
    rec_rates = []
    for chunk in range(num_chunks):
        start_chunk = chr_start + chunk * chunk_size
        end_chunk = min(chr_end, start_chunk + chunk_size)
        rates = [row['rate'] for row in rec_map_data if start_chunk <= row['start'] < end_chunk]
        avg_rate = sum(rates) / len(rates) if rates else 1.0
        rec_rates.append(avg_rate)
    return rec_rates
