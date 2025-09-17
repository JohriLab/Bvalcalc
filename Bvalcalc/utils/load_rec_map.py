import csv
import numpy as np
import sys

def load_rec_map(rec_map, calc_start, calc_end, chunk_size, chromosome):
    """
    Processes the recombination map file and returns
    the average recombination rate per chunk, weighted by the proportion of each 
    chunk that falls within a given recombination interval.
    
    The file format is CSV with each row having exactly three columns:
    chromosome (str), start (int), and rate (float). Only rows whose first
    column equals the provided `chromosome` are used; all others are skipped.
    
    Lines starting with '#' are treated as comments and automatically skipped.
    
    When the recombination map doesn't cover the full chromosome range:
    - Missing start coverage: Uses the first available rate from the map
    - Missing end coverage: Uses the last available rate from the map
    - Warnings are printed when coverage is incomplete
    
    Parameters:
        rec_map (str): Path to the recombination map file.
        calc_start (int): Starting position of the chromosome.
        calc_end (int): Ending position of the chromosome.
        chunk_size (int): Size of each chunk in base pairs.
        chromosome (str): Chromosome name to filter rows by.
    
    Returns:
        np.ndarray: Average recombination rate for each chunk.
    """
    rec_map_data = []
    with open(rec_map, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Skip header lines (starting with #)
            if row and row[0].startswith('#'):
                continue
            # must have exactly [str, int, float]
            if len(row) != 3:
                continue
            # only use rows matching this chromosome
            if row[0] != chromosome:
                continue
            try:
                start_val = int(row[1])
                rate_val  = float(row[2])
            except ValueError:
                continue
            rec_map_data.append({'start': start_val, 'rate': rate_val})
    
    rec_map_data.sort(key=lambda x: x['start'])

    intervals = []
    if rec_map_data and rec_map_data[0]['start'] > calc_start:
        # Use the first available rate for missing start coverage
        first_rate = rec_map_data[0]['rate']
        intervals.append({'start': calc_start, 'end': rec_map_data[0]['start'], 'rate': first_rate})
        print(f"WARNING: Recombination map for {chromosome} doesn't cover start of chromosome (positions {calc_start}-{rec_map_data[0]['start']-1}). Using first available rate ({first_rate}). Consider extending your map if this affects your analysis.")
    
    for i, entry in enumerate(rec_map_data):
        interval_start = entry['start']
        interval_end = rec_map_data[i+1]['start'] if i < len(rec_map_data)-1 else calc_end
        if interval_end > calc_start and interval_start < calc_end:
            intervals.append({
                'start': max(interval_start, calc_start),
                'end':   min(interval_end,   calc_end),
                'rate':  entry['rate']
            })
    
    if intervals:
        last_end = intervals[-1]['end']
        if last_end < calc_end:
            # Use the last available rate for missing end coverage
            last_rate = intervals[-1]['rate']
            intervals.append({'start': last_end, 'end': calc_end, 'rate': last_rate})
            # Only warn if gap is >1Mb
            gap_size = calc_end - last_end
            if gap_size > 1000000:
                print(f"WARNING: Recombination map for {chromosome} doesn't cover end of chromosome (positions {last_end}-{calc_end}, gap: {gap_size:,}bp). Using last available rate ({last_rate}). Consider extending your map if this affects your analysis.")
    else:
        # No recombination data found for this chromosome
        print(f"WARNING: No recombination map data found for chromosome {chromosome}. Using default rate (1.0) for entire chromosome.")
        intervals.append({'start': calc_start, 'end': calc_end, 'rate': 1.0})
    
    rec_rates = []
    num_chunks = (calc_end - calc_start + chunk_size - 1) // chunk_size
    for chunk in range(num_chunks):
        start_chunk = calc_start + chunk * chunk_size
        end_chunk   = min(calc_end, start_chunk + chunk_size)
        chunk_len   = end_chunk - start_chunk
        weighted_sum = 0.0
        
        for iv in intervals:
            o_start = max(start_chunk, iv['start'])
            o_end   = min(end_chunk,   iv['end'])
            if o_start < o_end:
                weighted_sum += iv['rate'] * (o_end - o_start)
        
        rec_rates.append(weighted_sum / chunk_len if chunk_len > 0 else 1.0)
    
    return np.array(rec_rates)
