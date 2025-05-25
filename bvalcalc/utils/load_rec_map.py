import csv
import numpy as np

def load_rec_map(rec_map, calc_start, calc_end, chunk_size, chromosome):
    """
    Processes the recombination map file (no header) and returns
    the average recombination rate per chunk, weighted by the proportion of each 
    chunk that falls within a given recombination interval.
    
    The file must be headerless, with each row having exactly three columns:
    chromosome (str), start (int), and rate (float). Only rows whose first
    column equals the provided `chromosome` are used; all others are skipped.
    
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
        intervals.append({'start': calc_start, 'end': rec_map_data[0]['start'], 'rate': 1.0})
    
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
            intervals.append({'start': last_end, 'end': calc_end, 'rate': 1.0})
    else:
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
