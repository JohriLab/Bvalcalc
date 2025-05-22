import csv
import numpy as np

def load_rec_map(rec_map, calc_start, calc_end, chunk_size):
    """
    Processes the recombination map CSV file using the csv module and returns
    the average recombination rate per chunk, weighted by the proportion of each 
    chunk that falls within a given recombination interval.
    
    The CSV file must have two columns with headers 'start' and 'rate'. 
    Each row defines the recombination rate for positions starting at the given 
    'start' until the next 'start' (or until calc_end for the last entry). If any 
    region in the chromosome is not covered by the map, a default rate of 1.0 is used.
    
    Parameters:
        rec_map (str): Path to the recombination map CSV file.
        calc_start (int): Starting position of the chromosome.
        calc_end (int): Ending position of the chromosome.
        chunk_size (int): Size of each chunk in base pairs.
    
    Returns:
        list: Average recombination rate for each chunk.
        
    Raises:
        ValueError: If the CSV file does not have 'start' and 'rate' as headers,
                    or if any row has an invalid value for these columns.
    """
    # Read and validate CSV data
    rec_map_data = []
    with open(rec_map, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
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
            rec_map_data.append({'start': start_val, 'rate': rate_val})
    
    rec_map_data.sort(key=lambda x: x['start']) # Ensure the data is sorted by start position.
    intervals = [] # Build a list of intervals covering the entire region [calc_start, calc_end)
    if rec_map_data and rec_map_data[0]['start'] > calc_start: # If the first map entry starts after calc_start, assign default rate 1 from calc_start up to that entry.
        intervals.append({'start': calc_start, 'end': rec_map_data[0]['start'], 'rate': 1.0})
    
    for i, entry in enumerate(rec_map_data): # Create intervals for each recombination map entry.
        interval_start = entry['start']
        if i < len(rec_map_data) - 1: # For non-last entries, the interval ends at the next entry's start.
            interval_end = rec_map_data[i+1]['start']
        else: # For the last entry, extend the interval to calc_end.
            interval_end = calc_end     
        if interval_end > calc_start and interval_start < calc_end: # Only add intervals that overlap the region of interest.
            intervals.append({
                'start': max(interval_start, calc_start),
                'end': min(interval_end, calc_end),
                'rate': entry['rate']})
    
    if intervals: # If the last interval doesn't reach calc_end, fill in with default rate 1.
        last_end = intervals[-1]['end']
        if last_end < calc_end:
            intervals.append({'start': last_end, 'end': calc_end, 'rate': 1.0})
    else: # If no intervals were added, cover the entire region with default rate 1.
        intervals.append({'start': calc_start, 'end': calc_end, 'rate': 1.0}) 
    
    rec_rates = []
    for chunk in range((calc_end - calc_start + chunk_size - 1) // chunk_size): # For each chunk, compute the weighted average rate.
        start_chunk = calc_start + chunk * chunk_size
        end_chunk = min(calc_end, start_chunk + chunk_size)
        chunk_length = end_chunk - start_chunk
        weighted_sum = 0.0
        for interval in intervals: # Calculate overlap between this chunk and each interval.
            overlap_start = max(start_chunk, interval['start'])
            overlap_end = min(end_chunk, interval['end'])
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                weighted_sum += interval['rate'] * overlap_length
        
        avg_rate = weighted_sum / chunk_length if chunk_length > 0 else 1.0 # The average is the weighted sum divided by the chunk length.
        rec_rates.append(avg_rate)
    
    return np.array(rec_rates) # Return rec rate per chunk
