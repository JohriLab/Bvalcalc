import csv
import numpy as np

def recmapHandler(rec_map, calc_start, calc_end, chunk_size):
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
    
    # Ensure the data is sorted by start position.
    rec_map_data.sort(key=lambda x: x['start'])
    
    # Build a list of intervals covering the entire region [calc_start, calc_end)
    intervals = []
    
    # If the first map entry starts after calc_start, assign default rate 1 from calc_start up to that entry.
    if rec_map_data and rec_map_data[0]['start'] > calc_start:
        intervals.append({'start': calc_start, 'end': rec_map_data[0]['start'], 'rate': 1.0})
    
    # Create intervals for each recombination map entry.
    for i, entry in enumerate(rec_map_data):
        interval_start = entry['start']
        # For non-last entries, the interval ends at the next entry's start.
        if i < len(rec_map_data) - 1:
            interval_end = rec_map_data[i+1]['start']
        else:
            # For the last entry, extend the interval to calc_end.
            interval_end = calc_end
        # Only add intervals that overlap the region of interest.
        if interval_end > calc_start and interval_start < calc_end:
            intervals.append({
                'start': max(interval_start, calc_start),
                'end': min(interval_end, calc_end),
                'rate': entry['rate']
            })
    
    # If the last interval doesn't reach calc_end, fill in with default rate 1.
    if intervals:
        last_end = intervals[-1]['end']
        if last_end < calc_end:
            intervals.append({'start': last_end, 'end': calc_end, 'rate': 1.0})
    else:
        # If no intervals were added, cover the entire region with default rate 1.
        intervals.append({'start': calc_start, 'end': calc_end, 'rate': 1.0})
    
    # Calculate the number of chunks over the region.
    num_chunks = (calc_end - calc_start + chunk_size - 1) // chunk_size
    rec_rates = []
    
    # For each chunk, compute the weighted average rate.
    for chunk in range(num_chunks):
        start_chunk = calc_start + chunk * chunk_size
        end_chunk = min(calc_end, start_chunk + chunk_size)
        chunk_length = end_chunk - start_chunk
        
        weighted_sum = 0.0
        
        # Calculate overlap between this chunk and each interval.
        for interval in intervals:
            overlap_start = max(start_chunk, interval['start'])
            overlap_end = min(end_chunk, interval['end'])
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                weighted_sum += interval['rate'] * overlap_length
        
        # The average is the weighted sum divided by the chunk length.
        avg_rate = weighted_sum / chunk_length if chunk_length > 0 else 1.0
        rec_rates.append(avg_rate)

    print("rates", avg_rate)
    
    return np.array(rec_rates)

def calcRLengths(blockstart, blockend, rec_rate_per_chunk, calc_start, calc_end, chunk_size, chunk_num):
    """
    Calculates the weighted lengths of each conserved block (gene), so that for example if the mean 
    recombination rate across the block is 0.5, this will return the length of the block multiplied by 0.5
    """
    num_chunks = (calc_end - calc_start) // chunk_size
    # Build chunk boundaries (note: length = num_chunks + 1)
    chunk_starts = calc_start + np.arange(0, num_chunks + 1) * chunk_size
    chunk_left  = chunk_starts           # shape: (num_chunks+1,)
    chunk_right = chunk_starts + chunk_size  # shape: (num_chunks+1,)

    # Ensure block boundaries are numpy arrays (and 1D)
    blockstart = np.asarray(blockstart).flatten()
    blockend   = np.asarray(blockend).flatten()

    # Compute the overlap between each block and each chunk interval:
    # For block i and chunk j, the overlap is:
    #   max(0, min(blockend[i], chunk_right[j]) - max(blockstart[i], chunk_left[j]))
    overlap = np.maximum(0, np.minimum(blockend[:, None], chunk_right[None, :]) -
                           np.maximum(blockstart[:, None], chunk_left[None, :]))
    weighted_overlap = overlap * rec_rate_per_chunk[None, :] # Multiply by the recombination rate for each chunk
    weighted_sum = np.sum(weighted_overlap, axis=1) # Sum over the chunk intervals for each block
    
    return weighted_sum

def calcRDistances(precise_blockstart, precise_blockend, precise_rates, precise_region_start, precise_region_end, chunk_size, pos_chunk_clean, chunk_num, chunk_start):

    num_chunks = (precise_region_end - precise_region_start) // chunk_size
    chunk_starts = precise_region_start + np.arange(0, num_chunks + 1) * chunk_size
    chunk_ends = np.minimum(chunk_starts + chunk_size, precise_region_end)
    this_chunk_idx = np.where(chunk_starts == chunk_start)[0][0] # The ID of this chunk in the chunk_starts array, e.g. if precise_chunks = 3, this will be [3] for chunk_num > 2
    chunk_end = chunk_ends[this_chunk_idx] # Focal chunk's end
    blockstart_chunks = (precise_blockstart - precise_region_start) // chunk_size
    blockend_chunks = (precise_blockend - precise_region_start) // chunk_size
    # blockend_rec_distances = np.empty(len(blockend_chunks))
    blockend_rec_distances = np.empty((len(blockend_chunks), len(pos_chunk_clean)))
    blockstart_rec_distances = np.empty((len(blockstart_chunks), len(pos_chunk_clean)))
    # blockstart_rec_distances = []

## FOR LOOP COULD BE IMPORVED BY USING NP ARRAY OPERATIONS INSTEAD!!!

    for block_idx in range(len(blockend_chunks)): # For blockends (i.e. upstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(pos_chunk_clean - precise_blockend[block_idx], 
                                       pos_chunk_clean - chunk_start) # To blockend if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockend_chunks[block_idx] < this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = chunk_ends[blockend_chunks[block_idx]] - precise_blockend[block_idx]  # Distances from the end of the block to the end of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockend_chunks[block_idx]] # Rec_distance to end of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(blockend_chunks[block_idx] + 1, this_chunk_idx)
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))
                    
        blockend_rec_distances[block_idx] = total_rec_distances

    for block_idx in range(len(blockstart_chunks)): # For blockstarts (i.e. downstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(precise_blockstart[block_idx] - pos_chunk_clean, 
                                       chunk_end - pos_chunk_clean) # To blockstart if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockstart_chunks[block_idx] > this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = precise_blockstart[block_idx] - chunk_starts[blockstart_chunks[block_idx]]   # Distances from the start of the block to the start of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockstart_chunks[block_idx]] # Rec_distance to start of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(this_chunk_idx + 1, blockstart_chunks[block_idx])
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))

        blockstart_rec_distances[block_idx] = total_rec_distances

    return blockend_rec_distances, blockstart_rec_distances

def calcRLengthsDistances_forchunks(upstream_indices, downstream_indices, rec_rate_per_chunk, relevant_upstream_psdc_lengths, relevant_downstream_psdc_lengths, chunk_index, chunk_size, relevant_upstream_pseudoblockends, relevant_downstream_pseudoblockstarts, chunk_starts, chunk_ends, chunk_rec_distances, num_chunks):

    ## Calculate relevant upstream and downstream rec lengths of pseudoblocks
    upstream_rec_rates = rec_rate_per_chunk[upstream_indices] # Relevant rec rates for pseudochunks upstream
    upstream_rec_lengths = upstream_rec_rates * relevant_upstream_psdc_lengths
    downstream_rec_rates = rec_rate_per_chunk[downstream_indices] # Relevant rec rates for pseudochunks downstream
    downstream_rec_lengths = downstream_rec_rates * relevant_downstream_psdc_lengths

    ## Calculate relevant upstream rec distances!
    ## CURRENTLY TAKING A LONG TIME TO CALCULATE 
    mean_rec_distance_focalchunk = rec_rate_per_chunk[chunk_index] * chunk_size / 2 - 0.5 # Note that this is distance to middle of focal chunk.

    if chunk_index == num_chunks -1: # If it's the final chunk (which may be not be full chunk_size length)
        end_focalchunk_distance = (chunk_ends[chunk_index] - chunk_starts[chunk_index])/2 # I think needs to be shifted by 1bp
        mean_rec_distance_focalchunk = end_focalchunk_distance

    upstream_distance_blockchunk = chunk_ends[upstream_indices] - relevant_upstream_pseudoblockends
    upstream_rec_distance_blockchunk = upstream_distance_blockchunk * rec_rate_per_chunk[upstream_indices] # This is rec distance from edge of pseudoblock to its chunk end

    upstream_overlapped_indices = [np.arange(u + 1, chunk_index) for u in upstream_indices]
    upstream_overlapped_rec_distances = np.array([chunk_rec_distances[idx].sum() for idx in upstream_overlapped_indices]) # This is rec distance spanned in fully overlapped chunks

    upstream_rec_distances = mean_rec_distance_focalchunk + upstream_rec_distance_blockchunk + upstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream

    ## Calculate downstream rec distances!
    ## CURRENTLY TAKING A LONG TIME TO CALCULATE 
    downstream_distance_blockchunk = relevant_downstream_pseudoblockstarts - chunk_starts[downstream_indices]
    downstream_rec_distance_blockchunk = downstream_distance_blockchunk * rec_rate_per_chunk[downstream_indices] # This is rec distance from edge of pseudoblock to its chunk start
    
    downstream_overlapped_indices = [np.arange(chunk_index + 1, d) for d in downstream_indices]
    downstream_overlapped_rec_distances = np.array([chunk_rec_distances[idx].sum() for idx in downstream_overlapped_indices]) # This is rec distance spanned in fully overlapped chunks

    downstream_rec_distances = mean_rec_distance_focalchunk + downstream_rec_distance_blockchunk + downstream_overlapped_rec_distances # Combined rec distance from middle of focal chunk to edge of pseudo"blocks" upstream


    return upstream_rec_lengths, downstream_rec_lengths, upstream_rec_distances, downstream_rec_distances
