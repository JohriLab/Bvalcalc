import csv
import numpy as np

def recmapHandler(rec_map, chr_start, chr_end, chunk_size):
    """
    Processes the recombination map CSV file using the csv module and returns
    the average recombination rate per chunk, weighted by the proportion of each 
    chunk that falls within a given recombination interval.
    
    The CSV file must have two columns with headers 'start' and 'rate'. 
    Each row defines the recombination rate for positions starting at the given 
    'start' until the next 'start' (or until chr_end for the last entry). If any 
    region in the chromosome is not covered by the map, a default rate of 1.0 is used.
    
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
    
    # Build a list of intervals covering the entire region [chr_start, chr_end)
    intervals = []
    
    # If the first map entry starts after chr_start, assign default rate 1 from chr_start up to that entry.
    if rec_map_data and rec_map_data[0]['start'] > chr_start:
        intervals.append({'start': chr_start, 'end': rec_map_data[0]['start'], 'rate': 1.0})
    
    # Create intervals for each recombination map entry.
    for i, entry in enumerate(rec_map_data):
        interval_start = entry['start']
        # For non-last entries, the interval ends at the next entry's start.
        if i < len(rec_map_data) - 1:
            interval_end = rec_map_data[i+1]['start']
        else:
            # For the last entry, extend the interval to chr_end.
            interval_end = chr_end
        # Only add intervals that overlap the region of interest.
        if interval_end > chr_start and interval_start < chr_end:
            intervals.append({
                'start': max(interval_start, chr_start),
                'end': min(interval_end, chr_end),
                'rate': entry['rate']
            })
    
    # If the last interval doesn't reach chr_end, fill in with default rate 1.
    if intervals:
        last_end = intervals[-1]['end']
        if last_end < chr_end:
            intervals.append({'start': last_end, 'end': chr_end, 'rate': 1.0})
    else:
        # If no intervals were added, cover the entire region with default rate 1.
        intervals.append({'start': chr_start, 'end': chr_end, 'rate': 1.0})
    
    # Calculate the number of chunks over the region.
    num_chunks = (chr_end - chr_start + chunk_size - 1) // chunk_size
    rec_rates = []
    
    # For each chunk, compute the weighted average rate.
    for chunk in range(num_chunks):
        start_chunk = chr_start + chunk * chunk_size
        end_chunk = min(chr_end, start_chunk + chunk_size)
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


# def calcRLengths(blockstart, blockend, rec_rate_per_chunk, chr_start, chr_end, chunk_size, chunk_num):
#     """
#     Calculates the weighted lengths of each conserved block (gene), so that for example if the mean 
#     recombination rate across the block is 0.5, this will return the length of the block multiplied by 0.5
#     """

#     # print("In calcRLengths: ", chr_end, chr_start, chunk_size)
#     # Number of chunks covering the chromosome
#     num_chunks = (chr_end - chr_start) // chunk_size
#     chunk_starts = chr_start + np.arange(0, num_chunks + 1) * chunk_size
    
#     # Determine chunk indices for blockstart and blockend
#     blockstart_chunks = (blockstart - chr_start) // chunk_size
#     blockend_chunks = (blockend - chr_start) // chunk_size
#     block_chunk_overlaps = []
#     block_chunk_lengths = []
#     rec_rate_weighted_sums = []

#     for block_idx in range(len(blockstart)):
#         chunks = np.arange(blockstart_chunks[block_idx],  blockend_chunks[block_idx] + 1)
#         block_chunk_overlaps.append(chunks)
#         chunk_lengths = np.minimum(blockend[block_idx], chunk_starts[chunks] + chunk_size) - np.maximum(blockstart[block_idx], chunk_starts[chunks])
#         block_chunk_lengths.append(chunk_lengths)
#         weighted_sum = np.sum(chunk_lengths * rec_rate_per_chunk[chunks])
#         rec_rate_weighted_sums.append(weighted_sum)
    
#     # print("Weighted recombinant length (rate * length) for each block using map:", rec_rate_weighted_sums)
    
#     return rec_rate_weighted_sums


def calcRLengths(blockstart, blockend, rec_rate_per_chunk, chr_start, chr_end, chunk_size, chunk_num):
    """
    Calculates the weighted lengths of each conserved block (gene), so that for example if the mean 
    recombination rate across the block is 0.5, this will return the length of the block multiplied by 0.5
    """
    num_chunks = (chr_end - chr_start) // chunk_size
    # Build chunk boundaries (note: length = num_chunks + 1)
    chunk_starts = chr_start + np.arange(0, num_chunks + 1) * chunk_size
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

    distance_chunk_start = pos_chunk_clean - chunk_start

    num_chunks = (precise_region_end - precise_region_start) // chunk_size
    chunk_starts = precise_region_start + np.arange(0, num_chunks + 1) * chunk_size
    chunk_ends = np.minimum(chunk_starts + chunk_size, precise_region_end)
    this_chunk_idx = np.where(chunk_starts == chunk_start)[0][0] # The ID of this chunk in the chunk_starts array, e.g. if precise_chunks = 3, this will be [3] for chunk_num > 2
    chunk_end = chunk_ends[this_chunk_idx] # Focal chunk's end
    blockstart_chunks = (precise_blockstart - precise_region_start) // chunk_size
    blockend_chunks = (precise_blockend - precise_region_start) // chunk_size
    blockend_rec_distances = []
    blockstart_rec_distances = []

## FOR LOOP COULD BE IMPORVED BY USING NP ARRAY OPERATIONS INSTEAD!!!

    for block_idx in range(len(blockend_chunks)): # For blockends (i.e. upstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(pos_chunk_clean - precise_blockend[block_idx], pos_chunk_clean - chunk_start) # To blockend if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockend_chunks[block_idx] < this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = chunk_ends[blockend_chunks[block_idx]] - precise_blockend[block_idx]  # Distances from the end of the block to the end of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockend_chunks[block_idx]] # Rec_distance to end of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(blockend_chunks[block_idx] + 1, this_chunk_idx)
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))
                    
        blockend_rec_distances.append(total_rec_distances)

    for block_idx in range(len(blockstart_chunks)): # For blockstarts (i.e. downstream genes)
        rec_distance_overlapped, rec_distance_blockchunk = 0, 0 # Set to 0 to allow for sums even when not relevant
        inchunk_distances = np.minimum(precise_blockstart[block_idx] - pos_chunk_clean, chunk_end - pos_chunk_clean) # To blockstart if block is within same chunk, else to chunk start
        rec_distance_focalchunk = inchunk_distances * precise_rates[this_chunk_idx]
        isin_diffchunk = (blockstart_chunks[block_idx] > this_chunk_idx) # 0 if in same chunk, 1 if in different chunk

        chunk_edge_distances = precise_blockstart[block_idx] - chunk_starts[blockstart_chunks[block_idx]]   # Distances from the start of the block to the start of its chunk
        rec_distance_blockchunk = chunk_edge_distances * precise_rates[blockstart_chunks[block_idx]] # Rec_distance to start of block's chunk NOT spanned chunks

        overlapped_chunks = np.arange(this_chunk_idx + 1, blockstart_chunks[block_idx])
        rec_distance_overlapped = np.sum(precise_rates[overlapped_chunks] * chunk_size) # Rec_distance in chunks that are overlapped

        total_rec_distances = np.array(rec_distance_focalchunk + isin_diffchunk * (rec_distance_blockchunk + rec_distance_overlapped))

        blockstart_rec_distances.append(total_rec_distances)

            # print("Here", chunk_starts, chunk_ends)#chunk_starts, this_chunk_idx)#, rec_distance_overlapped)



    # return all_rec_distances