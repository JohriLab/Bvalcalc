import numpy as np
from bisect import bisect_left, bisect_right

chunks = np.array([[1000, 2000], [2000, 3000]])  # Define chunk ranges: [start, end]

blockstart = np.array([500, 1500, 2500])  # Example blockstart positions
blockend = np.array([1200, 2800, 2700])   # Example blockend positions
flank_blockstart = blockstart - np.array([50, 100, 150])  # Example flank_blockstart adjustments
flank_blockend = blockend + np.array([50, 100, 150])      # Example flank_blockend adjustments

chunks = np.array([(10, 20), (50, 60)])
flank_blockstart = np.array([5, 15, 25, 35])
flank_blockend = np.array([60, 55, 45, 40])

def precompute_block_overlaps_optimized(chunks, flank_blockstart, flank_blockend):
    # Sort blocks by start positions, keeping end positions and indices tied
    sorted_blocks = sorted(zip(flank_blockstart, flank_blockend, range(len(flank_blockstart))))
    sorted_start = [b[0] for b in sorted_blocks]  # Sorted block starts
    sorted_end = [b[1] for b in sorted_blocks]    # Sorted block ends
    sorted_indices = [b[2] for b in sorted_blocks]  # Original indices

    # Create an empty list to store block indices for each chunk
    chunk_overlaps = []

    for chunk_start, chunk_end in chunks:
        # Find candidate blocks that might overlap the chunk
        start_idx = bisect_left(sorted_start, chunk_end)  # Blocks starting before chunk_end
        end_idx = bisect_right(sorted_end, chunk_start)   # Blocks ending after chunk_start
        
        # Ensure correct overlaps by explicitly checking bounds
        overlapping_blocks = [
            sorted_indices[i] for i in range(end_idx)
            if sorted_start[i] < chunk_end and sorted_end[i] > chunk_start
        ]
        chunk_overlaps.append(overlapping_blocks)

    return chunk_overlaps

chunk_overlaps = precompute_block_overlaps_optimized(chunks, flank_blockstart, flank_blockend)

# Output the relevant block indices for each chunk
# for i, overlap in enumerate(chunk_overlaps):
#     print(f"Chunk {i+1} overlaps with block indices: {overlap}")
##################################################################################################################################

def find_overlapping_ranges_for_chunks(chunks, ranges):
    overlaps = []
    
    # Iterate through each chunk and find overlapping ranges
    for chunk_start, chunk_end in chunks:
        chunk_overlaps = []
        for range_start, range_end in ranges:
            if range_start <= chunk_end and range_end >= chunk_start:
                chunk_overlaps.append((range_start, range_end))
        
        # Store the overlapping ranges for the current chunk
        overlaps.append(chunk_overlaps)
    
    return overlaps

# Test Case
chunks = [(100, 200), (250, 300)]
ranges = [(50, 150), (120, 180), (200, 250), (300, 350)]

overlapping_ranges = find_overlapping_ranges_for_chunks(chunks, ranges)
print("Overlapping ranges:", overlapping_ranges)

