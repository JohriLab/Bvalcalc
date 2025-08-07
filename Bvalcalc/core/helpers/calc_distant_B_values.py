from Bvalcalc.core.calculateB import calculateB_linear, calculateB_recmap
import numpy as np

def calc_distant_B_values(interfering_L_shape
    #     calc_chunks,
    #     l_per_chunk,
    #     chunks_with_low_rec,
    # chunk_size,
    # blockstart,
    # blockend,
    # chr_start,
    # chr_size,
    ):



    print("gotcha")

    distant_B = 0.90 # NEED to actually calculate this, same was as in process_single_chunk

    B_from_outside_local_interference_regime = np.full(interfering_L_shape, distant_B, dtype=float)

    return B_from_outside_local_interference_regime
        #chunk_idx, chunk_size, chr_start, chr_size, num_chunks, 
         #           precise_chunks, lperchunk, rec_rate_per_chunk, gc_rate_per_chunk):
