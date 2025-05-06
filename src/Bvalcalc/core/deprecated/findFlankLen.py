from core.calculateB import calculateB_linear

def findFlankLen(target_B, length_of_element, max_distance=1e6):
    """
    Find the minimum distance where B > target_B for a given element length,
    or return max_distance if no such distance exists.
    """
    low, high = 0, max_distance
    while low < high:
        mid = (low + high) // 2
        B = calculateB_linear(mid, length_of_element)
        if B > target_B:
            high = mid  # Narrow the search to smaller distances
        else:
            low = mid + 1  # Narrow the search to larger distances

    # After exiting the loop, check if the condition is met
    final_B = calculateB_linear(low, length_of_element)
    if final_B > target_B:
        return int(low)
    else:
        return int(max_distance)  # Return max_distance if no valid distance is found
    