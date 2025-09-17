import csv
from .header_utils import parse_headers, extract_header_info

def load_chr_sizes(file_path):
    chr_size_dict = {}
    valid_rows = 0

    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)

        for row in reader:
            # Skip header lines (starting with #)
            if row and row[0].startswith('#'):
                continue
                
            if len(row) < 2:
                continue

            chr_name = row[0].strip()
            size_str = row[1].strip()

            try:
                chr_size = int(size_str)
            except ValueError:
                continue

            chr_size_dict[chr_name] = chr_size
            valid_rows += 1

    if valid_rows == 0:
        raise ValueError(f"No valid chromosome size entries found in '{file_path}'. Expecting format: string,int per line.")

    return chr_size_dict


def get_chr_sizes_headers(file_path):
    """
    Get header lines from a chromosome sizes file.
    
    Args:
        file_path: Path to the chromosome sizes file
        
    Returns:
        List of header lines (with # prefix)
    """
    header_lines, _ = parse_headers(file_path)
    return header_lines


def get_chr_sizes_header_info(file_path):
    """
    Get structured header information from a chromosome sizes file.
    
    Args:
        file_path: Path to the chromosome sizes file
        
    Returns:
        HeaderInfo object with extracted information
    """
    header_lines = get_chr_sizes_headers(file_path)
    return extract_header_info(header_lines)
