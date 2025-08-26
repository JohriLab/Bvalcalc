"""
Sample data utilities for Bvalcalc.
"""
import os
import shutil
import sys
from pathlib import Path


def get_sample_data_dir():
    """
    Get the directory where sample data should be installed.
    
    Returns:
        str: Path to the sample data directory
    """
    # Use current working directory as default
    return os.getcwd()


def get_package_data_dir():
    """
    Get the directory where sample data is stored in the package.
    
    Returns:
        str: Path to the package data directory
    """
    try:
        # Try to use importlib.resources for modern Python packaging
        import importlib.resources as resources
        with resources.path('Bvalcalc', 'data') as data_path:
            return str(data_path)
    except (ImportError, FileNotFoundError):
        # Fallback for older Python versions or development environment
        # Get the directory where this module is located
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Go up to Bvalcalc root, then to data directory
        package_root = os.path.dirname(os.path.dirname(current_dir))
        return os.path.join(package_root, 'data')


def download_sample_data(force=False, quiet=False, target_dir=None):
    """
    Download Drosophila CDS sample data to user-accessible location.
    
    Args:
        force (bool): If True, overwrite existing data
        quiet (bool): If True, suppress output messages
        target_dir (str): Target directory (defaults to current working directory)
    
    Returns:
        bool: True if successful, False otherwise
    """
    source_dir = get_package_data_dir()
    if target_dir is None:
        target_dir = get_sample_data_dir()
    
    if not os.path.exists(source_dir):
        if not quiet:
            print(f"Error: Sample data source directory not found: {source_dir}")
        return False
    
    # Create target directory if it doesn't exist
    os.makedirs(target_dir, exist_ok=True)
    
    if not quiet:
        print(f"Downloading Drosophila CDS sample data to: {target_dir}")
    
    # Define which files to copy (Drosophila CDS only)
    files_to_copy = [
        'cds_noX.bed',           # CDS annotations
        'dmel_comeron_recmap.csv', # Drosophila recombination map
        'README.md'              # Documentation
    ]
    
    try:
        existing_files = []
        copied_files = []
        
        # Copy only the specified files
        for filename in files_to_copy:
            source_path = os.path.join(source_dir, filename)
            target_path = os.path.join(target_dir, filename)
            
            if not os.path.exists(source_path):
                if not quiet:
                    print(f"Warning: Source file not found: {filename}")
                continue
            
            if os.path.exists(target_path) and not force:
                existing_files.append(filename)
                continue
            
            shutil.copy2(source_path, target_path)
            copied_files.append(filename)
        
        # Show warnings for existing files
        if existing_files and not quiet:
            print(f"Warning: The following files already exist and were not overwritten:")
            for item in existing_files:
                print(f"  {item}")
            print()
        
        # Show copied files
        if copied_files and not quiet:
            for item in copied_files:
                print(f"  Copied {item}")
        
        if not quiet:
            print(f"\nDrosophila CDS sample data successfully downloaded!")
            print(f"Available files:")
            for item in sorted(os.listdir(target_dir)):
                item_path = os.path.join(target_dir, item)
                if os.path.isfile(item_path):
                    size = os.path.getsize(item_path)
                    print(f"  {item} ({size:,} bytes)")
                else:
                    print(f"  {item}/")
            print(f"\nYou can now use these files with Bvalcalc commands.")
            print(f"Example: bvalcalc --genome --params your_params.py --bedgff cds_noX.bed --rec_map dmel_comeron_recmap.csv")
        
        return True
        
    except Exception as e:
        if not quiet:
            print(f"Error downloading sample data: {e}")
        return False


def list_sample_data(quiet=False):
    """
    List available Drosophila CDS sample data files.
    
    Args:
        quiet (bool): If True, suppress output messages
    
    Returns:
        bool: True if successful, False otherwise
    """
    data_dir = get_sample_data_dir()
    
    if not os.path.exists(data_dir):
        if not quiet:
            print("Drosophila CDS sample data not found. Run 'bvalcalc --download_sample_data' to download sample files.")
        return False
    
    if not quiet:
        print(f"Drosophila CDS sample data available in: {data_dir}")
        print("Available files:")
        
        for item in sorted(os.listdir(data_dir)):
            item_path = os.path.join(data_dir, item)
            if os.path.isfile(item_path):
                size = os.path.getsize(item_path)
                print(f"  {item} ({size:,} bytes)")
            else:
                print(f"  {item}/")
    
    return True
