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
    # Get the directory where this module is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up to Bvalcalc root, then to data directory
    package_root = os.path.dirname(os.path.dirname(current_dir))
    return os.path.join(package_root, 'data')


def download_sample_data(force=False, quiet=False, target_dir=None):
    """
    Download sample data to user-accessible location.
    
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
        print(f"Downloading sample data to: {target_dir}")
    
    try:
        existing_files = []
        copied_files = []
        
        # Copy all files from source to target
        for item in os.listdir(source_dir):
            source_path = os.path.join(source_dir, item)
            target_path = os.path.join(target_dir, item)
            
            if os.path.isfile(source_path):
                if os.path.exists(target_path):
                    existing_files.append(item)
                    continue
                
                shutil.copy2(source_path, target_path)
                copied_files.append(item)
                    
            elif os.path.isdir(source_path):
                if os.path.exists(target_path):
                    existing_files.append(f"{item}/")
                    continue
                
                shutil.copytree(source_path, target_path, dirs_exist_ok=True)
                copied_files.append(f"{item}/")
        
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
            print(f"\nSample data successfully downloaded!")
            print(f"Available files:")
            for item in sorted(os.listdir(target_dir)):
                item_path = os.path.join(target_dir, item)
                if os.path.isfile(item_path):
                    size = os.path.getsize(item_path)
                    print(f"  {item} ({size:,} bytes)")
                else:
                    print(f"  {item}/")
            print(f"\nYou can now use these files with Bvalcalc commands.")
            print(f"Example: bvalcalc calculate -v {target_dir}/sample.vcf -p {target_dir}/sample_params.py")
        
        return True
        
    except Exception as e:
        if not quiet:
            print(f"Error downloading sample data: {e}")
        return False


def list_sample_data(quiet=False):
    """
    List available sample data files.
    
    Args:
        quiet (bool): If True, suppress output messages
    
    Returns:
        bool: True if successful, False otherwise
    """
    data_dir = get_sample_data_dir()
    
    if not os.path.exists(data_dir):
        if not quiet:
            print("Sample data not found. Run 'bvalcalc download-sample-data' to download sample files.")
        return False
    
    if not quiet:
        print(f"Sample data available in: {data_dir}")
        print("Available files:")
        
        for item in sorted(os.listdir(data_dir)):
            item_path = os.path.join(data_dir, item)
            if os.path.isfile(item_path):
                size = os.path.getsize(item_path)
                print(f"  {item} ({size:,} bytes)")
            else:
                print(f"  {item}/")
    
    return True
