import os
import sys

# keep your list of valid species here - updated for new naming convention
SPECIES = [
    'selfing',  # Special case - no element type
    'aratha_cds', 'aratha_phastcons',  # Arabidopsis thaliana
    'dromel_cds', 'dromel_utr', 'dromel_phastcons',  # Drosophila melanogaster
    'homsap_cds', 'homsap_phastcons',  # Homo sapiens
]

def check_generate_params_args(argv=None):
    """
    If '--generate_params' appears with no species or next flag,
    print a concise error and exit.
    """
    if argv is None:
        argv = sys.argv
    if '--generate_params' in argv:
        idx = argv.index('--generate_params')
        if idx == len(argv) - 1 or argv[idx+1].startswith('-'):
            print(f"Provide name of default template as an argument: {' '.join(SPECIES)}")
            sys.exit(1)

def generateParams(species, folder='.'):
    # Convert to lowercase for case-insensitive matching
    species_lower = species.lower()
    
    # Handle special case for selfing
    if species_lower == 'selfing':
        tpl_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '..', 'templates',
                'SelfingParams.py'
            )
        )
        dest_name = 'SelfingParams.py'
    else:
        # Handle new naming convention: GenSpe_ElementType_Params.py
        # Map from command format to file format
        species_mapping = {
            'aratha_cds': 'AraTha_Cds_Params.py',
            'aratha_phastcons': 'AraTha_Phastcons_Params.py',
            'dromel_cds': 'DroMel_Cds_Params.py',
            'dromel_utr': 'DroMel_Utr_Params.py',
            'dromel_phastcons': 'DroMel_Phastcons_Params.py',
            'homsap_cds': 'HomSap_Cds_Params.py',
            'homsap_phastcons': 'HomSap_Phastcons_Params.py',
        }
        
        if species_lower not in species_mapping:
            print(f"Error: Unknown species '{species}'.")
            print(f"Available options: {', '.join(SPECIES)}")
            print("Note: Species names are case-insensitive (e.g., 'ARATHA_CDS' works the same as 'aratha_cds')")
            sys.exit(1)
        
        filename = species_mapping[species_lower]
        tpl_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                '..', 'templates',
                filename
            )
        )
        dest_name = filename
    
    if not os.path.isfile(tpl_path):
        raise FileNotFoundError(f"Template for '{species}' not found at {tpl_path}")

    # Read the template
    with open(tpl_path, 'r') as tpl_file:
        content = tpl_file.read()
    print(f"Loaded template from:   {tpl_path}")

    # Ensure target folder exists
    os.makedirs(folder, exist_ok=True)

    # Write into <folder>/<dest_name>
    dest_path = os.path.join(folder, dest_name)
    with open(dest_path, 'w') as out_file:
        out_file.write(content)

    print(f"Wrote parameters to:     {dest_path}")
    print("Note that these are example parameters, please tailor to your population and analysis")