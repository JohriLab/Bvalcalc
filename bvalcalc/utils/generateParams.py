import os
import sys

# keep your list of valid species here
SPECIES = ['selfing', 'human', 'drosophila', 'arabidopsis', 'mouse', 'pfalciparum', 'celegans']

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
    species_cap = species.capitalize()
    tpl_path = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '..', 'templates',
            f'{species_cap}Params.py'
        )
    )
    if not os.path.isfile(tpl_path):
        raise FileNotFoundError(f"Template for '{species}' not found at {tpl_path}")

    # Read the template
    with open(tpl_path, 'r') as tpl_file:
        content = tpl_file.read()
    print(f"Loaded template from:   {tpl_path}")

    # Ensure target folder exists
    os.makedirs(folder, exist_ok=True)

    # Write into <folder>/<Species>Params.py
    dest_path = os.path.join(folder, f'{species_cap}Params.py')
    with open(dest_path, 'w') as out_file:
        out_file.write(content)

    print(f"Wrote parameters to:     {dest_path}")
    print("Note that these are example parameters, please tailor to your population and analysis")