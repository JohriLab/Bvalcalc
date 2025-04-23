import os
import shutil

def generateParams(species):
    species_cap = species.capitalize()
    tpl_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', 'templates', f'{species_cap}Params.py')
    )
    if not os.path.isfile(tpl_path):
        raise FileNotFoundError(f"Template for '{species}' not found at {tpl_path}")

    with open(tpl_path, 'r') as tpl_file:
        content = tpl_file.read()

    example_dir = os.path.join(os.getcwd(), "ExampleParams")
    os.makedirs(example_dir, exist_ok=True)

    example_path = os.path.join(example_dir, f"{species_cap}Params.py")
    with open(example_path, 'w') as out_file:
        out_file.write(content)

    dest_path = os.path.join(os.getcwd(), "newParams.py")
    shutil.copyfile(example_path, dest_path)

    print(f"✔ Loaded template from:   {tpl_path}")
    print(f"✔ Wrote ExampleParams to: {example_path}")
    print(f"✔ Duplicated as:          {dest_path}")
