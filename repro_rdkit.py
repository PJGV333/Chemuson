import sys
from rdkit import Chem
from chemuson.chemio.rdkit_io import smiles_to_molgraph, molgraph_to_rdkit_with_map

def repro():
    smiles = "O=C1C=CC(=O)C=C1"
    print(f"Testing SMILES: {smiles}")
    try:
        graph = smiles_to_molgraph(smiles)
        print("Graph built successfully.")
        print(f"Atoms: {len(graph.atoms)}")
        print(f"Bonds: {len(graph.bonds)}")
        
        # Check for self-loops
        for b in graph.bonds.values():
            if b.a1_id == b.a2_id:
                print(f"SELF-LOOP DETECTED in MolGraph: bond {b.id} ({b.a1_id} -> {b.a2_id})")

        # Try to convert back to RDKit
        print("Converting back to RDKit...")
        mol, id_map = molgraph_to_rdkit_with_map(graph)
        print("Conversion successful.")
        
    except Exception as e:
        print(f"CAUGHT EXCEPTION: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    repro()
