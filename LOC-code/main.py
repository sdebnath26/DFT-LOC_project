from utils.io_utils import open_xyz, open_log
from utils.mo_utils import get_top_atom_coeff
from utils.geometry_utils import distance, angle
from analysis.bond_analysis import bond_order1, valency_atom
from utils.print_utils import print_dict

def main():
    log_data = open_log("data/fb_pyrrole.log")
    xyz_data = open_xyz("data/pyrrole.xyz")

    print("Loaded log and xyz files.")
    print("Example bond order and valency checks:")

    # This is placeholder logic
    # Replace with actual function usage based on your flow
    # Example dummy atoms
    atom1, atom2 = "C1", "C2"
    print(f"Bond order for {atom1}-{atom2}: {bond_order1(atom1, atom2)}")
    print(f"Valency of {atom1}: {valency_atom(atom1)}")

if __name__ == "__main__":
    main()
