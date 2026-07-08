from gpuampy.io_tools import Crit
from gpuampy.bio_tools import Biomol

# === GPUAMPy Crit ===

# === Get the critical points === 
gpuam_data = Crit("PDB_moleculeCrit.log")
gpuam_data.read_data()

# === Get the bond critical points === 
bcp_data = gpuam_data.get_bcp()

# === Get the noncovalent critical points ===
bcp_data = gpuam_data.get_interaction_bcp(interaction="all")
#print(bcp_data)
bcp_data.save_csv("crit_df.csv")



# === GPUAMPy Bio ===

# === Parse PDB file ===
bio_data = Biomol("PDB_molecule.pdb")

# === Add PDB info ===

bcp_data = bio_data.biomol_info(bcp_data, interaction_cat="intermolecular")
#print(bcp_data)
bcp_data.save_csv("bio_df.csv")

