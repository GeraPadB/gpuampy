from gpuampy.io_tools import Crit

# === GPUAMPy Crit ===

# === Get the critical points === 
gpuam_data = Crit("FromPDB_moleculeCrit.log")
gpuam_data.read_data()

# === Get the bond critical points === 
bcp_data = gpuam_data.get_bcp()

# === Get the noncovalent critical points ===
bcp_data = gpuam_data.get_noncovalent_bcp()
bcp_data.save_csv("crit_df.csv")

