#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools for working with GPUAM Crit Calculations.
"""
import pandas as pd
import sys

class GPUAMdf(pd.DataFrame):

    @property
    def _constructor(self):
        return GPUAMdf

    def save_csv(self, filename, encoding="utf-8-sig"):
        self.to_csv(
            filename,
            index=False,
            encoding=encoding
        )




class Crit:
    """
    Class for recovering and storing data from a GPUAM Crit calculation.

    """

    def __init__(self, filename):
        self.filename = filename
        self.dataframes = [None,None,None,None]
        self.last_common_line = 51
        self.line_to_skip = 0
        self.n_cp_lines = [26,34,26,26]
        self.n_cp = [0,0,0,0]
        self.nuclei = None
        self.nuclei_columns = ['id','type','coord','density','laplacian']
        self.cp_columns = [
                      ['id','coord','density','normgrad','laplacian',
                       'rho/lap','G','K','virial','lagran_rho',
                       'G/rho','elf','lol','eigenvalues','shannon_rho',
                       'diseq_rho','complex_rho','shannon_sigma',
                       'diseq_sigma','complex_sigma'],
                      ['id','atom1_id','atom1_type','atom2_id',
                       'atom2_type','d_attract','d_att1','d_att2',
                       'angle','coord','density','normgrad','laplacian',
                       'rho/lap','G','K','virial','lagran_rho',
                       'G/rho','elf','lol','eigenvalues','ellip',
                       'eta','shannon_rho','diseq_rho','complex_rho',
                       'shannon_sigma','diseq_sigma','complex_sigma']
                     ]
        
    def check_run(self,id_line):
        """
        Check what type of run was performed, changing accordingly some
        default values and adapting the follow up processing.
        """
        run_type = id_line.split()[-1]
        if run_type=="GPUAM":
            print("GPUAM run to be processed.")
            self.run_type = "GPUAM"
        elif run_type=="cube3d":
            print("Cube3D run to be processed.")
            self.run_type = "CUBE3D"
            self.last_common_line = 39
            self.n_cp_lines = [21,25,21,21]
            self.cp_columns = [
                      ['id','coord','density','normgrad','laplacian',
                       'G','K','virial','H',#Hessian mat should be added here
                       'eigenvalues'],
                      ['id','atom1_id','atom1_type','atom2_id',
                       'atom2_type','coord','density','normgrad','laplacian',
                       'G','K','virial','H', #Hessian mat should be added here
                       'eigenvalues','ellip',
                       'eta']
                     ]

    def get_basic_data(self,data_lines):
        if self.run_type == "GPUAM":
            try:
                self.n_primitives = int(data_lines[34].split(':')[1])
                self.n_orbitals = int(data_lines[35].split(':')[1])
                self.n_nuclei = int(data_lines[36].split(':')[1])
                self.n_electrons = int(data_lines[37].split(':')[1].split()[0])
                self.n_cp[0] = int(data_lines[38].split(':')[1])
                self.n_cp[1] = int(data_lines[39].split(':')[1])
                self.n_cp[2] = int(data_lines[40].split(':')[1])
                self.n_cp[3] = int(data_lines[41].split(':')[1])
                self.n_totalcp = int(data_lines[42].split(':')[1])
                self.poincare_hopf = int(data_lines[45].split('=')[1])
            except:
                print("MOPAC results or your specific type of run analysis are not available for now.")
                sys.exit(1)
        elif self.run_type == "CUBE3D":
            try:
                self.n_nuclei = int(data_lines[25].split(':')[1])
                self.n_cp[0] = int(data_lines[26].split(':')[1])
                self.n_cp[1] = int(data_lines[27].split(':')[1])
                self.n_cp[2] = int(data_lines[28].split(':')[1])
                self.n_cp[3] = int(data_lines[29].split(':')[1])
                self.n_totalcp = int(data_lines[30].split(':')[1])
                self.poincare_hopf = int(data_lines[33].split('=')[1])
            except:
                print("A format error ocurred, check your input file or contact the development team.")
                sys.exit(1)

    def set_initial_data(self):
        self.n_dif_cp = ((self.n_cp[0] != 0) + (self.n_cp[1] != 0) 
          + (self.n_cp[2] != 0) + (self.n_cp[3] != 0))
        # 2 lines for the R,B,C CP header. 1 separator
        # line for each kind of CP (n_dif_cp). 2 more
        # lines from the first CP header.
        self.line_to_skip = ( self.last_common_line + self.n_nuclei 
                    + self.n_totalcp + 2 + self.n_dif_cp + 2 )

    def set_nuclei_data(self,data_lines):
        n_lines = self.n_nuclei
        #data_lines[line_to_skip + 5 * n_bcp_lines]
        nuclei_data = [[] for i in range(self.n_nuclei)]
        for i_nuclei in range(n_lines):
            line_data = data_lines[self.last_common_line+i_nuclei].split()
            nuclei_data[i_nuclei].append(int(line_data[0]))
            nuclei_data[i_nuclei].append(line_data[1])
            nuclei_data[i_nuclei].append(list(map(float,line_data[2:5])))
            nuclei_data[i_nuclei].append(float(line_data[5]))
            nuclei_data[i_nuclei].append(float(line_data[6]))
        self.nuclei = pd.DataFrame(nuclei_data,columns=self.nuclei_columns)


    def read_block(self,id,data_lines):
        n_lines = self.n_cp_lines[id]
        #data_lines[line_to_skip + 5 * n_bcp_lines]
        cp_data = [[] for i in range(self.n_cp[id])]
        hessian = False
        for i_cp in range(self.n_cp[id]):
            for n_line in range(n_lines):
                if n_line == 0: continue
                current_line = self.line_to_skip + n_line + (i_cp * n_lines)
                line_data = data_lines[current_line].split()
                if not line_data: continue
                elif line_data[0] == 'Critical':
                    cp_data[i_cp].append(int(line_data[2]))
                elif line_data[0] == 'Between':
                    if self.run_type == "GPUAM":
                        cp_data[i_cp].append(int(line_data[3]))
                        cp_data[i_cp].append(line_data[4])
                        cp_data[i_cp].append(int(line_data[6]))
                        cp_data[i_cp].append(line_data[7])
                    elif self.run_type == "CUBE3D":
                        cp_data[i_cp].append(int(line_data[4]))
                        cp_data[i_cp].append(line_data[5])
                        cp_data[i_cp].append(int(line_data[7]))
                        cp_data[i_cp].append(line_data[8])
                elif line_data[0] == 'Coordinates':
                    cp_data[i_cp].append(list(map(float,line_data[-3:])))
                elif line_data[0] == 'Eigenvalues':
                    cp_data[i_cp].append(list(map(float,line_data[-3:])))
                elif line_data[0] == '|': # Here the Hessian mat should be stored in CUBE runs
                    #n_line += 3
                    hessian = True
                    continue
                elif line_data[0] == 'Hessian':
                    continue
                elif line_data[0] == '|' and hessian:
                    hesian = False
                    continue
                else:
                    cp_data[i_cp].append(float(line_data[-1]))
        
        if id == 1:
            cp_columns = self.cp_columns[1]
        else:
            cp_columns = self.cp_columns[0]
        
        self.dataframes[id] = pd.DataFrame(cp_data,columns=cp_columns)
        self.line_to_skip += self.n_cp[id] * (self.n_cp_lines[id]) + 2
         

    def read_data(self):
        try:
            inp = open(self.filename,'r')
            data = inp.read()
            inp.close()
            data_lines = data.strip().split('\n')
            pass
        except FileNotFoundError:
            print(f"File not found: {self.filename}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

        self.check_run(data_lines[1])
        self.get_basic_data(data_lines)
        self.set_initial_data()
        self.set_nuclei_data(data_lines)

        for cp_id in range(4):
            if self.n_cp[cp_id] == 0: continue
            else: self.read_block(cp_id,data_lines)

    def get_nuclei(self):
        return self.nuclei
    
    def xyz_nuclei(self):
        xyz_filename = self.filename[:-8]
        df_nuclei = self.get_nuclei()
        lines = [f" {row['type']:<2} {'':5} {row['coord'][0]:8.5f} {row['coord'][1]:8.5f} {row['coord'][2]:8.5f}\n"
                 for index, row in df_nuclei.iterrows()]
        with open(f'{xyz_filename}.xyz', 'w') as xyz:
            xyz.write(f'{self.n_nuclei}\n\n')
            xyz.writelines(lines)
        return "Nuclei XYZ created."

    def get_data(self):
        return self.dataframes

    def get_nna(self):
        return self.dataframes[0]

    def get_bcp(self):
        return self.dataframes[1]

    def get_rcp(self):
        return self.dataframes[2]

    def get_ccp(self):
        return self.dataframes[3]
        
    def get_bcp_properties(self):
        """
        Return the BCP DataFrame including additional QTAIM descriptors.
        """
        
        bcp = self.get_bcp()
        
        if bcp is None:
            return None
        
        bcp = bcp.copy()
        
        # -------- H total energy --------
        h = bcp["G"] + bcp["virial"]
        idx = bcp.columns.get_loc("virial") + 1
        bcp.insert(idx, "H", h)
        
        # -------- Ratios virial/rho, H/rho and |virial| /G --------
        idx = bcp.columns.get_loc("G/rho") + 1
        
        bcp.insert(
            idx,
            "virial/rho",
            bcp["virial"] / bcp["density"]
        )
        
        bcp.insert(
            idx + 1,
            "H/rho",
            bcp["H"] / bcp["density"]
        )
        
        bcp.insert(
            idx + 2,
            "virial/G",
            bcp["virial"].abs() / bcp["G"]
        )
        
        return bcp

    def _get_covalent_bcp(self, rho_max=0.10, vg_max=2.0):
        """
        Return only covalent bond critical points.
    
        Parameters
        ----------
        rho_max : float, optional
            Electron density threshold for covalent interactions.
        vg_max : float, optional
            |V|/G threshold for covalent interactions.
    
        Returns
        -------
        pandas.DataFrame
            DataFrame containing only covalent bond critical points.
        """
    
        bcp = self.get_bcp_properties()
    
        if bcp is None:
            return None
    
        covalent = (
            (bcp["density"] >= rho_max) |
            (bcp["virial/G"] >= vg_max)
        )
    
        return bcp[covalent].reset_index(drop=True)
        
        
    def h_donor(self):
        """
        Return the atom bonded to each hydrogen atom from covalent BCPs.
    
        Returns
        -------
        dict
            Dictionary where hydrogen atom IDs are keys and values are
            tuples containing the bonded atom ID and atom type.
    
            Example:
            {
                64: (63, "C"),
                21: (20, "N")
            }
        """
    
        bonds = self._get_covalent_bcp()
    
        if bonds is None:
            return None
    
        h_donors = {}
    
        for _, row in bonds.iterrows():
    
            atom1_id = row["atom1_id"]
            atom1_type = row["atom1_type"]
    
            atom2_id = row["atom2_id"]
            atom2_type = row["atom2_type"]
    
            if atom1_type == "H":
                h_donors[atom1_id] = (
                    atom2_id,
                    atom2_type
                )
    
            elif atom2_type == "H":
                h_donors[atom2_id] = (
                    atom1_id,
                    atom1_type
                )
    
        return h_donors
        
        
    def _get_interaction_type(
        self,
        atom1_type,
        atom2_type,
        atom1_id,
        atom2_id,
        covalent,
        h_donors
    ):
        """
        Classify the chemical type of a BCP interaction.
    
        Returns
        -------
        str
            Contact classification.
        """
    
        # Covalent interaction
        if covalent:
            return "Covalent"
    
        # No hydrogen involved: heavy atom interaction
        if atom1_type != "H" and atom2_type != "H":
            return "Lewis int"
    
        # H...H interactions
        if atom1_type == "H" and atom2_type == "H":
    
            if atom1_id in h_donors and atom2_id in h_donors:
    
                _, donor1 = h_donors[atom1_id]
                _, donor2 = h_donors[atom2_id]
    
                if donor1 == donor2:
                    return "Dihydrogen bond"
    
                else:
                    return "H-H bond"
    
            return "H-H interaction"
    
        # H...X interactions
        if atom1_type == "H" and atom1_id in h_donors:
    
            _, donor_type = h_donors[atom1_id]
    
            if donor_type in ("N", "O", "S"):
                return "H bond"
    
            elif donor_type == "C":
                return "Nonconv H bond"
    
    
        if atom2_type == "H" and atom2_id in h_donors:
    
            _, donor_type = h_donors[atom2_id]
    
            if donor_type in ("N", "O", "S"):
                return "H bond"
    
            elif donor_type == "C":
                return "Nonconv H bond"
    
    
        return "H interaction"

    def get_interaction_bcp(self, interaction="all"):
        """
        Return bond critical points with interaction labels.
    
        Parameters
        ----------
        interaction : {"all", "covalent", "noncovalent"}, default="all"
            Type of interactions to return.
    
        Returns
        -------
        pandas.DataFrame
            Bond critical points with ``interaction`` and ``interaction_type`` columns.
        """
    
        bcp = self.get_bcp_properties()
    
        if bcp is None:
            return None
    
        bcp = bcp.copy()
    
        covalent_bcp = self._get_covalent_bcp()
    
        if covalent_bcp is None:
            return None
    
        # Identify covalent rows
        covalent_ids = set(covalent_bcp["id"])
        is_covalent = bcp["id"].isin(covalent_ids)
    
        h_donors = None
    
        if interaction in ("all", "noncovalent"):
            h_donors = self.h_donor()
    
        interaction_labels = []
        interaction_types = []
    
        for _, row in bcp.iterrows():
    
            atom1_id = row["atom1_id"]
            atom1_type = row["atom1_type"]
    
            atom2_id = row["atom2_id"]
            atom2_type = row["atom2_type"]
    
            cov = is_covalent.loc[row.name]
    
            # Generate interaction label
            if cov:
    
                label = f"{atom1_type}\u2014{atom2_type}"
    
            else:
    
                left = atom1_type
                right = atom2_type
    
                if h_donors is not None:
    
                    if atom1_type == "H" and atom1_id in h_donors:
                        _, donor_type = h_donors[atom1_id]
                        left = f"{donor_type}\u2014H"
    
                    if atom2_type == "H" and atom2_id in h_donors:
                        _, donor_type = h_donors[atom2_id]
                        right = f"H\u2014{donor_type}"
    
                label = f"{left}\u00b7\u00b7\u00b7{right}"
    
            interaction_labels.append(label)
    
            # Classify interaction type
            interaction_types.append(
                self._get_interaction_type(
                    atom1_type,
                    atom2_type,
                    atom1_id,
                    atom2_id,
                    cov,
                    h_donors
                )
            )
    
        # Insert columns after atom2_type
        position = bcp.columns.get_loc("atom2_type") + 1
    
        bcp.insert(position, "interaction", interaction_labels)
        bcp.insert(position + 1, "interaction_type", interaction_types)
    
        # Return selected interactions
        if interaction == "all":
    
            return GPUAMdf(
                bcp.reset_index(drop=True)
            )
    
        elif interaction == "covalent":
    
            return GPUAMdf(
                bcp[is_covalent].reset_index(drop=True)
            )
    
        elif interaction == "noncovalent":
    
            return GPUAMdf(
                bcp[~is_covalent].reset_index(drop=True)
            )
    
        else:
            raise ValueError(
                "interaction must be one of {'all', 'covalent', 'noncovalent'}."
            )
            
        
# Example usage for displaying critical points data:
# from gpuampy.io_tools import Crit

#gpuam_data = Crit("moleculeCrit.log")
#gpuam_data.read_data()
#bcp_data = gpuam_data.get_bcp()
#Obtener los no covalentes
#bcp_data = gpuam_data.get_interaction_bcp(interaction="all")

