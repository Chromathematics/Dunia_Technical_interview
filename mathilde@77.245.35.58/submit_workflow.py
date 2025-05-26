from fireworks import Firework, LaunchPad, Workflow
from src.firetasks.pa_np_packing import PANPPackingConfig, PANPPackingFireTask
from src.storage import LocalStorage
from pymatgen.core import Structure
import os
import subprocess
from pymatgen.analysis.local_env import CrystalNN
from collections import defaultdict
from fireworks import FiretaskBase, explicit_serialize

def submit_workflow(wf) -> None:
    launchpad = LaunchPad.from_file("my_launchpad.yaml")
    launchpad.add_wf(wf)

from pymatgen.io.cif import CifParser

def import_structure(cifdoc):
    parser = CifParser(cifdoc)
    structure = parser.get_structures()[0]
    return(structure)

def importation2():
    folder=r"temp_dir//"
    liste = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]
    most_recent_folder = max(liste, key=lambda f: os.path.getmtime(os.path.join(folder, f)))
    folder=folder+'//'+most_recent_folder
    liste=[f for f in os.listdir(folder) if f.endswith("cif")]
    structure=import_structure(folder+"//"+liste[0])
    return(folder+"//"+liste[0])


def default():
    structure=importation2()
    structure=Structure.from_file(structure)
    cnn=CrystalNN()
    coord_numbers= []
    species_to_coords = defaultdict(list)
    structure.add_oxidation_state_by_element({"Ag": 1.2,"O":-2,"Au":1.2})
    for i, site in enumerate(structure):
        try: 
            cn=cnn.get_cn(structure,i)
            species_to_coords[site.specie].append(cn)
            coord_numbers.append((i, site.specie, cn))
        except:
            coord_numbers.append((i, site.specie, 0))
    average={specie: sum(cns)/len(cns) for specie, cns in species_to_coords.items()}
    threshold = 0.8  # 
    undercoordinated_atoms = []

    for i, specie, cn in coord_numbers:
        mean=average[specie]
        if cn < threshold *mean:
            undercoordinated_atoms.append((i, specie.symbol, cn, round(mean, 2)))
    print(len(undercoordinated_atoms))
    print(len(coord_numbers))
    ratio=len(undercoordinated_atoms)/len(coord_numbers)
    print(f"estimated ratio : {ratio}")
    return(ratio)

@explicit_serialize
class UndercoordinationRatioTask(FiretaskBase):
    required_params = []
    def run_task(self, fw_spec):
        from pymatgen.core import Structure
        stored_files = fw_spec.get("stored_files", None)
        if not stored_files:
            raise RuntimeError("stored_files not found in fw_spec!")

        # Assuming the first file is the one you want
        structure_path = stored_files[0]  
        structure = Structure.from_file(structure_path)
        ratio = self.default2(structure)
        return {"undercoordination_ratio": ratio}

    def default2(self, structure):
        from pymatgen.analysis.local_env import CrystalNN
        from collections import defaultdict

        structure.add_oxidation_state_by_element({"Ag": 1.2, "O": -2, "Au": 1.2})
        cnn = CrystalNN()
        coord_numbers = []
        species_to_coords = defaultdict(list)

        for i, site in enumerate(structure):
            try:
                cn = cnn.get_cn(structure, i)
                species_to_coords[site.specie].append(cn)
                coord_numbers.append((i, site.specie, cn))
            except:
                coord_numbers.append((i, site.specie, 0))

        average = {specie: sum(cns) / len(cns) for specie, cns in species_to_coords.items()}
        threshold = 0.8
        undercoordinated_atoms = []

        for i, specie, cn in coord_numbers:
            mean = average[specie]
            if cn < threshold * mean:
                undercoordinated_atoms.append((i, specie.symbol, cn, round(mean, 2)))

        ratio = len(undercoordinated_atoms) / len(coord_numbers)
        print(f"The estimated undercoordinated/mean_coordination for every specie (O,Ag,Au) is {ratio}")
        return ratio

def tri_a_bulle(liste):
    n=len(liste)
    for i in range(n):
        bulle=True
        for j in range(n-i-1):
            if liste[j]>liste[j+1]:
                liste[j],liste[j+1]=liste[j+1],liste[j]
                bulle=False
        if bulle:
            break
    return(liste)
#print(tri_a_bulle(liste)) #OK






if __name__ == "__main__":
    storage = LocalStorage(local_storage_path="/home/mathilde/var/mmm/mathilde/temp_dir") # change this here

    fw1 = Firework(
        [
            PANPPackingFireTask(
                pa_np_packing_config=PANPPackingConfig(
                    oxidation_level=1.2,
                    random_seed=42,
                ).model_dump(),
                storage=storage.model_dump(),
                composition={"Au": 5, "Ag": 5, 'O': 10},
            ),
        ],
        name="Task 1",
    )
    fw2 = Firework(
        [
            UndercoordinationRatioTask()
        ],
        name="Test_ratio_1_NP",
        
    )

    workflow = Workflow(
        [fw1,fw2],
        {fw1: [fw2]},
        name="honestlyidk",
    )  # Add more Fireworks as needed here

    # Submit the workflow to the launchpad
    submit_workflow(workflow)

    # Now you need to run the Workflow with `rlaunch rapidfire` command
