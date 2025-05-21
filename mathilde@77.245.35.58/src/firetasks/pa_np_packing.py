from datetime import datetime
from typing import List

import ase
import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase.io import write

from src.storage import LocalStorage
from src.utils import clean_atoms

from fireworks import FWAction, explicit_serialize, FireTaskBase
from pydantic import BaseModel

# Constants
DISTANCE_BETWEEN_NPS = 10  # Distance between nanoparticles in Å
MINIMUM_OXIDATION_LEVEL = 0
MAXIMUM_OXIDATION_LEVEL = 1.5


### === COMPOSITION CLASSES === ###


class Space:
    """Base class for defining spatial constraints."""

    pass


class Box(Space):
    """Defines a cubic simulation box with a given length."""

    def __init__(self, length: float) -> None:
        self.length = length


class Composition:
    """Represents a chemical composition with elements and their atomic counts."""

    def __init__(self, atoms: List[str], numbers: List[int]) -> None:
        if len(atoms) != len(numbers):
            raise ValueError("Number of atoms and numbers must be equal")
        if not all(isinstance(n, int) for n in numbers):
            raise ValueError("Numbers must be integers")

        self.atoms = atoms
        self.numbers = numbers

    @property
    def string(self) -> str:
        """Returns a string representation of the composition (e.g., Cu10Ni20O45)."""
        return "".join(
            f"{a}{n}" for a, n in zip(self.atoms, self.numbers, strict=False)
        )

    def __eq__(self, value: object) -> bool:
        return isinstance(value, Composition) and self.string == value.string

    def __repr__(self) -> str:
        return self.string


### === ATOMIC PROPERTIES === ###


def get_atomic_radius(element: str) -> float:
    """Returns the covalent radius of the element in Angstrom."""
    return covalent_radii[atomic_numbers[element]]


def get_atomic_volume(element: str) -> float:
    """Estimates atomic volume using the atomic radius."""
    radius = get_atomic_radius(element)
    return (4 / 3) * np.pi * radius**3


### === NANOPARTICLE FUNCTIONS === ###


def create_oxidized_composition(
    composition: Composition, oxidation_level: float
) -> Composition:
    """
    Adds oxygen to the given composition based on the oxidation level.

    Args:
        composition (Composition): Initial composition.
        oxidation_level (float): Oxygen atoms added per total atoms.

    Returns:
        Composition: New oxidized composition.
    """
    oxygen_atoms = int(oxidation_level * sum(composition.numbers))
    return Composition(composition.atoms + ["O"], composition.numbers + [oxygen_atoms])


def estimate_nanoparticle_diameter(
    composition: Composition, packing_efficiency: float = 0.74
) -> float:
    """
    Estimates the nanoparticle diameter assuming spherical packing.

    Args:
        composition (Composition): Nanoparticle composition.
        packing_efficiency (float): Fraction of volume occupied by spheres.

    Returns:
        float: Estimated diameter in Angstrom.
    """
    total_volume = sum(
        get_atomic_volume(atom) * num
        for atom, num in zip(composition.atoms, composition.numbers, strict=False)
    )
    effective_volume = total_volume / packing_efficiency
    return 2 * (3 * effective_volume / (4 * np.pi)) ** (1 / 3)  # Diameter = 2 * radius


def compute_box_size(
    composition: Composition, distance_between_nps: float = DISTANCE_BETWEEN_NPS
) -> float:
    """
    Computes the simulation box length based on nanoparticle size and spacing.

    Args:
        composition (Composition): Nanoparticle composition.
        distance_between_nps (float): Desired spacing.

    Returns:
        float: Box length in Angstrom.
    """
    return estimate_nanoparticle_diameter(composition) + distance_between_nps


### === STRUCTURE GENERATION === ###


def create_initial_structure(
    composition: Composition,
    space: Box,
    packing_efficiency: float = 0.74,
    max_iterations: int = 1000,
    random_seed: int = 42,
) -> ase.Atoms:
    """
    Generates a packed atomic structure constrained in a sphere.

    Args:
        composition (Composition): Nanoparticle composition.
        space (Box): Simulation box.
        packing_efficiency (float): Packing efficiency.
        max_iterations (int): Max iterations to resolve overlaps.

    Returns:
        ase.Atoms: Packed atomic structure.
    """
    np.random.seed(random_seed)  # Set random seed for reproducibility

    nanoparticle_diameter = estimate_nanoparticle_diameter(
        composition, packing_efficiency
    )
    sphere_radius = nanoparticle_diameter / 2
    box_center = np.full(3, space.length / 2)

    atoms = np.repeat(composition.atoms, composition.numbers)
    radii = [get_atomic_radius(atom) for atom in atoms]

    def random_position_in_sphere(center, radius):
        while True:
            point = np.random.uniform(-radius, radius, 3) + center
            if np.linalg.norm(point - center) <= radius:
                return point

    positions = np.array(
        [random_position_in_sphere(box_center, sphere_radius) for _ in atoms]
    )

    # Iterative overlap resolution
    for iteration in range(max_iterations):
        overlaps_resolved = True
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                delta = positions[j] - positions[i]
                dist = np.linalg.norm(delta)
                min_dist = radii[i] + radii[j]

                if dist < min_dist:
                    overlaps_resolved = False
                    correction = (min_dist - dist) / 2
                    correction_vector = correction * delta / dist
                    positions[i] -= correction_vector
                    positions[j] += correction_vector

            # Ensure atoms remain within the sphere
            dist_to_center = np.linalg.norm(positions[i] - box_center)
            if dist_to_center + radii[i] > sphere_radius:
                positions[i] = (
                    box_center
                    + (positions[i] - box_center)
                    * (sphere_radius - radii[i])
                    / dist_to_center
                )

        if overlaps_resolved:
            print(f"All overlaps resolved at iteration {iteration}.")
            break

    return ase.Atoms(
        symbols=atoms.tolist(),
        positions=positions,
        cell=[(0, 0, space.length), (0, space.length, 0), (space.length, 0, 0)],
        pbc=True,
    )


class PANPPackingConfig(BaseModel):
    """
    Config model for IR spectrum postprocessing.
    """

    oxidation_level: float
    random_seed: int


@explicit_serialize
class PANPPackingFireTask(FireTaskBase):
    """
    FireWorks task for generating a packed nanoparticle structure.
    """

    @property
    def optional_params(self) -> list[str]:
        """
        List of optional parameters for the task.
        """
        return [
            "storage",
            "pa_np_packing_config",
            "composition",
        ]

    @property
    def task_config_name(self) -> str:
        return "pa_np_packing_config"

    @property
    def task_config(self) -> PANPPackingConfig:
        return PANPPackingConfig

    def run_task(self, fw_spec: dict) -> FWAction:
        # First, check if the parameters are explicitly passed
        config_data = self.get(self.task_config_name)
        storage_data = self.get("storage")

        # Ensure all required inputs are present
        if config_data is None:
            raise ValueError(
                f"Missing required input parameters: {self.task_config_name}."
            )
        if storage_data is None:
            raise ValueError("Missing required input parameters: 'storage'.")

        pa_np_packing_config = self.task_config(**config_data)

        if storage_data["type"] == "local":
            storage = LocalStorage(**storage_data)

        # First, check if the parameters are explicitly passed
        composition = self.get("composition") or fw_spec.get("composition")

        # Ensure all required inputs are present
        if composition is None:
            raise ValueError("Missing required input parameters: 'composition'.")

        composition = Composition(
            list(composition.keys()), list(composition.values())
        )
        oxidized_composition = create_oxidized_composition(
            composition, pa_np_packing_config.oxidation_level
        )

        print(
            f"Estimated diameters: {estimate_nanoparticle_diameter(oxidized_composition):.2f} Å (default), "
            f"{estimate_nanoparticle_diameter(oxidized_composition, 0.32):.2f} Å (loose)"
        )

        box = Box(length=compute_box_size(oxidized_composition))
        atoms = create_initial_structure(
            oxidized_composition,
            box,
            packing_efficiency=0.74,
            max_iterations=1000,
            random_seed=pa_np_packing_config.random_seed,
        )
        atoms = clean_atoms(atoms)

        write("packed_structure.cif", atoms)
        unique_id = f"{atoms.get_chemical_formula()}_PA_NP_packing_{datetime.now().strftime('%Y-%m-%d_%H:%M:%S')}"
        stored_files = [storage.upload(unique_id, "packed_structure.cif")]

        return FWAction(
            update_spec={"atoms": atoms.todict(), "stored_files": stored_files},
            stored_data={"atoms": atoms.todict(), "stored_files": stored_files},
        )