import ase

def clean_atoms(atoms: ase.Atoms) -> ase.Atoms:
    """Create a new ASE Atoms object keeping only the basic properties."""
    return ase.Atoms(
        symbols=atoms.get_chemical_symbols(),  # Element types
        positions=atoms.get_positions(),  # Atomic positions
        momenta=atoms.get_momenta(),  # Atomic momenta, we need those to continue MD trajectories e.g.
        cell=atoms.get_cell(),  # Unit cell
        pbc=atoms.get_pbc(),  # Periodic boundary conditions
    )