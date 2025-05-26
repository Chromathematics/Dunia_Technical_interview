# Dunia_Technical_interview
Technical Interview for Mathilde Roux candidate for internship

The idea was to propose a simple computational way to rank nanoparticles from their reactivity, but without using observables that need DFT nor energy calculations, which requiers complex background researches (functionals, pseudopotential, environmental explicit or implicit simulation), computational time, and hard to achieve on Python. 
 # installation of the virtual environment
```bash
pip install uv
uv venv .mathilde
source .mathilde/bin/activate
uv pip install fireworks numpy pydantic ase pymatgen
```


 # Running the Workflow
 in the terminal
 ```bash
python submit_workflow.py
rlaunch rapidfire
```

 # Machine answer
The code is a simple estimation of the S/V ratio of a nanoparticle, by estimating the number of undercoordinated atoms (which are more likely to be on the surface of the nanoparticle. Undercoordinated atoms are more reactive than bulk atoms, as the atoms "misses" electrons to reach stability. By estimating the number of binding of each atom specie, then the mean value, it is therefore possible to discriminate undercoordinated atoms. For a same number of atoms, nanoparticles with a lower ratio could tend to be more reactive.

