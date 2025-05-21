# Dunia Recruitment task
This code is using the FireWorks workflow framework to submit computational workflows.

## 1. Installation of the package
### 1.1 Install uv
```bash
pip install uv
```
or 
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
```

### 1.2 Create a venv with uv
```bash
uv venv .mathilde
source .mathilde/bin/activate
uv pip install fireworks numpy pydantic ase
```

### 1.3 Run env.sh
Make sure you run this everytime you open a new linux shell, because it sets the PYTHONPATH.
```bash
source env.sh
```

## 2. Your assignment
Nice! Now you are ready to add new FireTasks to the workflow (in `submit_workflow.py`) in order to refine the initial structure and/or compute some descriptors for OER.

### 2.2 How to run the workflow
1. Run `python submit_workflow.py`
2. Run `rlaunch rapidfire`