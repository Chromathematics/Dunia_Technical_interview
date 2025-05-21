from fireworks import Firework, LaunchPad, Workflow
from src.firetasks.pa_np_packing import PANPPackingConfig, PANPPackingFireTask
from src.storage import LocalStorage

def submit_workflow(wf) -> None:
    launchpad = LaunchPad.from_file("my_launchpad.yaml")
    launchpad.add_wf(wf)


if __name__ == "__main__":
    storage = LocalStorage(local_storage_path="/home/mathilde/temp_dir") # change this here

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

    workflow = Workflow(
        [fw1],
        {},
        name="Workflow",
    )  # Add more Fireworks as needed here

    # Submit the workflow to the launchpad
    submit_workflow(workflow)

    # Now you need to run the Workflow with `rlaunch rapidfire` command