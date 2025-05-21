import os
import shutil
from typing import Optional
from abc import abstractmethod

from pydantic import BaseModel


class Storage(BaseModel):
    type: str

    @abstractmethod
    def upload(self, job_id: str, file_path: str) -> str:
        """
        Uploads a file to the local storage.
        """
        raise NotImplementedError



class LocalStorage(Storage):
    local_storage_path: str
    type: str = "local"

    def upload(self, job_id: str, file_path: str) -> Optional[str]:
        """
        Uploads a file to the local storage.
        """
        if not os.path.exists(self.local_storage_path):
            os.makedirs(self.local_storage_path)

        raw_output_file_dir = os.path.join(self.local_storage_path, job_id)
        os.makedirs(raw_output_file_dir, exist_ok=True)

        raw_output_file_destination = os.path.join(
            raw_output_file_dir, os.path.basename(file_path)
        )

        if os.path.exists(file_path):
            shutil.move(file_path, raw_output_file_destination)
            return raw_output_file_destination
        else:
            return None
