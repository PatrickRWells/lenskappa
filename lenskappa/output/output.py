from abc import ABC, abstractmethod
from multiprocessing import Lock
from pathlib import Path
from typing import List
import pandas as pd

class outputHandler(ABC):

    def __init__(self, path: Path, *args, **kwargs):
        self._path = path
        self._lock = Lock()
    

    @abstractmethod
    def write_output(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def take_output(self, output, *args, **kwargs):
        pass



class csvOutputHandler(outputHandler):

    def __init__(self, path: Path, columns: List[str], *args, **kwargs):
        super().__init__(path, *args, **kwargs)
        self._columns = set(columns)
        self._df = pd.DataFrame(columns=columns)
    
    def write_output(self, *args, **kwargs):
        with self._lock:
            self._df.to_csv(self._path)

    def take_output(self, output: pd.DataFrame, *args, **kwargs):
        if set(output.columns) != self._columns:
            print("Error! Output handler was given an output with different columns!")
            return
        
        with self._lock:
            self._df = pd.concat([self._df, output])
