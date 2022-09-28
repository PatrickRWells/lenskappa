from abc import ABC, abstractmethod
from ast import parse
from multiprocessing import Lock
from pathlib import Path
from typing import List
import pandas as pd
from . import parser

class outputHandler(ABC):

    def __init__(self, path: Path, parser: parser.lenskappaOutputParser, *args, **kwargs):
        self._path = path
        self._parser = parser()
        self._lock = Lock()
    
    

    @abstractmethod
    def write_output(self, *args, **kwargs):
        pass
    
    def set_lock(self, lock):
        self._lock = lock

    def take_output(self, output, *args, **kwargs):
        parsed_output = self._parser(output)
        self._take_output(parsed_output)

    def _take_output(self, output, *args, **kwargs):
        pass


class csvOutputHandler(outputHandler):

    def __init__(self, path: Path, parser, columns: List[str], *args, **kwargs):
        super().__init__(path, parser, *args, **kwargs)
        self._columns = set(columns)
        self._df = pd.DataFrame(columns=columns)
    
    def __len__(self):
        return len(self._df)

    def write_output(self, *args, **kwargs):
        with self._lock:
            self._df.to_csv(self._path, *args, **kwargs)

    def _take_output(self, output: pd.DataFrame, *args, **kwargs):
        if set(output.columns) != self._columns:
            print("Error! Output handler was given an output with different columns!")
            return
        
        with self._lock:
            self._df = pd.concat([self._df, output], ignore_index=True)
