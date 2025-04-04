from abc import ABC, abstractmethod
from typing import Callable

from backend.structure.protein import Protein


class ES(ABC):
    @abstractmethod
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None] = None, callback_frequency: int = 1) -> Protein:
        ...