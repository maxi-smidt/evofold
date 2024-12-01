import numpy as np
import numpy.typing as npt

from dataclasses import dataclass


@dataclass
class Atom:
    atom_id: str
    x: float
    y: float
    z: float

    def get_position(self) -> npt.ArrayLike:
        return np.array([self.x, self.y, self.z])