"""Domain dataclasses for the MPCO post-processor."""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from .ops_elements import lookup_or_guess, ElemInfo


@dataclass
class NodeData:
    ids: np.ndarray          # (N,) int
    coords: np.ndarray       # (N, 3) float
    id_to_index: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.id_to_index = {int(nid): i for i, nid in enumerate(self.ids)}


@dataclass
class ElementGroup:
    elem_type: str           # e.g. "ASDShellQ4", "ElasticTimoshenkoBeam3d"
    hdf5_key: str            # e.g. "203-ASDShellQ4[201:0]"
    connectivity: np.ndarray # (M, cols) int — col 0 = elem_id, rest = node ids
    elem_ids: np.ndarray = field(default=None, repr=False)
    node_ids_per_elem: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        self.elem_ids = self.connectivity[:, 0]
        self.node_ids_per_elem = self.connectivity[:, 1:]


@dataclass
class SectionAssignment:
    section_key: str
    assignment: np.ndarray   # (P, 2) — [elem_id, fiber_index]
    fiber_data: np.ndarray   # (F, 3) — fiber coordinates
    num_fibers: int = 0

    def __post_init__(self):
        self.num_fibers = self.fiber_data.shape[0]

    @property
    def thickness(self) -> float:
        """Auto-detect section thickness from fiber z-coordinates."""
        z = self.fiber_data[:, 1]  # y-column holds through-thickness position
        if self.num_fibers <= 1:
            return abs(self.fiber_data[0, 2]) * 2 if self.num_fibers == 1 else 0.0
        spacing = (z.max() - z.min()) / (self.num_fibers - 1)
        return spacing * self.num_fibers

    @property
    def unique_elem_ids(self) -> np.ndarray:
        return np.unique(self.assignment[:, 0])


@dataclass
class BeamProfile:
    profile_id: int
    vertices: np.ndarray     # (N, 2) — 2D cross-section in local y-z

    @property
    def width(self) -> float:
        return self.vertices[:, 0].max() - self.vertices[:, 0].min()

    @property
    def height(self) -> float:
        return self.vertices[:, 1].max() - self.vertices[:, 1].min()


@dataclass
class ElementInfo:
    elem_id: int
    geometry_name: str
    section_name: str
    property_name: str


@dataclass
class ResultSubGroup:
    hdf5_group: str          # e.g. "203-ASDShellQ4[201:0:0]"
    elem_ids: np.ndarray
    num_gauss_pts: int
    num_fibers: int
    num_components: int
    component_names: list
    step_keys: list
    elem_type: str = ""      # e.g. "ASDShellQ4", parsed from hdf5_group key


@dataclass
class ResultDescriptor:
    category: str            # "ON_NODES" or "ON_ELEMENTS"
    result_name: str         # e.g. "DISPLACEMENT", "section.fiber.stress"
    sub_groups: list = field(default_factory=list)  # list of ResultSubGroup (element results only)
    node_ids: Optional[np.ndarray] = None           # for nodal results
    num_components: int = 0
    step_keys: list = field(default_factory=list)


@dataclass
class StageInfo:
    stage_index: int
    step_keys: list          # sorted step keys for this stage
    label: str = ""


@dataclass
class ModelState:
    nodes: NodeData
    element_groups: list                               # list of ElementGroup
    section_assignments: list                           # list of SectionAssignment
    beam_profiles: dict = field(default_factory=dict)   # profile_id -> BeamProfile
    beam_profile_assignments: dict = field(default_factory=dict)  # elem_id -> profile_id
    local_axes: dict = field(default_factory=dict)      # elem_id -> 3x3 rotation matrix
    element_info: dict = field(default_factory=dict)    # elem_id -> ElementInfo
    result_descriptors: list = field(default_factory=list)
    stages: list = field(default_factory=list)          # list of StageInfo
    total_steps: int = 0

    def get_elem_info(self, group: ElementGroup) -> Optional[ElemInfo]:
        """Look up element topology info from the registry."""
        return lookup_or_guess(group.elem_type, group.node_ids_per_elem.shape[1])

    def get_shell_groups(self) -> list:
        """Return all surface element groups (shells, 2D quads/tris)."""
        return [g for g in self.element_groups
                if (info := self.get_elem_info(g)) is not None and info.is_surface]

    def get_beam_groups(self) -> list:
        """Return all beam/truss element groups."""
        return [g for g in self.element_groups
                if (info := self.get_elem_info(g)) is not None and info.is_beam]

    def get_section_for_element(self, elem_id: int) -> Optional[SectionAssignment]:
        for sa in self.section_assignments:
            if elem_id in sa.unique_elem_ids:
                return sa
        return None
