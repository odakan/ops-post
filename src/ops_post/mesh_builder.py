"""Build PyVista meshes from the parsed model data."""

import numpy as np
import pyvista as pv

from .model import ModelState
from .ops_elements import lookup_or_guess
from .utils import detect_up_axis


class MeshBuilder:
    """Converts ModelState into PyVista meshes for 3D visualization.

    Shell/surface elements are visualized as:
      - Transparent extruded volumes (hex for quads, wedge for tris)
      - Mid-surface faces (quads or tris)
      - Thin lines at each fiber layer position
      - Colored spheres at Gauss-point x fiber locations
    """

    def __init__(self, model: ModelState):
        self.model = model
        self._shell_data = []
        self.shell_node_to_layer_points = {}
        self.shell_elem_to_cell = {}
        self.shell_node_to_hex_points = {}
        self.beam_node_to_hex_points = {}
        self.beam_elem_to_cell = {}

        # Detect up-axis from model geometry
        self._up_axis = detect_up_axis(model.nodes.coords)
        self._up_vector = np.zeros(3)
        self._up_vector[self._up_axis] = 1.0

        self._precompute_shell_geometry()

    def _precompute_shell_geometry(self):
        nodes = self.model.nodes
        for group in self.model.get_shell_groups():
            elem_info = self.model.get_elem_info(group)
            for row_idx in range(len(group.elem_ids)):
                eid = int(group.elem_ids[row_idx])
                nids = group.node_ids_per_elem[row_idx]
                sa = self.model.get_section_for_element(eid)
                thickness = sa.thickness if sa else 0.1

                node_coords = np.array([
                    nodes.coords[nodes.id_to_index[int(nid)]] for nid in nids
                ])

                rot = self.model.local_axes.get(eid)
                if rot is not None:
                    normal = rot[:, 2]
                else:
                    v1 = node_coords[1] - node_coords[0]
                    v2 = node_coords[-1] - node_coords[0]
                    normal = np.cross(v1, v2)
                    norm = np.linalg.norm(normal)
                    normal = normal / norm if norm > 1e-12 else self._up_vector.copy()

                self._shell_data.append({
                    "eid": eid,
                    "nids": nids,
                    "node_coords": node_coords,
                    "normal": normal,
                    "thickness": thickness,
                    "sa": sa,
                    "elem_info": elem_info,
                })

    def build_shell_extrusion_mesh(self) -> pv.UnstructuredGrid:
        """Build transparent extruded volumes showing shell thickness.

        Quads → hexahedra (8 pts), tris → wedges (6 pts).
        """
        if not self._shell_data:
            return pv.UnstructuredGrid()

        all_points = []
        all_cells = []
        all_cell_types = []
        all_elem_ids = []
        point_idx = 0
        self.shell_node_to_hex_points = {}

        for sd in self._shell_data:
            n_nodes = len(sd["node_coords"])
            offset = sd["normal"] * (sd["thickness"] / 2.0)
            bottom = sd["node_coords"] - offset
            top = sd["node_coords"] + offset
            vol_points = np.vstack([bottom, top])
            all_points.append(vol_points)

            n_total = 2 * n_nodes
            cell = [n_total] + list(range(point_idx, point_idx + n_total))
            all_cells.append(cell)
            if n_nodes == 3:
                all_cell_types.append(pv.CellType.WEDGE)
            else:
                all_cell_types.append(pv.CellType.HEXAHEDRON)
            all_elem_ids.append(sd["eid"])

            for local_i, nid in enumerate(sd["nids"]):
                nid = int(nid)
                if nid not in self.shell_node_to_hex_points:
                    self.shell_node_to_hex_points[nid] = []
                self.shell_node_to_hex_points[nid].append(point_idx + local_i)
                self.shell_node_to_hex_points[nid].append(point_idx + n_nodes + local_i)
            point_idx += n_total

        points = np.vstack(all_points)
        cells = np.hstack(all_cells)
        cell_types = np.array(all_cell_types, dtype=np.uint8)
        mesh = pv.UnstructuredGrid(cells, cell_types, points)
        mesh.cell_data["elem_id"] = np.array(all_elem_ids, dtype=np.int32)
        return mesh

    def build_shell_surface_mesh(self) -> pv.PolyData:
        """Build mid-surface faces (quads or tris)."""
        if not self._shell_data:
            return pv.PolyData()

        all_points = []
        all_faces = []
        all_elem_ids = []
        point_idx = 0
        self.shell_node_to_layer_points = {}
        self.shell_elem_to_cell = {}

        for cell_idx, sd in enumerate(self._shell_data):
            n_nodes = len(sd["node_coords"])
            all_points.append(sd["node_coords"])
            all_faces.extend([n_nodes] + list(range(point_idx, point_idx + n_nodes)))
            self.shell_elem_to_cell[sd["eid"]] = cell_idx
            all_elem_ids.append(sd["eid"])

            for local_i, nid in enumerate(sd["nids"]):
                nid = int(nid)
                if nid not in self.shell_node_to_layer_points:
                    self.shell_node_to_layer_points[nid] = []
                self.shell_node_to_layer_points[nid].append(point_idx + local_i)
            point_idx += n_nodes

        points = np.vstack(all_points)
        faces = np.array(all_faces, dtype=np.int64)
        mesh = pv.PolyData(points, faces)
        mesh.cell_data["elem_id"] = np.array(all_elem_ids, dtype=np.int32)
        return mesh

    def build_fiber_layer_edges(self) -> pv.PolyData:
        """Build edge lines for all fiber layers inside the extrusion."""
        if not self._shell_data:
            return pv.PolyData()

        all_points = []
        all_lines = []
        point_idx = 0

        for sd in self._shell_data:
            sa = sd["sa"]
            if sa is None:
                continue
            n_nodes = len(sd["node_coords"])
            for fi in range(sa.num_fibers):
                fiber_offset = sa.fiber_data[fi, 1]
                offset = sd["normal"] * fiber_offset
                layer_coords = sd["node_coords"] + offset
                all_points.append(layer_coords)
                for i in range(n_nodes):
                    j = (i + 1) % n_nodes
                    all_lines.extend([2, point_idx + i, point_idx + j])
                point_idx += n_nodes

        if not all_points:
            return pv.PolyData()

        points = np.vstack(all_points)
        lines = np.array(all_lines, dtype=np.int64)
        return pv.PolyData(points, lines=lines)

    def build_gauss_point_cloud(self, fiber_idx: int = None) -> tuple:
        """Build a point cloud at all GP x fiber positions inside shell elements.

        Uses per-element shape functions and GP coordinates from the registry.

        Args:
            fiber_idx: if None, build for ALL fibers; if int, only that fiber.

        Returns:
            (points, elem_ids, gp_indices, fiber_indices)
        """
        all_points = []
        all_elem_ids = []
        all_gp = []
        all_fiber = []
        gp_displ_info = []

        for sd in self._shell_data:
            sa = sd["sa"]
            if sa is None:
                continue
            elem_info = sd["elem_info"]
            if elem_info is None or elem_info.gp_coords is None:
                continue

            node_coords = sd["node_coords"]
            normal = sd["normal"]
            eid = sd["eid"]
            nids = sd["nids"]
            gp_coords = elem_info.gp_coords
            shape_fn = elem_info.shape_fn

            fibers = [fiber_idx] if fiber_idx is not None else range(sa.num_fibers)

            for fi in fibers:
                if fi >= sa.num_fibers:
                    continue
                fiber_offset = sa.fiber_data[fi, 1]
                thickness_offset = normal * fiber_offset

                for gpi in range(len(gp_coords)):
                    xi, eta = gp_coords[gpi]
                    N = shape_fn(xi, eta)
                    pos = N @ node_coords + thickness_offset
                    all_points.append(pos)
                    all_elem_ids.append(eid)
                    all_gp.append(gpi)
                    all_fiber.append(fi)
                    gp_displ_info.append((nids, N))

        if not all_points:
            self._gp_displ_info = []
            self._gp_node_indices = np.empty((0, 0), dtype=np.int32)
            self._gp_weights = np.empty((0, 0))
            return np.empty((0, 3)), np.array([]), np.array([]), np.array([])

        points = np.array(all_points)
        elem_ids = np.array(all_elem_ids, dtype=np.int32)
        gp_indices = np.array(all_gp, dtype=np.int32)
        fiber_indices = np.array(all_fiber, dtype=np.int32)

        self._gp_displ_info = gp_displ_info

        # Precompute vectorized arrays for fast displacement interpolation
        # Each GP has variable node count (3 for T3, 4 for Q4, 9 for Q9).
        # Pad to max width so we can use a single matrix multiply.
        max_nodes = max(len(nids) for nids, _ in gp_displ_info)
        n_gp = len(gp_displ_info)
        node_idx_arr = np.zeros((n_gp, max_nodes), dtype=np.int32)
        weight_arr = np.zeros((n_gp, max_nodes))
        id_to_index = self.model.nodes.id_to_index

        for i, (nids, weights) in enumerate(gp_displ_info):
            for j, nid in enumerate(nids):
                idx = id_to_index.get(int(nid), 0)
                node_idx_arr[i, j] = idx
                weight_arr[i, j] = weights[j]

        self._gp_node_indices = node_idx_arr  # (N_gp, max_nodes)
        self._gp_weights = weight_arr          # (N_gp, max_nodes)

        return points, elem_ids, gp_indices, fiber_indices

    def displace_gp_points(self, base_pts: np.ndarray,
                           displacements: np.ndarray,
                           scale_factor: float,
                           mask: np.ndarray = None) -> np.ndarray:
        """Apply interpolated displacement to GP cloud — vectorized.

        Args:
            base_pts: (N, 3) base GP positions (already masked if mask given)
            displacements: (n_nodes, 3) full displacement array
            scale_factor: displacement scale
            mask: boolean mask into the full GP arrays (None = all)

        Returns:
            (N, 3) displaced positions
        """
        if mask is not None:
            idx = self._gp_node_indices[mask]   # (N, max_nodes)
            wts = self._gp_weights[mask]        # (N, max_nodes)
        else:
            idx = self._gp_node_indices
            wts = self._gp_weights

        # disp_per_node: (N, max_nodes, 3)
        d = displacements[idx]
        # weighted sum: (N, 3)
        interp = np.einsum("ij,ijk->ik", wts, d)
        return base_pts + interp * scale_factor

    def build_beam_mesh(self) -> pv.UnstructuredGrid:
        """Build beam mesh (beams -> hexahedra using cross-section profile)."""
        nodes = self.model.nodes
        beam_groups = self.model.get_beam_groups()
        if not beam_groups:
            return pv.UnstructuredGrid()

        all_points = []
        all_cells = []
        all_cell_types = []
        all_elem_ids = []
        point_idx = 0

        for group in beam_groups:
            for row_idx in range(len(group.elem_ids)):
                eid = int(group.elem_ids[row_idx])
                nids = group.node_ids_per_elem[row_idx]
                n1_idx = nodes.id_to_index.get(int(nids[0]))
                n2_idx = nodes.id_to_index.get(int(nids[1]))
                if n1_idx is None or n2_idx is None:
                    continue

                p1 = nodes.coords[n1_idx]
                p2 = nodes.coords[n2_idx]

                prof_id = self.model.beam_profile_assignments.get(eid)
                if prof_id is not None and prof_id in self.model.beam_profiles:
                    verts_2d = self.model.beam_profiles[prof_id].vertices
                else:
                    s = 0.025
                    verts_2d = np.array([[-s, -s], [s, -s], [s, s], [-s, s]])

                rot = self.model.local_axes.get(eid)
                if rot is None:
                    dx = p2 - p1
                    dx_norm = np.linalg.norm(dx)
                    if dx_norm < 1e-12:
                        continue
                    local_x = dx / dx_norm
                    temp = self._up_vector if abs(np.dot(local_x, self._up_vector)) < 0.9 else np.array([1, 0, 0])
                    local_y = np.cross(local_x, temp)
                    local_y /= np.linalg.norm(local_y)
                    local_z = np.cross(local_x, local_y)
                    rot = np.column_stack([local_x, local_y, local_z])

                local_y = rot[:, 1]
                local_z = rot[:, 2]
                nverts = min(4, len(verts_2d))
                corners_1 = np.array([
                    p1 + verts_2d[j, 0] * local_y + verts_2d[j, 1] * local_z
                    for j in range(nverts)
                ])
                corners_2 = np.array([
                    p2 + verts_2d[j, 0] * local_y + verts_2d[j, 1] * local_z
                    for j in range(nverts)
                ])

                hex_points = np.vstack([corners_1, corners_2])
                all_points.append(hex_points)
                cell = [8] + list(range(point_idx, point_idx + 8))
                all_cells.append(cell)
                all_cell_types.append(pv.CellType.HEXAHEDRON)

                self.beam_elem_to_cell[eid] = len(all_elem_ids)
                all_elem_ids.append(eid)

                for local_i, nid in enumerate(nids):
                    nid = int(nid)
                    if nid not in self.beam_node_to_hex_points:
                        self.beam_node_to_hex_points[nid] = []
                    base = point_idx + local_i * 4
                    for j in range(4):
                        self.beam_node_to_hex_points[nid].append(base + j)
                point_idx += 8

        if not all_points:
            return pv.UnstructuredGrid()

        points = np.vstack(all_points)
        cells = np.hstack(all_cells)
        cell_types = np.array(all_cell_types, dtype=np.uint8)
        mesh = pv.UnstructuredGrid(cells, cell_types, points)
        mesh.cell_data["elem_id"] = np.array(all_elem_ids, dtype=np.int32)
        return mesh

    def apply_displacement(self, mesh, displacements: np.ndarray,
                           scale_factor: float, node_to_points: dict):
        """Apply scaled displacements to a mesh copy."""
        displaced = mesh.copy(deep=True)
        nodes = self.model.nodes
        for nid, pt_indices in node_to_points.items():
            node_idx = nodes.id_to_index.get(nid)
            if node_idx is None:
                continue
            disp = displacements[node_idx] * scale_factor
            for pi in pt_indices:
                if pi < displaced.n_points:
                    displaced.points[pi] += disp
        return displaced
