"""Result extraction, fiber/GP decoding, and mesh mapping."""

import re
import numpy as np

from .model import ModelState, ResultDescriptor, ResultSubGroup
from .mesh_builder import MeshBuilder
from .mpco_reader import read_step_data
from .ops_elements import lookup_or_guess
from .utils import compute_von_mises_shell, compute_magnitude


NODAL_COMPONENTS = {
    "DISPLACEMENT": ["Ux", "Uy", "Uz", "|U|"],
    "REACTION_FORCE": ["Fx", "Fy", "Fz", "|F|"],
}


class ResultProcessor:
    """Extracts results from HDF5 and maps them onto PyVista meshes."""

    def __init__(self, model: ModelState, mesh_builder: MeshBuilder):
        self.model = model
        self.builder = mesh_builder

    def get_available_results(self) -> dict:
        available = {"Nodal": {}, "Element": {}}
        for rd in self.model.result_descriptors:
            if rd.category == "ON_NODES":
                comps = NODAL_COMPONENTS.get(rd.result_name,
                    [f"C{i}" for i in range(rd.num_components)] + ["|mag|"])
                available["Nodal"][rd.result_name] = {"components": comps}
            else:
                if not rd.sub_groups:
                    continue
                sg0 = rd.sub_groups[0]
                comps = list(sg0.component_names)
                # Only offer von Mises if there are enough stress components
                if sg0.num_components >= 3:
                    comps.append("von Mises")
                max_fibers = max(sg.num_fibers for sg in rd.sub_groups)
                max_gp = max(sg.num_gauss_pts for sg in rd.sub_groups)
                gp_list = ["Average"] + [f"GP {i}" for i in range(max_gp)]
                available["Element"][rd.result_name] = {
                    "components": comps,
                    "gauss_points": gp_list,
                    "max_fibers": max_fibers,
                }
        return available

    def get_result_descriptor(self, category: str, result_name: str):
        cat_key = "ON_NODES" if category == "Nodal" else "ON_ELEMENTS"
        for rd in self.model.result_descriptors:
            if rd.category == cat_key and rd.result_name == result_name:
                return rd
        return None

    # --- Nodal results: mapped to mid-surface mesh ---

    def extract_nodal_result(self, mpco_path: str, stage_index: int,
                             step_key: str, result_name: str,
                             component_idx: int, h5file=None) -> np.ndarray:
        """Extract nodal result mapped to shell mid-surface mesh points."""
        data = read_step_data(mpco_path, stage_index, result_name, "ON_NODES", step_key, h5file=h5file)
        if "nodes" not in data:
            return None

        raw = data["nodes"]
        if component_idx < raw.shape[1]:
            values = raw[:, component_idx]
        else:
            values = compute_magnitude(raw)

        return self._map_nodal_to_surface(values)

    def _map_nodal_to_surface(self, values):
        mapping = self.builder.shell_node_to_layer_points
        if not mapping:
            return None
        n_pts = max(max(indices) for indices in mapping.values()) + 1
        result = np.full(n_pts, np.nan)
        for nid, pt_indices in mapping.items():
            node_idx = self.model.nodes.id_to_index.get(nid)
            if node_idx is None:
                continue
            val = values[node_idx]
            for pi in pt_indices:
                if pi < n_pts:
                    result[pi] = val
        return result

    # --- Element results as contour ---

    def extract_element_result_contour(self, mpco_path: str, stage_index: int,
                                        step_key: str, result_name: str,
                                        gauss_pt: int, fiber_idx: int,
                                        component_idx: int, h5file=None) -> np.ndarray:
        """Extract element result extrapolated to mid-surface mesh points."""
        data = read_step_data(mpco_path, stage_index, result_name, "ON_ELEMENTS", step_key, h5file=h5file)
        if not data:
            return None

        rd = self.get_result_descriptor("Element", result_name)
        if rd is None:
            return None

        node_values = {}

        for sg in rd.sub_groups:
            raw = data.get(sg.hdf5_group)
            if raw is None:
                continue

            n_elems = raw.shape[0]
            n_gp = sg.num_gauss_pts
            n_fib = sg.num_fibers
            n_comp = sg.num_components

            if fiber_idx >= n_fib:
                continue

            is_von_mises = component_idx >= n_comp

            # Get per-element-type extrapolation matrix
            elem_info = self._get_sg_elem_info(sg)

            if gauss_pt == -1:
                gp_values = np.zeros(n_elems)
                for gp in range(n_gp):
                    gp_values += self._extract_gp_fiber_component(
                        raw, gp, fiber_idx, component_idx if not is_von_mises else None,
                        n_fib, n_comp, is_von_mises
                    )
                gp_values /= n_gp
                self._assign_elem_values_to_nodes(sg.elem_ids, gp_values, node_values)
            else:
                all_gp_values = np.zeros((n_elems, n_gp))
                for gp in range(n_gp):
                    all_gp_values[:, gp] = self._extract_gp_fiber_component(
                        raw, gp, fiber_idx, component_idx if not is_von_mises else None,
                        n_fib, n_comp, is_von_mises
                    )
                self._extrapolate_to_nodes(sg.elem_ids, all_gp_values, node_values, elem_info)

        if not node_values:
            return None

        return self._node_values_to_surface(node_values)

    # --- Element fiber results as GP spheres ---

    def extract_element_result_for_gp_cloud(self, mpco_path: str, stage_index: int,
                                             step_key: str, result_name: str,
                                             component_idx: int,
                                             gp_elem_ids: np.ndarray,
                                             gp_indices: np.ndarray,
                                             fiber_indices: np.ndarray,
                                             h5file=None) -> np.ndarray:
        """Extract element fiber result as scalar per Gauss-point sphere."""
        data = read_step_data(mpco_path, stage_index, result_name, "ON_ELEMENTS", step_key, h5file=h5file)
        if not data:
            return None

        rd = self.get_result_descriptor("Element", result_name)
        if rd is None:
            return None

        n_pts = len(gp_elem_ids)
        result = np.full(n_pts, np.nan)

        for sg in rd.sub_groups:
            raw = data.get(sg.hdf5_group)
            if raw is None:
                continue

            n_fib = sg.num_fibers
            n_comp = sg.num_components
            is_von_mises = component_idx >= n_comp

            eid_to_row = {int(eid): i for i, eid in enumerate(sg.elem_ids)}

            for pt_i in range(n_pts):
                eid = int(gp_elem_ids[pt_i])
                row = eid_to_row.get(eid)
                if row is None:
                    continue

                fi = int(fiber_indices[pt_i])
                if fi >= n_fib:
                    continue

                gpi = int(gp_indices[pt_i])
                base_col = gpi * (n_fib * n_comp) + fi * n_comp

                if is_von_mises:
                    s11 = raw[row, base_col + 0]
                    s22 = raw[row, base_col + 1]
                    s12 = raw[row, base_col + 2]
                    if n_comp >= 5:
                        s33 = raw[row, base_col + 3]
                        s13 = raw[row, base_col + 4]
                        result[pt_i] = compute_von_mises_shell(s11, s22, s12, s33, s13)
                    else:
                        result[pt_i] = compute_von_mises_shell(s11, s22, s12)
                else:
                    result[pt_i] = raw[row, base_col + component_idx]

        return result

    # --- Internal helpers ---

    def _get_sg_elem_info(self, sg: ResultSubGroup):
        """Get ElemInfo for a result sub-group from its elem_type or HDF5 key."""
        if sg.elem_type:
            return lookup_or_guess(sg.elem_type, 0)
        # Fallback: parse from hdf5_group key
        match = re.match(r"\d+-(\w+)\[", sg.hdf5_group)
        if match:
            return lookup_or_guess(match.group(1), 0)
        return None

    def _extract_gp_fiber_component(self, raw, gp, fiber_idx, comp_idx,
                                     n_fib, n_comp, is_von_mises):
        base_col = gp * (n_fib * n_comp) + fiber_idx * n_comp
        if is_von_mises:
            s11 = raw[:, base_col + 0]
            s22 = raw[:, base_col + 1]
            s12 = raw[:, base_col + 2]
            if n_comp >= 5:
                s33 = raw[:, base_col + 3]
                s13 = raw[:, base_col + 4]
                return compute_von_mises_shell(s11, s22, s12, s33, s13)
            return compute_von_mises_shell(s11, s22, s12)
        else:
            return raw[:, base_col + comp_idx]

    def _assign_elem_values_to_nodes(self, elem_ids, elem_values, node_values):
        for group in self.model.get_shell_groups():
            eid_to_row = {int(eid): i for i, eid in enumerate(group.elem_ids)}
            for i, eid in enumerate(elem_ids):
                eid = int(eid)
                row = eid_to_row.get(eid)
                if row is None:
                    continue
                nids = group.node_ids_per_elem[row]
                val = elem_values[i]
                for nid in nids:
                    nid = int(nid)
                    if nid not in node_values:
                        node_values[nid] = []
                    node_values[nid].append(val)

    def _extrapolate_to_nodes(self, elem_ids, all_gp_values, node_values, elem_info):
        """Extrapolate GP values to corner nodes using per-element-type matrix."""
        if elem_info is None or elem_info.extrap_fn is None:
            # No extrapolation available — fall back to averaging
            avg = all_gp_values.mean(axis=1)
            self._assign_elem_values_to_nodes(elem_ids, avg, node_values)
            return

        E = elem_info.extrap_fn()
        nodal_extrap = all_gp_values @ E.T
        for group in self.model.get_shell_groups():
            eid_to_row = {int(eid): i for i, eid in enumerate(group.elem_ids)}
            for i, eid in enumerate(elem_ids):
                eid = int(eid)
                row = eid_to_row.get(eid)
                if row is None:
                    continue
                nids = group.node_ids_per_elem[row]
                for j, nid in enumerate(nids):
                    nid = int(nid)
                    if nid not in node_values:
                        node_values[nid] = []
                    node_values[nid].append(nodal_extrap[i, j])

    def _node_values_to_surface(self, node_values):
        """Average at shared nodes and expand to surface mesh points."""
        mapping = self.builder.shell_node_to_layer_points
        if not mapping:
            return None
        n_pts = max(max(indices) for indices in mapping.values()) + 1
        result = np.full(n_pts, np.nan)
        for nid, vals in node_values.items():
            avg = np.nanmean(vals)
            for pi in mapping.get(nid, []):
                if pi < n_pts:
                    result[pi] = avg
        return result

    def get_displacement_data(self, mpco_path: str, stage_index: int,
                              step_key: str, h5file=None) -> np.ndarray:
        data = read_step_data(mpco_path, stage_index, "DISPLACEMENT", "ON_NODES", step_key, h5file=h5file)
        return data.get("nodes")
