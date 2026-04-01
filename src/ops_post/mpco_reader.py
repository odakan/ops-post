"""MPCO (HDF5) and .cdata file reader."""

import re
import os
import numpy as np
import h5py

from .model import (
    NodeData, ElementGroup, SectionAssignment, BeamProfile, ElementInfo,
    ResultDescriptor, ResultSubGroup, StageInfo, ModelState,
)
from .utils import quaternion_to_rotation_matrix


def read_mpco(mpco_path: str) -> ModelState:
    """Read an MPCO file and its .postdata companion, returning the full ModelState.

    The .postdata file is generated from TCL model files by the tcl_parser
    module and contains local axes and element info.
    """
    mpco_path = os.path.abspath(mpco_path)
    postdata_path = mpco_path + ".postdata"

    with h5py.File(mpco_path, "r") as f:
        # Detect stages
        stage_keys = sorted(
            [k for k in f.keys() if k.startswith("MODEL_STAGE[")],
            key=lambda k: int(k.split("[")[1].rstrip("]"))
        )

        # Read model from first stage (geometry is same across stages)
        first_stage = f[stage_keys[0]]
        nodes = _read_nodes(first_stage)
        element_groups = _read_elements(first_stage)
        section_assignments = _read_section_assignments(first_stage)

        # Discover results across all stages
        stages = []
        result_descriptors = []
        all_step_keys = []

        for stage_key in stage_keys:
            stage_grp = f[stage_key]
            stage_idx = int(stage_key.split("[")[1].rstrip("]"))

            # Collect step keys from displacement (always present)
            disp_grp = stage_grp.get("RESULTS/ON_NODES/DISPLACEMENT/DATA")
            if disp_grp is not None:
                step_keys = sorted(disp_grp.keys(), key=lambda x: int(x.split("_")[1]))
            else:
                step_keys = []

            stages.append(StageInfo(
                stage_index=stage_idx,
                step_keys=step_keys,
                label=f"Stage {stage_idx}",
            ))
            all_step_keys.extend(step_keys)

        # Discover result descriptors from first stage that has results
        for stage_key in stage_keys:
            stage_grp = f[stage_key]
            result_descriptors = _discover_results(stage_grp, stages)
            if result_descriptors:
                break

    # Read companion .postdata file
    local_axes = {}
    beam_profiles = {}
    beam_profile_assignments = {}
    element_info = {}

    if os.path.exists(postdata_path):
        local_axes, beam_profiles, beam_profile_assignments, element_info = _read_postdata(postdata_path)

    return ModelState(
        nodes=nodes,
        element_groups=element_groups,
        section_assignments=section_assignments,
        beam_profiles=beam_profiles,
        beam_profile_assignments=beam_profile_assignments,
        local_axes=local_axes,
        element_info=element_info,
        result_descriptors=result_descriptors,
        stages=stages,
        total_steps=len(all_step_keys),
    )


def read_step_data(mpco_path: str, stage_index: int, result_name: str,
                   category: str, step_key: str, h5file=None) -> dict:
    """Read result data for a single step. Returns dict of sub_group_key -> np.ndarray.

    If h5file is provided (an open h5py.File), uses it directly instead of
    opening/closing the file.  This avoids repeated open/close overhead
    during animation.
    """
    def _read(f):
        result = {}
        stage_key = f"MODEL_STAGE[{stage_index}]"
        base = f"{stage_key}/RESULTS/{category}/{result_name}"
        grp = f.get(base)
        if grp is None:
            return result
        if category == "ON_NODES":
            data_grp = grp.get("DATA")
            if data_grp is not None and step_key in data_grp:
                result["nodes"] = data_grp[step_key][:]
        else:
            for sub_key in grp.keys():
                data_grp = grp.get(f"{sub_key}/DATA")
                if data_grp is not None and step_key in data_grp:
                    result[sub_key] = data_grp[step_key][:]
        return result

    if h5file is not None:
        return _read(h5file)

    with h5py.File(mpco_path, "r") as f:
        return _read(f)


def _read_nodes(stage_grp) -> NodeData:
    nodes_grp = stage_grp["MODEL/NODES"]
    ids = nodes_grp["ID"][:]
    coords = nodes_grp["COORDINATES"][:]
    return NodeData(ids=ids, coords=coords)


def _read_elements(stage_grp) -> list:
    groups = []
    elem_grp = stage_grp["MODEL/ELEMENTS"]
    for key in elem_grp.keys():
        # Parse element type from key like "203-ASDShellQ4[201:0]"
        match = re.match(r"\d+-(\w+)\[", key)
        elem_type = match.group(1) if match else key
        conn = elem_grp[key][:]
        groups.append(ElementGroup(
            elem_type=elem_type,
            hdf5_key=key,
            connectivity=conn,
        ))
    return groups


def _read_section_assignments(stage_grp) -> list:
    sa_grp = stage_grp.get("MODEL/SECTION_ASSIGNMENTS")
    if sa_grp is None:
        return []
    assignments = []
    for key in sa_grp.keys():
        sub = sa_grp[key]
        assignments.append(SectionAssignment(
            section_key=key,
            assignment=sub["ASSIGNMENT"][:],
            fiber_data=sub["FIBER_DATA"][:],
        ))
    return assignments


def _discover_results(stage_grp, stages: list) -> list:
    """Walk RESULTS/ to build ResultDescriptor objects."""
    descriptors = []
    results_grp = stage_grp.get("RESULTS")
    if results_grp is None:
        return descriptors

    # Collect all step keys across all stages
    all_step_keys = []
    for st in stages:
        all_step_keys.extend(st.step_keys)

    # Nodal results
    on_nodes = results_grp.get("ON_NODES")
    if on_nodes is not None:
        for result_name in on_nodes.keys():
            rgrp = on_nodes[result_name]
            node_ids = rgrp["ID"][:].flatten() if "ID" in rgrp else None
            data_grp = rgrp.get("DATA")
            if data_grp is None:
                continue
            sample_key = list(data_grp.keys())[0]
            num_components = data_grp[sample_key].shape[1] if len(data_grp[sample_key].shape) > 1 else 1
            descriptors.append(ResultDescriptor(
                category="ON_NODES",
                result_name=result_name,
                node_ids=node_ids,
                num_components=num_components,
                step_keys=all_step_keys,
            ))

    # Element results
    on_elements = results_grp.get("ON_ELEMENTS")
    if on_elements is not None:
        for result_name in on_elements.keys():
            rgrp = on_elements[result_name]
            sub_groups = []
            for sub_key in rgrp.keys():
                sub = rgrp[sub_key]
                meta = sub.get("META")
                if meta is None:
                    continue

                # Parse element type from sub_key (e.g. "203-ASDShellQ4[201:0:0]")
                sub_match = re.match(r"\d+-(\w+)\[", sub_key)
                sub_elem_type = sub_match.group(1) if sub_match else ""

                elem_ids = sub["ID"][:].flatten()
                gauss_ids = meta["GAUSS_IDS"][:].flatten()
                num_gp = len(gauss_ids)
                multiplicities = meta["MULTIPLICITY"][:].flatten()
                num_comp_arr = meta["NUM_COMPONENTS"][:].flatten()
                num_fibers = int(multiplicities[0])
                num_components = int(num_comp_arr[0])

                # Parse component names from COMPONENTS string
                comp_str = meta["COMPONENTS"][0]
                if isinstance(comp_str, bytes):
                    comp_str = comp_str.decode()
                component_names = _parse_component_names(comp_str, num_components)

                data_grp = sub.get("DATA")
                step_keys = sorted(data_grp.keys(), key=lambda x: int(x.split("_")[1])) if data_grp else []

                sub_groups.append(ResultSubGroup(
                    hdf5_group=sub_key,
                    elem_ids=elem_ids,
                    num_gauss_pts=num_gp,
                    num_fibers=num_fibers,
                    num_components=num_components,
                    component_names=component_names,
                    step_keys=step_keys,
                    elem_type=sub_elem_type,
                ))

            if sub_groups:
                descriptors.append(ResultDescriptor(
                    category="ON_ELEMENTS",
                    result_name=result_name,
                    sub_groups=sub_groups,
                    step_keys=all_step_keys,
                ))

    return descriptors


def _parse_component_names(comp_str: str, num_components: int) -> list:
    """Parse component names from the META/COMPONENTS string.

    Format: "0.1.2.3.4.Name0,Name1,Name2,Name3,Name4;..." (one block per GP)
    We only need the names from the first GP block.
    """
    first_gp = comp_str.split(";")[0]
    # Format: "0.1.2.3.4.Name0,Name1,Name2,Name3,Name4"
    parts = first_gp.split(".")
    # The last part contains comma-separated names
    if len(parts) > num_components:
        names_str = parts[num_components]
        names = names_str.split(",")
    else:
        names = [f"C{i}" for i in range(num_components)]

    # Clean up names - replace "Unknown" with generic labels
    cleaned = []
    for i, name in enumerate(names):
        name = name.strip()
        if not name or "Unknown" in name:
            cleaned.append(f"C{i}")
        else:
            cleaned.append(name)
    return cleaned


def _read_postdata(postdata_path: str):
    """Parse the .mpco.postdata companion file (generated from TCL).

    Returns (local_axes, beam_profiles, beam_profile_assignments, element_info).
    """
    local_axes = {}
    beam_profiles = {}
    beam_profile_assignments = {}
    element_info = {}

    with open(postdata_path, "r") as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == "*LOCAL_AXES":
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("*") or line.startswith("#"):
                    break
                if line:
                    parts = line.split()
                    if len(parts) >= 5:
                        eid = int(parts[0])
                        qw, qx, qy, qz = float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
                        local_axes[eid] = quaternion_to_rotation_matrix(qw, qx, qy, qz)
                i += 1
            continue

        if line == "*BEAM_PROFILE":
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("*") or line.startswith("#"):
                    break
                if line:
                    parts = line.split()
                    profile_id = int(parts[0])
                    n_verts = int(parts[1])
                    vertices = []
                    for _ in range(n_verts):
                        i += 1
                        xy = lines[i].strip().split()
                        vertices.append([float(xy[0]), float(xy[1])])
                    beam_profiles[profile_id] = BeamProfile(
                        profile_id=profile_id,
                        vertices=np.array(vertices),
                    )
                i += 1
            continue

        if line == "*BEAM_PROFILE_ASSIGNMENT":
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("*") or line.startswith("#"):
                    break
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        eid = int(parts[0])
                        prof_id = int(parts[1])
                        beam_profile_assignments[eid] = prof_id
                i += 1
            continue

        if line == "*ELEMENT_INFO":
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("*") or line.startswith("#"):
                    break
                if line:
                    parts = line.split(None, 2)
                    if len(parts) >= 2:
                        eid = int(parts[0])
                        elem_type = parts[1]
                        sec_name = parts[2] if len(parts) > 2 else ""
                        element_info[eid] = ElementInfo(
                            elem_id=eid,
                            geometry_name=elem_type,
                            section_name=sec_name,
                            property_name="",
                        )
                i += 1
            continue

        i += 1

    return local_axes, beam_profiles, beam_profile_assignments, element_info
