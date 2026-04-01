"""Parse OpenSees TCL model files into an in-memory model object.

Handles single-file models (everything in main.tcl) or multi-file models
(main.tcl sources nodes.tcl, elements.tcl, etc.).

Produces a TclModel with enough information to generate a .cdata-like
companion file for the post-processor.
"""

import os
import re
import numpy as np
from dataclasses import dataclass, field
from scipy.spatial.transform import Rotation


# ── Data classes ────────────────────────────────────────────────────

@dataclass
class TclNode:
    id: int
    x: float
    y: float
    z: float = 0.0


@dataclass
class TclElement:
    id: int
    elem_type: str          # e.g. "ASDShellQ4", "ElasticTimoshenkoBeam"
    node_ids: list           # [n1, n2, ...]
    section_id: int = None   # section or material tag used
    geom_transf_id: int = None
    local_dir: list = None   # -local vx vy vz (for shells)
    orient: list = None      # -orient ... (for zero-length)
    extra_args: dict = field(default_factory=dict)


@dataclass
class TclGeomTransf:
    id: int
    transf_type: str         # "Linear", "PDelta", "Corotational"
    vecxz: list              # [vx, vy, vz] local-x-z plane vector


@dataclass
class TclSection:
    id: int
    section_type: str        # "LayeredShell", "Elastic", "Fiber", etc.
    layers: list = field(default_factory=list)  # [(mat_id, thickness), ...]
    args: list = field(default_factory=list)


@dataclass
class TclMaterial:
    id: int
    mat_type: str
    args: list = field(default_factory=list)


@dataclass
class TclModel:
    """Complete parsed TCL model."""
    ndm: int = 3
    ndf: int = 6
    nodes: dict = field(default_factory=dict)          # id -> TclNode
    elements: dict = field(default_factory=dict)        # id -> TclElement
    geom_transfs: dict = field(default_factory=dict)    # id -> TclGeomTransf
    sections: dict = field(default_factory=dict)        # id -> TclSection
    materials: dict = field(default_factory=dict)       # id -> TclMaterial
    source_files: list = field(default_factory=list)    # files that were parsed

    def get_nodes_array(self) -> np.ndarray:
        """Return (N, 3) coordinate array sorted by node ID."""
        sorted_ids = sorted(self.nodes.keys())
        return np.array([[self.nodes[i].x, self.nodes[i].y, self.nodes[i].z]
                         for i in sorted_ids])

    def get_node_ids(self) -> np.ndarray:
        return np.array(sorted(self.nodes.keys()), dtype=int)

    def compute_local_axes(self, elem: TclElement) -> np.ndarray:
        """Compute 3x3 rotation matrix for an element's local axes.

        Returns rotation matrix R where columns are [local_x, local_y, local_z]
        in global coordinates, or None if not determinable.
        """
        if elem.local_dir is not None:
            # Shell with -local vx vy vz: this is the local-x direction
            return self._shell_local_axes(elem)
        if elem.geom_transf_id is not None:
            return self._beam_local_axes(elem)
        if elem.orient is not None:
            return self._zero_length_local_axes(elem)
        return None

    def _shell_local_axes(self, elem: TclElement) -> np.ndarray:
        """Compute rotation matrix for a shell with -local directive."""
        # -local gives the local x-axis direction
        lx = np.array(elem.local_dir[:3], dtype=float)
        lx_norm = np.linalg.norm(lx)
        if lx_norm < 1e-12:
            return None
        lx = lx / lx_norm

        # Compute element normal from node positions
        coords = [self._node_coords(nid) for nid in elem.node_ids[:4]]
        if any(c is None for c in coords):
            return None
        v1 = coords[1] - coords[0]
        v2 = coords[-1] - coords[0]
        normal = np.cross(v1, v2)
        n_norm = np.linalg.norm(normal)
        if n_norm < 1e-12:
            return None
        lz = normal / n_norm

        # local_y = local_z x local_x
        ly = np.cross(lz, lx)
        ly_norm = np.linalg.norm(ly)
        if ly_norm < 1e-12:
            return None
        ly = ly / ly_norm

        # Re-orthogonalize local_x
        lx = np.cross(ly, lz)
        return np.column_stack([lx, ly, lz])

    def _beam_local_axes(self, elem: TclElement) -> np.ndarray:
        """Compute rotation matrix for a beam element from geomTransf."""
        gt = self.geom_transfs.get(elem.geom_transf_id)
        if gt is None:
            return None

        p1 = self._node_coords(elem.node_ids[0])
        p2 = self._node_coords(elem.node_ids[1])
        if p1 is None or p2 is None:
            return None

        lx = p2 - p1
        lx_norm = np.linalg.norm(lx)
        if lx_norm < 1e-12:
            return None
        lx = lx / lx_norm

        # vecxz defines the local x-z plane
        vecxz = np.array(gt.vecxz, dtype=float)
        ly = np.cross(vecxz, lx)
        ly_norm = np.linalg.norm(ly)
        if ly_norm < 1e-12:
            # vecxz parallel to element axis — try fallback
            temp = np.array([0, 0, 1]) if abs(lx[2]) < 0.9 else np.array([1, 0, 0])
            ly = np.cross(temp, lx)
            ly_norm = np.linalg.norm(ly)
            if ly_norm < 1e-12:
                return None
        ly = ly / ly_norm

        lz = np.cross(lx, ly)
        return np.column_stack([lx, ly, lz])

    def _zero_length_local_axes(self, elem: TclElement) -> np.ndarray:
        """Compute rotation matrix for a zero-length element from -orient."""
        orient = elem.orient
        if len(orient) >= 6:
            lx = np.array(orient[0:3], dtype=float)
            ly = np.array(orient[3:6], dtype=float)
            lx_norm = np.linalg.norm(lx)
            ly_norm = np.linalg.norm(ly)
            if lx_norm < 1e-12 or ly_norm < 1e-12:
                return None
            lx = lx / lx_norm
            ly = ly / ly_norm
            lz = np.cross(lx, ly)
            return np.column_stack([lx, ly, lz])
        elif len(orient) >= 3:
            lx = np.array(orient[0:3], dtype=float)
            lx_norm = np.linalg.norm(lx)
            if lx_norm < 1e-12:
                return None
            lx = lx / lx_norm
            temp = np.array([0, 1, 0]) if abs(lx[1]) < 0.9 else np.array([0, 0, 1])
            ly = np.cross(lx, temp)
            ly = ly / np.linalg.norm(ly)
            lz = np.cross(lx, ly)
            return np.column_stack([lx, ly, lz])
        return None

    def _node_coords(self, nid: int) -> np.ndarray:
        n = self.nodes.get(nid)
        if n is None:
            return None
        return np.array([n.x, n.y, n.z])


# ── TCL parser ──────────────────────────────────────────────────────

def parse_tcl(path: str) -> TclModel:
    """Parse an OpenSees TCL model from a main file.

    Follows `source` commands to read included files.
    Handles backslash line continuation.

    Args:
        path: path to the main .tcl file

    Returns:
        TclModel with all parsed data
    """
    model = TclModel()
    base_dir = os.path.dirname(os.path.abspath(path))
    _parse_file(os.path.abspath(path), model, base_dir)
    return model


def _parse_file(path: str, model: TclModel, base_dir: str):
    """Parse a single TCL file, following source commands."""
    if not os.path.exists(path):
        return

    model.source_files.append(path)

    with open(path, "r") as f:
        raw_lines = f.readlines()

    # Join backslash-continued lines
    lines = _join_continued_lines(raw_lines)

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        tokens = _tokenize(line)
        if not tokens:
            continue

        cmd = tokens[0]

        if cmd == "source":
            if len(tokens) >= 2:
                src_path = os.path.join(base_dir, tokens[1])
                _parse_file(src_path, model, os.path.dirname(src_path))

        elif cmd == "model":
            _parse_model_cmd(tokens, model)

        elif cmd == "node":
            _parse_node(tokens, model)

        elif cmd == "element":
            _parse_element(tokens, model)

        elif cmd == "geomTransf":
            _parse_geom_transf(tokens, model)

        elif cmd == "section":
            _parse_section(tokens, model)

        elif cmd == "nDMaterial" or cmd == "uniaxialMaterial":
            _parse_material(tokens, model)


def _join_continued_lines(raw_lines: list) -> list:
    """Join lines ending with backslash."""
    result = []
    buf = ""
    for line in raw_lines:
        stripped = line.rstrip("\n\r")
        if stripped.endswith("\\"):
            buf += stripped[:-1] + " "
        else:
            buf += stripped
            result.append(buf)
            buf = ""
    if buf:
        result.append(buf)
    return result


def _tokenize(line: str) -> list:
    """Split a TCL line into tokens, respecting braces and quotes."""
    tokens = []
    i = 0
    n = len(line)
    while i < n:
        # Skip whitespace
        while i < n and line[i] in " \t":
            i += 1
        if i >= n:
            break
        # Skip inline comments
        if line[i] == "#" or (line[i] == ";" and i + 1 < n and line[i + 1] == "#"):
            break
        if line[i] == ";":
            i += 1
            continue
        # Brace-quoted
        if line[i] == "{":
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if line[j] == "{":
                    depth += 1
                elif line[j] == "}":
                    depth -= 1
                j += 1
            tokens.append(line[i + 1:j - 1])
            i = j
        # Double-quoted
        elif line[i] == '"':
            j = i + 1
            while j < n and line[j] != '"':
                j += 1
            tokens.append(line[i + 1:j])
            i = j + 1
        # Square bracket (TCL command substitution) — treat as one token
        elif line[i] == "[":
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if line[j] == "[":
                    depth += 1
                elif line[j] == "]":
                    depth -= 1
                j += 1
            tokens.append(line[i:j])
            i = j
        # Regular token
        else:
            j = i
            while j < n and line[j] not in " \t;#":
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


# ── Command parsers ─────────────────────────────────────────────────

def _parse_model_cmd(tokens, model):
    """Parse: model basic -ndm 3 -ndf 6"""
    for i, t in enumerate(tokens):
        if t == "-ndm" and i + 1 < len(tokens):
            model.ndm = int(tokens[i + 1])
        elif t == "-ndf" and i + 1 < len(tokens):
            model.ndf = int(tokens[i + 1])


def _parse_node(tokens, model):
    """Parse: node $tag $x $y [$z] [-mass ...]"""
    if len(tokens) < 4:
        return
    try:
        tag = int(tokens[1])
        x = float(tokens[2])
        y = float(tokens[3])
        z = float(tokens[4]) if len(tokens) > 4 and not tokens[4].startswith("-") else 0.0
        model.nodes[tag] = TclNode(id=tag, x=x, y=y, z=z)
    except (ValueError, IndexError):
        pass


def _parse_element(tokens, model):
    """Parse element commands using the ops_elements registry.

    Dispatches by ElemFamily to handle the different TCL argument layouts:
      - Surface (shell/quad/tri): nodes + section tag + optional -local
      - Beam/column: 2 nodes + geomTransf + optional inline section props
      - Truss: 2 nodes + A + material tag
      - Link: 2 nodes + options
      - Zero-length: 2 nodes + -mat/-dir + optional -orient
      - Unknown: fallback to consecutive integer node IDs
    """
    if len(tokens) < 4:
        return

    elem_type = tokens[1]
    try:
        eid = int(tokens[2])
    except ValueError:
        return

    from .ops_elements import lookup, ElemFamily

    info = lookup(elem_type)

    if info is not None:
        family = info.family
        if family in (ElemFamily.SHELL_QUAD, ElemFamily.SHELL_TRI,
                      ElemFamily.QUAD_2D, ElemFamily.TRI_2D):
            _parse_surface_element(tokens, elem_type, eid, info.num_nodes, model)
        elif family == ElemFamily.BEAM:
            _parse_beam(tokens, elem_type, eid, model)
        elif family == ElemFamily.TRUSS:
            _parse_two_node(tokens, elem_type, eid, model)
        elif family == ElemFamily.LINK:
            _parse_two_node(tokens, elem_type, eid, model)
        elif family == ElemFamily.ZERO_LENGTH:
            _parse_zero_length(tokens, elem_type, eid, model)
        else:
            _parse_generic_element(tokens, elem_type, eid, model)
    else:
        _parse_generic_element(tokens, elem_type, eid, model)


# ── Family-based parsers ────────────────────────────────────────────

def _parse_surface_element(tokens, elem_type, eid, n_nodes, model):
    """Parse shell, quad, and tri elements.

    Common layout: element $type $tag $n1 ... $nN $secTag [options]
    Options may include: -local vx vy vz
    """
    node_start = 3
    node_end = node_start + n_nodes
    sec_idx = node_end  # section/material tag follows the nodes

    if len(tokens) <= sec_idx:
        return
    try:
        nids = [int(tokens[i]) for i in range(node_start, node_end)]
        sec_id = int(tokens[sec_idx])
    except (ValueError, IndexError):
        return

    # Extract -local direction if present
    local_dir = _extract_flag_floats(tokens, "-local", 3)

    model.elements[eid] = TclElement(
        id=eid, elem_type=elem_type, node_ids=nids,
        section_id=sec_id, local_dir=local_dir,
    )


def _parse_beam(tokens, elem_type, eid, model):
    """Parse beam/column elements (all registered beam types).

    Handles multiple argument layouts by scanning for the geomTransf tag
    and optional inline section properties (A, Iy, Iz).
    """
    if len(tokens) < 5:
        return
    try:
        nids = [int(tokens[3]), int(tokens[4])]
    except (ValueError, IndexError):
        return

    transf_id = None
    sec_id = None

    # --- Locate geomTransf tag ---
    # Layout A (forceBeamColumn, dispBeamColumn): transfTag at tokens[5]
    # Layout B (Timoshenko, elasticBeamColumn): transfTag is the last integer
    #          before any dash-flags

    # Try layout A first: tokens[5] is small integer (typical transf IDs)
    if len(tokens) > 5:
        try:
            candidate = int(tokens[5])
            # If tokens[6] is not a number, this is likely layout A
            if len(tokens) <= 6 or tokens[6].startswith("-"):
                transf_id = candidate
        except ValueError:
            pass

    # If layout A didn't work, scan backwards for the last integer (layout B)
    if transf_id is None:
        for t in reversed(tokens[5:]):
            if t.startswith("-"):
                continue
            try:
                transf_id = int(t)
                break
            except ValueError:
                continue

    elem = TclElement(
        id=eid, elem_type=elem_type, node_ids=nids,
        section_id=sec_id, geom_transf_id=transf_id,
    )

    # --- Extract inline section properties for beam profile computation ---
    # Maps element type pattern -> (A_index, Iy_index, Iz_index) in tokens
    #
    # ElasticTimoshenkoBeam: element $type tag n1 n2 E  G  A  Jx Iy Iz Ay Az transfTag
    #   token indices:        [0]    [1]  [2] [3][4][5][6][7][8][9][10][11][12] [13]
    #
    # elasticBeamColumn 3D:  element $type tag n1 n2 A  E  G  Jx Iy Iz transfTag
    #   token indices:        [0]    [1]  [2] [3][4][5][6][7][8][9][10] [11]
    #
    _BEAM_INLINE_PROPS = {
        "Timoshenko":        (7, 9, 10),
        "elasticBeamColumn": (5, 9, 10),
    }
    for pattern, (a_idx, iy_idx, iz_idx) in _BEAM_INLINE_PROPS.items():
        if pattern in elem_type and len(tokens) > max(a_idx, iy_idx, iz_idx):
            try:
                elem._beam_A = float(tokens[a_idx])
                elem._beam_Iy = float(tokens[iy_idx])
                elem._beam_Iz = float(tokens[iz_idx])
            except (ValueError, IndexError):
                pass
            break

    model.elements[eid] = elem


def _parse_two_node(tokens, elem_type, eid, model):
    """Parse 2-node elements (truss, link) with no special options."""
    if len(tokens) < 5:
        return
    try:
        nids = [int(tokens[3]), int(tokens[4])]
    except ValueError:
        return
    model.elements[eid] = TclElement(
        id=eid, elem_type=elem_type, node_ids=nids,
    )


def _parse_zero_length(tokens, elem_type, eid, model):
    """Parse zero-length elements: 2 nodes + optional -orient."""
    if len(tokens) < 5:
        return
    try:
        nids = [int(tokens[3]), int(tokens[4])]
    except ValueError:
        return

    orient = _extract_flag_floats(tokens, "-orient", None)

    model.elements[eid] = TclElement(
        id=eid, elem_type=elem_type, node_ids=nids, orient=orient,
    )


def _parse_generic_element(tokens, elem_type, eid, model):
    """Fallback — extracts consecutive integer node IDs after the tag."""
    nids = []
    for t in tokens[3:]:
        try:
            nids.append(int(t))
        except ValueError:
            break
    if nids:
        model.elements[eid] = TclElement(
            id=eid, elem_type=elem_type, node_ids=nids,
        )


# ── Shared helpers ──────────────────────────────────────────────────

def _is_flag(token: str) -> bool:
    """Check if a token is a flag (e.g. -local, -orient) vs a negative number."""
    return (len(token) > 1 and token[0] == "-" and token[1].isalpha())


def _extract_flag_floats(tokens, flag, count):
    """Extract float values following a flag (e.g. -local vx vy vz).

    Args:
        tokens: full token list
        flag: flag string to search for (e.g. "-local")
        count: exact number of floats to read, or None to read all until next flag

    Returns list of floats, or None if flag not found.
    """
    for i, t in enumerate(tokens):
        if t == flag:
            vals = []
            for j in range(i + 1, len(tokens)):
                if _is_flag(tokens[j]):
                    break
                try:
                    vals.append(float(tokens[j]))
                except ValueError:
                    break
                if count is not None and len(vals) >= count:
                    break
            return vals if vals else None
    return None


def _parse_geom_transf(tokens, model):
    """Parse: geomTransf $type $tag $vecxzX $vecxzY $vecxzZ [options]"""
    if len(tokens) < 6:
        return
    try:
        transf_type = tokens[1]
        tag = int(tokens[2])
        vecxz = [float(tokens[3]), float(tokens[4]), float(tokens[5])]
    except (ValueError, IndexError):
        return

    model.geom_transfs[tag] = TclGeomTransf(
        id=tag, transf_type=transf_type, vecxz=vecxz,
    )


def _parse_section(tokens, model):
    """Parse section commands.

    LayeredShell: section LayeredShell $tag $nLayers $mat1 $t1 $mat2 $t2 ...
    Elastic:      section Elastic $tag $E $A $Iz ...
    """
    if len(tokens) < 3:
        return
    sec_type = tokens[1]
    try:
        tag = int(tokens[2])
    except ValueError:
        return

    layers = []
    if sec_type == "LayeredShell":
        if len(tokens) > 3:
            try:
                n_layers = int(tokens[3])
                for i in range(n_layers):
                    idx = 4 + i * 2
                    if idx + 1 < len(tokens):
                        mat_id = int(tokens[idx])
                        thickness = float(tokens[idx + 1])
                        layers.append((mat_id, thickness))
            except (ValueError, IndexError):
                pass

    model.sections[tag] = TclSection(
        id=tag, section_type=sec_type, layers=layers,
        args=tokens[3:],
    )


def _parse_material(tokens, model):
    """Parse: nDMaterial/uniaxialMaterial $type $tag ..."""
    if len(tokens) < 3:
        return
    mat_type = tokens[1]
    try:
        tag = int(tokens[2])
    except ValueError:
        return
    model.materials[tag] = TclMaterial(
        id=tag, mat_type=mat_type, args=tokens[3:],
    )


# ── .cdata generation ───────────────────────────────────────────────

def rotation_to_quaternion(R: np.ndarray) -> tuple:
    """Convert 3x3 rotation matrix to quaternion (qw, qx, qy, qz)."""
    r = Rotation.from_matrix(R)
    qx, qy, qz, qw = r.as_quat()  # scipy: scalar-last
    return qw, qx, qy, qz


def _compute_beam_profiles(model: TclModel) -> tuple:
    """Derive rectangular beam cross-section profiles from TCL data.

    Sources:
    - section Elastic (3D): section Elastic $tag $E $A $Iz $Iy $G $J ...
    - Inline Timoshenko: element ElasticTimoshenkoBeam $tag $n1 $n2 $E $G $A $Jx $Iy $Iz ...

    Assumes rectangular: A=b*h, Iz=b*h^3/12 => h=sqrt(12*Iz/A), b=A/h

    Returns (profiles, assignments):
      profiles: dict[profile_id -> (width, height)]
      assignments: dict[elem_id -> profile_id]
    """
    profiles = {}  # profile_id -> (b, h)
    assignments = {}  # elem_id -> profile_id

    # From section Elastic (3D): args = [E, A, Iz, Iy, G, J, ...]
    for sid, sec in model.sections.items():
        if sec.section_type == "Elastic" and len(sec.args) >= 4:
            try:
                A = float(sec.args[1])   # args[0]=E, args[1]=A
                Iz = float(sec.args[2])  # args[2]=Iz
                if A > 0 and Iz > 0:
                    h = np.sqrt(12.0 * Iz / A)
                    b = A / h if h > 1e-12 else np.sqrt(A)
                    profiles[sid] = (b, h)
            except (ValueError, IndexError):
                pass

    # Assign section-based profiles to beam elements
    for eid, elem in model.elements.items():
        if elem.section_id is not None and elem.section_id in profiles:
            assignments[eid] = elem.section_id

    # From inline Timoshenko beam element properties
    for eid, elem in model.elements.items():
        if eid not in assignments and elem.geom_transf_id is not None:
            _try_inline_beam_profile(elem, profiles, assignments)

    return profiles, assignments


def _try_inline_beam_profile(elem, profiles, assignments):
    """For ElasticTimoshenkoBeam, extract A/Iz/Iy from inline element args.

    element ElasticTimoshenkoBeam $tag $n1 $n2 $E $G $A $Jx $Iy $Iz $Ay $Az $transfTag
    """
    if "Timoshenko" not in elem.elem_type:
        return

    # The element args (after tag, n1, n2) contain: E G A Jx Iy Iz Ay Az transfTag
    # We stored the full token list... but we only have node_ids and geom_transf_id
    # We need to re-derive from the original parse. Let's use a synthetic profile_id.
    # Actually, the section data for Timoshenko beams is in the element command itself.
    # We need to extract it during parsing. Let me store it on the element.
    if hasattr(elem, '_beam_A') and elem._beam_A is not None:
        A, Iz, Iy = elem._beam_A, elem._beam_Iz, elem._beam_Iy
        if A > 0 and Iz > 0:
            h = np.sqrt(12.0 * Iz / A)
            b = A / h if h > 1e-12 else A
            # Use negative eid as synthetic profile ID to avoid collisions
            prof_id = -elem.id
            profiles[prof_id] = (b, h)
            assignments[elem.id] = prof_id


def write_cdata(model: TclModel, output_path: str):
    """Write a .postdata file from a parsed TCL model.

    Generates:
    - *LOCAL_AXES: quaternion orientation for each element
    - *BEAM_PROFILE: rectangular cross-section vertices
    - *BEAM_PROFILE_ASSIGNMENT: element -> profile mapping
    - *ELEMENT_INFO: element type and section info
    """
    # Compute beam profiles from section data
    beam_profiles, beam_assignments = _compute_beam_profiles(model)

    with open(output_path, "w") as f:
        # --- LOCAL_AXES ---
        f.write("#elem_id q.w q.x q.y q.z\n")
        f.write("*LOCAL_AXES\n")

        for eid in sorted(model.elements.keys()):
            elem = model.elements[eid]
            R = model.compute_local_axes(elem)
            if R is not None:
                qw, qx, qy, qz = rotation_to_quaternion(R)
                f.write(f"{eid} {qw} {qx} {qy} {qz}\n")

        # --- BEAM_PROFILE ---
        if beam_profiles:
            f.write("#profile_id n_vertices\\n  y z (per vertex)\n")
            f.write("*BEAM_PROFILE\n")

            for prof_id in sorted(beam_profiles.keys()):
                b, h = beam_profiles[prof_id]
                hb, hh = b / 2, h / 2
                verts = [(-hb, -hh), (hb, -hh), (hb, hh), (-hb, hh)]
                f.write(f"{prof_id} {len(verts)}\n")
                for y, z in verts:
                    f.write(f"{y} {z}\n")

        # --- BEAM_PROFILE_ASSIGNMENT ---
        if beam_assignments:
            f.write("#elem_id profile_id\n")
            f.write("*BEAM_PROFILE_ASSIGNMENT\n")

            for eid in sorted(beam_assignments.keys()):
                f.write(f"{eid} {beam_assignments[eid]}\n")

        # --- ELEMENT_INFO ---
        f.write("#elem_id elem_type section_name\n")
        f.write("*ELEMENT_INFO\n")

        for eid in sorted(model.elements.keys()):
            elem = model.elements[eid]
            sec_name = ""
            if elem.section_id is not None:
                sec = model.sections.get(elem.section_id)
                sec_name = f"{sec.section_type}_{elem.section_id}" if sec else f"Section_{elem.section_id}"
            f.write(f"{eid} {elem.elem_type} {sec_name}\n")

    print(f"Wrote {output_path}")


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    """CLI entry point: parse TCL files and print model summary."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python -m ops_post.tcl_parser <main.tcl> [--write-cdata output.cdata]")
        sys.exit(1)

    tcl_path = sys.argv[1]
    model = parse_tcl(tcl_path)

    print(f"Parsed {len(model.source_files)} file(s):")
    for f in model.source_files:
        print(f"  {f}")
    print(f"\nModel: ndm={model.ndm}, ndf={model.ndf}")
    print(f"Nodes: {len(model.nodes)}")
    print(f"Elements: {len(model.elements)}")

    # Count by type
    type_counts = {}
    for e in model.elements.values():
        type_counts[e.elem_type] = type_counts.get(e.elem_type, 0) + 1
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {t}: {c}")

    print(f"GeomTransf: {len(model.geom_transfs)}")
    print(f"Sections: {len(model.sections)}")
    for sid, sec in sorted(model.sections.items()):
        layers_info = f", {len(sec.layers)} layers" if sec.layers else ""
        print(f"  Section {sid}: {sec.section_type}{layers_info}")
    print(f"Materials: {len(model.materials)}")

    # Count elements with computable local axes
    n_axes = sum(1 for e in model.elements.values() if model.compute_local_axes(e) is not None)
    print(f"Local axes computable: {n_axes}/{len(model.elements)}")

    # Write cdata if requested
    if "--write-cdata" in sys.argv:
        idx = sys.argv.index("--write-cdata")
        if idx + 1 < len(sys.argv):
            write_cdata(model, sys.argv[idx + 1])
        else:
            print("Error: --write-cdata requires an output path")


if __name__ == "__main__":
    main()
