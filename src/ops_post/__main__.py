"""Entry point for ops-post.

Usage:
    ops-post model              # looks for model.tcl + model.mpco
    ops-post model.mpco         # same (extension stripped)
    ops-post model.tcl          # same
    ops-post                    # file dialog

Same-name convention: the TCL model and MPCO results share a base name.
If model.tcl exists and model.mpco.postdata is missing or stale,
the TCL is parsed automatically to generate it.
"""

import sys
import os
import ctypes


def _request_discrete_gpu():
    """Hint to the OS to use the discrete GPU for this process.

    On Windows, NVIDIA and AMD drivers check for exported symbols or
    global variables to decide which GPU to assign. We set them via
    ctypes before any OpenGL context is created.

    Falls back silently if the driver DLLs are not present.
    """
    if sys.platform != "win32":
        return

    # NVIDIA: NvOptimusEnablement = 1 in nvapi/nvcuda
    # AMD:    AmdPowerXpressRequestHighPerformance = 1
    try:
        ctypes.windll.kernel32.SetEnvironmentVariableW(
            "SHIM_MCCOMPAT", "0x800000001"  # undocumented Windows GPU hint
        )
    except Exception:
        pass

    # Set env vars that some drivers read before context creation
    os.environ.setdefault("MESA_GL_VERSION_OVERRIDE", "4.5")

    # The standard approach: export global symbols. Since we're in Python,
    # we load the DLLs and set the globals directly.
    for dll_name, var_name in [
        ("nvapi64", "NvOptimusEnablement"),
        ("nvapi", "NvOptimusEnablement"),
        ("atig6pxx", "AmdPowerXpressRequestHighPerformance"),
        ("atiadlxx", "AmdPowerXpressRequestHighPerformance"),
    ]:
        try:
            dll = ctypes.WinDLL(dll_name)
            var = ctypes.c_int.in_dll(dll, var_name)
            var.value = 1
        except (OSError, ValueError, AttributeError):
            continue


# Request discrete GPU before any Qt/VTK imports create an OpenGL context
_request_discrete_gpu()


def _resolve_paths(arg: str):
    """From a user argument, derive the base name and locate files.

    Returns (tcl_path_or_None, mpco_path, postdata_path).
    """
    # Strip known extensions to get the base name
    base = arg
    for ext in (".tcl", ".mpco", ".mpco.postdata", ".mpco.cdata"):
        if base.lower().endswith(ext):
            base = base[:-len(ext)]
            break

    base = os.path.abspath(base)
    tcl_path = base + ".tcl"
    mpco_path = base + ".mpco"
    postdata_path = mpco_path + ".postdata"

    tcl_path = tcl_path if os.path.exists(tcl_path) else None

    return tcl_path, mpco_path, postdata_path


def _needs_regeneration(tcl_path, postdata_path):
    """Check if .postdata needs to be (re)generated from .tcl."""
    if tcl_path is None:
        return False
    if not os.path.exists(postdata_path):
        return True
    # Regenerate if TCL is newer than postdata
    return os.path.getmtime(tcl_path) > os.path.getmtime(postdata_path)


def main():
    from PyQt5.QtWidgets import QApplication, QFileDialog

    app = QApplication(sys.argv)

    if len(sys.argv) >= 2:
        tcl_path, mpco_path, postdata_path = _resolve_paths(sys.argv[1])
    else:
        # File dialog — ask for the .mpco file
        mpco_path, _ = QFileDialog.getOpenFileName(
            None, "Open MPCO File", os.getcwd(), "MPCO Files (*.mpco)"
        )
        if not mpco_path:
            sys.exit(0)
        tcl_path, mpco_path, postdata_path = _resolve_paths(mpco_path)

    if not os.path.exists(mpco_path):
        print(f"Error: result file not found: {mpco_path}")
        sys.exit(1)

    # Auto-generate .postdata from TCL if needed
    if _needs_regeneration(tcl_path, postdata_path):
        print(f"Parsing: {tcl_path}")
        from .tcl_parser import parse_tcl, write_cdata
        tcl_model = parse_tcl(tcl_path)
        print(f"  {len(tcl_model.nodes)} nodes, {len(tcl_model.elements)} elements, "
              f"{len(tcl_model.sections)} sections")
        write_cdata(tcl_model, postdata_path)

    from .gui import MainWindow
    win = MainWindow(mpco_path)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
