
# 3rd party lib imports
import numpy as np

# gmap imports
import GMAP.src.tools.constants as GM_con
import GMAP.src.tools.file_handler as GM_fh
import GMAP.src.tools.parameter_parser as GM_pp
import GMAP.src.tools.print_tools as GM_pt

# own module imports
import TrEsp_code.TrEsp_clib as MC_tc

_ = GM_con.bohr  # to validify the import. The import is needed for exec.


def GM_change_coup_type(map_, system, oscix1, osc1, oscix2, osc2):
    return map_.name  # return TrEsp


def GM_post_init(map_, system):
    # load in the charges from all different kinds of oscillators
    map_.charges = {}
    map_charges = {}
    for osc in system.oscillators_ordered_coup[map_.name]:
        # if an issue occurs, no need to continue
        if not map_.success:
            return

        # If a map has a dedicated function, use that instead of interpreting
        # the provided file.

        if hasattr(osc.map.code, "CP_TrEsp_get_charges"):
            map_.charges[osc.oscix] = osc.map.code.CP_TrEsp_get_charges(
                    osc.map, system, osc)
        # only look for each type of singles once.
        elif osc.map.name in map_charges:
            map_.charges[osc.oscix] = map_charges[osc.map.name]
        else:
            map_charges[osc.map.name] = get_charges(map_, osc.map)
            map_.charges[osc.oscix] = map_charges[osc.map.name]

    MC_tc.init_map_for_clib(map_, system)


def GM_pre_run(map_, system):
    map_.osclens = [
        len(map_.charges.get(osc.oscix, [])) for osc in system.oscillators]
    # do this before converting osclens to array (faster)
    map_.all_used_ats = []
    for osclen, osc in zip(map_.osclens, system.oscillators):
        map_.all_used_ats.extend(osc.used_atoms[:osclen])

    map_.osclens = np.array(map_.osclens, dtype="int32")
    map_.osclens_c = np.ctypeslib.as_ctypes(map_.osclens)
    map_.oscstart = np.concatenate(
        (np.zeros(1, dtype="int32"), np.cumsum(map_.osclens)[:-1]),
        dtype="int32")  # set the dtype again, as linux changes it here.
    map_.oscstart_c = np.ctypeslib.as_ctypes(map_.oscstart)
    map_.all_used_ats = np.array(map_.all_used_ats, dtype="int32")
    map_.all_used_ats_c = np.ctypeslib.as_ctypes(map_.all_used_ats)

    map_.charge_array = []
    for osc in system.oscillators:
        map_.charge_array.extend(map_.charges.get(osc.oscix, []))
    map_.charge_array = np.array(map_.charge_array, dtype="float32")
    map_.charge_array_c = np.ctypeslib.as_ctypes(map_.charge_array)

    # map_.allpairs = np.array(map_.allpairs, dtype="int32")
    map_.allpairs_c = np.ctypeslib.as_ctypes(np.ravel(map_.allpairs))
    map_.n_allpairs = np.int32(map_.allpairs.shape[0])


def GM_calc_coupling(map_, system, hamiltonian):
    fpieps = np.float32(GM_con.e2i4pieps_angcm)
    hamiltonian_c = np.ctypeslib.as_ctypes(np.ravel(hamiltonian))
    map_.clib.calc_coupling(map_, system, fpieps, hamiltonian_c)


def get_charges(map_, oscmap):
    # This map contains the keyword for the TrEsp charges file, obtain
    # that file's name
    fname = gc_get_filename(map_, oscmap)
    if fname is None:
        return None

    contents = gc_get_file_contents(fname, map_, oscmap)
    if contents is None:
        return None

    keyword = f"{map_.name}.charges_multiply"
    # optional multiplication - if not needed, skip.
    if keyword not in oscmap.rawcore:
        return contents

    multiplier = gc_get_multiplier(keyword, map_, oscmap)
    if multiplier is None:
        return None

    return contents * multiplier


def gc_get_filename(map_, oscmap):
    """Helper function for get_charges. Gets filename for oscmap.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TrEsp map object
    oscmap: :class:`~GMAP.src.tools.map_reader.SingleMap`
        The map object of the map for which we'd like to obtain TrEsp
        charges.

    Returns
    -------
    fname : pathlib.Path or None
        The path to the file that stores the TrEsp charges. None is
        returned when something is wrong with the file.
    """

    fnameraw = oscmap.rawcore[f"{map_.name}.charges_filename"][0]
    fname = (oscmap.directory / fnameraw).resolve()
    if GM_fh.try_file(fname) is None:
        GM_pt.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file doesn't exist:\n"
            f"{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TrEsp_1", False
        )
        map_.success = False
        return None

    if not GM_fh.check_file_readability(fname, False, False):
        GM_pt.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file is of the wrong format:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TrEsp_2", False
        )
        map_.success = False
        return None
    return fname


def gc_get_file_contents(fname, map_, oscmap):
    """Helper function for get_charges. Gets contents of oscmaps'TrEsp
    file.

    Parameters
    ----------
    fname : pathlib.Path
        The path to the file that stores the TrEsp charges.
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TrEsp map object
    oscmap: :class:`~GMAP.src.tools.map_reader.SingleMap`
        The map object of the map for which we'd like to obtain TrEsp
        charges.

    Returns
    -------
    contents : `np.ndarray` or None
        A numpy array with the TrEsp charges from the file. Returns
        None if there was some issue with the contents of the file.
        The array has no certain datatype, that is enforced by the
        TrEsp map in a later stage.
    """

    contents = []
    with open(fname, encoding="utf-8") as fhand:
        for line in fhand:
            line = GM_pp.cleanline(line).strip()
            if line:
                contents.append(line)

    try:
        contents = [float(item) for item in contents]
    except Exception:
        GM_pt.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file has the wrong contents:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TrEsp_3", False
        )
        map_.success = False
        return None

    if len(contents) > len(oscmap.core.used_atoms):
        GM_pt.Printer.warning(
            f"\nThe map {oscmap.name} provided the following file to the "
            f"{map_.name} coupling map, but that file has too many contents:"
            f"\n{fname}\nPlease make sure the map is installed correctly. If "
            "the problem persists, please contact the author of the "
            f"{oscmap.name} map.",
            "map_TrEsp_4", False
        )
        map_.success = False
        return None

    contents = np.array(contents)
    return contents


def gc_get_multiplier(keyword, map_, oscmap):
    """Helper function for get_charges. Gets contents of oscmaps'TrEsp
    file.

    Parameters
    ----------
    keyword : str
        The name of the parameter used in oscmap's core.txt file to
        store the multiplication factor.
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TrEsp map object
    oscmap: :class:`~GMAP.src.tools.map_reader.SingleMap`
        The map object of the map for which we'd like to obtain TrEsp
        charges.

    Returns
    -------
    multiplier : float
        The number by which to multiply the TrEsp charges from the
        TrEsp charges file before use.
    """

    cmdstr = "multiplier = " + " ".join(oscmap.rawcore[keyword])
    mapdir = oscmap.directory
    pars = {}
    try:
        exec(cmdstr, globals(), pars)
    except Exception:
        GM_pt.Printer.warning(
            "\nCould not interpret the choice for the keyword "
            f"'{keyword}' in the file {mapdir / 'core.txt'}. "
            "Please make sure the choice only contains numbers (and "
            "optionally a single '.') that represent a decimal value. "
            "Alternatively, make sure it is a python-parsable string. ",
            "map_TrEsp_5", False
        )
        map_.success = False
        return None

    try:
        multiplier = float(pars["multiplier"])
    except Exception:
        GM_pt.Printer.warning(
            "\nCould not interpret the choice for the keyword "
            f"'multiply_freq' in the file {mapdir / 'core.txt'}. "
            "Please make sure the choice only contains numbers (and "
            "optionally a single '.') that represent a decimal value. "
            "Alternatively, make sure it is a python-parsable string. ",
            "map_TrEsp_5", False
        )
        map_.success = False
        return None

    return multiplier
