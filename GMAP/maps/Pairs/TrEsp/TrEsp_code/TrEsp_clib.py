
# standard lib imports
import ctypes as ct

# gmap imports
import GMAP.src.tools.coding_tools as GM_ct
import GMAP.src.tools.file_handler as GM_fh


class TrEsp_Clib(metaclass=GM_ct.Singleton):
    """Stores and manages all c functions for this map.

    Each (external) function in the library has it's own associated
    method on this class. The calls to C are ugly and convoluted due
    to c functions needing so many parameters (either single values
    or numpy arrays), so these methods make their calls more pythonic.
    They each require just the relevant classes, and unpack the required
    attributes themselves.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TrEsp map object in the code - this will be the object that
        gains the new c-library attribute.

    Notes
    -----
    C functions cannot return more than a single value. Therefore, most
    will write their output into an array provided as an input. This
    array must be of a c-friendly datatype, use
    ``np.ctypeslib.as_ctypes()``
    for creating these. This version must be saved along with the
    original, so, for example, you have both ``my_arr`` and
    ``my_arr_c``, where ``my_arr_c`` is defined as
    ``np.ctypeslib.as_ctypes(my_arr)``. This results in two views of the
    same array, meaning that any changes to any values made in
    ``my_arr`` will also apply to ``my_arr_c``, and vice versa.

    Attributes
    ----------
    clib : `ctypes.CDLL`
        The actual compiled c-code. Must be compiled to be a library,
        so a .dll (windows), .so (linux) or .dylib (macOS) file.
    """

    def __init__(self, map_):
        self.clib = ct.CDLL(str(map_.clibfile))

        self.clib.calc_coupling.restype = None
        self.clib.calc_coupling.argtypes = [
            ct.c_int,  # npairs
            ct.POINTER(ct.c_int),  # allpairs
            ct.POINTER(ct.c_int),  # noscats
            ct.POINTER(ct.c_int),  # oscstart
            ct.POINTER(ct.c_float),  # diff_charges
            ct.POINTER(ct.c_int),  # osc_used_ats
            ct.POINTER(ct.c_float),  # positions_box
            ct.POINTER(ct.c_float),  # boxvects
            ct.c_int,  # totosc
            ct.c_float,  # fpieps
            ct.POINTER(ct.c_float)  # hamiltonian
        ]

    def calc_coupling(self, map_, system, fpieps, hamiltonian_c):
        """Calculate all couplings for this map, this frame.

        This is done by calling the respective c function. It loops over
        all requested pairs and add them to the hamiltonian.

        Parameters
        ----------
        map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
            The TrEsp map object in the code - this will be the object
            that gains the new c-library attribute.
        system : :class:`~GMAP.src.tools.system_reader.System`
            The object that stores everything the program currently
            knows about the system being treated (names, numbers, types,
            masses, charges of all atoms, for example)
        fpieps : `np.float32`
            The value for 4 pi epsilon, as a float32.
        hamiltonian_c : `c_float_array_XX`
            The c-pointer to the hamiltonian. 'XX' in the type is
            variable, as it depends on the amount of singles in the
            calculation.
        """

        self.clib.calc_coupling(
            map_.n_allpairs,  # npairs
            map_.allpairs_c,  # allpairs
            map_.osclens_c,  # noscats
            map_.oscstart_c,  # oscstart
            map_.charge_array_c,  # diff_charges
            map_.all_used_ats_c,  # osc_used_ats
            system.positions_box_c,  # positions_box
            system.boxvects_c,  # boxvects
            system.nosc,  # totosc
            fpieps,
            hamiltonian_c  # hamiltonian
        )


def init_map_for_clib(map_, system):
    """Load the c library into the map.

    There are no warnings here - if the c library file is not present,
    the user will be presented with an 'ugly' python error.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TrEsp map object in the code - this will be the object that
        gains the new c-library attribute.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """

    map_.clibfile = map_.directory / "src"
    map_.clibfile /= "TrEsp_clib" + GM_fh.FileLocations.clib_extension
    map_.clib = TrEsp_Clib(map_)
