
# standard lib imports
import ctypes as ct

# gmap imports
import GMAP.src.tools.coding_tools as GM_ct
import GMAP.src.tools.file_handler as GM_fh


class TCC_Clib(metaclass=GM_ct.Singleton):
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
        The TCC map object in the code - this will be the object that
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

        self.clib.prep_coupling.restype = None
        self.clib.prep_coupling.argtypes = [
            ct.c_int,  # nosc
            ct.c_int,  # noscats
            ct.POINTER(ct.c_int),  # oscixarr
            ct.POINTER(ct.c_int),  # osc_used_ats,
            ct.POINTER(ct.c_float),  # positions_box
            ct.POINTER(ct.c_float),  # boxvects
            ct.POINTER(ct.c_int),  # dopro
            ct.c_float,  # alpha_gen
            ct.c_float,  # alpha_pro
            ct.POINTER(ct.c_float),  # v_gen
            ct.POINTER(ct.c_float),  # v_pro
            ct.POINTER(ct.c_float)  # tcc_v
        ]

        self.clib.calc_coupling.restype = None
        self.clib.calc_coupling.argtypes = [
            ct.c_int,  # npairs
            ct.POINTER(ct.c_int),  # allpairs
            ct.c_int,  # noscats
            ct.POINTER(ct.c_int),  # oscix_to_ix
            ct.POINTER(ct.c_int),  # osc_used_ats
            ct.POINTER(ct.c_float),  # positions_box
            ct.POINTER(ct.c_float),  # boxvects
            ct.POINTER(ct.c_int),  # dopro
            ct.POINTER(ct.c_float),  # q_gen
            ct.POINTER(ct.c_float),  # q_pro
            ct.POINTER(ct.c_float),  # dq_gen
            ct.POINTER(ct.c_float),  # dq_pro
            ct.c_float,  # fourPiEps
            ct.POINTER(ct.c_float),  # tcc_v
            ct.c_int,  # totosc
            ct.POINTER(ct.c_float)  # hamiltonian
        ]

    def prep_coupling(self, map_, system):
        """Prepare all couplings for this map, this frame.

        This is done by calling the respective c function. It loops over
        all singles that will be treated by this map, and calculate
        the properties of those singles that do not depend on the
        coupling partner.

        Parameters
        ----------
        map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
            The TrEsp map object in the code - this will be the object
            that gains the new c-library attribute.
        system : :class:`~GMAP.src.tools.system_reader.System`
            The object that stores everything the program currently
            knows about the system being treated (names, numbers, types,
            masses, charges of all atoms, for example)
        """

        self.clib.prep_coupling(
            map_.nosc,  # nosc
            map_.noscats,  # noscats
            map_.oscixlist_c,  # oscixarr
            map_.all_used_ats_c,  # osc_used_ats
            system.positions_box_c,  # positions_box
            system.boxvects_c,  # boxvects
            map_.dopro_c,  # dopro
            map_.alpha_gen,  # alpha_gen
            map_.alpha_pro,  # alpha_pro
            map_.v_gen_c,  # v_gen
            map_.v_pro_c,  # v_pro
            map_.map_tcc_v_c  # tcc_v
        )

    def calc_coupling(self, map_, system, hamiltonian_c):
        """Calculate all couplings for this map, this frame.

        This is done by calling the respective c function. It loops over
        all requested pairs and add them to the hamiltonian.

        Parameters
        ----------
        map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
            The TCC map object in the code - this will be the object
            that gains the new c-library attribute.
        system : :class:`~GMAP.src.tools.system_reader.System`
            The object that stores everything the program currently
            knows about the system being treated (names, numbers, types,
            masses, charges of all atoms, for example)
        hamiltonian_c : `c_float_array_XX`
            The c-pointer to the hamiltonian. 'XX' in the type is
            variable, as it depends on the amount of singles in the
            calculation.
        """

        self.clib.calc_coupling(
            map_.n_allpairs,  # npairs
            map_.allpairs_c,  # allpairs
            map_.noscats,  # noscats
            map_.oscix_to_ix_c,  # oscix_to_ix
            map_.all_used_ats_c,  # osc_used_ats
            system.positions_box_c,  # positions_box
            system.boxvects_c,  # boxvects
            map_.dopro_c,  # dopro
            map_.q_gen_c,  # q_gen
            map_.q_pro_c,  # q_pro
            map_.dq_gen_c,  # dq_gen
            map_.dq_pro_c,  # dq_pro
            map_.fourPiEps,  # fourPiEps
            map_.map_tcc_v_c,  # tcc_v
            system.nosc,  # totosc
            hamiltonian_c  # hamiltonian
        )


def init_map_for_clib(map_, system):
    """Load the c library into the map.

    There are no warnings here - if the c library file is not present,
    the user will be presented with an 'ugly' python error.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.PairMap`
        The TCC map object in the code - this will be the object that
        gains the new c-library attribute.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """

    map_.clibfile = map_.directory / "src"
    map_.clibfile /= "TCC_clib" + GM_fh.FileLocations.clib_extension
    map_.clib = TCC_Clib(map_)
