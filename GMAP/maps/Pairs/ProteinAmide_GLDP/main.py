
# 3rd party lib imports
# from numba import njit
import numpy as np

# gmap imports
import GMAP.src.tools.constants as GM_con
import GMAP.src.tools.math_functions as GM_mf
import GMAP.src.tools.print_tools as GM_pt


class NeighborMap:
    """Contains a single neighbormapping

    This class has been borrowed from the AmideBB map, which has similar
    maps in use (for C and N terminal neighbour influences), with one
    notable change: The application of phi and psi is reverted (as this
    map is saved into file transposed compared to the AmideBB map)

    This mapping is a grid. When provided with two values (one for each
    direction), a bilinearly interpolated value is obtained.
    After making an instance of this class, the main use is for the
    instance to have the function get_delta called (which is in practice
    the application of the map).

    Parameters
    ----------
    fname : str or `pathlib.Path`
        the name of the file that contains the map information
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.

    Attributes
    ----------
    dim : int
        The amount of rows and columns to expect from the file
    space : float or int
        The amount of degrees between two columns/rows. Should equal
        360/(dim-1)
    data : `np.ndarray`
        An array of shape (dim, dim) that contains the actual data
        points.
    """

    def __init__(self, fname, map_):
        self.dim = 13
        self.space = 30
        try:
            self.data = np.loadtxt(fname)
        except Exception as ex:
            GM_pt.Printer.warning(
                "\nCould not interpret the data "
                f" in the file {fname}. Please make sure "
                "the file was not changed since downloading, contains "
                "13 lines of 13 numbers, and no characters that cannot be "
                "interpreted as floats.",
                "map_ProteinAmide_GLDP_1", exception=ex
            )
            map_.success = False
            return
        if self.data.shape != (13, 13):
            GM_pt.Printer.warning(
                "\nCould not interpret the data "
                f" in the file {fname}. Please make sure "
                "the file was not changed since downloading, and contains "
                "13 lines of 13 numbers.",
                "map_ProteinAmide_GLDP_1"
            )
            map_.success = False

    def process_angle(self, angle):
        """Takes an angle in radians, and converts in to degrees, and
        also the column corresponding to the map.

        Parameters
        ----------
        angle : float
            The angle (in radians) to convert

        Returns
        -------
        angle : float
            The same angle, but in degrees
        angle_N : int
            The column that forms the 'left' boundary for interpolation.
            Found through int((angle + 180) // self.space)
        """

        angle *= GM_con.rad2deg
        angle_N = int((round(angle, 4) + 180) // self.space)
        if angle_N == (self.dim - 1):
            angle_N = self.dim - 2
        return angle, angle_N

    def get_coupling(self, Nosc, Cosc, system):
        """Calculate the coupling for this pair of oscillators.

        Note: the Tasumi map is transposed compared to the AmidBB maps.

        Parameters
        ----------
        Nosc, Cosc : :class:`~GMAP.src.tools.system_reader.Oscillator`
            The oscillators for which the shift should be calculated.
        system : :class:`~GMAP.src.tools.system_reader.System`
            The object that stores everything the program currently
            knows about the system being treated (names, numbers, types,
            masses, charges of all atoms, for example)

        Returns
        -------
        delta : float
            The amount that should be added to the frequency to
            compensate for interactions between these two oscillators.
        """

        # calculate the ramachandran angles
        Nbpos = Nosc.positions_box
        Cbpos = Cosc.positions_box
        phi_ang = GM_mf.dihedral_boxcoords(  # around NCA, angle of two Cs
            Nbpos[0], Nbpos[3], Nbpos[5], Cbpos[0], system.boxvects)
        phi_ang, phi_N = self.process_angle(phi_ang)
        psi_ang = GM_mf.dihedral_boxcoords(  # around CAC, angle of two Ns
            Nbpos[3], Nbpos[5], Cbpos[0], Cbpos[3], system.boxvects)
        psi_ang, psi_N = self.process_angle(psi_ang)

        if (0 <= phi_N < (self.dim - 1)) and (0 <= psi_N < (self.dim - 1)):
            y1 = self.data[psi_N, phi_N]
            y2 = self.data[psi_N+1, phi_N]
            y3 = self.data[psi_N+1, phi_N+1]
            y4 = self.data[psi_N, phi_N+1]

            u = (phi_ang % self.space) / self.space
            t = (psi_ang % self.space) / self.space

            # bilinear interpolation!
            J = (1-u)*(1-t)*y1 + (1-u)*t*y2 + u*t*y3 + u*(1-t)*y4

        else:
            GM_pt.Printer.warning(
                "Ill defined ramachandran angles found for residue "
                f"{Nosc.resnames[1]}{Nosc.resnums[1]}. The nearest "
                "neighbour shift will be set to zero.",
                "map_ProteinAmide_GLDP_2"
            )
            J = 0
        return J


def GM_calc_coupling(map_, system, hamiltonian):
    """Calculate all the couplings that should be determined by this map

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    hamiltonian : `np.ndarray`
        The hamiltonian of the full system. Consists of float32, has a
        column and a row for each oscillator.
    """

    for pair in map_.allpairs:
        J = map_.map_per_pair[tuple(pair)].get_coupling(
            system.oscillators[pair[0]],
            system.oscillators[pair[1]], system)
        hamiltonian[pair[0], pair[1]] = J
        hamiltonian[pair[1], pair[0]] = J


def GM_pre_run(map_, system):
    """determine what coupling needs which GLDP mapping.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everyting the program currently knows
        about the MD system.
    """

    map_.map_per_pair = {}
    for pair in map_.allpairs:
        map_.map_per_pair[tuple(pair)] = map_.neighbormaps[
            "GLDP" + determine_map(pair, map_, system)]


def determine_map(pair, map_, system):
    """
    For the provided pair, find out what map should be used to consider its
    N-term and C-term neighbor.

    The determined maps are saved as an attribute of the oscillators.

    Parameters
    ----------
    pair : Tuple of 2 int
        The oscillator indices of the two oscillators making up this
        pair.
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """
    osc1 = system.oscillators[pair[0]]
    osc2 = system.oscillators[pair[1]]
    if osc1.used_atoms[5] == osc2.used_atoms[2]:
        Nosc = osc1
        Cosc = osc2
    else:
        Nosc = osc2
        Cosc = osc1

    # no (pre-)prolines! easy!
    if "PRO" not in Cosc.resnames:
        return ""

    # for now, the pro-pro case is treated as if it is pro-gly. in the future,
    # a pro-pro map should be made!

    boxpos = Nosc.positions_box
    COvec = GM_mf.PBC_boxdiff_triclin(boxpos[1], boxpos[0], system.boxvects)
    NHvec = GM_mf.PBC_boxdiff_triclin(boxpos[4], boxpos[3], system.boxvects)

    # (both pro-pro(should for now be treated as pro-gly) and pro-gly)
    if Nosc.resnames[1] == "PRO":
        # In the original code, Pro-Pro is actually treated as Gly-Pro
        if (
            Cosc.resnames[1] == "PRO"
            and map_.run_pars.legacy_mode == "AmideImaps"
        ):
            bondtype = "GP"
        else:
            bondtype = "PG"
    else:
        bondtype = "GP"

    if bondtype == "GP":
        if GM_mf.dotprod(COvec, NHvec) < 0:
            return "_transGly_transPro"
        else:
            return "_cisGly_transPro"

    LorD = DLcheck(Nosc, Cosc, map_, system)  # -1 for D, 1 for L, 0 for nodir
    if GM_mf.dotprod(COvec, NHvec) < 0:
        if LorD < 0:
            return "_transDPro_transGly"
        else:
            return "_transPro_transGly"
    elif LorD < 0:
        return "_cisDPro_transGly"
    else:
        return "_cisPro_transGly"


def GM_post_init(map_, system):
    """Do some final initializations that need to happen before the
    calculation starts.

    Steps present:
    - Build the seven GLDP maps
    - Copy out the legacy_mode setting from master map.

    Parameters
    ----------
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)
    """

    # read in all parameters from the constants file
    mapdir = map_.directory / "maps"
    map_.neighbormaps = {
        fname.stem: NeighborMap(fname, map_) for fname in mapdir.iterdir()}

    # Copy the legacy_mode setting from the 'main' map to this one.
    rps = map_.run_pars
    main_runpars = rps.main_run_pars
    protam_rps = main_runpars.requested_pairmapdict["ProteinAmide"].run_pars
    rps.legacy_mode = protam_rps.legacy_mode


def DLcheck(osc1, osc2, map_, system):
    """Determine of an amino acid is of D or L chirality.

    Taken from AmideBB map.

    Parameters
    ----------
    osc1, osc2 : :class:`~GMAP.src.tools.system_reader.Oscillator`
        The oscillators surrounding the amino acid in question
    map_ : :class:`~GMAP.src.tools.map_reader.Map`
        The object that stores everything the program currently knows
        about this map.
    system : :class:`~GMAP.src.tools.system_reader.System`
        The object that stores everything the program currently knows
        about the system being treated (names, numbers, types, masses,
        charges of all atoms, for example)

    Returns
    -------
    chirality : int
        Whether the amino acid is in D or L configuration.
        -1 is returned for D (special), 1 is returned for L (default for
        most organisms), 0 is returned for groups where the CA atom
        is not chiral (like glycine for example)
    """

    # checks whether the amino acid between two oscillators is in L or D
    # configuration

    if osc2.resnames[0] in ("GLY", "FOR", "ETA", "GL2"):
        return 0

    atomCA = osc2.positions_box[2]
    atomC = osc2.positions_box[0]
    atomN = osc1.positions_box[3]

    first_at = system.residues.first_ix[osc2.resnums[0]]
    last_at = system.residues.last_ix[osc2.resnums[0]]
    for atix, atname in zip(
        system.atnums[first_at:last_at+1],
        system.atnames[first_at:last_at+1]
    ):
        if atname == "CB":
            atomCBix = atix
            break
    else:
        GM_pt.Printer.warning(
            "Warning! The residue between the following two oscillators "
            "does not have a CB atom, and thus its chirality cannot be "
            f"determined:\n{osc1}\n{osc2}\nPlease make sure you're applying "
            "the correct map to the correct system. If this map erraneously "
            "detects something it shouldn't, please contact the developers!",
            "map_ProteinAmide_GLDP_3"
        )
        map_.success = False
        return 0
    shift = atomC
    atomCB = system.positions[atomCBix] @ system.boxvects_inv - shift
    atomCB -= np.floor(atomCB + 0.5) - shift

    # used ats order:    res0{C O CA} res1{N H CA} ({N CD CA} for prepro)
    # osc1 is first, osc2 is last
    CACvec = GM_mf.PBC_boxdiff_triclin(atomC, atomCA, system.boxvects)
    CANvec = GM_mf.PBC_boxdiff_triclin(atomN, atomCA, system.boxvects)
    CACBvec = GM_mf.PBC_boxdiff_triclin(atomCB, atomCA, system.boxvects)

    CxN = GM_mf.crossprod(CACvec, CANvec)
    if GM_mf.dotprod(CxN, CACBvec) > 0:
        return -1
    else:
        return 1
