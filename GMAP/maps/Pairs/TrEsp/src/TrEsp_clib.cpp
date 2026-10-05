#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // calloc!

// These files are not from default libraries, or within this map folder.
// Instead, they should be manually included during installation:

// windows (my machine, edit path):
// cl.exe /LD /Fe: TrEsp_clib_Win64bit /I\github\GEMAIM-dev\GMAP\sourcefiles TrEsp_clib.cpp

// linux (Kai's cluster, edit path):
// g++ -fPIC -shared -o TrEsp_clib_Linux.so -I/scratch/p302934/GMAP_fin/GMAP/GMAP/sourcefiles TrEsp_clib.cpp
#include "vectormath.cpp"  // in GMAP sourcefiles directory

#ifdef _WIN32
    extern "C" {
        __declspec(dllexport) void calc_coupling(int npairs, int *allpairs, int *noscats, int *oscstart, float *diff_charges, int *osc_used_ats, float *positions_box, float *boxvects, int totosc,float fpieps,
        float *hamiltonian);
    }
#endif


extern "C" {
    
    void calc_coupling(
        int npairs, int *allpairs, int *noscats, int *oscstart,
        float *diff_charges, int *osc_used_ats, float *positions_box,
        float *boxvects, int totosc, float fpieps, float *hamiltonian
    ) {
        /*
        TrEsp multiplies the found J by 116141.70590152.
        Why, what is this number?

        1/4pieps = 8.9875517862(14) e9 Nm^2C^-2 (coulombs constant)
        e = 1.602176634 * e-19 C  (elementary charge)

        newton = kg m s-2
        1/4pieps = 8.9875517862(14) e9 kg m^3s-2C^-2 (coulombs constant)

        e^2/4pieps = 23.070775507894 e-29 kg m^3 s^-2
                   = 23.070775507894 e-29 Jm

        h = 6.62607015 e-34 Js
        c = 2.99792458 e8 m s-1
        hc = 19.864458571489287 e-26  Jm
        J = 1/hc = 0.05034116567542709302218551285258 e-26 m-1
        J = 1/hc = (1/1.98644586) * 10^25 1/m
        J = 5.03411656 * 10^22 1/cm (joule in wavenumbers)

        e^2/4pieps = 23.0707755078585 e-29 Jm
                   = 23.0707755078585 e-29 5.03411656 e22 m cm-1
                   = 116.14097303 e-7 m cm-1
        
        using ang: 116.140973210 e3 ang cm-1
                 = 116140.973210 ang cm-1
        */

        int pairix, oscix1, oscix2, ix1, ix2, osc1len, osc2len, TRix1, TRix2;
        float tempvec[3], diff[3], r2, ir, J;
        for (pairix = 0; pairix < npairs; pairix++) {
            // oscillator indices
            oscix1 = allpairs[pairix * 2];
            oscix2 = allpairs[pairix * 2 + 1];
            
            J = 0;

            // amount of atoms per oscilator
            osc1len = noscats[oscix1];
            osc2len = noscats[oscix2];

            // index in all-atom TrEsp arrays
            TRix1 = oscstart[oscix1];
            for (ix1 = 0; ix1 < osc1len; ix1++) {
                TRix2 = oscstart[oscix2];
                for (ix2 = 0; ix2 < osc2len; ix2++) {
                    // calculate difference vector between the two atoms
                    VM_PBC_diff_mod1(
                        &positions_box[osc_used_ats[TRix1] * 3],
                        &positions_box[osc_used_ats[TRix2] * 3],
                        // 3rd and 4th argument aren't used; last is.
                        diff, diff, tempvec);
                    VM_vect_at_matrix33(tempvec, boxvects, diff);
                    r2 = VM_veclen2(diff);
                    ir = r2 < 0.01 ? 1 : 1/sqrt(r2);

                    J += diff_charges[TRix1] * diff_charges[TRix2] * ir;
                    TRix2++;
                }
                TRix1++;
            }
            J *= fpieps;
            hamiltonian[oscix1 * totosc + oscix2] = J;
            hamiltonian[oscix2 * totosc + oscix1] = J;
        }
    }
}
