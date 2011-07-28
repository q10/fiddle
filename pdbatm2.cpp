/*
     #################################################################
     ##                                                             ##
     ##  subroutine pdbatm2  --  single atom internal to Cartesian  ##
     ##                                                             ##
     #################################################################

    Parameters
    4 atoms - the atom where we want to set the coords (atom), and the three atoms nearby (atom_a, atom_b, atom_c)
    bond_length - length between atom and atom_a
    angle1deg - desired atom-atom_a-atom_b angle.  Units in degrees
    angle2deg - can either be the dihedral angle or the atom_b-angle_a-angle_c angle.  Units in degrees
    chiral_flag - if 0, then we are trying to fit atom (labeled A in diagram) to
        A---aa---ab---ac

        Else, it is +1 or -1, in which case we are trying to fit atom to model, with left-hand or right-hand chirality
                 ab
                /
               /
      A------aa
              \
               \
                ac
*/

#include "common.h"

void pdbatm2(Atom * atom, Atom * atom_a, Atom * atom_b, Atom * atom_c, double bond_length, double angle1deg, double angle2deg, int chiral_flag) {
    // convert angles to radians, and get their sines and cosines
    double eps = 0.00000001;
    double angle1rad = DEG2RAD(angle1deg), angle2rad = DEG2RAD(angle2deg);
    double sin1 = sin(angle1rad), cos1 = cos(angle1rad), sin2 = sin(angle2rad), cos2 = cos(angle2rad);
    double * new_coords = new double[3];

    // if no second site given, place the atom at the origin
    if (atom_a == NULL) {
        for (int i = 0; i < 3; i++)
            new_coords[i] = 0.0;
        atom->set_coords(new_coords);
    }// if no third site given, place the atom along the z-axis
    else if (atom_b == NULL) {
        for (int i = 0; i < 3; i++)
            new_coords[i] = atom_a->coords[i];
        new_coords[2] += bond_length;
        atom->set_coords(new_coords);
    }//if no fourth site given, place the atom in the x,z-plane
    else if (atom_c == NULL) {
        double rab = atom_a->distance_from(atom_b);
        double xab = (atom_a->coords[0] - atom_b->coords[0]) / rab;
        double yab = (atom_a->coords[1] - atom_b->coords[1]) / rab;
        double zab = (atom_a->coords[2] - atom_b->coords[2]) / rab;
        double cosb = zab, sinb = sqrt(xab*xab + yab*yab), cosg, sing, xtmp, ztmp;

        if (sinb == 0.0) {
            cosg = 1.0; 
            sing = 0.0;            
        } else {
            cosg = yab / sinb;
            sing = xab / sinb;
        }
        xtmp = bond_length * sin1;
        ztmp = rab - bond_length * cos1;

        for (int i = 0; i < 3; i++)
            new_coords[i] = atom_b->coords[i];

        new_coords[0] += xtmp * cosg + ztmp * sing * sinb;
        new_coords[1] += -xtmp * sing + ztmp * cosg * sinb;
        new_coords[2] += ztmp * cosb;
        atom->set_coords(new_coords);
    }//general case where the second angle is a dihedral angle
    else if (chiral_flag == 0) {
        double rab = atom_a->distance_from(atom_b), rbc = atom_b->distance_from(atom_c);
        double xab = (atom_a->coords[0] - atom_b->coords[0]) / rab;
        double yab = (atom_a->coords[1] - atom_b->coords[1]) / rab;
        double zab = (atom_a->coords[2] - atom_b->coords[2]) / rab;
        
        double xbc = (atom_b->coords[0] - atom_c->coords[0]) / rbc;
        double ybc = (atom_b->coords[1] - atom_c->coords[1]) / rbc;
        double zbc = (atom_b->coords[2] - atom_c->coords[2]) / rbc;
        
        double cosine = xab * xbc + yab * ybc + zab * zbc;
        double sine = sqrt(max(1.0 - cosine*cosine, eps));

        if (abs(cosine) >= 1.0)
            cerr << "Warning: Undefined dihedral angle at " << atom << endl;

        double xt = zab * ybc - yab * zbc / sine;
        double yt = xab * zbc - zab * xbc / sine;
        double zt = yab * xbc - xab * ybc / sine;

        double xu = yt * zab - zt * yab;
        double yu = zt * xab - xt * zab;
        double zu = xt * yab - yt * xab;
        
        for (int i = 0; i < 3; i++)
            new_coords[i] = atom_a->coords[i];

        new_coords[0] += (xu * sin1 * cos2 + xt * sin1 * sin2 - xab * cos1) * bond_length;
        new_coords[1] += (yu * sin1 * cos2 + yt * sin1 * sin2 - yab * cos1) * bond_length;
        new_coords[2] += (zu * sin1 * cos2 + zt * sin1 * sin2 - zab * cos1) * bond_length;        
        atom->set_coords(new_coords);
    }//general case where the second angle is a bond angle
    else if (abs(chiral_flag) == 1) {
        double rba = atom_b->distance_from(atom_a), rac = atom_a->distance_from(atom_c);

        double xba = (atom_b->coords[0] - atom_a->coords[0]) / rba;
        double yba = (atom_b->coords[1] - atom_a->coords[1]) / rba;
        double zba = (atom_b->coords[2] - atom_a->coords[2]) / rba;

        double xac = (atom_a->coords[0] - atom_c->coords[0]) / rac;
        double yac = (atom_a->coords[1] - atom_c->coords[1]) / rac;
        double zac = (atom_a->coords[2] - atom_c->coords[2]) / rac;

        double xt = zba * yac - yba * zac;
        double yt = xba * zac - zba * xac;
        double zt = yba * xac - xba * yac;

        double cosine = xba * xac + yba * yac + zba * zac;
        double sine2 = max(1.0 - cosine*cosine, eps);

        if (abs(cosine) >= 1.0)
            cerr << "Warning: Defining atoms colinear at " << atom << endl;

        double a = (-cos2 - cosine * cos1) / sine2;
        double b = (cos1 + cosine * cos2) / sine2;
        double c = (1.0 + a * cos2 - b * cos1) / sine2;

        if (c > eps)
            c = chiral_flag * sqrt(c);
        else if (c < -eps) {
            c = sqrt(pow((a * xac + b * xba), 2.0) + pow((a * yac + b * yba), 2.0) + pow((a * zac + b * zba), 2.0));
            a /= c;
            b /= c;
            c = 0.0;
            // cerr << "Warning: Sum of bond angles too large at atom " << atom << endl;
        } else
            c = 0.0;

        for (int i = 0; i < 3; i++)
            new_coords[i] = atom_a->coords[i];
        new_coords[0] += (a * xac + b * xba + c * xt) * bond_length;
        new_coords[1] += (a * yac + b * yba + c * yt) * bond_length;
        new_coords[2] += (a * zac + b * zba + c * zt) * bond_length;
        atom->set_coords(new_coords);
    }
    delete new_coords;
    return;
}
