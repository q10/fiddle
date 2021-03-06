"""
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
"""
import warnings
from Utils.Exceptions import PDBATM2Warning
from numpy import array, radians, sin, cos, sqrt

def pdbatm2(atom, atom_a, atom_b, atom_c, bond_length, angle1deg, angle2deg, chiral_flag, debug_flag=True):
    assert(isinstance(debug_flag, bool))

#     convert angles to radians, and get their sines and cosines
    eps = 0.00000001
    angle1rad, angle2rad = radians(angle1deg), radians(angle2deg)
    sin1, cos1, sin2, cos2 = sin(angle1rad), cos(angle1rad), sin(angle2rad), cos(angle2rad)

#     if no second site given, place the atom at the origin
    if atom_a is None:
        atom.set_coords(array([0.0, 0.0, 0.0], 'f'))

#     if no third site given, place the atom along the z-axis
    elif atom_b is None:
        coords = atom_a.coords()
        coords[2] += bond_length
        atom.set_coords(coords)

#     if no fourth site given, place the atom in the x,z-plane
    elif atom_c is None:
        rab = atom_a.distance_from(atom_b)
        [xab, yab, zab] = (atom_a.coords() - atom_b.coords()) / rab
        cosb, sinb = zab, sqrt(xab**2 + yab**2)
        if sinb == 0.0:
            cosg, sing = 1.0, 0.0
        else:
            cosg = yab / sinb
            sing = xab / sinb
        xtmp, ztmp = bond_length*sin1, rab - bond_length*cos1
        coords = atom_b.coords() + array([xtmp*cosg + ztmp*sing*sinb, -xtmp*sing + ztmp*cosg*sinb, ztmp*cosb], 'f')
        atom.set_coords(coords)

#     general case where the second angle is a dihedral angle
    elif not chiral_flag:
        [xab, yab, zab] = (atom_a.coords() - atom_b.coords()) / atom_a.distance_from(atom_b)
        [xbc, ybc, zbc] = (atom_b.coords() - atom_c.coords()) / atom_b.distance_from(atom_c)
        cosine = xab*xbc + yab*ybc + zab*zbc
        sine = sqrt(max(1.0-cosine**2,eps))
        
        if abs(cosine) >= 1.0:
            warnings.warn("%s: Undefined dihedral angle at atom %i: %s" % PDBATM2Warning, atom.serial_number(), atom)

        [xt, yt, zt] = array([zab*ybc - yab*zbc, xab*zbc - zab*xbc, yab*xbc - xab*ybc], 'f') / sine
        [xu, yu, zu] = yt*zab - zt*yab, zt*xab - xt*zab, xt*yab - yt*xab

        coords = atom_a.coords()
        coords += array([xu*sin1*cos2+xt*sin1*sin2-xab*cos1, yu*sin1*cos2+yt*sin1*sin2-yab*cos1, zu*sin1*cos2+zt*sin1*sin2-zab*cos1], 'f') * bond_length
        atom.set_coords(coords)

#     general case where the second angle is a bond angle
    elif abs(chiral_flag) == 1:
        [xba, yba, zba] = (atom_b.coords() - atom_a.coords()) / atom_b.distance_from(atom_a)
        [xac, yac, zac] = (atom_a.coords() - atom_c.coords()) / atom_a.distance_from(atom_c)
        [xt, yt, zt] = zba*yac - yba*zac, xba*zac - zba*xac, yba*xac - xba*yac

        cosine = xba*xac + yba*yac + zba*zac
        sine2 = max(1.0-cosine**2, eps)

        if abs(cosine) >= 1.0:
            warnings.warn("%s: Defining atoms colinear at atom %i: %s" % PDBATM2Warning, atom.serial_number(), atom)

        [a, b] = -cos2 - cosine*cos1, cos1 + cosine*cos2
        [a, b, c] = array([a, b, 1.0 + a*cos2 - b*cos1], 'f') / sine2

        if c > eps:
            c = chiral_flag * sqrt(c)
        elif c < -eps:
            [a, b, c] = array([a, b, 0.0], 'f') / sqrt((a*xac+b*xba)**2 + (a*yac+b*yba)**2 + (a*zac+b*zba)**2)
            if debug_flag:
                warnings.warn("%s: Sum of bond angles too large at atom %i: %s" % PDBATM2Warning, atom.serial_number(), atom)
        else:
            c = 0.0

        coords = atom_a.coords()
        coords += array([a*xac + b*xba + c*xt, a*yac + b*yba + c*yt, a*zac + b*zba + c*zt], 'f') * bond_length
        atom.set_coords(coords)
    return