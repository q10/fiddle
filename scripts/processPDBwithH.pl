#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Math::Trig;
use Switch;

my $pdb = undef;
my $out = undef;
my $chain = "A";
my $altchain = "A";
my $ligand = undef;
my $nohydrogens = undef;
GetOptions("pdb=s" => \$pdb, "out=s" => \$out, "chain=s" => \$chain, "altchain=s" => \$altchain, "ligand=s" => \$ligand, "nohydrogens"  => \$nohydrogens);

die "Usage: $0 -pdb [input pdb] -out [output pdb]\n" unless (defined $pdb && defined $out);
$chain = uc($chain);
$altchain = uc($altchain);

my $natom = 0;
my @residue = qw( ALA ARG ASP ASN CYS GLN GLU GLY HIS ILE 
                  LEU LYS MET PHE PRO SER THR TRP TYR VAL 
                  HID HIE LYN 
                  MSE );
push (@residue, $ligand) if (defined $ligand);

my %resname = ();
my %resnumb = ();
my %crd = ();
my @ntres = my @ctres = ();

&readpdb($pdb);
#&checkgap;
&checkspecial;
open (OUT, "> $out") || die "$0: Unable to write to $out: $!\n";
&processpdb;
close OUT;

exit;


##  subroutine

####====####
sub readpdb {
    my ($pdb) = @_;
    my $atmnam = my $resnam = "";
    my $residx = -9999;
    my $resnum = 0;
    my $nt = my $ct = 1;
    open (PDB, $pdb) || die "$0: Unable to read from $pdb: $!\n";
    while (<PDB>) {
        my @ln = split;
        $ct = 1 if (/^TER/ || /^END/);
        if ($ct) {
            push (@ctres, $residx);
            $nt = 1;
            $ct = 0;
        }
        if (/^ATOM  / || /^HETATM/) {
            next unless (substr($_, 21, 1) eq $chain);
            next unless (substr($_, 16, 1) eq " " || substr($_, 16, 1) eq $altchain);
            $resnum ++ unless (&removespace(substr($_, 17, 3)) eq $resnam && 
                               &removespace(substr($_, 22, 4)) eq $residx);
            $atmnam = &removespace(substr($_, 12, 4));
            $resnam = &removespace(substr($_, 17, 3));
            $residx = &removespace(substr($_, 22, 4));
            my ($x, $y, $z) = (substr($_, 30, 8), substr($_, 38, 8), substr($_, 46, 8));
            my $match = 0;
            foreach my $name (@residue) {
                $match = 1 if ($resnam eq $name);
            }
            next unless ($match);
            $resname{$residx} = $resnam;
            $resnumb{$residx} = $resnum;
            $crd{$residx}->{$atmnam} = [$x, $y, $z];
            if ($nt) {
                push (@ntres, $residx);
                $nt = 0;
                $ct = 0;
            }
        }
    }
    close PDB;
}


####====####
sub checkgap {
    my $name = "CA";
    my $crdprev = undef;
    my $cutoff = 4.5;
    my $idxprev = undef;
    foreach my $idx (sort {$a <=> $b} keys %crd) {
        die "$name atom is not found in residue $idx.\n" unless (defined $crd{$idx}->{$name});
        if (defined $crdprev) {
            my $dist = 0;
            $dist += ($crd{$idx}->{$name}->[$_] - $crdprev->[$_]) ** 2 foreach (0 .. 2);
            die "$name-$name distance beyond cutoff $cutoff between residue $idxprev and $idx\n" if (sqrt($dist) > $cutoff);
        }
        $crdprev->[$_] = $crd{$idx}->{$name}->[$_] foreach (0 .. 2);
        $idxprev = $idx;
    }
}

####====####
sub checkspecial {
## disulfide bridge
    my $ncys = 0;
    my @cyscrd = ();
    my $cysname = "SG";
    my $cyscutoff = 2.5;
    foreach my $idx (sort {$a <=> $b} keys %crd) {
        if ($resname{$idx} eq "CYS") {
            $cyscrd[$ncys]->[$_] = $crd{$idx}->{$cysname}->[$_] foreach (0 .. 2);
            $cyscrd[$ncys]->[3] = $idx;
            $ncys ++;
        }
    }
    foreach my $n1 (0 .. $ncys-1) {
        foreach my $n2 (0 .. $ncys-1) {
            next if ($n1 >= $n2);
            my $dist = 0;
            $dist += ($cyscrd[$n1]->[$_] - $cyscrd[$n2]->[$_]) ** 2 foreach (0 .. 2);
            if (sqrt($dist) <= $cyscutoff) {
#               printf "%4d%4d\n", $cyscrd[$n1]->[3], $cyscrd[$n2]->[3], "\n";
                $resname{$cyscrd[$n1]->[3]} = "CYX";
                $resname{$cyscrd[$n2]->[3]} = "CYX";
            }
        }
    }
## histidine protonation state
=head
    my $nhis = 0;
    my @hiscrd1 = ();
    my @hiscrd2 = ();
    my $hisname1 = "ND1";
    my $hisname2 = "NE2";
    my $hiscutoff = 3.2;
    foreach my $idx (sort {$a <=> $b} keys %crd) {
        if ($resname{$idx} eq "HIS") {
            $hiscrd1[$nhis]->[$_] = $crd{$idx}->{$hisname1}->[$_] foreach (0 .. 2);
            $hiscrd1[$nhis]->[3] = $idx;
            $hiscrd2[$nhis]->[$_] = $crd{$idx}->{$hisname2}->[$_] foreach (0 .. 2);
            $hiscrd2[$nhis]->[3] = $idx;
            $nhis ++;
        }
    }
    foreach my $n1 (0 .. $nhis-1) {
        my $n3 = $hiscrd1[$n1]->[3];
        foreach my $n2 (sort {$a <=> $b} keys %crd) {
            next if ($hiscrd1[$n1]->[3] == $n2);
            my $dist11 = my $dist12 = my $dist21 = my $dist22 = 0;
            my $sat11 = my $sat12 = my $sat21 = my $sat22 = 0;
            $dist11 += ($hiscrd1[$n1]->[$_] - $crd{$n2}->{"N"}->[$_]) ** 2 foreach (0 .. 2);
            $dist12 += ($hiscrd1[$n1]->[$_] - $crd{$n2}->{"O"}->[$_]) ** 2 foreach (0 .. 2);
            $dist21 += ($hiscrd2[$n1]->[$_] - $crd{$n2}->{"N"}->[$_]) ** 2 foreach (0 .. 2);
            $dist22 += ($hiscrd2[$n1]->[$_] - $crd{$n2}->{"O"}->[$_]) ** 2 foreach (0 .. 2);
            $sat11 = 1 if (sqrt($dist11) <= $hiscutoff);
            $sat12 = 1 if (sqrt($dist12) <= $hiscutoff);
            $sat21 = 1 if (sqrt($dist21) <= $hiscutoff);
            $sat22 = 1 if (sqrt($dist22) <= $hiscutoff);
            next if (($sat11 + $sat12 + $sat21 + $sat22) == 0); # HIE            HID            HID            HIE
            printf "%4d%4d%8.2f%8.2f%8.2f%8.2f\n", $n3, $n2, sqrt($dist11), sqrt($dist12), sqrt($dist21), sqrt($dist22), "\n";
            die "Check His $n3 local environment.\n" if (($sat11 + $sat12 + $sat21 + $sat22) > 2);
            die "Check His $n3 local environment.\n" if (($sat11 && $sat12) || ($sat21 && $sat22));
            if ($sat11 || $sat22) {
                $resname{$n3} = "HIE";
            } elsif ($sat12 || $sat21) {
                $resname{$n3} = "HID";
            } elsif ($sat11 && $sat21) {
                $resname{$n3} = "HIS";
            } elsif ($sat12 && $sat22) {
                $resname{$n3} = "HIS";
            }
        }
    }
=cut
}

####====####
sub processpdb {
    my @bbatom = qw( N CA C O OXT );
    my @scatom = ();
    foreach my $idx (sort {$a <=> $b} keys %crd) {
        foreach my $name (@bbatom) {
            &printatom($idx, $name);
#           printf "%6d%6d%6s%6s\n", $idx, $resnumb{$idx}, $resname{$idx}, $name, "\n";
        }
        foreach my $n (@ntres) {
            if ($idx == $n) {
#               printf "%4s%6d\n", "NTER", $idx, "\n";
                if ($resname{$idx} eq "PRO") {
                    unless (defined $crd{$idx}->{"H1"} or $nohydrogens) {
                        &buildatom($idx,  "H1", $idx,   "N", $idx,  "CA",   $idx,   "C", 1.02, 109.5,   0.0,  0);
                        &printatom($idx, "H1");
                    }
                    unless (defined $crd{$idx}->{"H2"} or $nohydrogens) {
                        &buildatom($idx,  "H2", $idx,   "N", $idx,  "CA",   $idx,   "C", 1.02, 109.5, 240.0,  0);
                        &printatom($idx, "H2");
                    }
                } else {
                    unless (defined $crd{$idx}->{"H1"} or $nohydrogens) {
                        &buildatom($idx,  "H1", $idx,   "N", $idx,  "CA",   $idx,   "C", 1.02, 109.5, 180.0,  0);
                        &printatom($idx, "H1");
                    }
                    unless (defined $crd{$idx}->{"H2"} or $nohydrogens) {
                        &buildatom($idx,  "H2", $idx,   "N", $idx,  "CA",   $idx,  "H1", 1.02, 109.5, 108.0,  1);
                        &printatom($idx, "H2");
                    }
                    unless (defined $crd{$idx}->{"H3"} or $nohydrogens) {
                        &buildatom($idx,  "H3", $idx,   "N", $idx,  "CA",   $idx,  "H1", 1.02, 109.5, 108.0, -1);
                        &printatom($idx, "H3");
                    }
                }
            }
        }
        foreach my $n (@ctres) {
            if ($idx == $n) {
#               printf "%4s%6d\n", "CTER", $idx, "\n";
                unless (defined $crd{$idx}->{"OXT"}) {
                    &buildatom($idx, "OXT", $idx,   "C", $idx,  "CA",   $idx,   "O", 1.25, 117.0, 126.0,  1);
                    &printatom($idx, "OXT");
                }
            }
        }
        my $nt = 0;
        foreach my $n (@ntres) {
           $nt = 1 if ($idx == $n);
        }
        unless ($nt || $resname{$idx} eq "PRO") {
            unless (defined $crd{$idx}->{"H"} or $nohydrogens) {
                &buildatom($idx,   "H", $idx,   "N", $idx,  "CA", $idx-1,   "C", 1.02, 121.0, 118.0,  1);
                &printatom($idx, "H");
            }
        }
        if ($resname{$idx} eq "GLY") {
            unless (defined $crd{$idx}->{"HA3"} or $nohydrogens) {
                &buildatom($idx, "HA3", $idx,  "CA", $idx,   "N",   $idx,   "C", 1.11, 109.5, 107.9, -1);
                &printatom($idx, "HA3");
            }
        } else {
            unless (defined $crd{$idx}->{"HA"} or $nohydrogens) {
                &buildatom($idx,  "HA", $idx,  "CA", $idx,   "N",   $idx,   "C", 1.11, 109.5, 107.9, -1);
                &printatom($idx, "HA");
            }
        }
        switch ($resname{$idx}) {
            case ("ALA") { @scatom = $nohydrogens ? qw( CB ) : qw( CB HB1 HB2 HB3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"HB1"}) {
                               &buildatom($idx, "HB1", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.11, 109.4, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,   "HB1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,   "HB1", 1.11, 109.4, 109.4, -1);
                           }
                         }
            case ("ARG") { @scatom = $nohydrogens ? qw( CB CG CD NE CZ NH1 NH2 ) : qw( CB CG CD NE CZ NH1 NH2 HB2 HB3 HG2 HG3 HD2 HD3 HE 1HH1 2HH1 1HH2 2HH2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NE"}) {
                               &buildatom($idx,  "NE", $idx,  "CD", $idx,  "CG", $idx,    "CB", 1.45, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CZ"}) {
                               &buildatom($idx,  "CZ", $idx,  "NE", $idx,  "CD", $idx,    "CG", 1.35, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NH1"}) {
                               &buildatom($idx, "NH1", $idx,  "CZ", $idx,  "NE", $idx,    "CD", 1.35, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NH2"}) {
                               &buildatom($idx, "NH2", $idx,  "CZ", $idx,  "NE", $idx,   "NH1", 1.35, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx,  "CD", $idx,  "CG", $idx,    "NE", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HD3"}) {
                               &buildatom($idx, "HD3", $idx,  "CD", $idx,  "CG", $idx,    "NE", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HE"}) {
                               &buildatom($idx,  "HE", $idx,  "NE", $idx,  "CD", $idx,    "CZ", 1.02, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"1HH1"}) {
                               &buildatom($idx,"1HH1", $idx, "NH1", $idx,  "CZ", $idx,    "NE", 1.02, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HH1"}) {
                               &buildatom($idx,"2HH1", $idx, "NH1", $idx,  "CZ", $idx,  "1HH1", 1.02, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"1HH2"}) {
                               &buildatom($idx,"1HH2", $idx, "NH2", $idx,  "CZ", $idx,    "NE", 1.02, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HH2"}) {
                               &buildatom($idx,"2HH2", $idx, "NH2", $idx,  "CZ", $idx,  "1HH2", 1.02, 120.0, 120.0,  1);
                           }
                         }
            case ("ASP") { @scatom = $nohydrogens ? qw( CB CG OD1 OD2 ) : qw( CB CG OD1 OD2 HB2 HB3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,  "N",  $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.51, 107.8, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OD1"}) {
                               &buildatom($idx, "OD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.25, 117.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OD2"}) {
                               &buildatom($idx, "OD2", $idx,  "CG", $idx,  "CB", $idx,   "OD1", 1.25, 117.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 107.9,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 107.9, -1);
                           }
                         }
            case ("ASN") { @scatom = $nohydrogens ? qw( CB CG OD1 ND2 ) : qw( CB CG OD1 ND2 HB2 HB3 1HD2 2HD2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.51, 107.8, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OD1"}) {
                               &buildatom($idx, "OD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.22, 122.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"ND2"}) {
                               &buildatom($idx, "ND2", $idx,  "CG", $idx,  "CB", $idx,   "OD1", 1.34, 112.7, 124.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 107.9,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 107.9, -1);
                           }
                           unless (defined $crd{$idx}->{"1HD2"}) {
                               &buildatom($idx,"1HD2", $idx, "ND2", $idx,  "CG", $idx,    "CB", 1.02, 119.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HD2"}) {
                               &buildatom($idx,"2HD2", $idx, "ND2", $idx,  "CG", $idx,  "1HD2", 1.02, 119.0, 120.0,  1);
                           }
                         }
            case ("CYS") { @scatom = $nohydrogens ? qw( CB SG ) : qw( CB SG HB2 HB3 HG );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"SG"}) {
                               &buildatom($idx,  "SG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.82, 109.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "SG", 1.11, 109.4, 112.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "SG", 1.11, 109.4, 112.0, -1);
                           }
                           unless (defined $crd{$idx}->{"HG"}) {
                               &buildatom($idx,  "HG", $idx,  "SG", $idx,  "CB", $idx,    "CA", 1.34,  96.0, 180.0,  0);
                           }
                         }
            case ("GLN") { @scatom = $nohydrogens ? qw( CB CG CD OE1 NE2 ) : qw( CB CG CD OE1 NE2 HB2 HB3 HG2 HG3 1HE2 2HE2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.51, 107.8, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OE1"}) {
                               &buildatom($idx, "OE1", $idx,  "CD", $idx,  "CG", $idx,    "CB", 1.22, 122.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NE2"}) {
                               &buildatom($idx, "NE2", $idx,  "CD", $idx,  "CG", $idx,   "OE1", 1.34, 112.7, 124.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 107.9,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 107.9, -1);
                           }
                           unless (defined $crd{$idx}->{"1HE2"}) {
                               &buildatom($idx,"1HE2", $idx, "NE2", $idx,  "CD", $idx,    "CG", 1.02, 119.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HE2"}) {
                               &buildatom($idx,"2HE2", $idx, "NE2", $idx,  "CD", $idx,  "1HE2", 1.02, 119.0, 120.0,  1);
                           }
                         }
            case ("GLU") { @scatom = $nohydrogens ? qw( CB CG CD OE1 OE2 ) : qw( CB CG CD OE1 OE2 HB2 HB3 HG2 HG3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.51, 107.8, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OE1"}) {
                               &buildatom($idx, "OE1", $idx,  "CD", $idx,  "CG", $idx,    "CB", 1.25, 117.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OE2"}) {
                               &buildatom($idx, "OE2", $idx,  "CD", $idx,  "CG", $idx,   "OE1", 1.25, 117.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 107.9,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 107.9, -1);
                           }
                         }
            case ("GLY") { @scatom = $nohydrogens ? qw( ) : qw( HA2 );
                           unless (defined $crd{$idx}->{"HA2"}) {
                               &buildatom($idx, "HA2", $idx,  "CA", $idx,   "N",   $idx,   "C", 1.11, 109.5, 107.9,  1);
                           }
                         }
            case ("HIS") { @scatom = $nohydrogens ? qw( CB CG ND1 CD2 CE1 NE2 ) : qw( CB CG ND1 CD2 CE1 NE2 HB2 HB3 HD1 HD2 HE1 HE2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"ND1"}) {
                               &buildatom($idx, "ND1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.35, 126.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "ND1", 1.35, 126.0, 108.0,  1);
                           }
                           unless (defined $crd{$idx}->{"CE1"}) {
                               &buildatom($idx, "CE1", $idx, "ND1", $idx,  "CG", $idx,   "CD2", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NE2"}) {
                               &buildatom($idx, "NE2", $idx, "CD2", $idx,  "CG", $idx,   "ND1", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD1"}) {
                               &buildatom($idx, "HD1", $idx, "ND1", $idx,  "CG", $idx,    "CB", 1.02, 126.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx, "CD2", $idx,  "CG", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "CE1", $idx, "ND1", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx, "NE2", $idx, "CD2", $idx,   "CE1", 1.02, 126.0, 126.0,  1);
                           }
                         }
            case ("HID") { @scatom = $nohydrogens ? qw( CB CG ND1 CD2 CE1 NE2 ) : qw( CB CG ND1 CD2 CE1 NE2 HB2 HB3 HD1 HD2 HE1 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"ND1"}) {
                               &buildatom($idx, "ND1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.35, 126.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "ND1", 1.35, 126.0, 108.0,  1);
                           }
                           unless (defined $crd{$idx}->{"CE1"}) {
                               &buildatom($idx, "CE1", $idx, "ND1", $idx,  "CG", $idx,   "CD2", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NE2"}) {
                               &buildatom($idx, "NE2", $idx, "CD2", $idx,  "CG", $idx,   "ND1", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD1"}) {
                               &buildatom($idx, "HD1", $idx, "ND1", $idx,  "CG", $idx,    "CB", 1.02, 126.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx, "CD2", $idx,  "CG", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "CE1", $idx, "ND1", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                         }
            case ("HIE") { @scatom = $nohydrogens ? qw( CB CG ND1 CD2 CE1 NE2 ) : qw( CB CG ND1 CD2 CE1 NE2 HB2 HB3 HD2 HE1 HE2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"ND1"}) {
                               &buildatom($idx, "ND1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.35, 126.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "ND1", 1.35, 126.0, 108.0,  1);
                           }
                           unless (defined $crd{$idx}->{"CE1"}) {
                               &buildatom($idx, "CE1", $idx, "ND1", $idx,  "CG", $idx,   "CD2", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NE2"}) {
                               &buildatom($idx, "NE2", $idx, "CD2", $idx,  "CG", $idx,   "ND1", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx, "CD2", $idx,  "CG", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "CE1", $idx, "ND1", $idx,   "NE2", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx, "NE2", $idx, "CD2", $idx,   "CE1", 1.02, 126.0, 126.0,  1);
                           }
                         }
            case ("ILE") { @scatom = $nohydrogens ? qw( CB CG1 CG2 CD1 ) : qw( CB CG1 CG2 CD1 HB 2HG1 3HG1 1HG2 2HG2 3HG2 1HD1 2HD1 3HD1 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG1"}) {
                               &buildatom($idx, "CG1", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CG2"}) {
                               &buildatom($idx, "CG2", $idx,  "CB", $idx,  "CA", $idx,   "CG1", 1.54, 109.5, 109.5,  1);
                           }
                           unless (defined $crd{$idx}->{"CD1"}) {
                               &buildatom($idx, "CD1", $idx, "CG1", $idx,  "CB", $idx,    "CA", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB"}) {
                               &buildatom($idx,  "HB", $idx,  "CB", $idx,  "CA", $idx,   "CG1", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"2HG1"}) {
                               &buildatom($idx,"2HG1", $idx, "CG1", $idx,  "CB", $idx,   "CD1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"3HG1"}) {
                               &buildatom($idx,"3HG1", $idx, "CG1", $idx,  "CB", $idx,   "CD1", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"1HG2"}) {
                               &buildatom($idx,"1HG2", $idx, "CG2", $idx,  "CB", $idx,   "CG1", 1.11, 110.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HG2"}) {
                               &buildatom($idx,"2HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 110.0, 109.0,  1);
                           }
                           unless (defined $crd{$idx}->{"3HG2"}) {
                               &buildatom($idx,"3HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 110.0, 109.0, -1);
                           }
                           unless (defined $crd{$idx}->{"1HD1"}) {
                               &buildatom($idx,"1HD1", $idx, "CD1", $idx, "CG1", $idx,    "CB", 1.11, 110.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HD1"}) {
                               &buildatom($idx,"2HD1", $idx, "CD1", $idx, "CG1", $idx,  "1HD1", 1.11, 110.0, 109.0,  1);
                           }
                           unless (defined $crd{$idx}->{"3HD1"}) {
                               &buildatom($idx,"3HD1", $idx, "CD1", $idx, "CG1", $idx,  "1HD1", 1.11, 110.0, 109.0, -1);
                           }
                         }
            case ("LEU") { @scatom = $nohydrogens ? qw( CB CG CD1 CD2 HB2 HB3 ) : qw( CB CG CD1 CD2 HB2 HB3 HG 1HD1 2HD1 3HD1 1HD2 2HD2 3HD2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD1"}) {
                               &buildatom($idx, "CD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "CD1", 1.54, 109.5, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG"}) {
                               &buildatom($idx,  "HG", $idx,  "CG", $idx,  "CB", $idx,   "CD1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"1HD1"}) {
                               &buildatom($idx,"1HD1", $idx, "CD1", $idx,  "CG", $idx,    "CB", 1.11, 109.4, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HD1"}) {
                               &buildatom($idx,"2HD1", $idx, "CD1", $idx,  "CG", $idx,  "1HD1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"3HD1"}) {
                               &buildatom($idx,"3HD1", $idx, "CD1", $idx,  "CG", $idx,  "1HD1", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"1HD2"}) {
                               &buildatom($idx,"1HD2", $idx, "CD2", $idx,  "CG", $idx,    "CB", 1.11, 109.4, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HD2"}) {
                               &buildatom($idx,"2HD2", $idx, "CD2", $idx,  "CG", $idx,  "1HD2", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"3HD2"}) {
                               &buildatom($idx,"3HD2", $idx, "CD2", $idx,  "CG", $idx,  "1HD2", 1.11, 109.4, 109.4, -1);
                           }
                         }
            case ("LYS") { @scatom = $nohydrogens ? qw( CB CG CD CE NZ ) : qw( CB CG CD CE NZ HB2 HB3 HG2 HG3 HD2 HD3 HE2 HE3 HZ1 HZ2 HZ3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE"}) {
                               &buildatom($idx,  "CE", $idx,  "CD", $idx,  "CG", $idx,    "CB", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NZ"}) {
                               &buildatom($idx,  "NZ", $idx,  "CE", $idx,  "CD", $idx,    "CG", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx,  "CD", $idx,  "CG", $idx,    "CE", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HD3"}) {
                               &buildatom($idx, "HD3", $idx,  "CD", $idx,  "CG", $idx,    "CE", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx,  "CE", $idx,  "CD", $idx,    "NZ", 1.11, 109.4, 108.8,  1);
                           }
                           unless (defined $crd{$idx}->{"HE3"}) {
                               &buildatom($idx, "HE3", $idx,  "CE", $idx,  "CD", $idx,    "NZ", 1.11, 109.4, 108.8, -1);
                           }
                           unless (defined $crd{$idx}->{"HZ1"}) {
                               &buildatom($idx, "HZ1", $idx,  "NZ", $idx,  "CE", $idx,    "CD", 1.02, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HZ2"}) {
                               &buildatom($idx, "HZ2", $idx,  "NZ", $idx,  "CE", $idx,   "HZ1", 1.02, 109.5, 109.5,  1);
                           }
                           unless (defined $crd{$idx}->{"HZ3"}) {
                               &buildatom($idx, "HZ3", $idx,  "NZ", $idx,  "CE", $idx,   "HZ1", 1.02, 109.5, 109.5, -1);
                           }
                         }
            case ("LYN") { @scatom = $nohydrogens ? qw( CB CG CD CE NZ ) : qw( CB CG CD CE NZ HB2 HB3 HG2 HG3 HD2 HD3 HE2 HE3 HZ2 HZ3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE"}) {
                               &buildatom($idx,  "CE", $idx,  "CD", $idx,  "CG", $idx,    "CB", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"NZ"}) {
                               &buildatom($idx,  "NZ", $idx,  "CE", $idx,  "CD", $idx,    "CG", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx,  "CD", $idx,  "CG", $idx,    "CE", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HD3"}) {
                               &buildatom($idx, "HD3", $idx,  "CD", $idx,  "CG", $idx,    "CE", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx,  "CE", $idx,  "CD", $idx,    "NZ", 1.11, 109.4, 108.8,  1);
                           }
                           unless (defined $crd{$idx}->{"HE3"}) {
                               &buildatom($idx, "HE3", $idx,  "CE", $idx,  "CD", $idx,    "NZ", 1.11, 109.4, 108.8, -1);
                           }
                           unless (defined $crd{$idx}->{"HZ2"}) {
                               &buildatom($idx, "HZ2", $idx,  "NZ", $idx,  "CE", $idx,    "CD", 1.02, 109.5,  60.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HZ3"}) {
                               &buildatom($idx, "HZ3", $idx,  "NZ", $idx,  "CE", $idx,    "CD", 1.02, 109.5, 300.0,  0);
                           }
                         }
            case ("MET") { @scatom = $nohydrogens ? qw( CB CG SD CE ) : qw( CB CG SD CE HB2 HB3 HG2 HG3 HE1 HE2 HE3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"SD"}) {
                               &buildatom($idx,  "SD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.82, 109.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE"}) {
                               &buildatom($idx,  "CE", $idx,  "SD", $idx,  "CG", $idx,    "CB", 1.82,  96.3, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "SD", 1.11, 109.4, 112.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "SD", 1.11, 109.4, 112.0, -1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx,  "CE", $idx,  "SD", $idx,    "CG", 1.11, 112.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx,  "CE", $idx,  "SD", $idx,   "HE1", 1.11, 112.0, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HE3"}) {
                               &buildatom($idx, "HE3", $idx,  "CE", $idx,  "SD", $idx,   "HE1", 1.11, 112.0, 109.4, -1);
                           }
                         }
            case ("PHE") { @scatom = $nohydrogens ? qw( CB CG CD1 CD2 CE1 CE2 CZ ) : qw( CB CG CD1 CD2 CE1 CE2 CZ HB2 HB3 HD1 HD2 HE1 HE2 HZ );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD1"}) {
                               &buildatom($idx, "CD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "CD1", 1.39, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"CE1"}) {
                               &buildatom($idx, "CE1", $idx, "CD1", $idx,  "CG", $idx,    "CB", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE2"}) {
                               &buildatom($idx, "CE2", $idx, "CD2", $idx,  "CG", $idx,    "CB", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CZ"}) {
                               &buildatom($idx,  "CZ", $idx, "CE1", $idx, "CD1", $idx,    "CG", 1.39, 120.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD1"}) {
                               &buildatom($idx, "HD1", $idx, "CD1", $idx,  "CG", $idx,   "CE1", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx, "CD2", $idx,  "CG", $idx,   "CE2", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "CE1", $idx, "CD1", $idx,    "CZ", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx, "CE2", $idx, "CD2", $idx,    "CZ", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HZ"}) {
                               &buildatom($idx,  "HZ", $idx,  "CZ", $idx, "CE1", $idx,   "CE2", 1.10, 120.0, 120.0,  1);
                           }
                         }
            case ("PRO") { @scatom = $nohydrogens ? qw( CB CG CD ) : qw( CB CG CD HB2 HB3 HG2 HG3 HD2 HD3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 107.0, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 107.0,  27.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD"}) {
                               &buildatom($idx,  "CD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.54, 107.0, -34.6,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "CD", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx,  "CD", $idx,  "CG", $idx,     "N", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HD3"}) {
                               &buildatom($idx, "HD3", $idx,  "CD", $idx,  "CG", $idx,     "N", 1.11, 109.4, 109.4, -1);
                           }
                         }
            case ("SER") { @scatom = $nohydrogens ? qw( CB OG ) : qw( CB OG HB2 HB3 HG );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"OG"}) {
                               &buildatom($idx,  "OG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.41, 107.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "OG", 1.11, 109.4, 106.7,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "OG", 1.11, 109.4, 106.7, -1);
                           }
                           unless (defined $crd{$idx}->{"HG"}) {
                               &buildatom($idx,  "HG", $idx,  "OG", $idx,  "CB", $idx,    "CA", 0.94, 106.9, 180.0,  0);
                           }
                         }
            case ("THR") { @scatom = $nohydrogens ? qw( CB OG1 CG2 ) : qw( CB OG1 CG2 HB HG1 1HG2 2HG2 3HG2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"OG1"}) {
                               &buildatom($idx, "OG1", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.41, 107.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CG2"}) {
                               &buildatom($idx, "CG2", $idx,  "CB", $idx,  "CA", $idx,   "OG1", 1.54, 109.5, 107.7,  1);
                           }
                           unless (defined $crd{$idx}->{"HB"}) {
                               &buildatom($idx,  "HB", $idx,  "CB", $idx,  "CA", $idx,   "OG1", 1.11, 109.4, 106.7, -1);
                           }
                           unless (defined $crd{$idx}->{"HG1"}) {
                               &buildatom($idx, "HG1", $idx, "OG1", $idx,  "CB", $idx,    "CA", 0.94, 106.9, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"1HG2"}) {
                               &buildatom($idx,"1HG2", $idx, "CG2", $idx,  "CB", $idx,    "CA", 1.11, 110.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HG2"}) {
                               &buildatom($idx,"2HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 110.0, 109.0,  1);
                           }
                           unless (defined $crd{$idx}->{"3HG2"}) {
                               &buildatom($idx,"3HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 110.0, 109.0, -1);
                           }
                         }
            case ("TRP") { @scatom = $nohydrogens ? qw( CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 ) : qw( CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 HB2 HB3 HD1 HE1 HE3 HZ2 HZ3 HH2);
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD1"}) {
                               &buildatom($idx, "CD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.35, 126.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "CD1", 1.35, 126.0, 108.0,  1);
                           }
                           unless (defined $crd{$idx}->{"NE1"}) {
                               &buildatom($idx, "NE1", $idx, "CD1", $idx,  "CG", $idx,   "CD2", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE2"}) {
                               &buildatom($idx, "CE2", $idx, "NE1", $idx, "CD1", $idx,    "CG", 1.35, 108.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE3"}) {
                               &buildatom($idx, "CE3", $idx, "CD2", $idx, "CE2", $idx,   "NE1", 1.35, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CZ2"}) {
                               &buildatom($idx, "CZ2", $idx, "CE2", $idx, "CD2", $idx,   "CE3", 1.35, 120.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CZ3"}) {
                               &buildatom($idx, "CZ3", $idx, "CE3", $idx, "CD2", $idx,   "CE2", 1.35, 120.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CH2"}) {
                               &buildatom($idx, "CH2", $idx, "CZ2", $idx, "CE2", $idx,   "CD2", 1.35, 120.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD1"}) {
                               &buildatom($idx, "HD1", $idx, "CD1", $idx,  "CG", $idx,   "NE1", 1.10, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "NE1", $idx, "CD1", $idx,   "CE2", 1.05, 126.0, 126.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE3"}) {
                               &buildatom($idx, "HE3", $idx, "CE3", $idx, "CD2", $idx,   "CZ3", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HZ2"}) {
                               &buildatom($idx, "HZ2", $idx, "CZ2", $idx, "CE2", $idx,   "CH2", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HZ3"}) {
                               &buildatom($idx, "HZ3", $idx, "CZ3", $idx, "CE3", $idx,   "CH2", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HH2"}) {
                               &buildatom($idx, "HH2", $idx, "CH2", $idx, "CZ2", $idx,   "CZ3", 1.10, 120.0, 120.0,  1);
                           }
                         }
            case ("TYR") { @scatom = $nohydrogens ? qw( CB CG CD1 CD2 CE1 CE2 CZ OH ) : qw( CB CG CD1 CD2 CE1 CE2 CZ OH HB2 HB3 HD1 HD2 HE1 HE2 HH );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.50, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD1"}) {
                               &buildatom($idx, "CD1", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CD2"}) {
                               &buildatom($idx, "CD2", $idx,  "CG", $idx,  "CB", $idx,   "CD1", 1.39, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"CE1"}) {
                               &buildatom($idx, "CE1", $idx, "CD1", $idx,  "CG", $idx,    "CB", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CE2"}) {
                               &buildatom($idx, "CE2", $idx, "CD2", $idx,  "CG", $idx,    "CB", 1.39, 120.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CZ"}) {
                               &buildatom($idx,  "CZ", $idx, "CE1", $idx, "CD1", $idx,    "CG", 1.39, 120.0,   0.0,  0);
                           }
                           unless (defined $crd{$idx}->{"OH"}) {
                               &buildatom($idx,  "OH", $idx,  "CZ", $idx, "CE2", $idx,   "CE1", 1.36, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HD1"}) {
                               &buildatom($idx, "HD1", $idx, "CD1", $idx,  "CG", $idx,   "CE1", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HD2"}) {
                               &buildatom($idx, "HD2", $idx, "CD2", $idx,  "CG", $idx,   "CE2", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx, "CE1", $idx, "CD1", $idx,    "CZ", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx, "CE2", $idx, "CD2", $idx,    "CZ", 1.10, 120.0, 120.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HH"}) {
                               &buildatom($idx,  "HH", $idx,  "OH", $idx,  "CZ", $idx,   "CE2", 0.97, 108.0,   0.0,  0);
                           }
                         }
            case ("VAL") { @scatom = $nohydrogens ? qw( CB CG1 CG2 ) : qw( CB CG1 CG2 HB 1HG1 2HG1 3HG1 1HG2 2HG2 3HG2 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG1"}) {
                               &buildatom($idx, "CG1", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"CG2"}) {
                               &buildatom($idx, "CG2", $idx,  "CB", $idx,  "CA", $idx,   "CG1", 1.54, 109.5, 109.5, -1);
                           }
                           unless (defined $crd{$idx}->{"HB"}) {
                               &buildatom($idx,  "HB", $idx,  "CB", $idx,  "CA", $idx,   "CG1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"1HG1"}) {
                               &buildatom($idx,"1HG1", $idx, "CG1", $idx,  "CB", $idx,    "CA", 1.11, 109.4, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HG1"}) {
                               &buildatom($idx,"2HG1", $idx, "CG1", $idx,  "CB", $idx,  "1HG1", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"3HG1"}) {
                               &buildatom($idx,"3HG1", $idx, "CG1", $idx,  "CB", $idx,  "1HG1", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"1HG2"}) {
                               &buildatom($idx,"1HG2", $idx, "CG2", $idx,  "CB", $idx,    "CA", 1.11, 109.4, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"2HG2"}) {
                               &buildatom($idx,"2HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"3HG2"}) {
                               &buildatom($idx,"3HG2", $idx, "CG2", $idx,  "CB", $idx,  "1HG2", 1.11, 109.4, 109.4, -1);
                           }
                         }
            case ("CYX") { @scatom = $nohydrogens ? qw( CB SG ) : qw( CB SG HB2 HB3 );
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"SG"}) {
                               &buildatom($idx,  "SG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.82, 109.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "SG", 1.11, 109.4, 112.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "SG", 1.11, 109.4, 112.0, -1);
                           }
                         }
            case ("MSE") { @scatom = $nohydrogens ? qw( CB CG SD CE ) : qw( CB CG SD CE HB2 HB3 HG2 HG3 HE1 HE2 HE3 );
                           my $chi2 = my $chi3 = 180.0;
                           unless (defined $crd{$idx}->{"CB"}) {
                               &buildatom($idx,  "CB", $idx,  "CA", $idx,   "N", $idx-1,   "C", 1.54, 109.5, 107.8, -1);
                           }
                           unless (defined $crd{$idx}->{"CG"}) {
                               &buildatom($idx,  "CG", $idx,  "CB", $idx,  "CA", $idx,     "N", 1.54, 109.5, 180.0,  0);
                           }
                           if (defined $crd{$idx}->{"SE"}) {
                               $chi2 = &getgeometry($idx, "SE", $idx, "CG", $idx, "CB", $idx, "CA");
                           }
                           if (defined $crd{$idx}->{"CE"}) {
                               $chi3 = &getgeometry($idx, "CE", $idx, "SE", $idx, "CG", $idx, "CB");
                           }
                           &buildatom($idx,  "SD", $idx,  "CG", $idx,  "CB", $idx,    "CA", 1.82, 109.0, $chi2,  0);
                           &buildatom($idx,  "CE", $idx,  "SD", $idx,  "CG", $idx,    "CB", 1.82,  96.3, $chi3,  0);
                           unless (defined $crd{$idx}->{"HB2"}) {
                               &buildatom($idx, "HB2", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HB3"}) {
                               &buildatom($idx, "HB3", $idx,  "CB", $idx,  "CA", $idx,    "CG", 1.11, 109.4, 109.4, -1);
                           }
                           unless (defined $crd{$idx}->{"HG2"}) {
                               &buildatom($idx, "HG2", $idx,  "CG", $idx,  "CB", $idx,    "SD", 1.11, 109.4, 112.0,  1);
                           }
                           unless (defined $crd{$idx}->{"HG3"}) {
                               &buildatom($idx, "HG3", $idx,  "CG", $idx,  "CB", $idx,    "SD", 1.11, 109.4, 112.0, -1);
                           }
                           unless (defined $crd{$idx}->{"HE1"}) {
                               &buildatom($idx, "HE1", $idx,  "CE", $idx,  "SD", $idx,    "CG", 1.11, 112.0, 180.0,  0);
                           }
                           unless (defined $crd{$idx}->{"HE2"}) {
                               &buildatom($idx, "HE2", $idx,  "CE", $idx,  "SD", $idx,   "HE1", 1.11, 112.0, 109.4,  1);
                           }
                           unless (defined $crd{$idx}->{"HE3"}) {
                               &buildatom($idx, "HE3", $idx,  "CE", $idx,  "SD", $idx,   "HE1", 1.11, 112.0, 109.4, -1);
                           }
                         }
            else         { die "Residue $resname{$idx} is not defined.\n"; }
        }
        foreach my $name (@scatom) {
            &printatom($idx, $name);
        }
    }
}

####====####
sub printatom {
    my ($idx, $name) = @_;
    unless ($name eq "OXT") {
        die "Atom $name is not found in residue $resname{$idx} $idx.\n" unless (defined $crd{$idx}->{$name});
    }
    return unless (defined $crd{$idx}->{$name});
    $natom ++;
    my $name2 = undef;
    switch (length($name)) {
        case (1) { $name2 = " $name  "; }
        case (2) { $name2 = " $name "; }
        case (3) { $name2 = " $name"; }
        case (4) { $name2 = "$name"; }
        else     { die "Atom name length incorrect.\n"; }
    }
##  correct residue name
    switch ($resname{$idx}) {
        case ("MSE") { $resname{$idx} = "MET"; }
    }
    printf OUT "%6s%5d %4s %3s %1s%4d    ", "ATOM  ", $natom, $name2, $resname{$idx}, $chain, $idx;
    printf OUT "%8.3f", $crd{$idx}->{$name}->[$_] foreach (0 .. 2);
    print OUT "\n";
}

####====####
sub removespace {
    my ($char) = @_;
    my @store;
    foreach my $n (0 .. length($char)-1) {
        push (@store, substr($char, $n, 1)) unless (substr($char, $n, 1) eq " ");
    }
    $char = "";
    $char = $char.$_ foreach (@store);
    return $char;
}

####====####
sub buildatom {
    my ($idx0, $name0, $idx1, $name1, $idx2, $name2, $idx3, $name3, $bond, $angle1, $angle2, $chiral) = @_;
#   print "$idx0, $name0, $idx1, $name1, $idx2, $name2, $idx3, $name3, $bond, $angle1, $angle2, $chiral\n";
    my ($x0, $y0, $z0);

    my $eps = 0.00000001;
    my $radian = 57.29577951308232088;
    my $rad1 = $angle1/$radian;
    my $rad2 = $angle2/$radian;
    my $sin1 = sin($rad1);
    my $cos1 = cos($rad1);
    my $sin2 = sin($rad2);
    my $cos2 = cos($rad2);

    my $xa = $crd{$idx1}->{$name1}->[0];
    my $ya = $crd{$idx1}->{$name1}->[1];
    my $za = $crd{$idx1}->{$name1}->[2];
    my $xb = $crd{$idx2}->{$name2}->[0];
    my $yb = $crd{$idx2}->{$name2}->[1];
    my $zb = $crd{$idx2}->{$name2}->[2];
    my $xc = $crd{$idx3}->{$name3}->[0];
    my $yc = $crd{$idx3}->{$name3}->[1];
    my $zc = $crd{$idx3}->{$name3}->[2];

    if ($chiral == 0) {

        my $xab = $xa - $xb;
        my $yab = $ya - $yb;
        my $zab = $za - $zb;
        my $rab = sqrt($xab**2 + $yab**2 + $zab**2);
        $xab /= $rab;
        $yab /= $rab;
        $zab /= $rab;

        my $xbc = $xb - $xc;
        my $ybc = $yb - $yc;
        my $zbc = $zb - $zc;
        my $rbc = sqrt($xbc**2 + $ybc**2 + $zbc**2);
        $xbc /= $rbc;
        $ybc /= $rbc;
        $zbc /= $rbc;

        my $xt = $zab*$ybc - $yab*$zbc;
        my $yt = $xab*$zbc - $zab*$xbc;
        my $zt = $yab*$xbc - $xab*$ybc;
        my $cosine = $xab*$xbc + $yab*$ybc + $zab*$zbc;
        my $sine = sqrt(1.0-$cosine**2);

        $xt /= $sine;
        $yt /= $sine;
        $zt /= $sine;
        my $xu = $yt*$zab - $zt*$yab;
        my $yu = $zt*$xab - $xt*$zab;
        my $zu = $xt*$yab - $yt*$xab;

        $x0 = $xa + $bond*($xu*$sin1*$cos2 + $xt*$sin1*$sin2 - $xab*$cos1);
        $y0 = $ya + $bond*($yu*$sin1*$cos2 + $yt*$sin1*$sin2 - $yab*$cos1);
        $z0 = $za + $bond*($zu*$sin1*$cos2 + $zt*$sin1*$sin2 - $zab*$cos1);

    } elsif (abs($chiral) == 1) {

        my $xba = $xb - $xa;
        my $yba = $yb - $ya;
        my $zba = $zb - $za;
        my $rba = sqrt($xba**2 + $yba**2 + $zba**2);
        $xba /= $rba;
        $yba /= $rba;
        $zba /= $rba;

        print "##  ", $xa, "  ##  ", $xc, "  ##  ", $idx3, "  ##  ", $name3, "  ##", "\n" if (not defined $xa or not defined $xc);

        my $xac = $xa - $xc;
        my $yac = $ya - $yc;
        my $zac = $za - $zc;
        my $rac = sqrt($xac**2 + $yac**2 + $zac**2);
        $xac /= $rac;
        $yac /= $rac;
        $zac /= $rac;

        my $xt = $zba*$yac - $yba*$zac;
        my $yt = $xba*$zac - $zba*$xac;
        my $zt = $yba*$xac - $xba*$yac;
        my $cosine = $xba*$xac + $yba*$yac + $zba*$zac;
        my $sine2 = 1.0 - $cosine**2;

        my $a = (-$cos2 - $cosine*$cos1)/$sine2;
        my $b = ( $cos1 + $cosine*$cos2)/$sine2;
        my $c = (1.0+ $a*$cos2 - $b*$cos1)/$sine2;
        if ($c > $eps) {
            $c = $chiral*sqrt($c);
        } elsif ($c < -$eps) {
            $c = sqrt(($a*$xac + $b*$xba)**2 + ($a*$yac + $b*$yba)**2 + ($a*$zac + $b*$zba)**2);
            $a /= $c;
            $b /= $c;
            $c = 0.0;
         } else {
            $c = 0.0;
         }

         $x0 = $xa + $bond*($a*$xac + $b*$xba + $c*$xt);
         $y0 = $ya + $bond*($a*$yac + $b*$yba + $c*$yt);
         $z0 = $za + $bond*($a*$zac + $b*$zba + $c*$zt);

    }

    $crd{$idx0}->{$name0}->[0] = $x0;
    $crd{$idx0}->{$name0}->[1] = $y0;
    $crd{$idx0}->{$name0}->[2] = $z0;
}

####====####
sub getgeometry {
    my ($idx1, $name1, $idx2, $name2, $idx3, $name3, $idx4, $name4) = @_;
    my $radian = 57.29577951308232088;
    my $geometry = 999.999;

    foreach my $n (0 .. 2) {
        goto DONE unless (defined $crd{$idx1}->{$name1});
        goto DONE unless (defined $crd{$idx2}->{$name2});
        goto DONE unless (defined $crd{$idx3}->{$name3});
        goto DONE unless (defined $crd{$idx4}->{$name4});
    }

    my @x = my @y = my @z = ();
    $x[1] = $crd{$idx1}->{$name1}->[0];
    $y[1] = $crd{$idx1}->{$name1}->[1];
    $z[1] = $crd{$idx1}->{$name1}->[2];
    $x[2] = $crd{$idx2}->{$name2}->[0];
    $y[2] = $crd{$idx2}->{$name2}->[1];
    $z[2] = $crd{$idx2}->{$name2}->[2];
    $x[3] = $crd{$idx3}->{$name3}->[0];
    $y[3] = $crd{$idx3}->{$name3}->[1];
    $z[3] = $crd{$idx3}->{$name3}->[2];
    $x[4] = $crd{$idx4}->{$name4}->[0];
    $y[4] = $crd{$idx4}->{$name4}->[1];
    $z[4] = $crd{$idx4}->{$name4}->[2];

    if ($idx3 == -9999) {
        my $xab = $x[1] - $x[2];
        my $yab = $y[1] - $y[2];
        my $zab = $z[1] - $z[2];
        $geometry = sqrt($xab*$xab + $yab*$yab + $zab*$zab);
    } elsif ($idx4 == -9999) {
        my $xab = $x[1] - $x[2];
        my $yab = $y[1] - $y[2];
        my $zab = $z[1] - $z[2];
        my $xcb = $x[3] - $x[2];
        my $ycb = $y[3] - $y[2];
        my $zcb = $z[3] - $z[2];
        my $rab2 = $xab*$xab + $yab*$yab + $zab*$zab;
        my $rcb2 = $xcb*$xcb + $ycb*$ycb + $zcb*$zcb;
        my $rabc = sqrt($rab2*$rcb2);
        unless ($rabc == 0.0) {
            my $cosine = ($xab*$xcb + $yab*$ycb + $zab*$zcb) / $rabc;
            $cosine = &min(1.0, &max(-1.0, $cosine));
            $geometry = $radian * acos2($cosine);
        } else {
            die "Error found.\n";
        }
    } else {
        my $xba = $x[2] - $x[1];
        my $yba = $y[2] - $y[1];
        my $zba = $z[2] - $z[1];
        my $xcb = $x[3] - $x[2];
        my $ycb = $y[3] - $y[2];
        my $zcb = $z[3] - $z[2];
        my $xdc = $x[4] - $x[3];
        my $ydc = $y[4] - $y[3];
        my $zdc = $z[4] - $z[3];
        my $xt = $yba*$zcb - $ycb*$zba;
        my $yt = $xcb*$zba - $xba*$zcb;
        my $zt = $xba*$ycb - $xcb*$yba;
        my $xu = $ycb*$zdc - $ydc*$zcb;
        my $yu = $xdc*$zcb - $xcb*$zdc;
        my $zu = $xcb*$ydc - $xdc*$ycb;
        my $rt2 = $xt*$xt + $yt*$yt + $zt*$zt;
        my $ru2 = $xu*$xu + $yu*$yu + $zu*$zu;
        my $rtru = sqrt($rt2*$ru2);
        unless ($rtru == 0.0) {
           my $cosine = ($xt*$xu + $yt*$yu + $zt*$zu) / $rtru;
           $cosine = &min(1.0, &max(-1.0, $cosine));
           $geometry = $radian * acos($cosine);
           my $sign = $xba*$xu + $yba*$yu + $zba*$zu;
           $geometry = -$geometry if ($sign < 0.0);
        } else {
            die "Error found.\n";
        }
    }
DONE:
    return $geometry;
}

####====####
sub max {
    my ($a1, $a2) = @_;
    if ($a1 > $a2) {
        return $a1;
    } else {
        return $a2;
    }
}

####====####
sub min {
    my ($a1, $a2) = @_;
    if ($a1 < $a2) {
        return $a1;
    } else {
        return $a2;
    }
}
