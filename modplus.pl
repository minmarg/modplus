#!/usr/bin/perl

##
## (C)2007-2019 Mindaugas Margelevicius
## Institute of Biotechnology
## Vilnius University
##

use strict;
use FindBin;
use lib "$FindBin::Bin";
use nwalign;
use Getopt::Long;
use File::Basename;


my  $MYPROGNAME = basename( $0 );
#my  $MODELLER = '/usr/local/install/modeller-installed/modeller9v4/bin/mod9v4';
#my  $MODELLER = '/usr/local/install/modeller-installed/modeller9.13/bin/mod9.13';
#my  $MODELLER = '/share/data/install/modeller-installed/modeller9.13/bin/mod9.13';
my  $MODELLER = '/usr/local/install/modeller-installed/modeller9.15/bin/mod9.15';
my  $TMSCORE  = '/usr/local/install/TMscore/TMscore';
my  $WORKDIR = 'models';

my  $MODELSUFFIX  = 'B99990001';
my  $TRANSUFFIX = 'trans';
my  $PDBSUFFIX = 'pdb';
my  $ENTSUFFIX = 'ent';
my  $ORGSUFFIX = 'org';
my  $FASUFFIX = 'fa';
my  $PIREXT = 'pir';
my  $PYEXT = 'py';
my  $TMSUFFIX = 'tm';
my  $FRAGSUFFIX = 'frag';
my  $RASMOLSUFFIX = 'rsml';

my  $usage = <<EOIN;

Sequence-structure alignment and protein modelling helper tool.
(C)2007-2019 Mindaugas Margelevicius, 
Institute of Biotechnology, Vilnius University

Usage:
$MYPROGNAME <Parameters>

Parameters:

--in <file>        Input file of aligned sequences in FASTA
                   (see also --help).

--pdb <pdbdir>     Directory of template PDB files.

--dir <outdir>     Directory of output generated files.
           default=$WORKDIR

--all              Model all sequences in input file.
           By default, only the first sequence is modelled

--tm [<filename>]  Superimpose resulting models with TM-score.
                   Optionally, write the extracted TM-scores to
                   file <filename>.

--modeller <pathname> Full pathname to modeller executable.
           default=$MODELLER

--TMscore <pathname>  Full pathname to TMscore executable.
           default=$TMSCORE

--help             Description of the input format.

EOIN

my  $description = <<EODESC;

modplus is a helper tool to model the target protein sequence using
template structures and MODELLER (Sali and Blundell, J Mol Biol, 
234(3), 779-815, 1993). It resolves inconsistencies that may arise
between template sequences and structures. The correspondence
between a model and the native target structure is optionally 
calculated with TM-score (Zhang and Skolnick, Proteins, 57(4),
702-710, 2004).

Input format:

>targetname [chains];   [optional user information]
G-----VDILR-MDAVAF-------IWKQMGTSCE ...
>templ1name [chains];   [optional user information]
FANYDEHLIFEGMNEPRLVGHANEWWPELTNSDVV ...
>templ2name [chains];   [optional user information]
EMLAIAVETKPHFCCLVPEKRQEVTTEGGLDV--A ...
> ...
...

Description:

 The input file contains aligned sequences in FASTA. By default, 
the first sequence is the target to be modelled. Optionally, other 
sequences (templates) may be modelled as well (see below). One or 
more template sequences (two or more sequences in total) can be 
listed in the file. A target is modelled using all the template 
structures and alignments listed in the file.

 The names templ1name, templ2name,... (and targetname if options 
--all and/or --tm are specified) have to correspond to the 
filenames (without extension) of PDB structures. Required PDB files 
with extension .ent or .pdb added are looked up in the directory 
specified by the option --pdb.

 Template chains are extracted from pdb files. There may be several 
chain identifiers given by `chains' (e.g., SCOP domains compiled: 
ABC). If no `chains' are given, the first one is identified and 
used automatically.

 If option --all is specified, each sequence in the input file is 
modelled using the other sequence(s) as template(s).

 Option --tm implies a model quality check with TM-Score. The 
TM-Score output files will appear in the output directory specified 
by the option --dir. An optional filename following the --tm option 
specifies the following format for the resulting TM-scores written 
in this file: 

(TM: <TMScore1> <TMScore2>)

<TMScore1> is a score obtained with respect to the whole structure, 
and <TMScore2> is a score obtained along the alignment extent.

EODESC


my  $Hitfile;
my  $outputname;
my  $OUTPUT = \*STDOUT;

my  $PDBSDIR;

my  $ALLMOD = 0;
my  $RUNTM;
my  $CASPFORM = 0;
my  $MAXSBJCTS = 10;  ## maximum number of templates

my  $noargs = $#ARGV;
my  $result = GetOptions( 
               'in=s'      => \$Hitfile,
               'pdb=s'     => \$PDBSDIR,
               'dir=s'     => \$WORKDIR,
               'all'       => \$ALLMOD,
               'tm:s'      => \$RUNTM,
               'casp'      => \$CASPFORM,
               'modeller=s'=> \$MODELLER,
               'TMscore=s' => \$TMSCORE,
               'help|h'    => sub { print $usage; print $description; exit( 0 ); }
);

do { print $usage; exit( 1 ); }  if !$result || $noargs < 0;

die "ERROR: No input file.\n$usage" unless $Hitfile;
-f $Hitfile || die "ERROR: File $Hitfile does not exist.";

-f $MODELLER || die "ERROR: Modeller executable $MODELLER does not exist.";

if( defined $RUNTM ) {
    $outputname = $RUNTM if length( $RUNTM );
    $RUNTM = 1;
}
else {
    $RUNTM = 0;
}
( $RUNTM && ! -f $TMSCORE  ) && die "ERROR: TMscore executable $TMSCORE does not exist.";

die "ERROR: No directory of pdb files given." unless $PDBSDIR;

$PDBSDIR = basename( $PDBSDIR ) unless -d $PDBSDIR;

-d $PDBSDIR || die "ERROR: Directory $PDBSDIR does not exist.";


## ===================================================================
## modeller python script
##

my  $SCRIPT = <<EOF;

# Homology modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = './:$PDBSDIR/'

# Read in HETATM records from template PDBs
env.io.hetatm = True


a = automodel(env,
              alnfile  = '???????.pir',   # alignment filename
              knowns   = '???????',       # codes of the templates
              sequence = '???????')       # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)
##a.auto_align()                      # get an automatic alignment
a.make()                            # do the actual homology modeling

EOF

## ===================================================================

die "ERROR: Failed to mkdir $WORKDIR\n" unless( -d $WORKDIR || mkdir( $WORKDIR ));

my  $curdir = `pwd`;
chomp( $curdir );

$Hitfile = "$curdir/$Hitfile" unless( substr( $Hitfile, 0, 1 ) eq '/' );
$PDBSDIR = "$curdir/$PDBSDIR" unless( substr( $PDBSDIR, 0, 1 ) eq '/' );
$outputname = "$curdir/$outputname" if( $outputname && substr( $outputname, 0, 1 ) ne '/' );

chdir( $WORKDIR );

## -------------------------------------------------------------------
## main globals
##

my  %PDBHASH;
my  %FASHASH;

my  %RESDS;
InitData();


## -------------------------------------------------------------------
## Parse file for alignment and call processing method
##

my  $last;

my  $QUERY;
my  $SBJCT;
my  @SUBJECTS;
my  $querychain;
my  $sbjctchain;
my  $userinfo;

my  $QUERYFASTA;

my  $querytmscore;
my  $sbjcttmscore;
my  $queryfragtmscore;
my  $sbjctfragtmscore;

my  $s = -1;
my  $retcode;

## define field numbers for @SUBJECTS
##
my ($NAME, $CHNS, $STAR, $SEND, $INFO, $ALIG, $TMDM, $TMFR, $ENTR, $SEQF ) = 
        ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 );


open( FD, "<$Hitfile" ) or die "ERROR: Cannot open file $Hitfile";


while( <FD> ) {
    chomp;
    next if /^$/;
    $last = $_;

    if( $last =~ /^>/) {
        unless( $last =~ /^>([^\s;]+)\s*([^\s;]*)\s*;\s*(.*)$/ ) {
            die "ERROR: Invalid format of input file.";
        }
        unless( defined( $QUERY )) {
            $QUERY = $1;
            $querychain = $2;
            $userinfo = $3;
            $querychain =~ s/[^a-zA-Z0-9]/ /g;
        }
        else {
            $s = ++$#SUBJECTS;
            $SUBJECTS[$s]->[$NAME] = $SBJCT = $1;
            $SUBJECTS[$s]->[$CHNS] = $sbjctchain = $2;
            $SUBJECTS[$s]->[$INFO] = $3;
            die "ERROR: Maximum number of templates, $MAXSBJCTS, exceeded."
                if $MAXSBJCTS <= $s;

            $sbjctchain =~ s/[^a-zA-Z0-9]/ /g;
            $SUBJECTS[$s]->[$CHNS] = $sbjctchain;
        }
        next;
    }

    if( 0 <= $s ) { $SUBJECTS[$s]->[$ALIG] .= $last }
        else      { $QUERYFASTA .= $last }
}

close( FD );

die "ERROR: No query found in file." if !defined( $QUERY );
die "ERROR: No templates found in file." if $#SUBJECTS < 0;

$QUERYFASTA =~ s/\s//g;

foreach( @SUBJECTS ) {
    die "ERROR: No template name obtained. Check input file." unless defined( $_->[$NAME] );
    $_->[$ALIG] =~ s/\s//g;

    die "ERROR: Wrong length of sequences of $QUERY and $_->[$NAME]." 
        if length( $QUERYFASTA ) != length( $_->[$ALIG] ) || length( $_->[$ALIG] ) == 0;
}

printf( STDERR "Processing %s ...\n", $QUERY );
$retcode = 
    ProcessHit( $QUERY,
        $querychain,
        $userinfo,
        \$QUERYFASTA, \@SUBJECTS,
        \$querytmscore,
        \$queryfragtmscore
    );
if( $outputname && $retcode && $RUNTM ) {
    if( $outputname ) {
        open( $$OUTPUT, ">$outputname" ) or die "ERROR: Cannot open file $outputname for writing.";
    }
    WriteHit(   $OUTPUT,
                $QUERY,
                $querychain,
                $userinfo,
                \$QUERYFASTA, \@SUBJECTS,
                $querytmscore,
                $queryfragtmscore
    );
    close( $OUTPUT ) if $outputname;
}

chdir $curdir;
exit( 0 );

## -------------------------------------------------------------------
## Process alignment
##

sub ProcessHit
{
    my  $llquery = shift;
    my  $llquerychain = shift;
    my  $lluserinfo = shift;
    my  $rqufasta = shift;
    my  $refSUBJECTS = shift;
    my  $refquerytmscore = shift;
    my  $refqueryfragtmscore = shift;

    my  $s;
    my  $code = 1;
    my  $locdirname;

    my  $llsbjct;
    my  $llsbjctchain;
    my  $llsfasta;
    my  $llsbjctinfo;
    my  $pdbsbjctent;
    my  $sbjcttmscore;
    my  $sbjctfragtmscore;
    my  @copySUBJECTS;

    my  $tneeds = $ALLMOD || $RUNTM;



    my  $pdbqueryent  = "$PDBSDIR/$llquery.$ENTSUFFIX";

    if( $tneeds && !-f $pdbqueryent ) {
        $pdbqueryent = "$PDBSDIR/$llquery.$PDBSUFFIX";
        unless( -f $pdbqueryent ) {
            $pdbqueryent = "$llquery.$ENTSUFFIX";
            unless( -f $pdbqueryent ) {
                $pdbqueryent = "$llquery.$PDBSUFFIX";
                die "ERROR: PDB entry of $llquery not found." unless -f $pdbqueryent;
            }
        }
    }

    $locdirname = sprintf( "%s-%s", $llquery, $$refSUBJECTS[0]->[$NAME] );

    unless( -d $locdirname || mkdir( $locdirname )) {
        die "ERROR: Failed to create directory $locdirname";
    }

    if( $tneeds && !RunCommand( "cp $pdbqueryent $locdirname" )) {
        printf( STDERR "ERROR: Command `cp' failed.\n" );
        return 0;
    }


    foreach( @$refSUBJECTS )
    {
        $llsbjct = $_->[$NAME];
        $pdbsbjctent  = "$PDBSDIR/$llsbjct.$ENTSUFFIX";

        unless( -f $pdbsbjctent ) {
            $pdbsbjctent = "$PDBSDIR/$llsbjct.$PDBSUFFIX";
            unless( -f $pdbsbjctent ) {
                printf( STDERR "ERROR: PDB entry of %s not found.\n", $llsbjct );
                return 0;
            }
        }

        unless( RunCommand( "cp $pdbsbjctent $locdirname" )) {
            printf( STDERR "ERROR: Failed to copy %s.\n", $llsbjct );
            return 0;
        }

        $_->[$ENTR] = basename( $pdbsbjctent );
    }


    unless( chdir( $locdirname )) {
        printf( STDERR "ERROR: Failed to chdir %s\n", $locdirname );
        return 0;
    }

    $pdbqueryent  = basename( $pdbqueryent );

    if( $tneeds ) {
        unless( -f $pdbqueryent ) {
            chdir '..';
            printf( STDERR "ERROR: PDB entry %s not found.\n", $pdbqueryent );
            return 0;
        }
        unless( RunCommand( "cp $pdbqueryent $pdbqueryent.$ORGSUFFIX" )) {
            chdir '..';
            printf( STDERR "ERROR: Failed to make copy of target %s.\n", $pdbqueryent );
            return 0;
        }
    }

    foreach( @$refSUBJECTS )
    {
        unless( -f $_->[$ENTR] ) {
            chdir '..';
            printf( STDERR "ERROR: PDB entry %s not found.\n", $_->[$ENTR] );
            return 0;
        }
        unless( RunCommand( "cp $_->[$ENTR] $_->[$ENTR].$ORGSUFFIX" )) {
            chdir '..';
            printf( STDERR "ERROR: Failed to make copy of %s.\n", $_->[$ENTR] );
            return 0;
        }
        push @copySUBJECTS, [ @$_ ];
    }


    $code = 
    ProcessSymmetry(
        $locdirname,
        $llquery,
        $llquerychain,
        $lluserinfo,
        $$rqufasta,    \@copySUBJECTS, ## pass copies to enable edition
        $pdbqueryent,
        $refquerytmscore,
        $refqueryfragtmscore
    );

    if( !$code || !$ALLMOD ) {
        chdir '..';
        return $code;
    }

    ## model all the subjects next
    ##
    for( $s = 0; $s <= $#$refSUBJECTS; $s++ ) {
        splice( @copySUBJECTS );
        ## reinitialize all subjects again
        ##
        push @copySUBJECTS, [ @{$_} ] foreach( @{$refSUBJECTS} );

        ## exchange places of query and subject to be modelled
        ##
        $llsbjct      = $copySUBJECTS[$s]->[$NAME]; $copySUBJECTS[$s]->[$NAME] = $llquery;
        $llsbjctchain = $copySUBJECTS[$s]->[$CHNS]; $copySUBJECTS[$s]->[$CHNS] = $llquerychain;
        $llsbjctinfo  = $copySUBJECTS[$s]->[$INFO]; $copySUBJECTS[$s]->[$INFO] = $lluserinfo;
        $llsfasta     = $copySUBJECTS[$s]->[$ALIG]; $copySUBJECTS[$s]->[$ALIG] = $$rqufasta;
        $pdbsbjctent  = $copySUBJECTS[$s]->[$ENTR]; $copySUBJECTS[$s]->[$ENTR] = $pdbqueryent;

        $code = 
        ProcessSymmetry(
            $locdirname,
            $llsbjct,
            $llsbjctchain,
            $llsbjctinfo,
            $llsfasta,     \@copySUBJECTS,   ## pass copies to enable edition
            $pdbsbjctent,
            \$sbjcttmscore,
            \$sbjctfragtmscore
        );

        ## store results in the original structure of subjects
        ##
        $$refSUBJECTS[$s]->[$TMDM] = $sbjcttmscore;
        $$refSUBJECTS[$s]->[$TMFR] = $sbjctfragtmscore;

        unless( $code ) {
            chdir '..';
            return $code;
        }
    }

    chdir '..';
    return $code;
}

## -------------------------------------------------------------------
## write processed hit to file
##

sub WriteHit
{
    my  $reffile = shift;
    my  $llquery = shift;
    my  $llquerychain = shift;
    my  $lluserinfo = shift;
    my  $rqufasta = shift;
    my  $refSUBJECTS = shift;
    my  $llquerytmscore = shift;
    my  $llqueryfragtmscore = shift;
    my  $fail;

    $llquerytmscore = '--' unless defined( $llquerytmscore );
    $llqueryfragtmscore = '--' unless defined( $llqueryfragtmscore );

    printf( $reffile ">%-13s %-5s  %-24s (TM: %7s %7s )\n", 
            $llquery, $llquerychain, $lluserinfo, $llquerytmscore, $llqueryfragtmscore );
    $fail = 1 unless WrapFasta( $reffile, $rqufasta );

    foreach( @$refSUBJECTS )
    {
        my  $llsbjct      =  $_->[$NAME];
        my  $llsbjctchain =  $_->[$CHNS];
        my  $llsbjctinfo  =  $_->[$INFO];
        my  $rsbfasta     = \$_->[$ALIG];
        my  $llsbjcttmscore     = $_->[$TMDM];
        my  $llsbjctfragtmscore = $_->[$TMFR];

        $llsbjcttmscore = '--' unless defined( $llsbjcttmscore );
        $llsbjctfragtmscore = '--' unless defined( $llsbjctfragtmscore );

        printf( $reffile ">%-13s %-5s  %-24s (TM: %7s %7s )\n",
                $llsbjct, $llsbjctchain, $llsbjctinfo, $llsbjcttmscore, $llsbjctfragtmscore );
        $fail = 1 unless WrapFasta( $reffile, $rsbfasta );
    }
    print( $reffile "//\n" );
    return 0 if $fail;
    return 1;
}

## -------------------------------------------------------------------
## Process hit given query - subject
##

sub ProcessSymmetry
{
    my  $locdirname = shift;     ## local directory name
    my  $llquery = shift;        ## query code
    my  $llquerychain = shift;   ## query chain ids
    my  $lluserinfo = shift;     ## statistical significance of hit
    my  $queryfasta = shift;     ## query alignment sequence
    my  $copSUBJECTS = shift;    ## reference to all the subjects
    my  $pdbqueryent  = shift;   ## name of query pdb file
    my  $reftmscore = shift;     ## reference to TM score
    my  $reftmscorefrag = shift; ## reference to fragment TM score

    my  $modelname = "$WORKDIR/$locdirname/$llquery.$MODELSUFFIX.$PDBSUFFIX";
    my  $modelname2 = "$WORKDIR/$locdirname/$llquery.$MODELSUFFIX.$TRANSUFFIX.$PDBSUFFIX";

    my  $llsbjct;         ## subject code
    my  $llsbjctchain;    ## subject chain ids
    my  $llsbjctinfo;     ## user info (e-value, scop...) of subject
    my  $sbjctfasta;      ## reference ot subject alignment sequence
    my  $pdbsbjctent;     ## name of subject pdb file

    my  @AUXSBJCTS;       ## auxiliary vector of subjects

    my  $s;
    my  $code;

    my  $querydomseq;     ## reference to raw query fasta sequence
    my  $sbjctdomseq;     ## reference to raw subject fasta sequence

    my  @target;          ## target pdb structure to be modelled
    my  $qubegresnum;     ## beginning residue of target
    my ($querypos, $querychn ) = (1);  ## pdb residue number of the first aligned residue

    my  @template;        ## template pdb structure
    my  $sbbegresnum;     ## beginning residue of template
    my ($sbjctpos, $sbjctchn );        ## pdb residue number of the first aligned residue

    my  $copyqueryfasta = $queryfasta;


    if( $RUNTM ) {
        unless( exists $PDBHASH{$pdbqueryent} ) {
            return 0 unless ReadPDB( "$pdbqueryent.$ORGSUFFIX", \@{$PDBHASH{$pdbqueryent}} );
        }
        @target = @{$PDBHASH{$pdbqueryent}};
    }

    $$reftmscore = '--';
    $$reftmscorefrag = '--';

    return 1
        if( $RUNTM && 
           !ProcessFasta(
            $llquery,
            $llquerychain,
            \$queryfasta,
            1,             ## query !!
            \@target,
            \$querypos,
            \$querychn
        ));

    ## avoid such preallocation...
    ##$#AUXSBJCTS = $#{$copSUBJECTS};

    for( $s = 0; $s <= $#$copSUBJECTS; $s++ ) {
        $llsbjct       =  $$copSUBJECTS[$s]->[$NAME];
        $llsbjctchain  =  $$copSUBJECTS[$s]->[$CHNS];
        $llsbjctinfo   =  $$copSUBJECTS[$s]->[$INFO];
        $sbjctfasta    =  $$copSUBJECTS[$s]->[$ALIG];
        $pdbsbjctent   =  $$copSUBJECTS[$s]->[$ENTR];

        unless( exists $PDBHASH{$pdbsbjctent}) {
            return 0 unless ReadPDB( "$pdbsbjctent.$ORGSUFFIX", \@{$PDBHASH{$pdbsbjctent}});
        }

        $AUXSBJCTS[$s]->{TEMP} = [ @{$PDBHASH{$pdbsbjctent}} ]; ## assign value rather than reference
        $AUXSBJCTS[$s]->{NAME} = $llsbjct;     ## asign name
        $AUXSBJCTS[$s]->{CHNS} = $llsbjctchain;## asign chain ids
        $AUXSBJCTS[$s]->{INFO} = $llsbjctinfo; ## asign user info (e-value, scop...)
        $AUXSBJCTS[$s]->{ALIG} = $sbjctfasta;  ## asign alignment
        $AUXSBJCTS[$s]->{ALCO} = $sbjctfasta;  ## make a copy of alignment before processing!

        return 1
            if( !ProcessFasta(
                $llsbjct,
                $llsbjctchain,
                \$AUXSBJCTS[$s]->{ALIG},
                0,             ## not a query
                $AUXSBJCTS[$s]->{TEMP},
                \$sbjctpos,
                \$sbjctchn
            ));

        $AUXSBJCTS[$s]->{TPOS} = $sbjctpos;    ## beginning position in template structure
        $AUXSBJCTS[$s]->{TCHN} = $sbjctchn;    ## chain id of position in template structure
    }

    return 0 unless
    MakePIRFile( 
        $llquery,
        $querypos,
        $querychn,
        $llquerychain,
        $lluserinfo,
        \$queryfasta,
        \$copyqueryfasta,
        \@AUXSBJCTS
    );

    if( $queryfasta !~ /[a-zA-Z]/ ) {
        printf( STDERR "WARNING: %s sequence contains non residue symbols.\n", $llquery );
        return 0;
    }

    my  $nosbjcts = 1;
    do { if( $_ ) { $nosbjcts = 0; last } } foreach( @AUXSBJCTS );

    if( $nosbjcts ) {
        printf( STDERR "WARNING: No template sequences; target %s.\n", $llquery );
        return 0;
    }


    return 0 unless MakeScriptFile( $llquery, \@AUXSBJCTS );

    for( $s = 0; $s <= $#{$copSUBJECTS}; $s++ ) {
        next unless $AUXSBJCTS[$s];
        return 0 unless 
            AdjustTemplate(
                $AUXSBJCTS[$s]->{TEMP},
                $$copSUBJECTS[$s]->[$ENTR],
                $AUXSBJCTS[$s]->{TPOS},
                $AUXSBJCTS[$s]->{TCHN}
            );
    }

    if( RunModeller( $llquery ))
    {
        if( $RUNTM ) {
            if( AdjustModel( $llquery, $querypos, $querychn,
                             $llquerychain, 
                             $queryfasta,
                             $copyqueryfasta,
                             \@target, 
                             \@AUXSBJCTS ))
            {
                printf( STDERR "FINAL MODEL: %s\n", $modelname2 );
                AdjustTarget( \@target,  $pdbqueryent, $querypos, $querychn );
                ObtainTMscore( $llquery, $pdbqueryent, $reftmscore, $reftmscorefrag );
            }
        }
        else {
            printf( STDERR "FINAL MODEL: %s\n", $modelname );
        }
    }

    return 1;
}

## -------------------------------------------------------------------

sub InitData 
{
    %RESDS = (
        'ALA' => 'A',     'ARG' => 'R',     'ASN' => 'N',     'ASP' => 'D',
        'CYS' => 'C',     'GLU' => 'E',     'GLN' => 'Q',     'GLY' => 'G',
        'HIS' => 'H',     'ILE' => 'I',     'LEU' => 'L',     'LYS' => 'K',
        'MET' => 'M',     'PHE' => 'F',     'PRO' => 'P',     'SER' => 'S',
        'THR' => 'T',     'TRP' => 'W',     'TYR' => 'Y',     'VAL' => 'V',
    );
}

## -------------------------------------------------------------------
## get next residue in pdb file 
##

sub NextRes
{
    my  $refpdb = shift;  ## reference tp pdb data
    my  $chains = shift;  ## chain ids
    my  $nnn = shift;     ## reference to current record
    my  $resname = shift; ## reference to next residue
    my  $chainid = shift; ## ref. to chain id of next residue
    my  $resnum = shift;  ## reference to next residue's number
    my  $reshet = shift;  ## reference to flag of whether next residue is HETATM

    my  $prevnum = $$resnum;
    my  $first = $$nnn == 0;
    my  $nechain = 0;
    my  $atom = 0;


    for( ; $$nnn <= $#$refpdb && 
       ( !defined( $$resnum ) || ( $$resnum eq $prevnum )); )
    {
        $$nnn++;
        last if $#$refpdb < $$nnn;
        return 0 if $refpdb->[$$nnn] =~ /^(?:END|TER)/;
        next if $refpdb->[$$nnn] =~ /^(?:SIGATM|ANISOU|SIGUIJ)/;
        $atom = $refpdb->[$$nnn] =~ /^(?:ATOM|HETATM)/;
        next unless $atom;
        ## 18-20 Res name, 23-26 Res No., 27 Ins. code
        ##
        $$resname = substr( $refpdb->[$$nnn], 17, 3 );
        $$chainid = substr( $refpdb->[$$nnn], 21, 1 );
        $$resnum = substr( $refpdb->[$$nnn], 22, 5 );
        $$reshet = 0;
        $$reshet = 1 if $refpdb->[$$nnn] =~ /^HETATM/;
        $nechain = $chains &&( $chains !~ /$$chainid/ );
        $prevnum = $$resnum if $nechain;
    }
    return 0 if( $#$refpdb < $$nnn || $nechain );
    return $atom;
}

## -------------------------------------------------------------------
## extract one-letter aa sequence from pdb data
##

sub GetPDBSequence
{
    my  $refpdb = shift;  ## reference tp pdb data
    my  $sequen = shift;  ## reference to sequence to be constructed
    my  $chains = shift;  ## chain ids

    my  $reshash = \%RESDS;

    my  $resname = '';
    my  $chainid = '';
    my  $resnum = '';
    my  $reshet = '';
    my  $nrec = 0;

    unless( $refpdb && ref( $refpdb )) {
        printf( STDERR "ERROR: GetPDBSequence: Reference expected.\n");
        return 0;
    }
    unless( $sequen && ref( $sequen )) {
        printf( STDERR "ERROR: GetPDBSequence: Reference expected.\n");
        return 0;
    }
    $$sequen = '';

    while( NextRes( $refpdb, $chains, \$nrec, \$resname, \$chainid, \$resnum, \$reshet )) {
        if( $reshet || !exists $reshash->{$resname} ) {
            $$sequen .= 'X';
            next;
        }
        $$sequen .= $reshash->{$resname};
    }
    return 1;
}

## -------------------------------------------------------------------
## remove residue from pdb data given its pdb number
##

sub RemovePDBResidue
{
    my  $refpdb = shift;  ## reference tp pdb data
    my  $pdbnum = shift;  ## pdb-style number of residue to be removed
    my  $chains = shift;  ## chain ids

    my  $nextnum;
    my  $resname;
    my  $chainid;
    my  $resnum;
    my  $reshet;
    my  $nrec = 0;
    my  $coords = qr/^(?:ATOM|HETATM|SIGATM|ANISOU|SIGUIJ)/;

    while( NextRes( $refpdb, $chains, \$nrec, \$resname, \$chainid, \$resnum, \$reshet )) {
        if( $resnum eq $pdbnum && ( $chains? $chains =~ /$chainid/: 1 ))
        {
            for( ; $nrec <= $#$refpdb &&
                   $refpdb->[$nrec] =~ /$coords/ &&
                 ( $nextnum = substr( $refpdb->[$nrec], 22, 5 )) eq $resnum; )
            {
                splice( @$refpdb, $nrec, 1 );
            }
            last;
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## Process fasta sequence and return the beginning index in PDB style 
## entry
##

sub ProcessFasta
{
    my  $domain       = shift; ## domain identifier
    my  $chains       = shift; ## chain ids
    my  $reffasta     = shift; ## reference to alignment string in fasta
    my  $isquery      = shift; ## whether a query is about to be processed
    my  $refpdb       = shift; ## reference to whole pdb structure
    my  $refpos       = shift; ## to be computed
    my  $refchn       = shift; ## ref. to chain id of resnum

    my  $resname;
    my  $chainid;
    my  $resnum;
    my  $reshet;
    my  $nrec = 0;
    my ($m, $f, $p );
    my ($maa, $faa, $paa );
    my  @positrem;  ## serial numbers of residues to be removed from structure 

    undef $$refpos;

    my  $algnfasta;
    my  $algnpdbsq;
    my  $pdbseqn;
    my  $nwalign = nwalign->new();

    my  $rawseq = $$reffasta;
    $rawseq =~ s/[\-]//g;

    return 0 unless GetPDBSequence( $refpdb, \$pdbseqn, $chains );

    ## make alignment between sequence and corresponding structure
    ##
    $nwalign->Query( $rawseq );
    $nwalign->Sbjct( $pdbseqn );
    $nwalign->Align( \$algnfasta, \$algnpdbsq );
##
##print( STDERR "$rawseq\n$pdbseqn\naln:\n$algnfasta\n$algnpdbsq\n--\n\n");## *** *** ***
##
    if( length( $algnfasta ) != length( $algnpdbsq )) {
        printf( STDERR "ERROR: Wrong length of alignment beween sequence and structure: %s.\n", $domain );
        return 0;
    }

    ## process alignment sequence and structure
    ##
    for( $m = 0, $f = 0, $p = 0; $f < length( $algnfasta ); $f++, $p++ )
    {
        $faa = substr( $algnfasta, $f, 1 );
        $paa = substr( $algnpdbsq, $p, 1 );

        if( $paa ne '-' ) {
            unless( NextRes( $refpdb, $chains, \$nrec, \$resname, \$chainid, \$resnum, \$reshet )) {
                printf( STDERR "WARNING: Unable to find residue %s in structure; entry %s.\n", $paa, $domain );
                return 0;
            }
        }

        if( $faa ne '-' ) {
            $m++ while $m < length( $$reffasta ) && ( $maa = substr( $$reffasta, $m, 1 )) eq '-';

            if( $maa ne $faa ) {
                printf( STDERR "ERROR: Sequences do not match each other while processing %s: %s(%d)<>%s(%d).\n",
                        $domain, $maa, $m, $faa, $f );
                return 0;
            }
            if( $faa ne 'X' && $paa ne '-' && $paa ne 'X' && $faa ne $paa ) {
                printf( STDERR "WARNING: Sequence-structure inconsistency for %s: pos %d: %s <=> %s; ".
                               "Adjusted automatically.\n", $domain, $m, $faa, $paa );
                $faa = 'X';
            }

            unless(( defined( $$refpos ) || $paa eq '-' )) {
                $$refpos = $resnum;
                $$refchn = $chainid;
            }

            ## NOTE:
            ## when modeller models BLK residues . it assigns them to the 
            ## corresponding residues from template
            ## when it models X residues they are assigned to UNK residues
            ## with backbone atoms only
            ##
            SWITCH: {
                ## missing structural residue a sequence has
                if( $paa eq '-' )                            { substr( $$reffasta, $m, 1 ) = '-'; last SWITCH }
                if( $paa eq 'X' && $faa eq 'X' &&  $reshet ) { substr( $$reffasta, $m, 1 ) = '.'; last SWITCH }
                if( $paa eq 'X' && $faa eq 'X' && !$reshet ) { substr( $$reffasta, $m, 1 ) = 'X'; last SWITCH }
                if( $paa eq 'X' &&  $reshet )                { substr( $$reffasta, $m, 1 ) = '.'; last SWITCH }
                ##if( $paa eq 'X' && !$reshet )                { substr( $$reffasta, $m, 1 ) = 'X'; last SWITCH } ##NOTE:PREVIOUS CODE
                if( $paa eq 'X' && !$reshet )                { substr( $$reffasta, $m, 1 ) = '.'; last SWITCH }
                if( $faa eq 'X' ) { substr( $$reffasta, $m, 1 ) = substr( $algnpdbsq, $p, 1 ); last SWITCH }
                if( $faa ne $paa ){
                    printf( STDERR "WARNING: Wrong alignment beween sequence and structure; entry %s.\n", $domain );
                    return 0;
                }
            }

            $m++ if $m < length( $$reffasta );
        }
        else {
            if( $paa eq '-' ) {
                print( STDERR "WARNING: Symmetric deletions in sequence and structure alignment; entry %s.\n", $domain );
                next;
            }
            ## missing sequence residue a structure has.
            ## if not the beginning and trailing gaps...
            push @positrem, $resnum
                if( defined( $maa ) && $m < length( $$reffasta )&&( $chains? $chains =~ /$chainid/: 1 ));
        }
    }

    ## remove residues from structure, that sequence lacks
    ##
    for( $p = $#positrem; 0 <= $p; $p-- ) {
        RemovePDBResidue( $refpdb, $positrem[$p], $chains );
    }

    return 1;
}

## -------------------------------------------------------------------
## make PIR alignment file
##

sub MakePIRFile
{
    my  $llquery = shift;
    my  $querypos = shift;
    my  $querychn = shift;   ## chain id found in structure file
    my  $querychain = shift; ## chain ids given in description line of alignments
    my  $lluserinfo    = shift;
    my  $refqueryfasta = shift;
    my  $refcopyqueryfasta = shift;
    my  $refAUXSBJCTS = shift;

    my  $reshash = \%RESDS;

    my  $llsbjct;
    my  $sbjctpos;
    my  $sbjctchain;
    my  $refsbjctfasta;
    my  $refcopysbjctfasta;
    my ($qchain, $schain ) = (' ',' ');

    my ($s, $n, $qr, $sr );
    my  $filename = "$llquery.$PIREXT";

    if( length( $$refqueryfasta ) != length( $$refcopyqueryfasta )) {
        printf( STDERR "ERROR: Lengths of fasta sequences are inconsistent.\n" );
        return 0;
    }

    foreach( @{$refAUXSBJCTS} ) {
        if( length( $_->{ALIG} ) != length( $_->{ALCO} ) ||
            length( $_->{ALIG} ) != length( $$refqueryfasta ) ) {
            printf( STDERR "ERROR: Lengths of subject fasta sequences are inconsistent.\n" );
            return 0;
        }
    }

    my  $owncopyqueryfasta = $$refqueryfasta;
    my  @owncopysbjctfasta;

    push @owncopysbjctfasta, $_->{ALIG} foreach( @$refAUXSBJCTS );


    ## make here one more adjustment of fasta sequences in order modeller
    ## be happy; remove alignment positions in which query BLK 
    ## residues ('.') are aligned with gaps of template(s)
    ##

    for( $n = 0; $n < length( $owncopyqueryfasta ); $n++ ) 
    {
        $qr = substr( $owncopyqueryfasta, $n, 1 );
        ##NOTE:if inserted
        if( $qr =~ /[A-Z]/ && ( scalar(grep{/$qr/}values(%$reshash)) < 1 )) {
            ##assign X to non-standard residues 
            substr( $owncopyqueryfasta, $n, 1 ) = 'X';
        }
        next unless $qr eq '.';

        my $gaps = 1;

        foreach( @owncopysbjctfasta ) {
            $sr = substr( $_, $n, 1 );
            ## do NOT uncomment next line since modeller fails to find
            ## sequence in structure if HETATM is deleted
            do { $gaps = 0; last } if( $sr ne '-' ); ##&& ( $sr ne '.' );
        }
        next unless $gaps;

        ## alternatively, do not remove alignment position,
        ## assign X instead
        substr( $owncopyqueryfasta, $n, 1 ) = 'X';
        ##substr( $owncopyqueryfasta, $n, 1 ) = '';
        ##substr( $_, $n, 1 ) = '' foreach( @owncopysbjctfasta );
        ##$n--;
    }


    unless( open( PIR, ">$filename" )) {
        printf( STDERR "ERROR: MakePIRFile: Failed to open %s\n", $filename );
        return 0;
    }

    $qchain = substr( $querychain, 0, 1 ) if( $querychain && length( $querychain ));
    $qchain = $querychn if $querychn;

    printf( PIR ">P1;%s\n", $llquery );
    printf( PIR "sequence:%s:%-4s:%1s:%4s:%1s: : : : %s\n",
                $llquery,
                $querypos,
                $qchain, ' ',
                $qchain, "$lluserinfo"
    );
    printf( PIR "%s*\n", $owncopyqueryfasta );

    for( $s = 0; $s <= $#$refAUXSBJCTS; $s++ )
    {
        if( $owncopysbjctfasta[$s] !~ /[a-zA-Z]/ ) {
            printf( STDERR "WARNING: Sequence of template %s contains not residue symbols; omited.\n",
                    $$refAUXSBJCTS[$s]->{NAME} );
            ## delete that subject
            ##
            $$refAUXSBJCTS[$s] = 0;
            next;
        }
        $schain = ' ';
        $schain = substr( $$refAUXSBJCTS[$s]->{CHNS}, 0, 1 )
                      if( $$refAUXSBJCTS[$s]->{CHNS} && length( $$refAUXSBJCTS[$s]->{CHNS}));
        ## chain id found in pdb structure
        $schain = $$refAUXSBJCTS[$s]->{TCHN} if $$refAUXSBJCTS[$s]->{TCHN};

        printf( PIR ">P1;%s\n", $$refAUXSBJCTS[$s]->{NAME} );
        printf( PIR "structure:%s:%-4s:%1s:%4s:%1s: : : : %s\n",
                    $$refAUXSBJCTS[$s]->{NAME}, 
                    $$refAUXSBJCTS[$s]->{TPOS},
                    $schain, ' ',
                    $schain,
                   "$$refAUXSBJCTS[$s]->{INFO}"
        );
        printf( PIR "%s*\n", $owncopysbjctfasta[$s] );
    }

    close( PIR );
    return 1;
}

## -------------------------------------------------------------------
## make modeller python script file
##

sub MakeScriptFile
{
    my  $llquery = shift;
    my  $refAUXSBJCTS = shift;

    my  $locscript = $SCRIPT;
    my  $filename = "$llquery.$PYEXT";
    my  $llsubjects;

    foreach( @$refAUXSBJCTS ) {
        next unless $_;
        $llsubjects .= ', ' if $llsubjects;
        $llsubjects .= "'$_->{NAME}'";
    }

    $locscript =~ s/(\s+alnfile\s*=\s*)'([^\s']+)'/$1'$llquery.$PIREXT'/;
    $locscript =~ s/(\s+knowns\s*=\s*)('[^\s']+')/$1( $llsubjects )/;
    $locscript =~ s/(\s+sequence\s*=\s*)'([^\s']+)'/$1'$llquery'/;

    unless( open( PY, ">$filename" )) {
        printf( STDERR "ERROR: Failed to open %s\n", $filename );
        return 0;
    }
    print( PY $locscript );
    close( PY );
    return 1;
}

## -------------------------------------------------------------------
## run modeller
##

sub RunModeller
{
    my  $llquery = shift;
    my  $filename = "$llquery.$PYEXT";

    unless( -f $filename ) {
        printf( STDERR "ERROR: Missing modeller script file %s\n", $filename );
        return 0;
    }
    unless( RunCommand( "$MODELLER $filename" )) {
        printf( STDERR "ERROR: Modeller failed.\n" );
        return 0;
    }
    printf( STDERR "Modeller succeeded.\n" );
    return 1;
}

## -------------------------------------------------------------------
## adjust template structure by modifying alternation symbol so that 
## modeller'd understand
##

sub AdjustTemplate
{
    my  $template = shift;
    my  $pdbsbjctent = shift;
    my  $pdbsbjctstart = shift; ## residue number
    my  $pdbsbjctchain = shift; ## chain id

    my  @newtemplate = @$template;
    my ($rec, $begrec );
    my ($alt, $altnxt, $altprv );
    my ($ins, $insnxt, $insprv );
    my ($chn, $chnnxt, $chnprv );
    my  @altsyms;
    my ($resnum, $resnxt, $resprv, $largest );
    my  $startenumcheck = 0;
    my  $space = ' ';
    my  $coords = qr/^(?:ATOM|HETATM|SIGATM|ANISOU|SIGUIJ)/;


    ## process for alternate location indicators
    ##
    for( $rec = 0; $rec <= $#newtemplate; ) {
        last if( $startenumcheck && $newtemplate[$rec] !~ /$coords/ );
        do { $rec++; next; } if( $newtemplate[$rec] !~ /$coords/ );

        undef @altsyms;
        $altprv = $alt;
        $chnprv = $chn;
        $insprv = $ins;
        $resprv = $resnum;

        ## 17 Alternate location indicator
        ## 23-26 Res No., 27 Ins. code
        ##
        $altnxt = $alt = substr( $newtemplate[$rec], 16, 1 );
        $chnnxt = $chn = substr( $newtemplate[$rec], 21, 1 );
        $resnxt = $resnum = substr( $newtemplate[$rec], 22, 5 );
        $insnxt = $ins = substr( $resnum, 4, 1 );

        $begrec = $rec;
        $startenumcheck = 1 if( $resnum eq $pdbsbjctstart && $chn eq $pdbsbjctchain );

        for( ; $rec <= $#newtemplate && 
            $resnum eq $resnxt && $ins eq $insnxt && $chn eq $chnnxt; )
        {
            push @altsyms, $altnxt  if( $altnxt ne $space && !grep {/$altnxt/} @altsyms );

            $rec++;
            last unless( $rec <= $#newtemplate && $newtemplate[$rec] =~ /$coords/ );
            $altnxt = substr( $newtemplate[$rec], 16, 1 );
            $chnnxt = substr( $newtemplate[$rec], 21, 1 );
            $resnxt = substr( $newtemplate[$rec], 22, 5 );
            $insnxt = substr( $resnxt, 4, 1 );
        }

        if( 0 <= $#altsyms ) {
            for( my $altrec = $begrec; $altrec < $rec; $altrec++ ) {
                substr( $newtemplate[$altrec], 16, 1 ) = $space if 
                    substr( $newtemplate[$altrec], 16, 1 ) eq $altsyms[0];
            }
        }

        if( $startenumcheck ) {
            ## adjust numeration of residues if the next residue
            ## number is less than the previous one
            ##
            if( defined( $resprv ) && $resnum <= $resprv ) {
                $largest = $resprv;
                $resnum  = $largest + 1;
                for( my $numrec = $begrec; $numrec < $rec; $numrec++ ) {
                    substr( $newtemplate[$numrec], 22, 5 ) = ##along with insertion code
                        sprintf( "%4d%s", $resnum, substr( $newtemplate[$numrec], 26, 1 ));
                }
            }
        }
    }

    ## write adjusted template to file
    ##
    return 0 unless WritePDB( $pdbsbjctent, \@newtemplate );
    return 1;
}

## -------------------------------------------------------------------
## adjust target model with respect to native structure by changing 
## residue enumeration
##

sub AdjustModel
{
    my  $llquery = shift;
    my  $querypos = shift;
    my  $querychn = shift;   ## chain id from pdb structure
    my  $querychain = shift; ## chain ids from alignment file
    my  $queryfasta = shift;
    my  $copyqueryfasta = shift;
    my  $reftarget = shift;
    my  $refAUXSBJCTS = shift;

    my ($modseqn, $natseqn );
    my ($algnmod, $algnnat );
    my ($mresname, $mreschn, $mresnum, $mreshet );
    my ($sresname, $sreschn, $sresnum, $sreshet );
    my ($mprename, $mprechn, $mprenum, $mprehet );
    my ($sprename, $sprechn, $sprenum, $sprehet );
    my ($rename, $rechn, $renum, $rehet, $proced );
    my ($runchn, $runnum );
    my (@model, $modchn );
    my (@transmodel );
    my  $modelname = "$llquery.$MODELSUFFIX.$PDBSUFFIX";
    my  $outputname = "$llquery.$MODELSUFFIX.$TRANSUFFIX.$PDBSUFFIX";
    my  $coords = qr/^(?:ATOM|HETATM|SIGATM|ANISOU|SIGUIJ)/;
    my  $HETATM = 'HETATM';
    my  $HETLEN = length( $HETATM );

    my ($m, $s, $maa, $saa );
    my ($mrec, $srec, $rrec, $trec ) = (0,0,0,0);

    my  $nwalign = nwalign->new();

    unless( -f $modelname ) {
        printf( STDERR "WARNING: AdjustModel: No model %s found: Check modeller log file.\n", $modelname );
        return 0;
    }
    return 0 unless ReadPDB( $modelname, \@model, $modchn );

    return 0 unless GetPDBSequence( \@model, \$modseqn );
    return 0 unless GetPDBSequence( $reftarget, \$natseqn, $querychain ); ##may be compiled of several chains

    ## make alignment between model and native structure
    ##
    $nwalign->Query( $modseqn );
    $nwalign->Sbjct( $natseqn );
    $nwalign->Align( \$algnmod, \$algnnat );
##
##print( STDERR "AdjustModel:\n\n$modseqn\n$natseqn\naln:\n$algnmod\n$algnnat\n--\n\n");## *** ***
##
    if( length( $algnmod ) != length( $algnnat )) {
        printf( STDERR "ERROR: AdjustModel: Wrong length of alignment beween model and structure: %s.\n", $llquery );
        return 0;
    }

    ## process alignment of model and structure
    ##
    for( $m = 0, $s = 0; $m < length( $algnmod ); $m++, $s++ )
    {
        $maa = substr( $algnmod, $m, 1 );
        $saa = substr( $algnnat, $s, 1 );

        if( $saa ne '-' ) {
            $sprename = $sresname; $sprechn = $sreschn; $sprenum = $sresnum; $sprehet = $sreshet;
            unless( NextRes( $reftarget, $querychain, \$srec, \$sresname, \$sreschn, \$sresnum, \$sreshet )) {
                printf( STDERR "WARNING: AdjustModel: Unable to locate residue %s in structure; %s.\n", $saa, $llquery );
                return 0;
            }
        }

        if( $maa ne '-' ) {
            $mprename = $mresname; $mprechn = $mreschn; $mprenum = $mresnum; $mprehet = $mreshet;
            unless( NextRes( \@model, $modchn, \$mrec, \$mresname, \$mreschn, \$mresnum, \$mreshet )) {
                printf( STDERR "WARNING: AdjustModel: Unable to locate residue %s in model; %s.\n", $maa, $llquery );
                return 0;
            }
        }

        next if $maa eq '-';

        $proced = 0;
        $rename = '';
        $rechn = $sreschn;
        $renum = $sresnum;
        $rehet = $sreshet;

        if( $mresname eq 'UNK' && $saa eq '-' ) {
            next; ##if do nothing
            ##remove residue X aligned with gap
            for( $rrec = $mrec; $rrec <= $#model; $rrec++ ) {
                last if $model[$rrec] =~ /^(?:END|TER)/;
                next unless $model[$rrec] =~ /$coords/;

                $runchn = substr( $model[$rrec], 21, 1 );
                $runnum = substr( $model[$rrec], 22, 5 );
                last if $runnum ne $mresnum;

                splice( @model, $rrec, 1 );
                $proced = 1;
            }
            ##restore previous residue
            if( $proced ) {
                $mresname = $mprename; $mreschn = $mprechn; $mresnum = $mprenum; $mreshet = $mprehet;
                $mrec--;
            }
            next;
        }
        elsif( $maa ne $saa || $saa eq '-' ) {
            ##increase serial number of previous residue
            next unless $mprenum;
            $renum = sprintf( "%4d ", $mprenum + 1 );
            $rechn = $mprechn;
            $rehet = 0 if $saa eq '-' || $maa ne 'X';
        }
        elsif( $maa eq 'X' && $saa eq 'X' && $mresname ne $sresname ) {
            ##change residue name too
            $rename = $sresname;
        }

        for( $rrec = $mrec; $rrec <= $#model; $rrec++ )
        {
            if( $model[$rrec] =~ /$coords/ ) {
                $runchn = substr( $model[$rrec], 21, 1 );
                $runnum = substr( $model[$rrec], 22, 5 );
                last if $runnum ne $mresnum
            }
            $transmodel[$trec] = $model[$rrec];
            if( $model[$rrec] =~ /$coords/ )
            {
                substr( $transmodel[$trec], 17, 3 ) = $rename if $rename;
                substr( $transmodel[$trec], 21, 1 ) = $rechn  if $rechn;
                substr( $transmodel[$trec], 22, 5 ) = $renum  if $renum;
                substr( $transmodel[$trec], 0, $HETLEN ) = $HETATM if $rehet;
            }
            $trec++;
        }

# COMMENTED BLOCK BELOW TO BE REMOVED!
#        $proced = 0;
#        for( $rrec = $mrec; $rrec <= $#model; $rrec++ )
#        {
#            last if $model[$rrec] =~ /^(?:END|TER)/;
#            next unless $model[$rrec] =~ /$coords/;
#
#            $runchn = substr( $model[$rrec], 21, 1 );
#            $runnum = substr( $model[$rrec], 22, 5 );
#            last if $runnum ne $mresnum;
#
#            substr( $model[$rrec], 17, 3 ) = $rename if $rename;
#            substr( $model[$rrec], 21, 1 ) = $rechn  if $rechn;
#            substr( $model[$rrec], 22, 5 ) = $renum  if $renum;
#            substr( $model[$rrec], 0, $HETLEN ) = $HETATM if $rehet;
#            $proced = 1;
#        }
#        if( $proced ) {
#            $mresname = $rename if $rename;
#            $mreschn = $rechn   if $rechn;
#            $mresnum = $renum   if $renum;
#            $mreshet = $rehet   if $rehet;
#        }
    }

    return 0 unless WritePDB( $outputname, \@transmodel );
    return 1;
}

## -------------------------------------------------------------------
## adjust target model with respect to native structure by changing 
## residue enumeration
##

sub AdjustModelObs
{
    my  $llquery = shift;
    my  $querypos = shift;
    my  $querychn = shift;   ## chain id from pdb structure
    my  $querychain = shift; ## chain ids from alignment file
    my  $queryfasta = shift;
    my  $copyqueryfasta = shift;
    my  $reftarget = shift;
    my  $refAUXSBJCTS = shift;

    my  @sbjctfasta;
    my  @copysbjctfasta;

    my ($trec, $mrec, $s, $n ); 
    my ($dot, $msk, $gap );
    my ($tarreschn, $tarresnum, $tarresnxt );
    my ($modreschn, $modresnum, $modresnxt );
    my  @model;
    my  $modelname = "$llquery.$MODELSUFFIX.$PDBSUFFIX";
    my  $outputname = "$llquery.$MODELSUFFIX.$TRANSUFFIX.$PDBSUFFIX";
    my  $coords = qr/^(?:ATOM|HETATM|SIGATM|ANISOU|SIGUIJ)/;


    foreach( @$refAUXSBJCTS ) {
        next unless $_;
        push @sbjctfasta, $_->{ALIG};
        push @copysbjctfasta, $_->{ALCO};
    }

    unless( -f $modelname ) {
        printf( STDERR "WARNING: No model %s found: Check modeller log file.\n", $modelname );
        return 0;
    }
    return 0 unless ReadPDB( $modelname, \@model );

    ## get to the specified position of the target structure
    ##
    for( $trec = 0; $trec <= $#$reftarget; $trec++ ) {
        next unless $reftarget->[$trec] =~ /$coords/;
        ## 23-26 Res No., 27 Ins. code
        $tarreschn = substr( $reftarget->[$trec], 21, 1 );
        $tarresnum = substr( $reftarget->[$trec], 22, 5 );
        last if( $tarresnum eq $querypos && $tarreschn eq $querychn );
    }
    unless( $tarresnum ) {
        printf( STDERR "ERROR: Unable to find corresponding beginning residue of %s\n",  $llquery );
        return 0;
    }

    ## get to the first coordinate section of the model
    ##
    for( $mrec = 0; $mrec <= $#model; $mrec++ ) {
        next unless $model[$mrec] =~ /$coords/;
        $modreschn = substr( $model[$mrec], 21, 1 );
        $modresnum = substr( $model[$mrec], 22, 5 );
        last;
    }
    unless( $modresnum ) {
        printf( STDERR "ERROR: Unable to find beginning residue of %s\n", $modelname );
        return 0;
    }

    for( $n = 0; $n < length( $queryfasta ); ) {
        do { $n++; next; } if substr( $queryfasta, $n, 1 ) ne '-';
        substr( $queryfasta, $n, 1 ) = '';
        substr( $copyqueryfasta, $n, 1 ) = '';
        ## make fasta copies of originals compliant with modified ones
        for( $s = 0; $s <= $#sbjctfasta; $s++ ) {
            substr( $sbjctfasta[$s],     $n, 1 ) = '';
            substr( $copysbjctfasta[$s], $n, 1 ) = '';
        }
    }
    if( length( $queryfasta ) != length( $copyqueryfasta )) {
        printf( STDERR "ERROR: Lengths of modified and original sequences do not match.\n" );
        return 0;
    }
    for( $s = 0; $s <= $#sbjctfasta; $s++ ) {
        if( length( $sbjctfasta[$s] ) != length( $copysbjctfasta[$s] ) ||
            length( $sbjctfasta[$s] ) != length( $queryfasta )) {
            printf( STDERR "ERROR: Lengths of template sequences are not consistent.\n" );
        }
    }

    ## reenumerate model
    ##
    for( $n = 0; $mrec <= $#model && $model[$mrec] =~ /$coords/ && 
            $trec <= $#$reftarget && $reftarget->[$trec] =~ /$coords/; $n++ ) 
    {
        ## do NOT uncomment line for $gap; see above
        ## $msk is used if there's no need to model positions
        ## which were masked with Xs in given alignment hit sequences
        ##
        $dot =  substr( $queryfasta, $n, 1 ) eq '.';
        $gap =  1;
        do { if( substr( $_, $n, 1 ) ne '-' ) { $gap = 0; last } } foreach( @sbjctfasta );

        $msk = substr( $queryfasta, $n, 1 ) ne 'X' &&
               substr( $copyqueryfasta, $n, 1 ) eq 'X';

        for( $s = 0; $s <= $#sbjctfasta && !$msk; $s++ ) {
            $msk = substr( $sbjctfasta[$s], $n, 1 ) ne 'X' &&
                   substr( $copysbjctfasta[$s], $n, 1 ) eq 'X';
        }

        ## if target has HETATM aligned with gaps, leave model 
        ## position unchanged
        ##
        unless( $dot && $gap ) {
            ## process model enumeration
            ##
            for(; $mrec <= $#model && $model[$mrec] =~ /$coords/ &&
                ( $modresnxt = substr( $model[$mrec], 22, 5 )) eq $modresnum; )
            {
                ## if we have a missing part in the native pdb structure OR
                ## the original sequence contains masked X residue
                if( $dot || $msk ) { 
                    splice( @model, $mrec, 1 );
                }else {
                    substr( $model[$mrec], 21, 1 ) = $tarreschn;
                    substr( $model[$mrec], 22, 5 ) = $tarresnum;
                    $mrec++;
                }
            }
            $modresnum = $modresnxt;
        }

        ## move to the next target residue
        ##
        for(; $trec <= $#$reftarget && $reftarget->[$trec] =~ /$coords/ &&
            ( $tarreschn = substr( $reftarget->[$trec], 21, 1 )) &&
            ( $tarresnxt = substr( $reftarget->[$trec], 22, 5 )) eq $tarresnum;
              $trec++ )
        {
        } 
        $tarresnum = $tarresnxt;
    }

    $mrec++ while( $mrec <= $#model && $model[$mrec] !~ /$coords/ );
    if( $mrec <= $#model ) {
        printf( STDERR "ERROR: Not all residues processed from model %s\n", $modelname );
        return 0;
    }
    return 0 unless WritePDB( $outputname, \@model );
    return 1;
}

## -------------------------------------------------------------------
## adjust target so that its residue enumeration is sequentially 
## ordered disregarding insertion codes that TMscore does not 
## understand
##

sub AdjustTarget
{
    my  $target = shift;
    my  $pdbqueryent = shift;
    my  $pdbquerystart = shift; ## residue number from pdb structure
    my  $pdbquerychain = shift; ## chain id from pdb structure

    ## write adjusted target structure to file
    ##
    return 0 unless WritePDB( $pdbqueryent, $target );
    return 1;
}

## -------------------------------------------------------------------
## run TMscore given model and native structure; 
## do reverse superposition as well to obtain fragment-based TM-score
##

sub ObtainTMscore
{
    my  $llquery = shift;
    my  $pdbqueryent = shift;
    my  $reftmscore = shift;
    my  $reftmscorefrag = shift;

    my  $pdbqueryorg = "$pdbqueryent.$ORGSUFFIX";
    my  $modelname = "$llquery.$MODELSUFFIX.$TRANSUFFIX.$PDBSUFFIX";
    my  $supnameout = "$llquery.$RASMOLSUFFIX";
    my  $tmscoreout = "$llquery.$TMSUFFIX";
    my  $supnamefragout = "$llquery.$FRAGSUFFIX.$RASMOLSUFFIX";
    my  $tmscorefragout = "$llquery.$FRAGSUFFIX.$TMSUFFIX";
    my  $verify = 1;

    unless( -f $modelname ) {
        printf( STDERR "ERROR: No model found: %s\n", $modelname );
        return 0;
    }
    unless( -f $pdbqueryent ) {
        printf( STDERR "ERROR: No native structure file found: %s\n", $pdbqueryent );
        return 0;
    }

    return 0 unless RunTMscore( $modelname, $pdbqueryorg, $supnameout, $tmscoreout );
    return 0 unless ProcessTMscoreOutput( $tmscoreout, $reftmscore, $verify );


    return 0 unless RunTMscore( $pdbqueryorg, $modelname, $supnamefragout, $tmscorefragout );
    return 0 unless ProcessTMscoreOutput( $tmscorefragout, $reftmscorefrag, $verify );
    return 1;
}

## -------------------------------------------------------------------
## run TMscore
##

sub RunTMscore
{
    my  $modelname = shift;
    my  $pdbqueryent = shift;
    my  $superpname = shift;
    my  $tmscoreout = shift;

    unless( -f $modelname ) {
        printf( STDERR "ERROR: RunTMscore: File not found: %s\n", $modelname );
        return 0;
    }
    unless( -f $pdbqueryent ) {
        printf( STDERR "ERROR: RunTMscore: File not found: %s\n", $pdbqueryent );
        return 0;
    }
    unless( RunCommand( "$TMSCORE $modelname $pdbqueryent -o $superpname >$tmscoreout" )) {
        printf( STDERR "ERROR: TMscore failed.\n" );
        return 0;
    }
    printf( STDERR "TMscore succeeded.\n" );
    return 1;
}

## -------------------------------------------------------------------
## process output from TMscore and return the TM-score
## TODO: make verification more thorough with respect to
## sequence from PIR file
##

sub ProcessTMscoreOutput
{
    my  $filename = shift;
    my  $refscore = shift;
    my  $verify = shift;

    my  $modseq;
    my  $natseq;
    my  $begin = 0;
    my ($n, $mr, $nr, $pmr, $pnr );

    unless( -f $filename ) {
        printf( STDERR "ERROR: No TMscore output found: %s\n", $filename );
        return 0;
    }
    unless( open( TFD, $filename )) {
        printf( STDERR "ERROR: Unable to open %s\n", $filename );
        return 0;
    }
    while( <TFD> ) {
        chomp;
        if( /^TM\-score\s+=\s+([0-9\.eE\-\+]+)\s+/ ) {
            $$refscore = $1;
            last unless $verify;
            next;
        }
        do { $begin = 1; next;} if /^\(\":\"\s+denotes\s+the\s+residue\s+pairs/;
        next unless $begin;
        do { $modseq = $_  }      if $begin == 1;
        do { $natseq = $_; last;} if $begin == 3;
        $begin++;
    }
    close( TFD );
    return 1 unless $verify;

    if( length( $modseq ) != length( $natseq )) {
        printf( STDERR "ERROR: Lengths of model (%s) and native structures are inconsistent.\n",
                $filename );
        return 0;
    }
    unless( length( $modseq )) {
        printf( STDERR "WARNING: Length of model %s is 0 ( %s %s ).\n", $filename );
        $$refscore = '--';
        return 1;
    }

    for( $n = 0; $n <= length( $modseq ); $n++ ) {
        $pmr = $mr if defined( $mr );
        $pnr = $nr if defined( $nr );
        $mr = substr( $modseq, $n, 1 ) if $n < length( $modseq );
        $nr = substr( $natseq, $n, 1 ) if $n < length( $natseq );

        next unless( defined( $pmr ) && defined( $pnr ));
        next if( $pmr eq '-' || $pnr eq '-' );
        ## TMscore misaligns then model begins with
        ## residue with insertion code
        if( $pmr ne $pnr && ( $mr eq '-' || $nr eq '-' )) {
            printf( STDERR "WARNING: Misaligned residue No. %d in %s.\n", $n, $filename );
            next;
        }
        if( $pmr ne $pnr ) {
            printf( STDERR "WARNING: TMScore alignment %s discrepancies: pos %d: %s <> %s\n", 
                $filename, $n-1, $pmr, $pnr );
##            $$refscore = '--';
##            last;
        }
    }
    return 1;
}




## ===================================================================
## read PDB file into list given by reference
##

sub ReadPDB {
    my  $filename = shift;
    my  $reflines = shift;

    unless( $reflines && ref( $reflines )) {
        printf( STDERR "ERROR: ReadPDB: Reference expected.\n" );
        return 0;
    }
    unless( open( PFD, $filename )) {
        printf( STDERR "ERROR: ReadPDB: Failed to open: %s\n", $filename );
        return 0;
    }
    while( <PFD> ) {
        chomp;
        push @{$reflines}, $_;
    }
    close( PFD );
    return 1;
}

## -------------------------------------------------------------------
## write PDB information vectors to file
##

sub WritePDB {
    my  $filename = shift;
    my  $reflines = shift;

    unless( $reflines && ref( $reflines )) {
        printf( STDERR "ERROR: WritePDB: Reference expected.\n" );
        return 0;
    }
    unless( open( PFD, ">$filename" )) {
        printf( STDERR "ERROR: WritePDB: Failed to open %s\n", $filename );
        return 0;
    }
    print( PFD "$_\n" ) foreach( @{$reflines} );
    close( PFD );
    return 1;
}

## -------------------------------------------------------------------
## read Fasta sequence from file
##

sub ReadFasta {
    my  $filename = shift;
    my  $refseque = shift;
    my  $start = 0;

    unless( $refseque && ref( $refseque )) {
        printf( STDERR "ERROR: ReadFasta: Reference expected.\n" );
        return 0;
    }
    unless( open( FFD, $filename )) {
        printf( STDERR "ERROR: ReadFasta: Failed to open %s\n", $filename );
        return 0;
    }
    while( <FFD> ) {
        chomp;
        next if /^$/;
        last if /^>/ && $start;
        do { $start = 1; $$refseque = ''; next; } if /^>/;
        next unless $start;
        $$refseque .= uc( $_ );
    }
    close( FFD );
    $$refseque =~ s/[\n\r\s\-]//g;
    return 1;
}

## -------------------------------------------------------------------
## wrap Fasta sequence and write it to file
##

sub WrapFasta {
    my  $reffile = shift;
    my  $reffasta = shift;
    my  $padding = shift; ## padding at the beginning of each line
    my  $width = 60;
    my  $line;

    unless( $reffile && ref( $reffile )) {
        printf( STDERR "ERROR: WrapFasta: Reference expected.\n" );
        return 0;
    }
    unless( ref( $reffile ) eq 'GLOB' || ref( $reffile ) eq 'SCALAR' ) {
        printf( STDERR "ERROR: WrapFasta: Wrong reference.\n" );
        return 0;
    }
    unless( $reffasta && ref( $reffasta )) {
        printf( STDERR "ERROR: WrapFasta: Reference expected.\n" );
        return 0;
    }
    for( my $n = 0; $n < length( $$reffasta ); $n += $width ) {
        if( $n && $padding ) {
            $line = sprintf( "%${padding}s%s\n", ' ', substr( $$reffasta, $n, $width ));
        } else {
            $line = sprintf( "%s\n", substr( $$reffasta, $n, $width ));
        }
        if( ref( $reffile ) eq 'SCALAR' ) {
                 $$reffile .= $line;
        } else { print( $reffile $line );
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## run system command
##

sub CheckStatus
{
    RunCommand();
}

sub RunCommand
{
    my  $cmdline = shift;
    my  $retstatus = shift;

    system( "$cmdline" ) if $cmdline;

    if( $? == -1 ) {
        printf( STDERR "ERROR: Failed to execute command: $!\n" );
        return 0;
    }
    if( $? & 127 ) {
        printf( STDERR "ERROR: Command ran died with signal %d, %s coredump.\n",
                       ( $? & 127 ), ( $? & 128 )? 'with' : 'without' );
        return 0;
    }
    else {
        if(( $? >> 8 ) != 0 ) {
            unless( $retstatus ) {
                printf( STDERR "ERROR: Command failed and exited with status %d\n", $? >> 8 );
                return 0
            }
            return( $? >> 8 );
        }
    }
    return 1;
}

## <<>>

