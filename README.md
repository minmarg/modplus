modplus, a sequence-structure alignment and protein modelling 
helper tool

(C)2007-2019 Mindaugas Margelevicius, 
Institute of Biotechnology, Vilnius University

# Description

  modplus is a helper tool to model the target protein sequence using
  template structures and MODELLER (Sali and Blundell, J Mol Biol, 
  234(3), 779-815, 1993). It resolves inconsistencies that may arise
  between template sequences and structures. The correspondence
  between a model and the native target structure is optionally 
  calculated with TM-score (Zhang and Skolnick, Proteins, 57(4),
  702-710, 2004).

# Input

   The input file contains aligned sequences in FASTA:

```
>targetname [chains];   [optional user information]
G-----VDILR-MDAVAF-------IWKQMGTSCE ...
>templ1name [chains];   [optional user information]
FANYDEHLIFEGMNEPRLVGHANEWWPELTNSDVV ...
>templ2name [chains];   [optional user information]
EMLAIAVETKPHFCCLVPEKRQEVTTEGGLDV--A ...
> ...
...
```

By default, the first sequence is the target to be modelled. 
Optionally, other sequences (templates) may be modelled as well 
(see below). One or more template sequences (two or more sequences 
in total) can be listed in the file. A target is modelled using all 
the template structures and alignments listed in the file.

 The names targetname, templ1name, templ2name, ... have to 
correspond to the filenames (without extension) of PDB structures.
Required PDB files with extension `.ent` or `.pdb` added are looked 
up in the directory specified by the option `--pdb`.

 Template chains are extracted from pdb files. There may be several 
chain identifiers given by `chains` (e.g., SCOP domains compiled: 
ABC). If no `chains` are given, the first one is identified and 
used automatically.

 If option `--all` is specified, each sequence in the input file is 
modelled using the other sequence(s) as template(s).

 Option `--tm` implies a model quality check with TM-Score. The 
TM-Score output files will appear in the output directory specified 
by the option `--dir`. An optional filename following the `--tm` 
option specifies the following format for the resulting TM-scores 
written in this file: 

(TM: `<TMScore1>` `<TMScore2>`)

`<TMScore1>` is a score obtained with respect to the whole 
structure, and `<TMScore2>` is a score obtained along the alignment 
extent.

Contact: <mindaugas.margelevicius@bti.vu.lt>

