(* ::Package:: *)

BeginPackage["PhononIrep`",{"SpaceGroupIrep`"}]

getcharandsymm ::usage = ""
phonopyIrepIO ::usage = ""
calcPhononIrep ::usage = "calcPhononIrep[{{\"supercell\"\[Rule]{a,b,c},
	\"unitcell\" \[Rule] path of unitcell cell strcuture file,
	\"force\" ->  path of force constants file,
	\"kset\"\[Rule]{{\!\(\*SubscriptBox[\(kx\), \(1\)]\),\!\(\*SubscriptBox[\(ky\), \(1\)]\),\!\(\*SubscriptBox[\(kz\), \(1\)]\)},{\!\(\*SubscriptBox[\(kx\), \(2\)]\),\!\(\*SubscriptBox[\(ky\), \(2\)]\),\!\(\*SubscriptBox[\(kz\), \(2\)]\)},...},
	\"showRep\"\[Rule] True or False,
	\"degeneracytolerance\" -> tolerance of determining the degeneracy in phonopy,
	\"symprec\" -> tolerance of determining the symmetry in phonopy] 
	generates the irreducible representations for phonon of the k points in \"kset\" option, where for opinions \"supercell\", \"unitcell\", \"force\", \"kset\" are shoule be set almost every calculations, but for \"degeneracytolerance\" and \"symprec\" the default values are OK.\"showRep\" is a opinion that tell PhononIrep whether to show the results in a Grid form.";


Begin["`Private`"];





ss=StartExternalSession["Python"];
ExternalEvaluate[ss,"
import numpy as np
import phonopy
from phonopy.structure.symmetry import Symmetry
from phonopy.interface.vasp import write_vasp
def getcharandsymm(supercell,unitcell_filename,force_constants_filename,kset,symprec,degeneracy_tolerance):
    chardict={}
    ph = phonopy.load(supercell_matrix=supercell,
            unitcell_filename=unitcell_filename,
            force_constants_filename=force_constants_filename,
            symprec=symprec)
    write_vasp("<>"'POSCARasdasd'"<>",ph._primitive)
    symm=Symmetry(ph._primitive, symprec=symprec).get_dataset()
    symm=dict((k,symm[k]) for k in ('number','rotations','translations','transformation_matrix') if k in symm)
    symm['rotations']=np.array(symm['rotations'],dtype='int64')
    chardict['symmetry']=symm
    for k in kset:
        ph.set_irreps(k,
                is_little_cogroup=False,
                nac_q_direction=None,
                degeneracy_tolerance=degeneracy_tolerance)
        char={}
        char['rotations']=np.array( ph._irreps._conventional_rotations,dtype='int64')
        char['q-position']=ph._irreps._q
        char['normal_modes']=[]
        for i, deg_set in enumerate(ph._irreps._degenerate_sets):
            char['normal_modes'].append({
                'band_indices': [i+1 for i in deg_set],
                'frequency': ph._irreps._freqs[deg_set[0]],
                'characters': ph._irreps._characters[i]})
        chardict[k]=char
    return chardict

def getcharandsymmfs(supercell,unitcell_filename,force_constants_filename,kset,symprec,degeneracy_tolerance):
    chardict={}
    ph = phonopy.load(supercell_matrix=supercell,
            unitcell_filename=unitcell_filename,
            primitive_matrix=[-1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2],
            force_sets_filename=force_constants_filename,
            symprec=symprec)
    symm=Symmetry(ph._primitive, symprec=symprec).get_dataset()
    symm=dict((k,symm[k]) for k in ('number','rotations','translations','transformation_matrix') if k in symm)
    symm['rotations']=np.array(symm['rotations'],dtype='int64')
    chardict['symmetry']=symm
    for k in kset:
        ph.set_irreps(k,
                is_little_cogroup=False,
                nac_q_direction=None,
                degeneracy_tolerance=degeneracy_tolerance)
        char={}
        char['rotations']=np.array( ph._irreps._conventional_rotations,dtype='int64')
        char['q-position']=ph._irreps._q
        char['normal_modes']=[]
        for i, deg_set in enumerate(ph._irreps._degenerate_sets):
            char['normal_modes'].append({
                'band_indices': [i+1 for i in deg_set],
                'frequency': ph._irreps._freqs[deg_set[0]],
                'characters': ph._irreps._characters[i]})
        chardict[k]=char
    return chardict

"];

    
Options[getcharandsymm] = {
   "supercell" -> {1, 1, 1},
   "unitcell" ->  "POSCAR-unitcell",
   "force" ->  "FORCE_CONSTANTS",
   "kset" -> {{0, 0, 0}},
   "symprec" -> 0.0001,
   "degeneracytolerance" -> 0.0001};
getcharandsymm[OptionsPattern[]] := Module[{
   supercell = OptionValue["supercell"],
   unitcell = OptionValue["unitcell"],
   force = OptionValue["force"],
   kset = N@OptionValue["kset"],
   symprec = N@OptionValue["symprec"],
   degeneracytolerance = N@OptionValue["degeneracytolerance"],output},
output=ExternalEvaluate[
   ss, <|"Command" -> "getcharandsymm", 
    "Arguments" -> {supercell, unitcell, force, kset, symprec, 
      degeneracytolerance}|>];
 output["symmetry"]= Rationalize/@Normal/@output["symmetry"];
(* Print[output];*)
output
];


phonopyIrepIO[char_] := Module[{nele, SymmData, CharData,
   nsymm, SOC = 0, SymmetryOperation, NK, kvector,
   PositionofPointGroupElement, BandChar,
   OutputStr,
   qpoints,
   tm
   },
  SymmData = char["symmetry"];
  qpoints = Drop[Keys[char], 1];
  (*Print[qpoints];*)
  tm=Rationalize[SymmData["transformation_matrix"],0.000001];
  nsymm = Length@SymmData["rotations"];
  SymmetryOperation = 
   Flatten/@Transpose@{#rotations,Mod[Round[N@#translations,0.00001],1] }&[SymmData];
   (*Print[SymmetryOperation];*)
  kvector = <||>;
  BandChar = <||>;
  PositionofPointGroupElement = <||>;
  Do[
   
   CharData = char[q];
   nele = Last[Last[#"band_indices" & /@ CharData["normal_modes"]]];
   
   kvector[q] = Normal@CharData["q-position"];
   
   (*Print[SymmData["rotations"],CharData["rotations"]];*)
   PositionofPointGroupElement[q] = 
    FirstPosition[# & /@ 
        SymmData["rotations"], #] & /@
     (Simplify[Inverse[tm] . # . tm] & /@ Normal@CharData["rotations"]);
    (*Print[PositionofPointGroupElement,"\n",Normal@CharData["rotations"],"\n",SymmData["rotations"],"\n",
    tm
    ];*)
   PositionofPointGroupElement[q] = 
    Flatten[PositionofPointGroupElement[q]];
   BandChar[
     q] = {#"band_indices", #frequency, {N[Re[#], 10], N[Im[#], 10]}&/@ Chop@Normal[#characters]} & /@ 
     CharData["normal_modes"];
     (*Print[BandChar];*)
   BandChar[q] = MapAt[{#[[1]], Length[#]} &, BandChar[q], {;; , 1}];
   BandChar[q] = Flatten /@ BandChar[q];
   , {q, qpoints}];
  (*Print[BandChar];*)
  OutputStr = StringRiffle[{nele}, " "];
  OutputStr = 
   OutputStr <> "\n" <> ToString[SOC] <> "\n" <> ToString[nsymm] <> 
    "\n";
  NK = Length[qpoints];
  Do[
   OutputStr = OutputStr <> str <> " 1 0 0 0 0 0 1 0" <> "\n",
   {str, StringRiffle[#, " "] & /@ SymmetryOperation}];
  OutputStr = OutputStr <> ToString[NK] <> "\n";
  Do[
   OutputStr = OutputStr <> StringRiffle[kvector[q], " "] <> "\n";
   , {q, qpoints}];
  Do[
   OutputStr = 
    OutputStr <> ToString[Length[PositionofPointGroupElement[q]]] <> 
     "\n";
   OutputStr = 
    OutputStr <> StringRiffle[PositionofPointGroupElement[q], " "] <> 
     "\n";
     (*Print[PositionofPointGroupElement[q]];*)
   Do[
    OutputStr = OutputStr <> str <> "\n",
    {str, StringRiffle[#, " "] & /@ BandChar[q]}];
   , {q, qpoints}];
  OutputStr;
  Export[FileNameJoin[{$TemporaryDirectory, "trace.txt"}], OutputStr];
  ];


Options[calcPhononIrep] = {
   "supercell" -> {1, 1, 1},
   "unitcell" ->  "POSCAR-unitcell",
   "force" ->  "FORCE_CONSTANTS",
   "kset" -> {{0, 0, 0}},
   "symprec" -> 0.0001,
   "degeneracytolerance" -> 0.0001,
   "showRep"->False};
calcPhononIrep[OptionsPattern[]] := Module[{
   supercell = OptionValue["supercell"],
   unitcell = OptionValue["unitcell"],
   force = OptionValue["force"],
   kset = N@OptionValue["kset"],
   symprec = N@OptionValue["symprec"],
   degeneracytolerance = N@OptionValue["degeneracytolerance"],
   showRep = OptionValue["showRep"],
   char, tracedata, tr, rep1
   },
  char = getcharandsymm[
    "supercell" -> supercell,
    "unitcell" -> unitcell,
    "force" -> force,
    "kset" -> kset,
    "symprec" -> symprec,
    "degeneracytolerance" -> degeneracytolerance];
  phonopyIrepIO[char];
  tracedata = 
   readVasp2trace[FileNameJoin[{$TemporaryDirectory, "trace.txt"}]];
   (*tr = autoConvTraceToBC["C:\\Users\\zhang\\Documents\\POSCARasdasd", tracedata];*)
  tr = autoConvTraceToBC[unitcell, tracedata];
  rep1 = getBandRep[char["symmetry"]["number"], "", 
    tr["trace"]];(*;
  MatrixForm /@ rep1["rep"]*)
  Print["SG No:",char["symmetry"]["number"]];
  If[showRep,showIrep[rep1];rep1,rep1]
  ];


showIrep[repoutput_]:=Module[{kpath,rep,kinfo,head},
{kpath,rep,kinfo}=repoutput/@Keys[repoutput];
(*Print[kinfo];*)
Do[
head=PadRight[{kinfo[[n]][[;;2]]},4,SpanFromLeft];
Print@Grid[Join[{head},{{"Band","Frequency","dim","Irep"}},rep[[n]]],Frame->All];
,{n,Length[kinfo]}]
]


End[];
EndPackage[];



