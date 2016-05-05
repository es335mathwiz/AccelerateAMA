(* Mathematica Package *)
BeginPackage["AccelerateAMA`", {"ProtectedSymbols`","JLink`","SymbolicAMA`", "NumericAMA`", "AMAModel`","Experimental`", "Format`","AMAModelDefinition`"}]


getModelDims::usage="getModelDims[modDir_String,modName_String] returns {Length[vars],getLags[eqns],getLeads[eqns],Length[params],linearQ[eqns]}"

getParams::usage="getParams[modName_String]"

parseMod::usage="parseMod[srcDir_String,fName_String,targDir_String]\n" <>
"parse dynare model to xml to mathematia  AMAModelDefinition";



sameEqns::usage = "sameEqns  "

linearQ::usage = "linearQ  "

mkNewDir::usage = "mkNewDir  "

firstOnPath::usage = "firstOnPath  "

collectData::usage = "collectData  "

getLags::usage = "getLags  "

getLeads::usage = "getLeads  "

allLinear::usage = "allLinear  "

allNonLinear::usage = "allNonLinear  "

trySolveSS::usage = "trySolveSS  "

makeSomeSSSubs::usage = "makeSomeSSSubs  "

makeSSValSubs::usage = "makeSSValSubs  "

compEigSpace::usage = "compEigSpace  "

ssSubs::usage = "ssSubs  "

preEvals::usage = "preEvals  "

trySolve::usage = "trySolve  "

tryFindRoot::usage = "tryFindRoot  "


tryCompEvals::usage = "tryCompEvals  "

allSSSolvedQ::usage = "allSSSolvedQ  "

uniqueSSSoln::usage = "uniqueSSSoln  "

tryCompEvecs::usage = "tryCompEvecs  "

genSSEqns::usage = "genSSEqns  "


(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

makeSSValSubs[vars_List]:=#[t_]->ToExpression[ToString[#]<>"SSVal"]&/@vars;



defineDouble[name_String]:=
"double "<>name<>";\n"
$defaultExamplesDir="dynareExamples/uniqueExamples/";
$defaultResDir="theLinRes/";




preEvals[modName_String]:=preEvals[$defaultExamplesDir,modName,$defaultResDir]

preEvals[theDir_String,modName_String,targDir_String]:=
Module[{parseTime,vars,ig,params,eqns,notSubs,paramSubs,
hmatTime,hmat,
arTime,zf,hf,
amatTime,amat,
lilTime,lilMat,cols},
Print["analyzing ",modName];
{parseTime,{vars,ig,params,ig,{ig,eqns},notSubs,ig}}=
Timing[parseMod[theDir,modName,targDir]];
Global`getAMAEqns[modName]=eqns;
Global`getAMAVars[modName]=vars;
paramSubs=#[[1]]->#[[2]]&/@params;
Global`getParamSubs[modName]=paramSubs;
Global`getParamNames[modName]=First/@params;
Global`getExampleParams[modName]=(First /@params)//.paramSubs;
{hmatTime,hmat}=Timing[equationsToMatrix[eqns,vars]/.makeSSValSubs[vars]];
Global`getHmat[modName]=hmat;
{arTime,{zf, hf}} = Timing[symbolicAR[hmat]];
Global`getZf[modName]=zf;
{amatTime,amat} = Timing[symbolicTransitionMatrix[hf]];
{lilTime,{lilMat,cols}}=
Timing[symbolicEliminateInessentialLags[{amat,Range[Length[amat]]}]];
Global`getLilMat[modName]=lilMat;
Global`getCols[modName]=cols;
(*
{more,lilMat}=Timing[FullSimplify[lilMat,TimeConstraint->Global`$tConst]];lilTime=lilTime+more;*)
{parseTime,hmatTime,arTime,amatTime,lilTime,paramSubs,eqns,hmat,vars}
]



tryCompEvals[modName_String]:=tryCompEvals[$defaultExamplesDir,modName,$defaultResDir]


tryCompEvals[theDir_String,modName_String,targDir_String]:=
Module[{},	Print["working on ", modName];
{evalsTime,evals}= Timing[TimeConstrained[
	Eigenvalues[Transpose[Global`getLilMat[modName]//.With[{subs=Global`getSolveSS[modName]},
		If[Or[subs==={},subs===$Aborted],{},subs[[1]]]]]],Global`$tConst]];
	Global`getEvals[modName]=evals;giveUpQ[evals,"after evals"];
	Print["done evals"];
{evalsTime,evals}]


tryCompEvecs[modName_String]:=tryCompEvecs[$defaultExamplesDir,modName,$defaultResDir]

tryCompEvecs[theDir_String,modName_String,targDir_String]:=
Module[{hmat=Global`getHmat[modName]},	Print["working on ", modName];
	{lilevecsTime,lilevecs}=Timing[TimeConstrained[compEigSpace[Global`getLilMat[modName],Global`getEvals[modName],
		Join[Global`getParamSubs[modName],Global`getFRootSS[modName]]],Global`$tConst]];giveUpQ[lilevecs,"after evecs"];
{evecsTime,evecs}=Timing[toLarge[lilevecs,Global`getCols[modName],qRowLength[modName]]];
Global`getEvecs[modName]=evecs;
	Print["done evecs"];qmat=compQ[Global`getZf[modName],evecs];
{bmatTime,bmat}=Timing[compB[qmat]];
	Print["done bmat"];
getB[modName]=bmat;
	{sTime,theS}=Timing[obStruct[hmat,bmat]];
	Global`getS[modName]=theS;
		Print["done obstruct"];
{evecsTime,bmatTime,sTime,evecs,bmat,theS}]

qRowLength[modName_String]:=With[{hm=Global`getHmat[modName]},Length[hm[[1]]]-Length[hm]]

allSSSolvedQ[vars_List,soln_List]:=Length[vars]===Length[soln];
uniqueSSSoln[solnS_List]:=Length[solnS]==1


trySolve[modName_String]:=trySolve[$defaultExamplesDir,modName,$defaultResDir]


trySolve[theDir_String,modName_String,targDir_String]:=
Module[{},	Print["working on ", modName];
{solveTime,solveSoln}=Timing[trySolveSS[Global`getAMAEqns[modName]/.Global`AMAModelDefinition[modName][[-2]],Global`getAMAVars[modName]]];
Global`getSolveSS[modName]=solveSoln;
{solveTime,solveSoln}]


tryFindRoot[modName_String]:=tryFindRoot[$defaultExamplesDir,modName,$defaultResDir]

tryFindRoot[theDir_String,modName_String,targDir_String]:=
Module[{},
Print["working on ", modName];
{fRootTime,fRootSoln}=Timing[makeSomeSSSubs[modName]];
Global`getFRootSS[modName]=fRootSoln;
{fRootTime,fRootSoln}]

giveUpQ[expr_,msg_String]:=If[Not[FreeQ[expr,$Aborted]],Throw[expr,msg]]

zapInnovSubs=eps[_][_]->0;

trySolveSS[modName_String]:=trySolveSS[Global`getAMAEqns[modName],Global`getAMAVars[modName]]
trySolveSS[eqns_List,vars_List]:=
Module[{},
ssEqns=eqns//.ssSubs;
someSSSubs=TimeConstrained[Solve[Thread[ssEqns==0],vars],Global`$tConst];
If[FreeQ[someSSSubs,Solve],someSSSubs/.
(#->ToExpression[ToString[#]<>"SSVal"]&/@vars),{}
]]


ssSubs={eps[_][_]->0,xx_[t+_.]:>xx};


genSSEqns[modName_String]:=With[{theVars=Sort[Global`getAMAVars[modName]]},
	With[{ssValVarSubs=makeSSValSubs[theVars]},
		With[{ssEqns=(Global`getAMAEqns[modName]/.ssValVarSubs)/.ssSubs,
			ssEqnVars=Through[theVars[t]]/.ssValVarSubs},
		{ssEqnVars,ssEqns,diffRow[#,ssEqnVars]&/@ssEqns}]]]

diffRow[aRow_,vars_List]:=D[aRow,#]&/@vars

makeTry[name_String]:=
With[{vars=Global`getAMAVars[name]},
	With[{ssVarNames=Last/@makeSSValSubs[vars],
	solveSS=(Global`getSolveSS[name])//.Global`getParamSubs[name]/.Global`AMAModelDefinition[name][[-2]]/.zapInnovSubs},
	If[solveSS==={}||solveSS===$Aborted,
{ToExpression[ToString[#]<>"SSVal"],Random[]}&/@vars,
prepSolvesForFRoot[solveSS[[1]],ssVarNames]]]]

prepSolvesForFRoot[theSubs_List,allVars_List]:=With[{subbed=theSubs/.Rule->List},
	With[{unAssn=Complement[allVars,First/@subbed]},
		With[{randAssn=#->Random[]& /@unAssn},Join[subbed/.randAssn,randAssn/.Rule->List]]]]

makeSomeSSSubs[name_String]:=
With[{ssEqns=makeSubbedEqns[name],
try=makeTry[name]},
FindRoot @@ {ssEqns,try}]




doNullSpace[lilMat_?MatrixQ,eval_]:=
If[eval==={},{},With[{theRes=
NullSpace[Transpose[lilMat]-eval*IdentityMatrix[Length[lilMat]]]},
	If[theRes==={},Throw[{},"after evals"],theRes]]]

compEigSpace[lilMat_?MatrixQ,evals_List,paramSubs_List] := 
   Join @@ (doNullSpace[lilMat, evals[[#1]]]) &  /@ 
largeLocs[evals//.paramSubs]


makeSSValSubs[vars_List]:=#[t_]->ToExpression[ToString[#]<>"SSVal"]&/@vars;

largeLocs[theVals_List]:=
With[{mags=(#>1)&/@Abs[theVals]},Flatten[Position[mags,True]]]


obStruct[hmat_?MatrixQ,bigB_?MatrixQ]:=
With[{neq=Length[hmat],lTau=Length[bigB[[1]]],lTheta=Length[bigB]},
With[{hMinus=hmat[[All,Range[lTau+neq]]],
hPlus=hmat[[All,lTau+neq+Range[lTheta]]]},
hMinus+blockMatrix[{{zeroMatrix[neq],hPlus . bigB}}]]]


compQ[zf_?MatrixQ,{}]:=zf

compQ[zf_?MatrixQ,evs_?MatrixQ]:=
With[{qmat=Join[zf,evs]},qmat]

compB[qmat_?MatrixQ]:=
With[{qcols=Length[qmat[[1]]],qrows=Length[qmat]},
With[{qr=qmat[[All,qcols-qrows+Range[qrows]]],
ql=qmat[[All,Range[qcols-qrows]]]},
(-Inverse[qr].ql)]]







toLarge[{},lilCols_List,cols_Integer]:={}

toLarge[lil_?MatrixQ,lilCols_List,cols_Integer]:=
Module[{bigEvecs},
bigEvecs=ConstantArray[0,{Length[lil],cols}];
bigEvecs[[All,lilCols]]=lil;
bigEvecs]

getLags[eqns_List] := With[{}, 
    With[{lgld = Union[Cases[eqns, (xx_)[yy:t + (zz_.)] -> zz, Infinity]]}, 
     -Min[lgld]]]
getLeads[eqns_List] := With[{}, 
    With[{lgld = Union[Cases[eqns, (xx_)[yy:t + (zz_.)] -> zz, Infinity]]}, 
     Max[lgld]]]



parseMod[srcDir_String,fName_String,targDir_String]:=
Module[{cmd},
JavaNew["gov.frb.ma.msu.DynareToAMAModel",
srcDir<>fName<>".mod",targDir<>fName<>".xml",fName];
cmd=If[$OperatingSystem=="Unix",
       StringForm[
"java " <> "-cp "<>$jarDir<>"  org.apache.xalan.xslt.Process -IN `3``2`.xml  -XSL /msu/home/m1gsa00/RES2/mathAMA/AndersonMooreAlgorithm/AndersonMooreAlgorithm/AMAModel2Mma.xsl -OUT `3``2`.mth",srcDir,fName,targDir],
              StringForm[
"java " <> "-cp "<>$jarDir<>"  org.apache.xalan.xslt.Process -IN `3``2`.xml  -XSL g:/RES2/mathAMA/AndersonMooreAlgorithm/AndersonMooreAlgorithm/AMAModel2Mma.xsl -OUT `3``2`.mth",srcDir,fName,targDir]];
Run[cmd];
Get[targDir<>fName<>".mth"];
Global`AMAModelDefinition[fName]]


mkNewDir[dirName_String]:=If[Not[FileExistsQ[dirName]],CreateDirectory[dirName]]
firstOnPath[dirName_String]:=If[System`$Path[[1]]=!=dirName,PrependTo[System`$Path,dirName]]


$jarDir=If[$OperatingSystem=="Unix",
           "/msu/res1/Software/xalan-j_2_7_1/xalan.jar",
           "r:/Software/xalan-j_2_7_1/xalan.jar"]

$tmpDir=$TemporaryDirectory <> "/GaryModDims/";
If[Not[FileExistsQ[$tmpDir]], CreateDirectory[$tmpDir]]
getModelDims[modDir_String,modName_String]:=
Module[{vars,ig,params,eqns,notSubs,tDir=$tmpDir},
mkNewDir[tDir];firstOnPath[tDir];
{vars,ig,params,ig,{ig,eqns},notSubs,ig}=parseMod[modDir,modName,tDir];
System`$Path=Drop[System`$Path,1];
{Length[vars],getLags[eqns],getLeads[eqns],Length[params],linearQ[eqns]}]




doMod[preDo_String,fName_String]:=doMod[preDo,fName,preDo]


doMod[preDo_String,fName_String,targDir_String]:=
Module[{cmd},
System`$Path=PrependTo[System`$Path,targDir];
JavaNew["gov.frb.ma.msu.DynareToAMAModel",preDo<>fName<>".mod",targDir<>fName<>".xml",fName];
homeDir=If[Global`windowsQ[],"g:","/msu/home/m1gsa00"];
System.out.println("osname="+nameOS);
cmd=StringForm[
"java " <> "-cp "<>$jarDir<>" org.apache.xalan.xslt.Process -IN `3``2`.xml  -XSL `1`/RES2/mathAMA/AndersonMooreAlgorithm/AndersonMooreAlgorithm/AMAModel2Mma.xsl -OUT `3``2`.mth",homeDir,fName,targDir];
Run[cmd];
Get[fName<>".mth"];
Global`AMAModelDefinition[fName]]


mkNewDir[dirName_String]:=If[Not[FileExistsQ[dirName]],CreateDirectory[dirName]]
firstOnPath[dirName_String]:=If[System`$Path[[1]]=!=dirName,PrependTo[System`$Path,dirName]]


collectData[modDir_String,modName_String]:=
Module[{vars,ig,params,eqns,notSubs,tDir=$tmpDir},
mkNewDir[tDir];firstOnPath[tDir];
{vars,ig,params,ig,{ig,eqns},notSubs,ig}=parseMod[modDir,modName,tDir];
System`$Path=Drop[System`$Path,1];
{Length[vars],getLags[eqns],getLeads[eqns],Length[params],linearQ[eqns]}]


linearQ[eqns_List]:=FreeQ[equationsToMatrix[eqns],t]




sameEqns[modDir_String,modNameA_String,modNameB_String]:=
Module[{vars,ig,params,eqns,notSubs,tDir=$tmpDir},
mkNewDir[tDir];firstOnPath[tDir];
{vars,ig,params,ig,{ig,eqnsA},notSubs,ig}=parseMod[modDir,modNameA,tDir];
{vars,ig,params,ig,{ig,eqnsB},notSubs,ig}=parseMod[modDir,modNameB,tDir];
System`$Path=Drop[System`$Path,1];
eqnsA === eqnsB]
 

allLinear[theDir_String,modName_String,targDir_String]:=
Module[{parseTime,vars,ig,params,eqns,notSubs,paramSubs,
hmatTime,hmat,
arTime,zf,hf,
amatTime,amat,
lilTime,lilMat,cols,
evalsTime,evals,
lilevecsTime,lilevecs,
evecsTime,evecs,
bmatTime,bmat,
sTime,theS},
Print["parsing " <>modName];
{parseTime,{vars,ig,params,ig,{ig,eqns},notSubs,ig}}=
Timing[parseMod[theDir,modName,targDir]];
paramSubs=#[[1]]->#[[2]]&/@params;
Print["gen hmat"];
{hmatTime,hmat}=Timing[equationsToMatrix[eqns,vars]];
Print["gen symbolicAR"];
{arTime,{zf, hf}} = Timing[symbolicAR[hmat]];
Print["gen amat"];
{amatTime,amat} = Timing[symbolicTransitionMatrix[hf]];
Print["shrink hmat"];
{lilTime,{lilMat,cols}}=
Timing[symbolicEliminateInessentialLags[{amat,Range[Length[amat]]}]];
Print["gen eigenvalues"];
{evalsTime,evals}= Timing[Eigenvalues[Transpose[lilMat]]];
Print["gen evecs"];
{lilevecsTime,lilevecs}=Timing[compEigSpace[lilMat,evals,paramSubs]];
{evecsTime,evecs}=Timing[toLarge[lilevecs,cols,Length[zf[[1]]]]];
Print["compute bmat"];qmat=compQ[zf,evecs];
{bmatTime,bmat}=Timing[compB[qmat]];
{phifmatTime,{bagain,phimat,fmat}}=Timing[symbolicComputeBPhiF[hmat,qmat]];
Print["gen smat"];
{sTime,theS}=Timing[obStruct[hmat,bmat]];
{parseTime,hmatTime,arTime,amatTime,lilTime,evalsTime,lilevecsTime,evecsTime,bmatTime,sTime,phifmatTime,paramSubs,eqns,bmat,theS,hmat,vars,phimat,fmat}
]

allNonLinear[theDir_String,modName_String,targDir_String]:=
Module[{parseTime,vars,ig,params,eqns,notSubs,paramSubs,
hmatTime,hmat,
arTime,zf,hf,
amatTime,amat,
lilTime,lilMat,cols,
evalsTime,evals,
lilevecsTime,lilevecs,
evecsTime,evecs,
bmatTime,bmat,
sTime,theS},
{parseTime,{vars,ig,params,ig,{ig,eqns},notSubs,ig}}=
Timing[parseMod[theDir,modName,targDir]];
paramSubs=#[[1]]->#[[2]]&/@params;
{solveTime,solveSoln}=Timing[trySolveSS[eqns,vars]];
{fRootTime,fRootSoln}=Timing[makeSomeSSSubs[modName]];
{hmatTime,hmat}=Timing[equationsToMatrix[eqns,vars]/.makeSSValSubs[vars]];
{arTime,{zf, hf}} = Timing[symbolicAR[hmat]];Print["done ar"];
{amatTime,amat} = Timing[symbolicTransitionMatrix[hf]];
{lilTime,{lilMat,cols}}=
Timing[symbolicEliminateInessentialLags[{amat,Range[Length[amat]]}]];Print["done inessential"];
{more,lilMat}=Timing[FullSimplify[lilMat,TimeConstraint->Global`$tConst]];lilTime=lilTime+more;
{evalsTime,evals}= Timing[TimeConstrained[Eigenvalues[Transpose[lilMat]],Global`$tConst]];Print["done evals"];
{lilevecsTime,lilevecs}=Timing[TimeConstrained[compEigSpace[lilMat,evals,Join[paramSubs,fRootSoln]],Global`$tConst]];
{evecsTime,evecs}=Timing[toLarge[lilevecs,cols,Length[zf[[1]]]]];
qmat=compQ[zf,evecs];
{bmatTime,bmat}=Timing[compB[qmat]];
{sTime,theS}=Timing[obStruct[hmat,bmat]];
{parseTime,hmatTime,arTime,amatTime,lilTime,evalsTime,lilevecsTime,evecsTime,bmatTime,sTime,solveTime,fRootTime,paramSubs,eqns,bmat,theS,solveSoln,fRootSoln,hmat,vars}
]




End[] (* End Private Context *)

EndPackage[]  
Print["done reading AccelerateAMA"]
