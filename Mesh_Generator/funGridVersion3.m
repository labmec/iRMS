(* ::Package:: *)

FilterClosePoints[points_, distanceperm_]:=Block[{fpoints=points, 
fdistanceperm=distanceperm,i,npoints,distance,j,mindistance,newlist},
newlist=fpoints;
npoints=Length[newlist];
For[i=1,i<= npoints,i++,

For[j=i+1,j<= npoints,j++,
distance= Norm[newlist[[i]]-newlist[[j]]];
If[distance< fdistanceperm,
If[i!=j,
newlist= Delete[newlist,j];
npoints=npoints-1;
j--;
];
];
];

];
Return[newlist];
];



GenQuadGrid[data_,hx_,hy_]:=Block[{fdata=data, fhx=hx, fhy=hy,minx,miny,maxx,maxy,
deltax, deltay,xpoints,ypoints,genPoints,i,j,valx,valy,genIndextope,genIndexBase,conttope,
contbase},
 minx=Min[fdata[[1;;All,1]]];
 miny=Min[fdata[[1;;All,2]]];
 maxx=Max[fdata[[1;;All,1]]];
 maxy=Max[fdata[[1;;All,2]]];
 
 deltax = (maxx - minx)/(fhx);
 deltay = (maxy - miny)/(fhy);

 genPoints=ConstantArray[0,{fhx+1,fhy+1}];
 genIndextope=ConstantArray[0,{fhx+1,fhy+1}];
 genIndexBase=ConstantArray[0,{fhx+1,fhy+1}];
 conttope=0;
 contbase=(fhx+1)*(fhy+1);
 For[j=1, j<= fhy+1, j++,
    valy = maxy - (j-1)*deltay;
    For[i=1, i<= fhx+1, i++,
    conttope++;
    contbase++;
    valx= minx + (i-1)*deltax;  
      genPoints[[j,i]]={valx,valy,2020};
      genIndextope[[j,i]]= conttope;
      genIndexBase[[j,i]]= contbase;
    ];
 ];
 
 Return[{genPoints,genIndextope,genIndexBase}];
];




GenMultivariateInterForm[data_]:=Block[{fdata=data,nvals, newtable, pointform,i},
nvals = Length[fdata];
newtable={};
For[i=1,i<= nvals,i++,
	pointform={{fdata[[i,1]],fdata[[i,2]]},fdata[[i,3]]};
	AppendTo[newtable,pointform];
];
Return[newtable];
];



FilterInternalPoints[Region_,data_,interpoltope_,interpolbase_]:=Block[{fRegion=Region, fdata=data,
npoints,newlist,fhx,fhy,boundaries,posx,posj,index,linesindex,i,j,CordinatesMatrixtope,CordinatesMatrixbase,
pointsBaseIndex,isInternal,pointsTopeIndex,comparationPoint,cordinatesPoints,
indexquadtope,indexquadbase,finterpoltope=interpoltope,finterpolbase=interpolbase,point,valInterpoltope,valInterpolbase},
fhx= Dimensions[fdata[[1]]][[1]];
fhy= Dimensions[fdata[[1]]][[2]];
CordinatesMatrixtope = fdata[[1]];
CordinatesMatrixbase = fdata[[1]];

pointsTopeIndex = fdata[[2]];
pointsBaseIndex=fdata[[3]];
 newlist = {}; 
 cordinatesPoints={};
 indexquadtope={};
 indexquadbase={};
 For[j=1, j<= fhy, j++,
    For[i=1, i<= fhx, i++,
      comparationPoint=CordinatesMatrixtope[[j,i]][[1;;2]];
      isInternal = RegionMember[fRegion,comparationPoint];
      If[!isInternal,
         pointsBaseIndex[[j,i]]= 0;
         pointsTopeIndex[[j,i]]= 0;,
         (*calculo de la coordenada z para el tope y base*)
         point=CordinatesMatrixtope[[j,i]];
         valInterpoltope= finterpoltope[point[[1]],point[[2]]];
      
        If[!(valInterpoltope>1550 && valInterpoltope < 2000),
         
             valInterpoltope= 1750;
         ];
        
         
         valInterpolbase= finterpolbase[point[[1]],point[[2]]];
     
        If[valInterpolbase>2120 || valInterpolbase < 1500,
        
             valInterpolbase= 1810;
         ];
          
          
         CordinatesMatrixtope[[j,i,3]]=valInterpoltope;
         CordinatesMatrixbase[[j,i,3]]=valInterpolbase;
         
     
      ];
      
    ];
 ];

 Return[{CordinatesMatrixtope,CordinatesMatrixbase,pointsTopeIndex,pointsBaseIndex} ];
];



LineExist[lines_,points_]:=Block[{flines=lines, fpoints=points,npoints,i},
 npoints = Length[flines];
 If[npoints==0, Return[0]];
 For[i=1, i<= npoints, i++,
 If[flines[[i,2]][[1]] == fpoints[[1]] && flines[[i,2]][[2]] == fpoints[[2]],
    Return[flines[[i,1]]];
    Break[];
 ];
 If[ flines[[i,2]][[2]] == fpoints[[1]] && flines[[i,2]][[1]] == fpoints[[2]],
     Return[-flines[[i,1]]];
    Break[];
 ];
 
 ];
 Return[0];
];



GenerateElementLines[lines_, points_,IniLineIndex_]:=Block[{index,fIniLineIndex=IniLineIndex,nlines,fpoints=points,flines=lines,
i,newLen,candidates,ElementLineIndex,max},
max=4;

  nlines = Length[flines];
  For[i=1, i<= max, i++,
    If[fpoints[[i]]==0,
    fpoints=Delete[fpoints,i];
    i--;
    max--;
    ];
  ];
  AppendTo[fpoints,fpoints[[1]]];
  
  newLen= Length[fpoints];
  ElementLineIndex={};
  For[i=1, i<newLen,i++,
    candidates = {fpoints[[i]],fpoints[[i+1]]};
    index = LineExist[flines, candidates ];
    If[index==0 && nlines!=0,
    AppendTo[flines,{flines[[-1,1]]+1, candidates}];
    index=flines[[-1,1]];
    
    nlines++;
    ];
    If[index==0 && nlines==0,
    AppendTo[flines,{fIniLineIndex, candidates}];
    index=fIniLineIndex;
    nlines++;
    ];
    AppendTo[ElementLineIndex,index];
  ];
  Return[{ElementLineIndex,flines}];
];



SetAttributes[CreateSurfaceConnectivities, HoldAll];
CreateSurfaceConnectivities[matrixnodeindex_,iniSurfaceIndex_, IniLineIndex_, CordinatesMatrix_]:=Block[{fmatrixnodeindex=matrixnodeindex,nnodex,nnodey,surfacematrix
,finiSurfaceIndex=iniSurfaceIndex, surfacecorners,fIniLineIndex=IniLineIndex,linestest,countsurface,i,j,CurrentIndex,RIndex, RIIndex,
IIndex,vecIndex,count,indexLines,indexandcandidate,ElementLineIndexVec},

   nnodex = Dimensions[fmatrixnodeindex][[1]];
   nnodey = Dimensions[fmatrixnodeindex][[2]];
   surfacematrix=ConstantArray[0,{nnodex,nnodey}];
   surfacecorners=ConstantArray[0,{nnodex,nnodey}];
   linestest={};
   countsurface=finiSurfaceIndex-1;
   ElementLineIndexVec={};
   For[j=1, j<=nnodex,j++,
    For[i=1, i<=nnodey,i++,
     If[fmatrixnodeindex[[j,i]]!=0,
        CurrentIndex=fmatrixnodeindex[[j,i]];
        RIndex = fmatrixnodeindex[[j,i+1]];
        RIIndex = fmatrixnodeindex[[j+1,i+1]];
        IIndex = fmatrixnodeindex[[j+1,i]];
        vecIndex={CurrentIndex,RIndex,RIIndex,IIndex};
        count=Count[vecIndex,0];
        If[count >= 2,
        Continue[];
        ];
 
	(*	If[count == 0,
		  
          CordinatesMatrix[[j+1,i+1,3]] = ForceCoplanarity[{j,i},CordinatesMatrix];
          
        ];*)
        
        indexandcandidate = GenerateElementLines[linestest,vecIndex,fIniLineIndex ];
        
        indexLines = indexandcandidate[[1]];
        linestest=indexandcandidate[[2]];
       
        countsurface++;
        surfacematrix[[j,i]]={countsurface,indexLines};
        surfacecorners[[j,i]]=vecIndex;
        AppendTo[ElementLineIndexVec,{countsurface,indexLines}];
     ];
   ];
   
   ];
  Return[{surfacematrix,linestest,ElementLineIndexVec,fmatrixnodeindex,surfacecorners}];
];




ForceCoplanarity[currentPoint_,CordinatesMatrix_]:=Block[{fcurrentPoint=currentPoint,j,i,CurrentCordPoint, RCordPoint,
 RICordPoint, ICordPoint, AA, BB, CC, den, Correction,valor},
        j=fcurrentPoint[[1]];
        i=fcurrentPoint[[2]];
        CurrentCordPoint=CordinatesMatrix[[j,i]];
      
        RCordPoint = CordinatesMatrix[[j,i+1]];
        RICordPoint = CordinatesMatrix[[j+1,i+1]];
        ICordPoint = CordinatesMatrix[[j+1,i]];
        
        AA =  RCordPoint - CurrentCordPoint;
        BB =  ICordPoint - CurrentCordPoint;
        CC =  RICordPoint - CurrentCordPoint;
   
        den= (AA[[1]]*BB[[2]] - AA[[2]]*BB[[1]]);
        Correction = ((AA[[3]]*BB[[2]] - AA[[2]]*BB[[3]])*CC[[1]] + (AA[[1]]*BB[[3]] - AA[[3]]*BB[[1]])*CC[[2]] )/(den);
        
        valor = CordinatesMatrix[[j+1,i+1]][[3]] + Correction;
        
        Return[valor];
        
];



CalcPointsAndIndex[PointsCordinatesMatrix_, IndexCordinatesMatrix_]:=Block[{dim, nx,ny, PointsCordinatesVec,PointsIndexVec,i,j },
dim = Dimensions[IndexCordinatesMatrix];
nx = dim[[1]];
ny= dim[[2]];
PointsCordinatesVec={};
PointsIndexVec={};
For[i=1, i<= nx, i++,
 For[j=1, j<= ny, j++,
    If[IndexCordinatesMatrix[[i,j]]!=0,
    AppendTo[PointsCordinatesVec, PointsCordinatesMatrix[[i,j]]];
    AppendTo[PointsIndexVec, IndexCordinatesMatrix[[i,j]]];
    ];
 ];
];


Return[{PointsCordinatesVec,PointsIndexVec }];
];



CreateInternalSurfConnectAndVolumes[topeinformation_, baseinformation_,wellVector_]:= Block[{ftopeinformation=topeinformation,
fbaseinformation=baseinformation,dim,nelx,nely,VolSurfaceConnectivities, SurfaceTOP, SurfaceBASE,LinesTope,LinesBase,
 AllLinesTopeAndBase,AllSurfaces,VolVec,countvol,TopeMatrixNodeIndex, BaseMatrixNodeIndex, countSurfaceIndex, i,j,
  CurrentIndexTOP, CurrentIndexBASE, RIndexTOP,RIIndexTOP,IIndexTOP, vecIndexTOP, RIndexBASE,RIIndexBASE,IIndexBASE,vecIndexBASE,
  count,VolIndexAndSurface,numfaces,numlines,vecIndexNodesTOP,vecIndexNodesBASE,TopeMatrixNodeIndexREAL,BaseMatrixNodeIndexREAL,
  comparacion,neigLinesSup, neigLinesLef, fwellVector=wellVector},
  TopeMatrixNodeIndex = ftopeinformation[[1]];
  BaseMatrixNodeIndex = fbaseinformation[[1]];

  TopeMatrixNodeIndexREAL = ftopeinformation[[5]];
  BaseMatrixNodeIndexREAL = fbaseinformation[[5]];
  (*SetTopeSurface and Lines*)
 (* SetBaseSurface and Lines*)
  
  dim = Dimensions[TopeMatrixNodeIndex];
  nelx = dim[[1]];
  nely = dim[[2]];
  VolSurfaceConnectivities=ConstantArray[0,{nelx,nely}];
  SurfaceTOP=ftopeinformation[[3]];
  SurfaceBASE=fbaseinformation[[3]];
  LinesTope=ftopeinformation[[2]];
  LinesBase=fbaseinformation[[2]];
  AllLinesTopeAndBase=Join[LinesTope,LinesBase];
  AllSurfaces=Join[SurfaceTOP,SurfaceBASE];

  numfaces=Length[AllSurfaces];
  numlines=Length[AllLinesTopeAndBase];
  
  VolVec={};
  countvol=0;
  countSurfaceIndex = SurfaceBASE[[-1,1]];
  For[j=1, j<= nely, j++,
  
    For[i=1, i<= nelx, i++,
 
        If[Length[TopeMatrixNodeIndex[[j,i]]] !=0,
        
        vecIndexTOP =TopeMatrixNodeIndex[[j,i,2]];
        vecIndexBASE =BaseMatrixNodeIndex[[j,i,2]];
        vecIndexNodesTOP= TopeMatrixNodeIndexREAL[[j,i]];
        vecIndexNodesBASE=BaseMatrixNodeIndexREAL[[j,i]];
        
       neigLinesSup={};
        neigLinesLef={};
     
        comparacion={};
        If[Length[VolSurfaceConnectivities[[j-1,i]]]!=0,
         neigLinesSup = VolSurfaceConnectivities[[j-1,i]];
       
        ];
        If[Length[VolSurfaceConnectivities[[j,i-1]]]!=0,
         neigLinesLef = VolSurfaceConnectivities[[j,i-1]];
          
        ];
        
        If[Length[neigLinesSup]!=0 && Length[neigLinesLef] != 0,
          comparacion = Join[neigLinesSup,neigLinesLef];
        ];
        If[Length[neigLinesSup]==0 && Length[neigLinesLef] != 0,
          comparacion = neigLinesLef;
        ];
         If[Length[neigLinesSup]!=0 && Length[neigLinesLef] == 0,
          comparacion = neigLinesSup;
        ];
        
      
        VolIndexAndSurface = GenerateElementSurfaces[AllSurfaces, AllLinesTopeAndBase,vecIndexTOP, 
        vecIndexBASE,vecIndexNodesTOP, vecIndexNodesBASE, comparacion, comparacionsup ];
        AllSurfaces =VolIndexAndSurface[[2]];
        AllLinesTopeAndBase =VolIndexAndSurface[[3]];
        countvol++;
        VolSurfaceConnectivities[[j,i]] = {countvol, VolIndexAndSurface[[1]]};
        AppendTo[VolVec,{countvol, VolIndexAndSurface[[1]]}];
       ];
       ];  
  
  ];
  nlines=Length[AllLinesTopeAndBase];
  nfaces=Length[AllSurfaces];
  Return[{VolSurfaceConnectivities,VolVec, AllLinesTopeAndBase[[numlines+1;;nlines]],AllSurfaces[[numfaces+1;;nfaces]] }];
];



SurfaceExist[Allsurfaces_,SurfaceIndexLines_]:=Block[{fAllsurfaces=Allsurfaces, fSurfaceIndexLines=SurfaceIndexLines,nSidescomp,
nSurfaces, i, count, candidateSides,nsidesCandidate, j,k, comSide },
   nSidescomp = Length[fSurfaceIndexLines];
   nSurfaces = Length[fAllsurfaces];
   
   For[i=1, i<= nSurfaces, i++,
        count =0;
        candidateSides = fAllsurfaces[[i,2]];
        nsidesCandidate = Length[candidateSides];
        
        If[nsidesCandidate!=nSidescomp,
        Continue[];
        ];
        For[j=1, j<= nsidesCandidate, j++,
            comSide = candidateSides[[j]];
            For[k=1, k<=nSidescomp, k++,
                If[Abs[comSide]==Abs[fSurfaceIndexLines[[k]]], (*no importa la orientacion*)
                count++;   
                ];
              ];
              If[count==0,
              Continue[];
              ];
        ];
        
        If[count==nSidescomp,
          Return[fAllsurfaces[[i,1]]];
          Break[];
        ];
   ];
   
   Return[0];
];




GenerateElementSurfaces[AllSurfacesVec_, AllLinesVec_, TopeElemenPointsVec_, 
BaseElemenPointsVec_, TOPPOINTS_, BASEPOINTS_, comparacion_, comparsup_]:=Block[{fAllSurfacesVec=AllSurfacesVec,
fAllLinesVec=AllLinesVec, fTopeElemenPointsVec=TopeElemenPointsVec, fBaseElemenPointsVec=BaseElemenPointsVec,  vecTope,vecTopeVerif, 
i, vecBase, volsurfaces,  countsurface, indexandcandidate,indexLines, surExistQ, vecIndex, 
fIniLineIndex, surfindex,vecPointsTop, vecPointsBase, line, point,FTOPPOINTS=TOPPOINTS, FBASEPOINTS=BASEPOINTS,
fcomparacion=comparacion, fcomparsup=comparsup,lines},
   
  
   volsurfaces={SurfaceExist[fAllSurfacesVec, fTopeElemenPointsVec], SurfaceExist[fAllSurfacesVec, fBaseElemenPointsVec]};
   countsurface=fAllSurfacesVec[[-1,1]];
   fIniLineIndex =fAllLinesVec[[-1,1]]+1;
   (*TOMA EL TOPE Y LA BASE*)
   vecPointsTop=Abs[FTOPPOINTS];
   vecPointsBase=Abs[FBASEPOINTS];
   vecPointsTop=DeleteCases[vecPointsTop, 0];
   vecPointsBase=DeleteCases[vecPointsBase, 0];
 
   
   AppendTo[vecPointsBase,vecPointsBase[[1]]];
   AppendTo[vecPointsTop,vecPointsTop[[1]]];
   
 
   For[i=1,i< Length[vecPointsTop],i++,
   
        vecIndex={vecPointsTop[[i]],vecPointsTop[[i+1]],vecPointsBase[[i+1]],vecPointsBase[[i]]};
        (*lines = ConvertLineIndexToPointIndex[fAllSurfacesVec,fAllLinesVec, Length[fAllLinesVec]+1]*);
        indexandcandidate = GenerateElementLines[fAllLinesVec,vecIndex, 1000]; 
        indexLines = indexandcandidate[[1]];
        surfindex = SurfaceExist[fAllSurfacesVec,indexLines];
        If[surfindex==0,
         fAllLinesVec=indexandcandidate[[2]];
         countsurface++;
         surfindex=countsurface;
         AppendTo[fAllSurfacesVec,{countsurface,indexLines}];
        ];
        AppendTo[volsurfaces,surfindex];
   ];
   Return[{volsurfaces,fAllSurfacesVec,fAllLinesVec }];
   
];



GenGmshPoint[point_,i_,h_]:=Block[{fpoint=point,fi=i,fh=h,val},
val=5;
StringJoin["Point(",ToString[fi],")= {",ToString[DecimalForm[fpoint[[1]]]],", ",
ToString[DecimalForm[fpoint[[2]]]],", ",ToString[DecimalForm[fpoint[[3]]]],", ",
ToString[fh],"};"]];




GenListGmshPoints[points_,indexes_,hs_]:=Block[{fpoints = points, findexes= indexes,
fhs=hs,fnpoints,i,pointsTransforms},
pointsTransforms={};
fnpoints = Length[fpoints];
For[i=1,i<=  fnpoints,i++,
AppendTo[pointsTransforms,GenGmshPoint[fpoints[[i]],findexes[[i]],fhs[[i]]]];
];
Return[pointsTransforms];
];




GenListGmshLine[indexpointtoLine_,indexline_]:=Block[{findexpointtoLine=indexpointtoLine,findexline=indexline, nPoints,Lines,i,
line,index2,index1,indexl},
nPoints = Length[findexpointtoLine];  
Lines={};
	For[i=1,i<= nPoints,i++,
	index1=findexpointtoLine[[i,1]];
	index2=findexpointtoLine[[i,2]];
	indexl = findexline[[i]];
	line =GenGmshLine[index1,index2,indexl];
	AppendTo[Lines,line];
	];
	Return[Lines];
];





GenListGmshSurface[indexLineSurface_,indexline_, wells_, surfacesWithWells_]:=Block[{findexLineSurface=indexLineSurface,findexline=indexline, 
nPoints,Surfaces,i,indexs,surface},
nPoints = Length[findexLineSurface];  
count=1;
Surfaces={};
	For[i=1,i<= nPoints,i++,
	    indexs = findexline[[i]];
	    Print[indexs];
	    If[indexs==surfacesWithWells[[count]],
	        Print["pase"];
	        surface =GenGmshSurfaceWithWell[findexLineSurface[[i]],indexs, wells[[count]]];
	        AppendTo[Surfaces,surface[[1]]];
		    AppendTo[Surfaces,surface[[2]]];
		    AppendTo[Surfaces,surface[[3]]];
		    If[count<Length[surfacesWithWells],
		    count++;
		    ];
		    
	    ,
	    surface =GenGmshSurface[findexLineSurface[[i]],indexs];
		AppendTo[Surfaces,surface[[1]]];
		
		AppendTo[Surfaces,surface[[2]]];
	    ];
		
	];
	Return[Surfaces];
];




GenGmshSurface[indexlines_, indexl_]:=Block[{findexlines=indexlines, findexl=indexl,line,surface,
npoints,lineloop},
 npoints= Length[findexlines];
 If[npoints==3,
 lineloop = StringJoin["Line Loop(",ToString[findexl],")= {",ToString[findexlines[[1]]],", ",
 ToString[findexlines[[2]]],",", ToString[findexlines[[3]]],"};" ];,
 lineloop = StringJoin["Line Loop(",ToString[findexl],")= {",ToString[findexlines[[1]]],", ",
 ToString[findexlines[[2]]],",", ToString[findexlines[[3]]],",",ToString[findexlines[[4]]],"};" ];
 ];
 surface = StringJoin["Surface(",ToString[findexl],")= {",ToString[findexl],"};"];
 

 Return[{lineloop,surface}];
];



GenGmshSurfaceWithWell[indexlines_, indexl_, well_]:=Block[{findexlines=indexlines, findexl=indexl,line,surface,
npoints,lineloop},
 npoints= Length[findexlines];
 
 lineloop1 = StringJoin["Line Loop(",ToString[findexl],")= {",ToString[findexlines[[1]]],", ",
 ToString[findexlines[[2]]],",", ToString[findexlines[[3]]],",",ToString[findexlines[[4]]],"};" ];
 
  lineloop2 = StringJoin["Curve Loop(",ToString[findexl+1500000],")= {",ToString[well[[1]]],", ",
 ToString[well[[2]]],",", ToString[well[[3]]],",",ToString[well[[4]]],"};" ];
 
 surface = StringJoin["Surface(",ToString[findexl],")= {",ToString[findexl], ", ",  ToString[findexl + 1500000],"};"];
 

 Return[{lineloop1,lineloop2,surface}];
];


GenGmshVolume[indexlines_, indexl_]:=Block[{findexlines=indexlines, findexl=indexl,line,surface,
npoints,lineloop},
 npoints= Length[findexlines];
 If[npoints==5,
 lineloop = StringJoin["Surface Loop(",ToString[findexl],")= {",ToString[findexlines[[1]]],", ",
 ToString[findexlines[[2]]],",", ToString[findexlines[[3]]],",",ToString[findexlines[[4]]],",",ToString[findexlines[[5]]],"};" ];,
 
 lineloop = StringJoin["Surface Loop(",ToString[findexl],")= {",ToString[findexlines[[1]]],", ",
 ToString[findexlines[[2]]],",", ToString[findexlines[[3]]],",",ToString[findexlines[[4]]],",",ToString[findexlines[[5]]],
 ",",ToString[findexlines[[6]]],"};" ];
 ];
 surface = StringJoin["Volume(",ToString[findexl],")= {",ToString[findexl],"};"];
 Return[{lineloop,surface}];
];




GenGmshLine[index1_,index2_,indexl_]:=Block[{findex1=index1,findex2=index2, findexl=indexl,line},
 line = StringJoin["Line(",ToString[findexl],")= {",ToString[findex1],", ",ToString[findex2],"};"];
 Return[line];
];



ConfigurateCurveLoop[data_,indexlineloop_]:=Block[{fdata=data,findexlineloop=indexlineloop,out},
   out=StringJoin["Curve Loop(",ToString[findexlineloop],")=",ToString[fdata],";"];
   Return[out];
];



GenListGmshVolume[indexSurfacesVol_,indexVol_]:=Block[{findexSurfacesVol=indexSurfacesVol,findexVol=indexVol, 
nPoints,Surfaces,i,indexs,volume},
nPoints = Length[findexSurfacesVol];  
Surfaces={};
	For[i=1,i<= nPoints,i++,
	    indexs = 1000+findexVol[[i]];
		volume =GenGmshVolume[findexSurfacesVol[[i]],indexs];
		
		AppendTo[Surfaces,volume[[1]]];
		AppendTo[Surfaces,volume[[2]]];
	];
	Return[Surfaces];
];



TakeSurfacesIndexBoundaries[volumesMatrix_]:=Block[{dim, nx,ny,i,j, CurrentVol, VolSup, VolInf, VolR, VolI,indexSurfaces },
dim = Dimensions[volumesMatrix];
nx= dim[[1]];
ny= dim[[2]];

indexSurfaces={};
For[i=1, i<= nx, i++,
For[j=1, j<= ny, j++,

   If[Length[volumesMatrix[[i,j]]]!=0,
      CurrentVol = volumesMatrix[[i,j]];
      VolSup = volumesMatrix[[i-1,j]];
      VolInf = volumesMatrix[[i+1,j]];
      VolR = volumesMatrix[[i,j+1]];
      VolI = volumesMatrix[[i,j-1]];
      
      If[Length[volumesMatrix[[i,j]][[2]]]==6,
      If[ Length[VolSup]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,3]]];
      ];
       If[ Length[VolR]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,4]]];
      ];
        If[ Length[VolInf]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,5]]];
      ];
      If[ Length[VolI]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,6]]];
      ];,
     
      If[ (Length[VolSup]==0 && Length[VolR]==0) || (Length[VolSup]==0 && Length[VolI]==0),
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,3]]];
      ];
        If[ Length[VolInf]==0 && Length[VolR]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,4]]];
      ];
         If[ Length[VolInf]==0 && Length[VolI]==0,
      AppendTo[indexSurfaces, volumesMatrix[[i,j]][[2,5]]];
      ];
      
      
      ];
   ];
];

];
Return[indexSurfaces];

];



WellsDefinition[tope_, base_, wellsvector_]:=Block[{ftope=tope, fbase=base, fwellsvector=wellsvector, nfilas, ncols, matrixwelss,count,i,j,pair, pairs},
nfilas = Dimensions[ftope][[1]];
ncols = Dimensions[ftope][[2]];
matrixwelss=ConstantArray[0,{nfilas, ncols}];
count=1;
pairs={};
For[i=1, i<= ncols, i++,
For[j=1, j<= ncols, j++,
	If[ftope[[i,j]]!=0,
		matrixwelss[[i,j]]=1;
		If[ftope[[i,j]]==fwellsvector[[count,1]],
			matrixwelss[[i,j]]=fwellsvector[[count,2]];
			If[count<Length[fwellsvector],
			count++;
			];
			pair={i,j};
			AppendTo[pairs,pair];
		];
	];
];
];
Return[{matrixwelss, pairs}];
];


WellsPointsCreation[topeinformation_,baseinformation_, wellsinformation_, pointsIndexTop_, pointsCoordinatesTop_,pointsIndexBas_, pointsCoordinatesBase_]:=
Block[{ftopeinformation=topeinformation, fbaseinformation=baseinformation, 
fwellsinformation= wellsinformation},
    shift=100000;
	nwells = Length[wellsinformation[[2]]];
	matrixsurfacetope=ftopeinformation[[1]];
	matrixLinestope=ftopeinformation[[2]];
	matrixsurfacebase=fbaseinformation[[1]];
	matrixLinesbase=fbaseinformation[[2]];
	points={};
	v=0;
	IndexSurfaceTop={};
	IndexSurfaceBase={};
	For[iwell=1, iwell<=nwells, iwell++,
		fila = fwellsinformation[[2]][[iwell,1]];
		cols = fwellsinformation[[2]][[iwell,2]];
		
		
		matrixsurfacetope=ftopeinformation[[1]];
		matrixLinestope=ftopeinformation[[2]];
		pointsIndex=pointsIndexTop;
		pointsCoordinates=pointsCoordinatesTop;
		For[count=1, count<=2, count++,
		
		If[count==2,
			matrixsurfacetope=matrixsurfacebase;
			matrixLinestope=matrixLinesbase;
			pointsIndex=pointsIndexBas;
		    pointsCoordinates=pointsCoordinatesBase;
		];
		
		line1topeIndex = Abs[matrixsurfacetope[[fila, cols]][[2,1]]];
		line2topeIndex = Abs[matrixsurfacetope[[fila, cols]][[2,2]]];
		points1Index = Position[matrixLinestope[[All,1]], line1topeIndex][[1,1]];
		(*primer ponto*)
		points1IndexPosition= matrixLinestope[[points1Index,2]][[1]];
		cordinatesPosition = Position[pointsIndex,points1IndexPosition][[1,1]];
		cordinates1Real1 = pointsCoordinates[[cordinatesPosition]];
		(*Segundo Ponto*)
		points1IndexPosition= matrixLinestope[[points1Index,2]][[2]];
		cordinatesPosition = Position[pointsIndex,points1IndexPosition][[1,1]];
		cordinates1Real2 = pointsCoordinates[[cordinatesPosition]];
		(*Tercer Ponto*)
		points1Index = Position[matrixLinestope[[All,1]], line2topeIndex][[1,1]];
		points1IndexPosition= matrixLinestope[[points1Index,2]][[2]];
		cordinatesPosition = Position[pointsIndex,points1IndexPosition][[1,1]];
		cordinates1Real3 = pointsCoordinates[[cordinatesPosition]];
		centralPoint ={ (cordinates1Real1[[1]]+ cordinates1Real2[[1]])/2 ,(cordinates1Real2[[2]] + cordinates1Real3[[2]])/2 , (cordinates1Real2[[3]] + cordinates1Real3[[3]])/2  };
		r=10;
		AppendTo[points,centralPoint];
		ponint1 ={centralPoint[[1]]+ r ,centralPoint[[2]],centralPoint[[3]]};
		AppendTo[points,ponint1];
		ponint2 ={centralPoint[[1]] ,centralPoint[[2]]+ r,centralPoint[[3]]};
		AppendTo[points,ponint2];
		ponint3 = {centralPoint[[1]]-r ,centralPoint[[2]],centralPoint[[3]]};
		AppendTo[points,ponint3];
		ponint4 = {centralPoint[[1]] ,centralPoint[[2]]-r,centralPoint[[3]]};
		AppendTo[points,ponint4];
		];
	];
	index = Table[indexi,{indexi, shift, shift + Length[points],1}];
	indexTope={};
	indexBase={};
	point={};
	v=0;
	For[i=1,i<= Length[index],i++,
	AppendTo[point,index[[i]]];
	  If[Mod[i,5]==0,
	    If[v==0,
	    AppendTo[indexTope,point];
	    v=1;,
	     AppendTo[indexBase,point];
	    v=0;
	    ];
	    point={};
	  ];
	
	];
	
	Return[{points,index, indexTope, indexBase}];
	

];


WellLines[WellsInformation_]:=Block[{i},
pointsTop = WellsInformation[[3]];
pointsBase = WellsInformation[[4]];
nwells = Dimensions[pointsTop][[1]];
LinesTOP={};
LinesBASE={};
IndexLineTop={};
IndexLineBase={};
shift=100000;
	For[i=1, i<= nwells, i++,
	point ={pointsTop[[i,2]],pointsTop[[i,1]],pointsTop[[i,3]]};  
	AppendTo[LinesTOP,point];
	point ={pointsTop[[i,3]],pointsTop[[i,1]],pointsTop[[i,4]]};
	AppendTo[LinesTOP,point];  
	point ={pointsTop[[i,4]],pointsTop[[i,1]],pointsTop[[i,5]]}; 
	AppendTo[LinesTOP,point];
	point ={pointsTop[[i,5]],pointsTop[[i,1]],pointsTop[[i,2]]}; 
	AppendTo[LinesTOP,point];

	];
	For[i=1, i<= nwells, i++,
	
	point ={pointsBase[[i,2]],pointsBase[[i,1]],pointsBase[[i,3]]};  
	AppendTo[LinesBASE,point];
	point ={pointsBase[[i,3]],pointsBase[[i,1]],pointsBase[[i,4]]};
	AppendTo[LinesBASE,point];  
	point ={pointsBase[[i,4]],pointsBase[[i,1]],pointsBase[[i,5]]}; 
	AppendTo[LinesBASE,point];
	point ={pointsBase[[i,5]],pointsBase[[i,1]],pointsBase[[i,2]]}; 
	AppendTo[LinesBASE,point];
	];
	IndexLineTop = Table[j,{j,shift, shift+nwells*4-1, 1}];
	IndexLineBase = Table[j,{j,shift+nwells*4 , shift+(2 nwells)*4 -1, 1}];
	PintsSeparateT=Table[{IndexLineTop[[i]],IndexLineTop[[i+1]],IndexLineTop[[i+2]],IndexLineTop[[i+3]] },{i,1,Length[IndexLineTop],4}];
	PintsSeparateB=Table[{IndexLineBase[[i]],IndexLineBase[[i+1]],IndexLineBase[[i+2]],IndexLineBase[[i+3]] },{i,1,Length[IndexLineBase],4}];
    PintsSeparateLateral=Table[{i,i+1,i+2,i+3 },{i,shift+(2 nwells)*4 ,shift+(3 nwells)*4 -1,4}];

    Return[{LinesTOP, LinesBASE, IndexLineTop, IndexLineBase,PintsSeparateT, PintsSeparateB,PintsSeparateLateral}];	
];



CreateCircGmshFormat[points_, lines_]:=Block[{i},
nCirc = Length[points];
circles={};
For[i=1, i<= nCirc, i++,
circ = StringJoin["Circle( ", ToString[lines[[i]]], ")= {" , ToString[points[[i,1]]], ", ",ToString[points[[i,2]]], ", ", ToString[points[[i,3]]],"};" ];
AppendTo[circles,circ];
];
Return[circles];
];



WellsLateralLines[wellsinformation_]:=Block[{i},
   TopePoints = wellsinformation[[1]];
   BasePoints = wellsinformation[[2]];
   indexLateralLines = wellsinformation[[7]]//Flatten;
   topemod ={};
   basemod={};
   Print[Length[TopePoints]];
   For[i=1,i<= Length[TopePoints],i++,
   pointop = TopePoints[[i,1]];
      AppendTo[topemod,pointop];

       If[i==Length[TopePoints],
      pointop = TopePoints[[i,3]];
      AppendTo[topemod,pointop];
     ];
      
      pointop = BasePoints[[i,1]];
      AppendTo[basemod,pointop];
  
     If[i==Length[TopePoints],
      pointop = BasePoints[[i,3]];
      AppendTo[basemod,pointop];
     ];
   
   ];
   
   
   Lines={};
   For[i=1, i<= Length[indexLateralLines], i++,
       line = {indexLateralLines[[i]],{topemod[[i]],basemod[[i]]}};
        AppendTo[Lines,line];
   ];
   Return[Lines];
];



TakePerfSurface[matrixsurface_, vectorCord_]:=Block[{},
    nsurfaces= Length[vectorCord];
    Index={};
    For[i=1, i<= nsurfaces,i++,
    ind = matrixsurface[[vectorCord[[i,1]],vectorCord[[i,2]] ]][[1]];
    AppendTo[Index,ind];
    ];
    Return[Index];
];
