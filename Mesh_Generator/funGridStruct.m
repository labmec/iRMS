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





FilterInternalPoints[Region_,data_]:=Block[{fRegion=Region, fdata=data,
npoints,newlist,fhx,fhy,boundaries,posx,posj,index,linesindex,i,j,CordinatesMatrixtope,CordinatesMatrixbase,
pointsBaseIndex,isInternal,pointsTopeIndex,comparationPoint,cordinatesPoints,
indexquadtope,indexquadbase,point,valInterpoltope,valInterpolbase},
fhx= Dimensions[fdata[[1]]][[1]];
fhy= Dimensions[fdata[[1]]][[2]];

CordinatesMatrixbase = fdata[[1]];
pointsBaseIndex=fdata[[2]];
 newlist = {}; 
 cordinatesPoints={};
 indexquadtope={};
 indexquadbase={};
 For[j=1, j<= fhy, j++,
    For[i=1, i<= fhx, i++,
      comparationPoint=CordinatesMatrixbase[[j,i]][[1;;2]];
      isInternal = RegionMember[fRegion,comparationPoint];
      If[!isInternal,
         pointsBaseIndex[[j,i]]= 0;
    ];
 ];
 ];

 Return[{CordinatesMatrixbase,pointsBaseIndex} ];
];



ConvertIndexesToPoints[vectorIndex_, MatrixCordinates_]:=Block[{i, index, nindex, icord,
jcord, cordPoint, Cordinates,nx},
     nx=Dimensions[MatrixCordinates][[1]];
     Cordinates ={};
     nindex = Length[vectorIndex];
     For[i=1, i<= nindex, i++,
        index = vectorIndex[[i]];
        icord = IntegerPart[(index-1)/nx]+1;
        jcord = index - (icord-1)*nx;
        cordPoint = MatrixCordinates[[icord, jcord]];
        cordPoint[[3]]=0.0;
        AppendTo[Cordinates, cordPoint];
     ];
     
     Return[Cordinates];

];



ConvertIndexesListToPointsFaultsList[vectorIndex_, MatrixCordinates_,inipoint_]:=Block[{fvectorIndex=vectorIndex, 
fMatrixCordinates=MatrixCordinates,faultsPoints, nFaults, i, finipoint=inipoint, faultsIndexPoints, indexPooint},
	nFaults = Length[fvectorIndex];
	faultsPoints={};
	faultsIndexPoints={};
	For[i=1, i<= nFaults, i++,
	   AppendTo[faultsPoints,ConvertIndexesToPoints[fvectorIndex[[i]],fMatrixCordinates]];
	   indexPooint = Table[finipoint+w,{w,0,Length[fvectorIndex[[i]]]-1,1}];
	   AppendTo[faultsIndexPoints,indexPooint];
	   finipoint = finipoint+Length[fvectorIndex[[i]]];
	];
	Return[{faultsPoints,faultsIndexPoints }];
];




FaultsLinesConnect [FaultsIndex_, iniline_]:=Block[{fFaultsIndex=FaultsIndex, finiline=iniline, index, AllLines,linesIndexes, nFaults, i ,  npoints, lines},
     nFaults = Length[fFaultsIndex];
     AllLines ={};
     linesIndexes={};
     For[i=1, i<= nFaults, i++,
       npoints = Length[fFaultsIndex[[i]]];
       lines = Table[{fFaultsIndex[[i,w]],fFaultsIndex[[i,w+1]]},{w,1,npoints-1}];
       index = Table[finiline + q, {q, 0, Length[lines]-1,1 }];
       finiline += Length[lines]+1;
       AppendTo[AllLines, lines];
       AppendTo[linesIndexes, index];
     ];
     Return[{AllLines, linesIndexes}];

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






GenGmshLine[index1_,index2_,indexl_]:=Block[{findex1=index1,findex2=index2, findexl=indexl,line},
 line = StringJoin["Line(",ToString[findexl],")= {",ToString[findex1],", ",ToString[findex2],"};"];
 Return[line];
];




WellsDefinitionPoints[vecCordinates_, vecRaios_, vecFactors_, iniPoint_]:=Block[{i},
    nwells = Length[vecCordinates];
    val = iniPoint;
    AllPointsCirc={}; 
    AllConnectPointsIndex={};   
    AllConnectPointsIndexExtern={};   
    PointsWithOutDef={};
    PointsWithOutDefExtern={};
    iniCircExtern= iniPoint + (nwells*5);
    AllPointsCircExter={};
    
    For[i=1, i<= nwells, i++,
     cord = vecCordinates[[i]];
     rad = vecRaios[[i]];
     fac=vecFactors[[i]];
     pointsCirc = CreatePointsCirc[cord, rad,0];
     pointsCircExt = CreatePointsCirc[cord, rad*fac,0];
     pointsCircExt = pointsCircExt[[2;;5]];
     
     AppendTo[AllPointsCirc,pointsCirc];
     AppendTo[AllPointsCircExter,pointsCircExt];
     
     indexPointsCirc = Table[j,{j,val,val+4,1}];
     indexPointsCircExtern = Table[j,{j,iniCircExtern,iniCircExtern+3,1}];
     
     
     val = val+5;
     iniCircExtern = iniCircExtern+4;
     AppendTo[AllConnectPointsIndex, indexPointsCirc];
     AppendTo[AllConnectPointsIndexExtern, indexPointsCircExtern];
     
     For[j=1, j<= 5, j++,
     AppendTo[PointsWithOutDef,pointsCirc[[j]]];
     ];
      For[j=1, j<= 4, j++,
     AppendTo[PointsWithOutDefExtern,pointsCircExt[[j]]];
     ];
    ];
    
    
    
    
	Return[{AllPointsCirc, AllConnectPointsIndex,PointsWithOutDef, AllPointsCircExter, AllConnectPointsIndexExtern,PointsWithOutDefExtern}];
];


CreatePointsCirc[Cordinates_, Raio_, val_]:=Block[{pointsVec},
    point1 = {Cordinates[[1]]+Raio,Cordinates[[2]],val};
    point2 = {Cordinates[[1]],Cordinates[[2]]+Raio,val};
    point3 = {Cordinates[[1]]-Raio,Cordinates[[2]],val};
    point4 = {Cordinates[[1]],Cordinates[[2]]-Raio,val};
    pointsVec={Cordinates,point1,point2,point3,point4 };
    Return[pointsVec];
];




WellLines[CircInformationIndex_, CircInformationIndexExtern_, IniLineIndex_]:=Block[{i},
nwells = Length[CircInformationIndex];


LinesInteral={};
LinesExtern={};
	Print[nwells];
	For[i=1, i<= nwells, i++,
	(*INTERNAL CIRC LINES*)
	point ={CircInformationIndex[[i,2]],CircInformationIndex[[i,1]],CircInformationIndex[[i,3]]};  
	AppendTo[LinesInteral,point];
	point ={CircInformationIndex[[i,3]],CircInformationIndex[[i,1]],CircInformationIndex[[i,4]]};
	AppendTo[LinesInteral,point];  
	point ={CircInformationIndex[[i,4]],CircInformationIndex[[i,1]],CircInformationIndex[[i,5]]}; 
	AppendTo[LinesInteral,point];
	point ={CircInformationIndex[[i,5]],CircInformationIndex[[i,1]],CircInformationIndex[[i,2]]}; 
	AppendTo[LinesInteral,point];
	
	(*EXTERNAL CIRC LINES*)
	point ={CircInformationIndexExtern[[i,1]],CircInformationIndex[[i,1]],CircInformationIndexExtern[[i,2]]};  
	AppendTo[LinesExtern,point];
	point ={CircInformationIndexExtern[[i,2]],CircInformationIndex[[i,1]],CircInformationIndexExtern[[i,3]]};
	AppendTo[LinesExtern,point];  
	point ={CircInformationIndexExtern[[i,3]],CircInformationIndex[[i,1]],CircInformationIndexExtern[[i,4]]}; 
	AppendTo[LinesExtern,point];
	point ={CircInformationIndexExtern[[i,4]],CircInformationIndex[[i,1]],CircInformationIndexExtern[[i,1]]}; 
	AppendTo[LinesExtern,point];

	];
	
	indexLinesCirc = Table[{i,i+1,i+2,i+3}, {i, IniLineIndex,IniLineIndex + nwells*(4)-1 , 4}];
	iniExtCirc = indexLinesCirc[[-1,-1]]+1;
	indexLinesCircExtern = Table[{i,i+1,i+2,i+3}, {i, iniExtCirc,  iniExtCirc + nwells*(4)-1, 4}];

    Return[{LinesInteral, indexLinesCirc, LinesExtern, indexLinesCircExtern}];	
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




GenGmshSurfaceWellList2[AllindexlinesInternals_, AllindexlinesExternals_]:=Block[{fAllindexlinesInternals=AllindexlinesInternals,fAllindexlinesExternals=AllindexlinesExternals,i,nwells,surface, iniLineLoopExternal, iniLineLoopInternal,
indexSurface, AllSurfaces, LineLoopsExt, indexInternal, indexExternal,LineLoopsInt, SurfacesIndex},
    nwells = Length[fAllindexlinesInternals];
    iniLineLoopInternal=1000;
    iniLineLoopExternal = iniLineLoopInternal+nwells+1;
    indexSurface=1;
    AllSurfaces={};
    LineLoopsExt={};
    LineLoopsInt={};
    SurfacesIndex={};
    For[i=1, i<= nwells, i++,
        indexInternal = fAllindexlinesInternals[[i]];
        indexExternal = fAllindexlinesExternals[[i]];
        surface = GenGmshSurfaceWell[indexInternal, indexExternal,iniLineLoopInternal, iniLineLoopExternal,indexSurface ];
        AppendTo[AllSurfaces, surface];
        AppendTo[LineLoopsExt,iniLineLoopExternal];
         AppendTo[LineLoopsInt,iniLineLoopInternal];
          AppendTo[SurfacesIndex,indexSurface];
        iniLineLoopInternal++;
        iniLineLoopExternal++;
        indexSurface++;
        
    ];
 
    Return[{AllSurfaces, LineLoopsExt, LineLoopsInt, SurfacesIndex}];   
];




GenGmshSurfaceWell[indexlinesInternals_, indexlinesExternals_, indInter_, indExtern_,indexSurface_]:=Block[{findexlinesInternals=indexlinesInternals,
findexlinesExternals=indexlinesExternals, findInter=indInter, findExtern=indExtern, findexSurface=indexSurface, lineloop1, lineloop2,surface},
	
  lineloop1 = StringJoin["Curve Loop(",ToString[findInter],")= {",ToString[findexlinesInternals[[1]]],", ",
  ToString[findexlinesInternals[[2]]],",", ToString[findexlinesInternals[[3]]],",",ToString[findexlinesInternals[[4]]],"};" ];
 
  lineloop2 = StringJoin["Curve Loop(",ToString[findExtern],")= {",ToString[findexlinesExternals[[1]]],", ",
  ToString[findexlinesExternals[[2]]],",", ToString[findexlinesExternals[[3]]],",",ToString[findexlinesExternals[[4]]],"};" ];
 
  surface= StringJoin["Plane Surface(",ToString[findexSurface],")= {",ToString[findInter], ", ",  ToString[findExtern],"};"];
   
   Return[{lineloop1,lineloop2,surface}];
 
];


ReservoirSurface[linesReservoir_, linesExternWells_, loopReservoirIndex_, surfIndex_]:=Block[{flinesReservoir= linesReservoir, 
flinesExternWells=linesExternWells, floopReservoirIndex=loopReservoirIndex, fsurfIndex=surfIndex},

	lineloop = StringJoin["Line Loop(", ToString[floopReservoirIndex], ")=", ToString[linesReservoir],";" ];
	resSurface = ToString[Join[{floopReservoirIndex}, flinesExternWells]];
	resSurfaceGmshFormat = StringJoin["Plane Surface(", ToString[fsurfIndex],")= ", resSurface,";"];
	Return[{lineloop,resSurfaceGmshFormat}];

];



TakeIndexProdAndInj[AllWellLines_, AllWellsTypes_]:=Block[{fAllWellLines=AllWellLines, fAllWellsTypes=AllWellsTypes,
productors, injectors,nwellss, i},
     injectors={};
     productors={};
     nwellss = Length[fAllWellsTypes];
     For[i=1, i<=nwellss,i++,
        If[fAllWellsTypes[[i]]==1,
        AppendTo[productors,fAllWellLines[[i]]];,
        
        AppendTo[injectors,fAllWellLines[[i]]];
        ];
     
       ];
       productors = productors //Flatten;
       injectors= injectors //Flatten;
       Return[{productors, injectors}];
];




TransfiniteLines[LinesIndex_, transval_]:=Block[{fLinesIndex=LinesIndex, ftransval=transval},
   trans = StringJoin["Transfinite Curve",ToString[fLinesIndex]," = ", ToString[ftransval],";"];
   Return[trans];
];


RecombineSurface[LinesIndex_]:=Block[{fLinesIndex=LinesIndex,trans},
   trans = StringJoin["Recombine Surface",ToString[fLinesIndex],";"];
   Return[trans];
];


PhysicalTagSurface[name_,listSurfaces_ ]:=Block[{fname=name, flistSurfaces=listSurfaces,pysurf},
      pysurf=StringJoin["Physical Surface(\"",fname,"\"",")=" ,ToString[flistSurfaces],";"];
      Return[pysurf];
];




PhysicalTagCurve[name_,listLines_ ]:=Block[{fname=name, flistLines=listLines,pylines},
      pylines=StringJoin["Physical Curve(","\"",fname,"\"",")=", ToString[flistLines],";"];
      Return[pylines];
];



JoinFaults[faultsIndexPointsVec_, faultsPointsVec_,joinTypeVec_, pointsResVec_, indexPointRes_]:=Block[{ffaultsIndexPointsVec=faultsIndexPointsVec,
ffaultsPointsVec=faultsPointsVec, fjoinTypeVec=joinTypeVec, fpointsResVec=pointsResVec, findexPointRes=indexPointRes ,i,j,
nfaults,quadContorn,quadFaults, fakeVectorRes, fakeVectorFaults, condType,PointAn,currentIndex, closePoint, fakeIndex, index    },
     nfaults = Length[ffaultsIndexPointsVec];
     quadContorn = Nearest[fpointsResVec];
     quadFaults = Nearest[ffaultsPointsVec];
     fakeVectorRes = findexPointRes//Flatten;
     fakeVectorFaults = ffaultsIndexPointsVec//Flatten;
     test={};
     pointsTest= {};
     For[i=1, i<= nfaults, i++,
      
         For[j=1, j<= 2, j++,
             condType = fjoinTypeVec[[i,j]];
               If[j==1,
                If[Length[ffaultsIndexPointsVec[[i]]]>2 ,
                   If[condType!= 1 && condType !=2,
                      AppendTo[test, ffaultsIndexPointsVec[[i]]];
                      AppendTo[pointsTest, ffaultsPointsVec[[i]]];
                   ];,
                  AppendTo[test, ffaultsIndexPointsVec[[i]][[2;;-2]]];
                  AppendTo[pointsTest, ffaultsPointsVec[[i]][[2;;-2]]];
               ];
              ];
         
             If[j==1,
               PointAn = ffaultsPointsVec[[i,1]];
               currentIndex = ffaultsIndexPointsVec[[i,1]];,
               PointAn = ffaultsPointsVec[[i,-1]];
               currentIndex = ffaultsIndexPointsVec[[i,-1]];
               
             ];
             If[condType==1, 
                closePoint = quadContorn[PointAn];
                fakeIndex = Position[closePoint, fpointsResVec];
                index = fakeVectorRes[[fakeIndex]];
                
             ];
             If[condType==2, 
                closePoint = quadFaults[PointAn];
                fakeIndex = Position[closePoint, ffaultsPointsVec];
                index = fakeVectorFaults[[fakeIndex]];
             ];
             
              If[j==1,
               
               ffaultsIndexPointsVec[[i,1]]=index;,
               ffaultsIndexPointsVec[[i,-1]]=index;
               
             ];
             
         ];
     ];
];



ModifyBoundaryPoints[fakePoints_, realPoints_]:=Block[{ffakePoints=fakePoints, frealPoints=realPoints,point,modPoint,npoints,
function, i,position,newpoints},
     npoints = Length[ffakePoints];
     function = Nearest[frealPoints];
     newpoints={};
     For[i=1, i<=npoints, i++,
         
         point = ffakePoints[[i]];
         modPoint = Nearest[frealPoints,point];
         position = Position[frealPoints,modPoint];
         frealPoints=Delete[frealPoints,position];
         modPoint[[1,3]]=0.0;
         If[Length[modPoint]==2,
            modPoint = (modPoint[[1]]+modPoint[[2]])/2;
         ];
         AppendTo[newpoints,modPoint[[1]]];
 
     ];
     Return[newpoints];
];



CheckExistPoints[IndexFaultsPoints_, PointsFaultsCordinates_, ResIndex_]:=Block[{fIndexFaultsPoints=IndexFaultsPoints, fPointsFaultsCordinates=PointsFaultsCordinates,
,fResIndex=ResIndex, fakeIndexPoints,nPoints, indexes,pointsCoordin,indexAn, findDuplicates, existInReservoir, i, nduplicates,
 deleteIndex, fakePos },
    fIndexFaultsPoints = fIndexFaultsPoints //Flatten;
    nPoints = Length[fIndexFaultsPoints];
    indexes = {};
    pointsCoordin={};
    For[i =1, i<=nPoints, i++,
         indexAn = fIndexFaultsPoints[[i]];
         findDuplicates = Position[fIndexFaultsPoints, indexAn] // Flatten;
         existInReservoir = Position[fResIndex, indexAn]//Flatten;
         If[Length[existInReservoir]==0,
             nduplicates = Length[findDuplicates];
             
             If[nduplicates>1 ,
                  
                   For[j=2, j<=nduplicates, j++,
                   deleteVal = findDuplicates[[j]];
                   fIndexFaultsPoints=Delete[fIndexFaultsPoints,deleteVal];
                   fPointsFaultsCordinates=Delete[fPointsFaultsCordinates,deleteVal];
                   nPoints=Length[fIndexFaultsPoints];
                  
                 ];
             ];,
              
              fIndexFaultsPoints=Delete[fIndexFaultsPoints, i];
              fPointsFaultsCordinates=Delete[fPointsFaultsCordinates,i];
              nPoints=Length[fIndexFaultsPoints];
              i--;
         ];
    ];
     Return[{fIndexFaultsPoints, fPointsFaultsCordinates}];
];



FindClosePoint[point_,listPoints_, listIndex_]:=Block[{fpoint, flistPoints=listPoints, flistIndex=listIndex},
   nPoints = Length[listPoints];
   men = flistIndex[[1]];
   distMen = 1.0*10^15;
   For[i=1, i<= nPoints, i++,
       Distance = Norm[fpoint - flistPoints[[i]]];
       If[Distance<distMen,
           distMen = Distance;
           men=flistIndex[[i]];
       ];
   ];
   Return[men];
];



CalcContornPoints[data_]:=Block[{fdata=data,nx,ny,newtable,i,
valr,vall,vals,vali,sum,j,newmatrix,testtable,nneigh,ListBoundaryPoints,valsr,valsl,valir
,valil,auxTable},

testtable = ExpandGrid[fdata];
auxTable = ExpandGrid[fdata];
nx = Dimensions[testtable][[1]];
ny = Dimensions[testtable][[2]];

ListBoundaryPoints={};

For[i=2, i< nx, i++,
	For[j=2,j< ny, j++,
	If[auxTable[[i,j]]== 0, 
		Continue[];
	];
	nneigh=0;
	vals = auxTable[[i-1,j]];
	vali = auxTable[[i+1,j]];
	valr = auxTable[[i,j+1]];
	vall = auxTable[[i,j-1]];
	valsr = auxTable[[i-1,j+1]];
	valsl = auxTable[[i-1,j-1]];
	valir = auxTable[[i+1,j+1]];
	valil = auxTable[[i+1,j-1]];
	If[vals!=0,
	nneigh++;
	];
	If[vali!=0,
	nneigh++;
	];
	If[valr!=0,
	nneigh++;
	];
	If[vall!=0,
	nneigh++;
	];
	If[valsr!=0,
	nneigh++;
	];
	If[valsl!=0,
	nneigh++;
	];
	If[valir!=0,
	nneigh++;
	];
	If[valil!=0,
	nneigh++;
	];
	If[nneigh==8  ,
	  testtable[[i,j]]=0;
	];
	If[(nneigh==7 && valsr==0) ||(nneigh==7 && valsl==0) ||(nneigh==7 && valil==0)
	||(nneigh==7 && valir==0)  ,
	  testtable[[i,j]]=0;
	];
	
	If[((nneigh==6 && valsr==0) &&(nneigh==6 && valir==0)) ||((nneigh==6 && valsl==0)
	&&(nneigh==6 && valil==0))  ,
	  testtable[[i,j]]=0;
	];
	
	];
];

Return[testtable];
];




ExpandGrid[data_]:=Block[{fdata=data,nx,ny,newgrid,i,j,globali,globalj},
nx=Dimensions[fdata][[1]];
ny=Dimensions[fdata][[1]];
newgrid=ConstantArray[0,{nx+2,ny+2}];
For[i=1,i<=nx,i++,
For[j=1,j<=ny,j++,
   globali=i+1;
   globalj=j+1;
   newgrid[[globali,globalj]]=fdata[[i,j]];
   
];
];
Return[newgrid];
];




CalcBoundariesPointsIndex[data_]:=Block[{fdata=data,
dim, nx, ny, valores, i, j, firstx,firsty,condition,cadena,currentPosX, currentPosY,
lastPosX, lastPosY,Sval, Ival, Lval,Rval, SRval, SLval, IRval,ILval,pase, bandera,newpoint,index,pointindex,currentPointIndex,
auxPosX,auxPosY},

	dim= Dimensions[fdata];
	nx= dim[[1]];
	ny= dim[[2]];
	{firstx,firsty}=CalcFirstSecond[fdata][[1]];
	{auxPosX,auxPosY}=CalcFirstSecond[fdata][[2]];
	
	condition=0;
	cadena ={fdata[[firstx,firsty]],fdata[[auxPosX,auxPosY]]};

	pase=0;
	While[condition ==0,
	    pase++;
	    bandera=0;
		currentPosX = auxPosX;
		currentPosY = auxPosY;
		Sval = fdata[[currentPosX+1,currentPosY]];
		Ival = fdata[[currentPosX-1,currentPosY]];
		Lval = fdata[[currentPosX,currentPosY-1]];
		Rval = fdata[[currentPosX,currentPosY+1]];
		SRval= fdata[[currentPosX+1,currentPosY+1]];
		SLval= fdata[[currentPosX-1,currentPosY-1]];
		IRval= fdata[[currentPosX-1,currentPosY+1]];
		ILval= fdata[[currentPosX+ 1,currentPosY-1]];
		newpoint={firstx,firsty};
		
		If[Sval!=0 && (IsChainMember[cadena,Sval]==False) &&bandera==0,
				auxPosX++;
				newpoint={currentPosX+1,currentPosY};
				bandera=1;
				
		];
		
		If[Ival!=0 && (IsChainMember[cadena,Ival]==False) &&bandera==0,
			 auxPosX--;
			  newpoint={currentPosX-1,currentPosY};
			  bandera=1;

		];
		
		
		If[Lval!=0 && (IsChainMember[cadena,Lval]==False)&&bandera==0,
			auxPosY--;
			   newpoint={currentPosX,currentPosY-1};
			   bandera=1;
			    
		];
	
	
		If[Rval!=0 && (IsChainMember[cadena,Rval]==False)&&bandera==0,
		   auxPosY++;
		   newpoint={currentPosX,currentPosY+1};
		   bandera=1;
	    
			   
		];
		
		If[SRval!=0 && (IsChainMember[cadena,SRval]==False)&&bandera==0,
				   newpoint={currentPosX+1,currentPosY+1};
		      auxPosY++;
		      auxPosX++;
			  bandera=1;
			
		];
		If[SLval!=0 && (IsChainMember[cadena,SLval]==False)&&bandera==0,
		auxPosY--;
	    auxPosX--;
		newpoint={currentPosX-1,currentPosY-1};
		bandera=1;
		];
		
		If[IRval!=0 &&(IsChainMember[cadena,IRval]==False)&&bandera==0,
		auxPosX--;
		auxPosY++;
	    
	    newpoint={currentPosX-1,currentPosY+1};
		bandera=1;
		];
		If[ILval!=0 && (IsChainMember[cadena,ILval]==False)&&bandera==0,
        auxPosX++;
		auxPosY--;
				  newpoint={currentPosX+1,currentPosY-1};
				  bandera=1;
		];
		If[newpoint!={firstx,firsty},
		AppendTo[cadena,fdata[[newpoint[[1]],newpoint[[2]]]]];
		];
		If[newpoint=={firstx,firsty},
		AppendTo[cadena,fdata[[newpoint[[1]],newpoint[[2]]]]];
		condition=1;
		];
		
		If[pase>10000,
		condition=1;
		];
			
	];
	
	Return[cadena];
];


IsChainMember[data_,point_]:=Block[{fdata=data,fpoint=point,i,n},
	n = Length[fdata];
	For[i=1,i<= n, i++,
		If[fdata[[i]]==fpoint ,
		Return[True];
		Break[];
		];
	];

	Return[False];
	
];





CalcFirstSecond[data_]:=Block[{fdata=data,fval=val,dim,nx,ny,condition, i, j,ncandidates,FirstSecond},
	dim= Dimensions[fdata];
	nx= dim[[1]];
	ny= dim[[2]];
	condition =0;
	ncandidates =0;
	FirstSecond = {};
	For[i=1,i<= nx, i++,
		For[j=1,j<= ny, j++,
			If[fdata[[i,j]]!=0,
				AppendTo[FirstSecond,{i,j}];
				ncandidates++;
				If[ncandidates==2,
				Return[FirstSecond];
				condition=1;
				Break[];
				];
				];
		     ];
		If[condition==1,
		Break[];
		];
	];
];




