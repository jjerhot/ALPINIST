(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["BDplotSettings`"]
baseStyle::usage="Basic font setup"
baseStyleBoldIt::usage="Basic font setup with bold and italic option"
baseStyleLeg90CL::usage="Basic font setup for 90CL plot legend"
LogTick::usage="List of ticks for x and y axis"
LogTickLegend::usage="List of ticks for z axis"
reportColorRange::usage="Returns the colors and range for given plot. Use as {plot,colors,range} = reportColorRange['original contour plot']"
trimPoint::usage="display number n with given number of sig.digits,trim trailing decimal point"
colorLegend::usage="Makes log z axis legend according to the input color and range"
colorLegendLinear::usage="Makes z axis legend according to the input color and range"
at::usage="convenience function to position objects in Graphics"
display::usage="Same as Graphics,but with fixed PlotRange"
cropVector::usage="Crop for vector graphics"
CleanContourPlot::usage="Remove reduntant data in contour plots (effectively decreases size for contour plots)"
expCol::usage="Colors for exp contours (0 - shading, 1 - style)"
expLeg::usage="Colors for exp legend"
expCont::usage="90% CL contours"

baseLegend90CL::usage="Basic setup for 90CL legend"
dashedLegend90CL::usage="90CL legend for dashed contours (2\[Gamma]s)"
regionsLegendDS90CL::usage="Setup for region plots legend"
gridLines3mesons::usage="Mask mass region around \[Pi],\[Eta],\[Eta]'"
gridLines2mesons::usage="Mask mass region around \[Eta],\[Eta]'"
gridLines2mesonsThin::usage="Mask mass region around \[Eta],\[Eta]'"
yLabelxMass::usage="Setup for axis labels. ylabel as input, xlabel fixed as ALP mass."
yLabelNoDimxMassIndex::usage="Setup for axis labels. ylabel as input, xlabel fixed as ALP mass."
plotSettings90CL::usage="Basic setup for 90CL plots"
plotReducedSettings90CL::usage="Basic setup for 90CL plots for ALP mass > 500 MeV"

ALPPlotStd::usage="Generate standard ALP exclusion plot for given dataset and experiment"
ALPPlotStdMask::usage="Generate standard ALP exclusion plot for given dataset and experiment, masking regions around pi0/eta/eta'"
ALPPlotReduced::usage="Generate standard ALP exclusion plot for given dataset and experiment, zoomed to higher masses, masking regions around pi0/eta/eta'"
ALPPlotDashed::usage="Generate standard ALP exclusion plot with dashed contours for given dataset and experiment, zoomed to higher masses, masking regions around pi0/eta/eta'"
ALPPlotKOTO::usage="Generate standard ALP exclusion plot for given dataset and experiment together with KOTO"
ALPPlotFFMask::usage="Generate fermion-coupled ALP exclusion plot for given dataset and experiment"
DSPlot::usage="Generate standard DS exclusion plot for given dataset and experiment"
DPPlot::usage="Generate standard DP exclusion plot for given dataset and experiment"

Begin["`Private`"]

baseStyle={FontSize-> 17,FontFamily->"CMU Serif",SingleLetterItalics-> False,FontColor-> Black};
baseStyleBoldIt={FontSize-> 17,FontFamily->"CMU Serif",SingleLetterItalics-> False,FontColor-> Black,Bold,Italic};
baseStyleLeg90CL={FontSize-> 13,FontFamily->"CMU Serif",SingleLetterItalics-> False,FontColor-> Black};
LogTick[min_,max_]:=Flatten[Table[If[j==1,{j*10^i,Superscript[10,IntegerPart[i]],{.02,0}},{j*10^i,Null,{.01,0}}],{i,Floor[min],Ceiling[max],1},{j,1,9}],1];
LogTickLegend[min_,max_]:=Flatten[Table[If[j==1,{Log10[j*10^i],Superscript[10,IntegerPart[i]],{.3,0}},{Log10[j*10^i],Null,{.1,0}}],{i,Floor[min],Ceiling[max],1},{j,1,9}],1];

reportColorRange[plotFunction_]:=
	Module[
		{p,min,max,plotHead,plotBody,colFunc,colScale,h,b,cf,cfs},
		{plotHead,plotBody}=First@Cases[Hold[plotFunction],h_[b__]->{h,Hold[b]},1];(*plotBody is kept inside Hold because it will later be wrapped in plotHead which may have attribute HoldAll*)
		colFunc=Replace[(*Replace wraps string colorfunction names in ColorData[...]*)
			First@Join[
				Cases[plotBody,HoldPattern[ColorFunction->cf_]->cf],{ColorData["AlpineColors"]} (*Extract color function or use default*)
			],
			s_String:>ColorData[s]
		];
		colScale=First@Append[
			Cases[plotBody,HoldPattern[ColorFunctionScaling->cfs_]->cfs],
			True (*ColorFunctionScaling is True by default*)
		];
		(*Turn off ColorFunction and scaling:*)
		plotBody=plotBody/.HoldPattern[ColorFunction->_]|HoldPattern[ColorFunctionScaling->_]->Sequence[];
		(*Make plot with Hue because it takes a single argument that's linear in the heigh value.*)
		{min,max}={Min[#],Max[#]}&@Flatten@Last@Reap[
			p=Apply[
				plotHead,
				Join[(*Join creates a single Hold[...] expression,and Apply replaces hold with plotHead:*)
					plotBody,
					Hold[
						(*Collect the function values-color function scaling must be turned off:*)
						ColorFunction->{(Sow[#];Hue[#])&},ColorFunctionScaling->False
					]
					(*Hue now could have arguments outside the interval[0,1].We'll recover them when restoring the original ColorFunction.*)
				]
			]
		];(*Reap collects the unscaled height values,and we keep the extremal values max,min to rescale if desired:*)
		{
			If[(*If p consists of ploygons colored by Hue,convert it back to the original color function:*)
				Cases[p,Hue[_],Infinity]=!={},If[ (*Recover function values from Hue and plot them with desired color,scaled/unscaled:*) colScale,p/. Hue[x_]:>colFunc[(x-min)/(max-min)],p/. Hue[x_]:>colFunc[x]],
				(*If plotHead is one of my private custom plot functions,then p may be a raster image where Hue no longer appears explicitly.Then we have to re-do the entire plot with the original color function:*)
				plotFunction
			],
			colFunc,
			{min,max}
		}(*The auxiliary Hue is replaced by the original ColorFunction stored in colFunc.*)
	];

SetAttributes[reportColorRange,HoldAll];

trimPoint[n_,digits_]:=(*display number n with given number of sig.digits,trim trailing decimal point*)
	NumberForm[n,digits,NumberFormat->(DisplayForm@RowBox[Join[{StringTrim[#1,RegularExpression["\\.$"]]},If[#3!="",{"\[Times]",SuperscriptBox[#2,#3]},{}]]]&)];

Options[colorLegend]=
	{LabelStyle->Directive[baseStyle],FrameLabel-> None,Background->Transparent,FrameStyle->None,RoundingRadius->10,"ColorSwathes"->None,"LeftLabel"->False,"Digits"->5,Contours->None,BoxFrame->0,"ColorBarFrameStyle"->Black,ImageSize->Automatic};

colorLegend[cFunc_,rawRange_,OptionsPattern[]]:=
	Module[
		{frameticks,tickPositions,nColor,nTick,range=N@Round[rawRange,10^Round[Log10[Abs@First@Differences[{-1.5,.5}]]]/1000],colors,contours=OptionValue[Contours],colorBarLabelStyle=OptionValue[LabelStyle],colorBarFrameStyle=OptionValue["ColorBarFrameStyle"],outerFrameStyle=OptionValue[FrameStyle],colorSwathes=OptionValue["ColorSwathes"],filling,origDim},
		(*Here we decide how many color gradations to diplay-either a given number,equally spaced,or "continuous," i.e.256 steps:*)
		Switch[
			colorSwathes,_?NumericQ,nColor=colorSwathes;
			colors=(Range[nColor]-1/2)/nColor;
			nTick=nColor,_,nColor=256;
			colors=(Range[nColor]-1)/(nColor-1);
			nTick=1
		];
		(*Number of labels is nTick+1,unless changed by numerical Contours setting below:*)
		Switch[
			contours,_?NumericQ,tickPositions=(range[[1]]+(range[[-1]]-range[[1]]) (Range[contours+1]-1)/contours);
			,List[Repeated[_?NumericQ]],tickPositions=contours,_,tickPositions=(range[[1]]+(range[[-1]]-range[[1]]) (Range[nTick+1]-1)/nTick);
		];
		frameticks=
			{If[TrueQ[OptionValue["LeftLabel"]],Reverse[#],#]&@{None,LogTickLegend[range[[1]],range[[-1]]]},{None,None}};
			(*{If[TrueQ[OptionValue["LeftLabel"]],Reverse[#],#]&@{None,Function[{min,max},{#,trimPoint[#,OptionValue["Digits"]],{0,.1}}&/@tickPositions]},{None,None}};*)

		filling=
		Graphics[(*Create strip of colored,translated unit squares.If colorSwathes are selected,colors don't vary inside squares.Otherwise,colors vary linearly in each of 256 squares to get smooth gradient using VertexColors:*)
				MapIndexed[
					{Translate[Polygon[{{0,0},{1,0},{1,1},{0,1}},VertexColors->{cFunc[#[[1]]],cFunc[#[[1]]],cFunc[#[[2]]],cFunc[#[[2]]]}],{0,#2[[1]]-1}]}&,
					Transpose[
						If[
							colorSwathes===None,
							{Most[colors],Rest[colors]}(*Offset top versus bottom colors of polygons to create linear VertexColors*),
							{colors,colors} (*Top and bottom colors are same when uniform colorSwathes are desired*)
						]
					]
				], (**End MapIndexed**)
						(*Options for inset Graphics:*)
						ImagePadding->0,PlotRangePadding->0,AspectRatio->Full(*AspectRatio\[Rule]Full allows colored squares to strecth with resizing in the following.*)
			];
		origDim=ImageDimensions[#]&@filling;
		(*DisplayForm@FrameBox replaces Framed because it allows additional BoxFrame option to specify THICKNESS of frame:*)
		DisplayForm@FrameBox[ 
			(*Wrapped in Pane to allow unlimited resizing:*)
			Pane@Graphics[
				Inset[
					ImageResize[#,{Last[origDim]/8,Last[origDim]}]&@filling//Rasterize,
					(*Options for Inset:*)
					{0,First[range]},{0,0},{(range[[-1]]-range[[1]])/8,range[[-1]]-range[[1]]},ContentSelectable->True (*this sets the size of the inset in the enclosing Graphics whose PlotRange is given next:*)
				],
				(*Options for enclosing Graphics:*)
				PlotRange->{{0,(range[[-1]]-range[[1]])/8},range[[{1,-1}]]},Frame->True,FrameLabel->OptionValue[FrameLabel],LabelStyle-> colorBarLabelStyle,FrameTicks->frameticks,FrameTicksStyle->{colorBarLabelStyle,Black},FrameStyle->colorBarFrameStyle,ImageSize->OptionValue[ImageSize]],
		(*Options for FrameBox:*)
		Background->OptionValue[Background],FrameStyle->outerFrameStyle,RoundingRadius->OptionValue[RoundingRadius],BoxFrame->OptionValue[BoxFrame]
		]
	];
Options[colorLegendLinear]=
	{LabelStyle->Directive[baseStyle],FrameLabel-> None,Background->Transparent,FrameStyle->None,RoundingRadius->10,"ColorSwathes"->None,"LeftLabel"->False,"Digits"->3,Contours->None,BoxFrame->0,"ColorBarFrameStyle"->Black,ImageSize->Automatic};

colorLegendLinear[cFunc_,rawRange_,OptionsPattern[]]:=
	Module[
		{frameticks,tickPositions,nColor,nTick,range=If[rawRange[[-1]]<1,{Floor[rawRange[[1]],0.01],Floor[rawRange[[-1]],0.01]},{Floor[rawRange[[1]],10],Floor[rawRange[[-1]],10]}],colors,contours=OptionValue[Contours],colorBarLabelStyle=OptionValue[LabelStyle],colorBarFrameStyle=OptionValue["ColorBarFrameStyle"],outerFrameStyle=OptionValue[FrameStyle],colorSwathes=OptionValue["ColorSwathes"],filling,origDim},
		(*Here we decide how many color gradations to diplay-either a given number,equally spaced,or "continuous," i.e.256 steps:*)
		Switch[
			colorSwathes,_?NumericQ,nColor=colorSwathes;
			colors=(Range[nColor]-1/2)/nColor;
			nTick=nColor,_,nColor=256;
			colors=(Range[nColor]-1)/(nColor-1);
			nTick=1
		];
		(*Number of labels is nTick+1,unless changed by numerical Contours setting below:*)
		Switch[
			contours,_?NumericQ,tickPositions=(range[[1]]+(range[[-1]]-range[[1]]) (Range[contours+1]-1)/contours);
			,List[Repeated[_?NumericQ]],tickPositions=contours,_,tickPositions=(range[[1]]+(range[[-1]]-range[[1]]) (Range[nTick+1]-1)/nTick);
		];
		frameticks=
			(*{If[TrueQ[OptionValue["LeftLabel"]],Reverse[#],#]&@{None,LogTickLegend[range[[1]],range[[-1]]]},{None,None}};*)
			{If[TrueQ[OptionValue["LeftLabel"]],Reverse[#],#]&@{None,Function[{min,max},{#,trimPoint[#,OptionValue["Digits"]],{-0.3,0}}&/@tickPositions]},{None,None}};

		filling=
		Graphics[(*Create strip of colored,translated unit squares.If colorSwathes are selected,colors don't vary inside squares.Otherwise,colors vary linearly in each of 256 squares to get smooth gradient using VertexColors:*)
				MapIndexed[
					{Translate[Polygon[{{0,0},{1,0},{1,1},{0,1}},VertexColors->{cFunc[#[[1]]],cFunc[#[[1]]],cFunc[#[[2]]],cFunc[#[[2]]]}],{0,#2[[1]]-1}]}&,
					Transpose[
						If[
							colorSwathes===None,
							{Most[colors],Rest[colors]}(*Offset top versus bottom colors of polygons to create linear VertexColors*),
							{colors,colors} (*Top and bottom colors are same when uniform colorSwathes are desired*)
						]
					]
				], (**End MapIndexed**)
						(*Options for inset Graphics:*)
						ImagePadding->0,PlotRangePadding->0,AspectRatio->Full(*AspectRatio\[Rule]Full allows colored squares to strecth with resizing in the following.*)
			];
		origDim=ImageDimensions[#]&@filling;
		(*DisplayForm@FrameBox replaces Framed because it allows additional BoxFrame option to specify THICKNESS of frame:*)
		DisplayForm@FrameBox[ 
			(*Wrapped in Pane to allow unlimited resizing:*)
			Pane@Graphics[
				Inset[
					ImageResize[#,{Last[origDim]/8,Last[origDim]}]&@filling//Rasterize,
					(*Options for Inset:*)
					{0,First[range]},{0,0},{(range[[-1]]-range[[1]])/8,range[[-1]]-range[[1]]},ContentSelectable->True (*this sets the size of the inset in the enclosing Graphics whose PlotRange is given next:*)
				],
				(*Options for enclosing Graphics:*)
				PlotRange->{{0,(range[[-1]]-range[[1]])/8},range[[{1,-1}]]},Frame->True,FrameLabel->OptionValue[FrameLabel],LabelStyle-> colorBarLabelStyle,FrameTicks->frameticks,FrameTicksStyle->{colorBarLabelStyle,Black},FrameStyle->colorBarFrameStyle,ImageSize->OptionValue[ImageSize]],
		(*Options for FrameBox:*)
		Background->OptionValue[Background],FrameStyle->outerFrameStyle,RoundingRadius->OptionValue[RoundingRadius],BoxFrame->OptionValue[BoxFrame]
		]
	];

at[position_,scale_:Automatic][obj_]:=
	(*convenience function to position objects in Graphics*)
	Inset[obj,position,{Left,Bottom},scale];

display[g_,opts:OptionsPattern[]]:=
	Module[	
		{frameOptions=FilterRules[{opts},Options[Graphics]]},
		(*Same as Graphics,but with fixed PlotRange*)
		Graphics[g,PlotRange->{{0,1},{0,1}},Evaluate@Apply[Sequence,frameOptions]]
	];

cropVector[g_,x_,y_,w_,h_]:=
	Graphics[
		Inset[g,{x,y},{0,0}],
		PlotRange->{{0,1},{0,1}},
		ImageSize->{w,h},
		AspectRatio->Full
	];

CleanContourPlot[cp_]:=Module[{points,groups,regions,lines},groups=Cases[cp,{style__,g_GraphicsGroup}:>{{style},g},Infinity];
points=First@Cases[cp,GraphicsComplex[pts_,___]:>pts,Infinity];
regions=Table[Module[{group,style,polys,edges,cover,graph},{style,group}=g;
polys=Join@@Cases[group,Polygon[pt_,___]:>pt,Infinity];
edges=Join@@(Partition[#,2,1,1]&/@polys);
cover=Cases[Tally[Sort/@edges],{e_,1}:>e];
graph=Graph[UndirectedEdge@@@cover];
{Sequence@@style,FilledCurve[List/@Line/@First/@Map[First,FindEulerianCycle/@(Subgraph[graph,#]&)/@ConnectedComponents[graph],{3}]]}],{g,groups}];
lines=Cases[cp,_Tooltip,Infinity];
Graphics[GraphicsComplex[points,{regions,lines}],Sequence@@Options[cp]]];
(*setup for 90CL plots*)
expCol=<|
"BEBC"->{Opacity[0.5,Gray],None},
"CHARM"->{Opacity[0.5,Darker[Gray,0.2]],None},
"NuCal"->{Opacity[0.6,Lighter[Gray,0.3]],None},
"DarkQuest"->{Transparent,Opacity[0.5,Darker[Green,0.5]]},
"DarkQuestPhase2"->{Transparent,Opacity[0.5,Darker[Green,0.7]]},
"DUNE"->{Transparent,Opacity[0.5,Red]},
"ORCA"->{Transparent,Opacity[0.5,Darker[Yellow,0.4]]},
"NuTeV"->{Opacity[0.7,Darker[Gray,0.5]],None},
"NA62"->{Transparent,Opacity[0.5,Darker[Red,0.6]]},
"HIKE"->{Transparent,Opacity[0.5,Darker[Red,0.6]]},
"SHADOWS"->{Transparent,Opacity[0.5,Black]},
"KLEVER"->{Transparent,Opacity[0.5,Darker[Yellow,0.7]]},
"KLEVERext"->{Transparent,Opacity[0.5,Darker[Yellow,0.4]]},
"SHiP"->{Transparent,Opacity[0.5,RGBColor[0.05,0.2,0.75]]},
"SHiPecn4"->{Transparent,Opacity[0.5,RGBColor[0.05,0.2,0.75]]},
"E137"->{Lighter[Blue,0.8],None},
"E141"->{Lighter[Blue,0.8],None},
"KOTOpnn"->{Opacity[0.7,Orange],None},
"KOTOexclPnn"->{Transparent,{Thickness[0.008],Opacity[0.7,Orange]}},
"KOTOdump"->{Transparent,Opacity[0.7,Green]},
"KOTO2pnn"->{Transparent,{Thickness[0.008],Opacity[0.7,Darker[Orange,0.5]]}},
"KOTO2dump"->{Transparent,Opacity[0.7,Darker[Green,0.5]]}
|>;
expLeg=<|
"BEBC"->{"BEBC",{Thickness[0.05],Opacity[0.5,Gray]}},
"CHARM"->{"CHARM",{Thickness[0.05],Opacity[0.5,Darker[Gray,0.2]]}},
"NuCal"->{"NuCal",{Thickness[0.05],Opacity[0.6,Lighter[Gray,0.3]]}},
"DarkQuest"->{"DarkQuest",Opacity[0.5,Darker[Green,0.5]]},
"DarkQuestPhase2"->{"\!\(\*SubscriptBox[\(DarkQuest\), \(Phase2\)]\)",Opacity[0.5,Darker[Green,0.5]]},
"DUNE"->{"DUNE",Opacity[0.5,Red]},
"ORCA"->{"ORCA",Opacity[0.5,Darker[Yellow,0.4]]},
"NuTeV"->{"NuTeV",{Thickness[0.05],Opacity[0.7,Darker[Gray,0.5]]}},
"NA62"->{"NA62 (\!\(\*SuperscriptBox[\(10\), \(18\)]\))",Opacity[0.5,Darker[Red,0.6]]},
"HIKE"->{"HIKE (5\!\(\*SuperscriptBox[\(\[Times]10\), \(19\)]\))",Opacity[0.5,Darker[Red,0.6]]},
"SHADOWS"->{"SHADOWS (5\!\(\*SuperscriptBox[\(\[Times]10\), \(19\)]\))",Opacity[0.5,Black]},
"KLEVER"->{"KLEVER (\!\(\*SuperscriptBox[\(10\), \(18\)]\))",{Transparent,Opacity[0.5,Darker[Yellow,0.7]]}},
"KLEVERext"->{"\!\(\*SubscriptBox[\(KLEVER\), \(EXT\)]\) (\!\(\*SuperscriptBox[\(10\), \(18\)]\))",{Transparent,Opacity[0.5,Darker[Yellow,0.4]]}},
"SHiP"->{"SHiP",Opacity[0.5,RGBColor[0.05,0.2,0.75]]},
"SHiPecn4"->{"\!\(\*SubscriptBox[\(SHiP\), \(ECN4\)]\) (\!\(\*SuperscriptBox[\(10\), \(20\)]\))",Opacity[0.5,RGBColor[0.05,0.2,0.75]]},
"E137"->{"E137, E141",{Thickness[0.05],Lighter[Blue,0.8]}},
"E141"->{"E141",{Thickness[0.05],Lighter[Blue,0.8]}},
"KOTOpnn"->{"KOTO (2015)",{Thickness[0.05],Opacity[0.7,Orange]}},
"KOTOexclPnn"->{"KOTO",{Thickness[0.008],Opacity[0.7,Orange]}},
"KOTOdump"->{"\!\(\*SubscriptBox[\(KOTO\), \(DUMP\)]\)",Opacity[0.7,Green]},
"KOTO2pnn"->{"KOTO2",{Thickness[0.008],Opacity[0.7,Darker[Orange,0.5]]}},
"KOTO2dump"->{"\!\(\*SubscriptBox[\(KOTO\), \(DUMP\)]\)2",Opacity[0.7,Darker[Green,0.5]]}
|>;
expCont=<|
"BEBC"->2.3,
"CHARM"->2.3,
"NuCal"->3.6,
"DarkQuest"->10.,(*1.44E18 + background*)
"DarkQuestPhase2"->10.,(*1E20 + background*)
"DUNE"->0.23,(*1E21 x 10 (yr of data-taking)*)
"ORCA"->2.3,(*8E19*)
"NuTeV"->2.3,
"NA62"->2.3,(*1E18*)
"HIKE"->0.046,(*1E18 x 50*)
"SHADOWS"->0.46,(*1E19 x 5*)
"KLEVER"->2.3,(*1E18*)
"KLEVERext"->2.3,(*1E18*)
"SHiP"->2.3,(*6E20*)
"SHiPecn4"->2.3,(*1E20*)
"E137"->2.3,
"E141"->2.3,
"KOTOpnn"->2.3,(*already scaled*)
"KOTOexclPnn"->4.2,
"KOTOdump"->0.23,(*10x statistics*)
"KOTO2pnn"->3.14,
"KOTO2dump"->0.23(*10x statistics*)
|>;
baseLegend90CL[expList_]:=LineLegend[{Map[expLeg[#][[2]]&,expList]},{Map[expLeg[#][[1]]&,expList]},Spacings-> 0.15,LabelStyle->baseStyleLeg90CL];
dashedLegend90CL=LineLegend[{Dashed,Darker[Gray,0.9]},{"2\[Gamma] only"},Spacings-> 0.2,LabelStyle->baseStyleLeg90CL];
regionsLegendDS90CL=LineLegend[{Directive[Lighter[Orange,0.8],Thickness[0.05]],Directive[Lighter[Blue,0.8],Thickness[0.05]]},{"LHCb","NA62 Run1 (\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\)\[RightArrow]\[Pi]\[Nu]\[Nu], \!\(\*SuperscriptBox[
StyleBox[\"\[Pi]\",\nFontSlant->\"Italic\"], \(0\)]\)\[RightArrow]inv)"},Spacings-> 0.15,LabelStyle->baseStyleLeg90CL];

gridLines3mesons={Method->{"GridLinesInFront"->True},GridLines->{{0.135,0.547,0.957},{}},GridLinesStyle->{{Darker[Gray,0.4],Thickness->0.01},{}}};
gridLines2mesons={Method->{"GridLinesInFront"->True},GridLines->{{0.547,0.94},{}},GridLinesStyle->{{Darker[Gray,0.4],Thickness->0.025},{}}};
gridLines2mesonsThin={Method->{"GridLinesInFront"->True},GridLines->{{0.547,0.94},{}},GridLinesStyle->{{Darker[Gray,0.4],Thickness->0.01},{}}};
gridLines3vectors={Method->{"GridLinesInFront"->True},GridLines->{{0.77526,0.78265,1.019461},{}},GridLinesStyle->{{Darker[Gray,0.4],Thickness->0.01},{}}};

yLabelxMass[ylabel_]:=FrameLabel->{StringForm["`` [GeV]",Subscript[Style["m",Italic],Style["a",Italic]]],StringForm["`` [``]",ylabel,Superscript["GeV",-1]]};
yLabelNoDimxMassIndex[ylabel_,xindex_]:=FrameLabel->{StringForm["`` [GeV]",Subscript[Style["m",Italic],Style[xindex,Italic]]],StringForm["``",ylabel]};
plotSettings90CLALPff={PlotRange->{{0.01,3.},{1*10^-10,1*10^-3},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thickness[0.005],FrameTicks->{{LogTick[-10,-3],None},{LogTick[-2,0],None}},LabelStyle->Directive[baseStyle]};
plotSettings90CL={PlotRange->{{0.0001,2.9},{1*10^-9,1*10^-1},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thickness[0.005],FrameTicks->{{LogTick[-11,-1],None},{LogTick[-4,0],None}},LabelStyle->Directive[baseStyle]};
plotReducedSettings90CL={PlotRange->{{0.4,2.9},{1*10^-10,1*10^-5},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thick,FrameTicks->{{LogTick[-11,-1],None},{LogTick[-4,0],None}},LabelStyle->Directive[baseStyle]};
plotKOTOSettings90CL={PlotRange->{{0.04,1.},{1*10^-9,1*10^-4},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thickness[0.005],FrameTicks->{{LogTick[-10,-1],None},{LogTick[-2,0],None}},LabelStyle->Directive[baseStyle]};

plotSettings90CLDS={PlotRange->{{0.01,5.},{1*10^-13,1*10^-3},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thickness[0.005],FrameTicks->{{LogTick[-13,-3],None},{LogTick[-2,0],None}},LabelStyle->Directive[baseStyle]};

plotSettings90CLDP={PlotRange->{{0.01,5.},{1*10^-10,1*10^-3},{10^-8,Full}},PlotRangePadding-> None,ScalingFunctions->{"Log10","Log10"},ImageSize->Medium,BaseStyle-> Thickness[0.005],FrameTicks->{{LogTick[-10,-3],None},{LogTick[-3,0],None}},LabelStyle->Directive[baseStyle]};

ALPPlotStd[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotSettings90CL]]];

ALPPlotStdMask[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotSettings90CL],Evaluate[gridLines3mesons]]];

ALPPlotReduced[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotReducedSettings90CL],Evaluate[gridLines2mesons]]];

ALPPlotDashed[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,Transparent},ContourStyle->{Dashed,If[expCol[experiment][[2]]==None,expCol[experiment][[1]],expCol[experiment][[2]]]},Evaluate[plotReducedSettings90CL],Evaluate[gridLines2mesons]]];

ALPPlotKOTO[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotKOTOSettings90CL]]];

ALPPlotFFMask[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotSettings90CLALPff],Evaluate[gridLines3mesons]]];

DSPlot[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotSettings90CLDS]]];

DPPlot[data_,experiment_]:=CleanContourPlot[ListContourPlot[data,Contours->{expCont[experiment]},ContourShading->{Transparent,expCol[experiment][[1]]},ContourStyle->expCol[experiment][[2]],Evaluate[plotSettings90CLDP],Evaluate[gridLines3vectors]]];

End[]

EndPackage[]
