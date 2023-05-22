Get["PacletManager`PacletManager`"]
SetOptions[{Plot, LogPlot, LogLogPlot, LogLinearPlot}, {Frame -> True, FrameStyle -> Black,
        ImageSize -> Large,
        LabelStyle -> Directive[Large, Black, FontFamily -> "Palatino"],
        PlotStyle -> AbsoluteThickness[3]}]
SetOptions[{ListPlot, ListLogPlot, ListLogLogPlot, ListLogLinearPlot,
    ListLinePlot}, {Frame -> True, FrameStyle -> Black,
        ImageSize -> Large,
        LabelStyle -> Directive[Large, Black, FontFamily -> "Palatino"],
        PlotStyle -> AbsoluteThickness[3], Joined -> True}]
SetOptions[{ContourPlot,ListContourPlot,ParametricPlot}, {Frame -> True, FrameStyle -> Black,
    LabelStyle -> Directive[Large, Black, FontFamily -> "Palatino"],
    ImageSize -> Large,
    AspectRatio -> 1/GoldenRatio
}]
SetOptions[{ParametricPlot}, {Frame -> True, FrameStyle -> Black,
    LabelStyle -> Directive[Large, Black, FontFamily -> "Palatino"],
    ImageSize -> Large,
    AspectRatio -> 1/GoldenRatio,
    PlotStyle -> AbsoluteThickness[3]
}]
SetOptions[{Plot3D,ParametricPlot3D},
	   {ImageSize -> Large,
	       BoxStyle -> Black,
	       LabelStyle -> Directive[Large, Black, FontFamily -> "Palatino"],
	       PlotStyle -> AbsoluteThickness[3]}]
SetOptions[{ListDensityPlot},
	   {ImageSize -> Large,
	       FrameStyle -> Black,
	       LabelStyle -> Directive[Large, Black, FontFamily ->"Palatino"]}]

BeginPackage["Color`"]
Color::usage = "Mathematica's default color set"
  Begin["`Private`"]
  RGBData = {"#5E81B5","#E19C24","#8FB032","#EB6235","#8778B3","#C56E1A","#5D9EC7","#FFBF00","#A5609D","#929600","#E95536","#6685D9","#F89F13","#BC5B80","#47B66D"};
Table[Color[i] = RGBColor[RGBData[[i]]],{i,Length[RGBData]}]
End[]
EndPackage[]
  
