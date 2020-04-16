#!/usr/bin/env wolframscript
(* ::Package:: *)

If[Dimensions[$ScriptCommandLine][[1]]<4, 
  Print["ERROR: ./X mHNLi mHNLf deltamHNL"];
  Exit[];
  ]

mhnli=ToExpression[$ScriptCommandLine[[2]]];
mhnlf=ToExpression[$ScriptCommandLine[[3]]];
paso=ToExpression[$ScriptCommandLine[[4]]];

hnlrangetest=Mod[IntegerPart[(mhnlf-mhnli)*1000+0.00001],IntegerPart[paso*1000+0.00001]];
If [hnlrangetest!=0., 
  Print["ERROR: Rango de mHNL incorrecto."];
  Exit[];
  ]

Print["Config OK"]

