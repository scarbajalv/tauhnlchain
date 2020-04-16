#!/usr/bin/env wolframscript
(* ::Package:: *)

#!/usr/bin/env wolframscript

Print["Dirac meson2HNL production & decay"];

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


(*HNL DECAYS TOO MESONS/LEPTONS*)(*ESTE SCRIPT GENERA LOS ARCHIVOS \
PCONFIG QUE UTILIZA PYTHIA8.*)(*NOTAR QUE SIEMPRE TOMA COMO ARGUMENTO \
LA MASA DEL HNL*)
e = 0.5109989461*10^-3;
mu = 105.6583745*10^-3;
tau = 1.77686;
leptons = {e, mu, tau};
tfactor = 1.5193*10^24(*Para pasar de s a GeV^-1*); 
mmc =  3*10^11(*para pasar de s a mm/c*); 
ckm = ({{0.9742, 	0.2243, 0.0365}, 
				{0.218, 	0.997, 	0.0422}, 
				{0.0081, 	0.0394, 1.019}});
gf = 1.166378*10^-5;
xw = 0.2223; (*sin^2(theta_w)*)
gll = -(1/2) + xw;
glr = xw;

(* Extras Semileptonic Bondarenko *)
d421 = {410.1*10^-15*tfactor, 1.86483, Null, Null, 421};
k311 = {Null, 0.497611, Null, Null, 311};
k313 = {Null, 0.89555, Null, Null, 313};

(*charged pseudoscalar*)
d431 = {504*10^-15*tfactor, 1.96834, 249.0*10^-3, ckm[[2, 2]], 431, 0.0555};
d431bar = {504*10^-15*tfactor, 1.96834, 249.0*10^-3, ckm[[2, 2]], -431};
d411 = {1040*10^-15*tfactor, 1.86958, 211.9*10^-3, ckm[[2, 1]], 411, 0.0012};
d411bar = {1040*10^-15*tfactor, 1.86958, 211.9*10^-3, ckm[[2, 1]], -411};
k321 = {1.238*10^-8*tfactor, 493.677*10^-3, 155.6*10^-3, ckm[[1, 2]], 321};
pi211 = {2.6*10^-8*tfactor, 139.57061*10^-3, 130.2*10^-3, ckm[[1, 1]], 211};

(*Charged vector*)
k323 = {Null, 895.55*10^-3, 0.1827, ckm[[1, 2]], 323}(*K*(892)*);
rho770c = {Null, 775.26*10^-3, 0.162, ckm[[1, 1]]};

(*neutral pseudoscalar*)
pi111 = {Null, 134.977*10^-3, 130.2*10^-3, Null, 111};
eta221 = {Null, 547.862*10^-3, 81.7*10^-3, Null, 221}(*eta*); 
eta331 = {Null, 957.78*10^-3, -94.7*10^-3, Null, 331}(*eta prime*);

(*Neutral vector*)
rho770n = {Null, 775.26*10^-3, 0.162, 1 - 2 xw};
omega782 = {Null, 782.65*10^-3, 0.153, 4/3 xw};
phi1020 = {Null, 1.019461, 0.234, 4/3 xw - 1, 333};

lambda[a_, b_, c_] := a^2 + b^2 + c^2 - 2 a*b - 2 a*c - 2 b*c;

i1[xu_, xd_, xl_] := 
  Re[12*NIntegrate[
     Re[1/x (x - xl^2 - xd^2) (1 + xu^2 - x) Sqrt[
       lambda[x, xl^2, xd^2]*lambda[1, x, xu^2]]], {x, (xd + 
        xl)^2, (1 - xu)^2}]];

datae=Import["./mainconfig/mixe.dat"];
datamu=Join[Table[{i,1},{i,0.01,0.1,0.01}],Import["./mainconfig/mixmu.dat"]];
datatau=Import["./mainconfig/mixtau.dat"];
fmixe=Interpolation[datae,InterpolationOrder->1];
fmixmu=Interpolation[datamu,InterpolationOrder->1];
fmixtau=Interpolation[datatau,InterpolationOrder->1];

(*****************************************************************)
(***************** HNL DECAY W (BONDARENKO) **********************)

(* 1. (3.1) l1 nu2 l2 *)
width1[l1_, l2_, n_, u_] := 
  If[n > l1 + l2, 
  	Re[(gf^2 n^5)/(192 Pi^3)*u^2*i1[0, (l2/n), (l1/n)]],
    0];
enumumu[n_, u_] := width1[e, mu, n, u];
enutautau[n_, u_] := Re[width1[e, tau, n, u] UnitStep[n - e - tau]];
munuee[n_, u_] := Re[width1[mu, e, n, u] UnitStep[n - mu - e]];
munutautau[n_, u_] := Re[width1[mu, tau, n, u] UnitStep[n - mu - tau]];
taunuee[n_, u_] := Re[width1[tau, e, n, u] UnitStep[n - tau - e]];
taunumumu[n_, u_] := 
  Re[width1[tau, mu, n, u] UnitStep[n - tau - mu]];
totalw1[n_, mixe_, mixmu_, mixtau_] := 
  enumumu[n, mixe] + enutautau[n, mixe] + munuee[n, mixmu] + 
   munutautau[n, mixmu] + taunuee[n, mixtau] + taunumumu[n, mixtau];

(* 2a. (3.4) nu1 f2 f2bar *)
c1fa = 1/4 (1 - 4 xw + 8 xw^2);
c2fa = 1/2 xw*(2*xw - 1);
width2a[l_, n_, u_] := If[n > 2*l,
   Module[{x, laux, g},
    x = l/n;
    g[t_] := 1 - 3 t^2 - (1 - t^2) Sqrt[1 - 4 t^2];
    If[l == e, (* La masa del electrón es muy pequeña asi que usamos una serie de \
								Taylor para g. *)     
     laux = Log[((D[g[t], {t, 6}] /. t -> 0)/6!*(x)^6 + (
           D[g[t], {t, 8}] /. t -> 0)/8!*(x)^8 + (
           D[g[t], {t, 10}] /. t -> 0)/
           10!*(x)^10)/(x^2 (1 + Sqrt[1 - 4 x^2]))],
     laux = 
       Log[(1 - 3 x^2 - (1 - x^2) Sqrt[1 - 4 x^2])/(
        x^2 (1 + Sqrt[1 - 4 x^2]))];];
    Re[(gf^2 n^5)/(192 Pi^3)*
      u^2 (c1fa ((1 - 14 x^2 - 2 x^4 - 12 x^6) Sqrt[1 - 4 x^2] + 
           12 x^4 (x^4 - 1) laux) + 
        4 c2fa (x^2 (2 + 10 x^2 - 12 x^4) Sqrt[1 - 4 x^2] + 
           6 x^4 (1 - 2 x^2 + 2 x^4) laux))]],
   0];
nuemumu[n_, u_] := width2a[mu, n, u];
nuetautau[n_, u_] := width2a[tau, n, u];
numuee[n_, u_] := width2a[e, n, u];
numutautau[n_, u_] := width2a[tau, n, u];
nutauee[n_, u_] := width2a[e, n, u];
nutaumumu[n_, u_] := width2a[mu, n, u];
totalw2a[n_, mixe_, mixmu_, mixtau_] := 
  nuemumu[n, mixe] + nuetautau[n, mixe] + numuee[n, mixmu] + 
   numutautau[n, mixmu] + nutauee[n, mixtau] + nutaumumu[n, mixtau];

(* 2b. (3.4) nu1 f1 f1bar *)
c1fb = 1/4 (1 + 4 xw + 8 xw^2);
c2fb = 1/2 xw*(2*xw + 1);
width2b[l_, n_, u_] := If[n > 2*l, Module[{x, laux, g},
    x = l/n;
    g[t_] := 1 - 3 t^2 - (1 - t^2) Sqrt[1 - 4 t^2];
    If[l == e,
     laux = 
      Log[((D[g[t], {t, 6}] /. t -> 0)/6!*(x)^6 + (
           D[g[t], {t, 8}] /. t -> 0)/8!*(x)^8 + (
           D[g[t], {t, 10}] /. t -> 0)/
           10!*(x)^10)/(x^2 (1 + Sqrt[1 - 4 x^2]))],
     laux = 
       Log[(1 - 3 x^2 - (1 - x^2) Sqrt[1 - 4 x^2])/(
        x^2 (1 + Sqrt[1 - 4 x^2]))];];
    Re[(gf^2 n^5)/(192 Pi^3)*
      u^2 (c1fb ((1 - 14 x^2 - 2 x^4 - 12 x^6) Sqrt[1 - 4 x^2] + 
           12 x^4 (x^4 - 1) laux) + 
        4 c2fb (x^2 (2 + 10 x^2 - 12 x^4) Sqrt[1 - 4 x^2] + 
           6 x^4 (1 - 2 x^2 + 2 x^4) laux))]], 0];
nueee[n_, u_] := width2b[e, n, u];
numumumu[n_, u_] := width2b[mu, n, u];
nutautautau[n_, u_] := width2b[tau, n, u];
totalw2b[n_, mixe_, mixmu_, mixtau_] := 
  nueee[n, mixe] + numumumu[n, mixmu] + nutautautau[n, mixtau];

(* 3a. (3.5) nu1 nu2 nu2bar *)
width3a[n_, u_] := Re[(gf^2 n^5)/(768 Pi^3) u^2];
totalw3a[n_, mixe_, mixmu_, mixtau_] := 
  2*width3a[n, mixe] + 2*width3a[n, mixmu] + 2*width3a[n, mixtau];

(* 3b. (3.5) nu1 nu1 nu1bar *)
width3b[n_, u_] := Re[(2*gf^2 n^5)/(768 Pi^3) u^2];
totalw3b[n_, mixe_, mixmu_, mixtau_] := 
  width3b[n, mixe] + width3b[n, mixmu] + width3b[n, mixtau];

(* 4. (3.6) l- P+ *)
width4[meson_, l_, n_, u_] := If[n > meson[[2]] + l, Module[{xl, xh},
    xl = l/n;
    xh = meson[[2]]/n;
    Re[(gf^2 meson[[3]]^2 meson[[4]]^2 u^2 n^3)/(
      16 Pi) ((1 - xl^2)^2 - xh^2 (1 + xl^2)) Sqrt[
      lambda[1, xh^2, xl^2]]]], 0];
ePi[n_, u_] := width4[pi211, e, n, u];
eK[n_, u_] := width4[k321, e, n, u];
eD[n_, u_] := width4[d411, e, n, u];
muPi[n_, u_] := width4[pi211, mu, n, u];
muK[n_, u_] := width4[k321, mu, n, u];
muD[n_, u_] := width4[d411, mu, n, u];
tauPi[n_, u_] := width4[pi211, tau, n, u];
tauK[n_, u_] := width4[k321, tau, n, u];
tauD[n_, u_] := width4[d411, tau, n, u];
totalw4[n_, mixe_, mixmu_, mixtau_] := 
  ePi[n, mixe] + eK[n, mixe] + eD[n, mixe] + muPi[n, mixmu] + 
   muK[n, mixmu] + muD[n, mixmu] + tauPi[n, mixtau] + 
   tauK[n, mixtau] + tauD[n, mixtau];

(* 5. (3.7) nu P0 *)
width5[meson_, n_, u_] := If[n > meson[[2]], Module[{xh},
    xh = meson[[2]]/n;
    Re[(gf^2 meson[[3]]^2 n^3)/(32 Pi)*u^2*(1 - xh^2)^2]], 0];
nupi0[n_, u_] := width5[pi111, n, u];
nueta[n_, u_] := width5[eta221, n, u];
nuetaprime[n_, u_] := width5[eta331, n, u];
totalw5[n_, mixe_, mixmu_, mixtau_] := 
  nupi0[n, mixe] + nupi0[n, mixmu] + nupi0[n, mixtau] + 
   nueta[n, mixe] + nueta[n, mixmu] + nueta[n, mixtau] + 
   nuetaprime[n, mixe] + nuetaprime[n, mixmu] + 
   nuetaprime[n, mixtau];

(* 6. (3.8) l- V+ *)
width6[meson_, l_, n_, u_] := If[n > meson[[2]] + l, 
	Module[{xl, xh},
    xl = l/n;
    xh = meson[[2]]/n;
    Re[(gf^2*meson[[3]]^2 meson[[4]]^2 u^2 n^3)/(
      16 Pi) ((1 - xl^2)^2 + xh^2 (1 + xl^2) - 2 xh^4) Sqrt[
      lambda[1, xh^2, xl^2]]]], 
   0];
eKstar[n_, u_] := width6[k323, e, n, u];
erhoc[n_, u_] := width6[rho770c, e, n, u];
muKstar[n_, u_] := width6[k323, mu, n, u];
murhoc[n_, u_] := width6[rho770c, mu, n, u];
tauKstar[n_, u_] := width6[k323, tau, n, u];
taurhoc[n_, u_] := width6[rho770c, tau, n, u];
totalw6[n_, mixe_, mixmu_, mixtau_] := 
  eKstar[n, mixe] + erhoc[n, mixe] + muKstar[n, mixmu] + 
   murhoc[n, mixmu] + tauKstar[n, mixtau] + taurhoc[n, mixtau];

(* 7. (3.9) nu V0 *)
width7[meson_, n_, u_] := If[n > meson[[2]], Module[{xh},
    xh = meson[[2]]/n;
    Re[(gf^2 meson[[4]]^2 meson[[3]]^2 u^2 n^3)/(
      32 Pi) (1 + 2 xh^2) (1 - xh^2)^2]], 0];
nurho[n_, u_] := width7[rho770n, n, u];
nuomega[n_, u_] := width7[omega782, n, u];
nuphi[n_, u_] := width7[phi1020, n, u];
totalw7[n_, mixe_, mixmu_, mixtau_] := 
  nurho[n, mixe] + nurho[n, mixmu] + nurho[n, mixtau] + 
   nuomega[n, mixe] + nuomega[n, mixmu] + nuomega[n, mixtau] + 
   nuphi[n, mixe] + nuphi[n, mixmu] + nuphi[n, mixtau];

(* TOTAL W *)
totalw[n_, mixe_, mixmu_, mixtau_] := 
 	totalw1[n, mixe, mixmu, mixtau] + 
 	totalw2a[n, mixe, mixmu, mixtau] + totalw2b[n, mixe, mixmu, mixtau] + 
  totalw3a[n, mixe, mixmu, mixtau] + totalw3b[n, mixe, mixmu, mixtau] + 
  totalw4[n, mixe, mixmu, mixtau] +
  totalw5[n, mixe, mixmu, mixtau] + 
  totalw6[n, mixe, mixmu, mixtau] +
  totalw7[n, mixe, mixmu, mixtau]

	(****** HNL ******)
	(* CC Dirac HNL  / 16 canales *)
expBrAll[n_,mixe_,mixmu_,mixtau_]:={
	(* 1. (3.1) l1 nu2 l2 *)
	"2000000001:oneChannel = 1 "<>ToString[AccountingForm[enumumu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 14 -13",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[enutautau[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 16 -15", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munuee[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 12 -11",  
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munutautau[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 16 -15",  (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunuee[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 12 -11",  (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 14 -13",  (* tau *)
	(* 2a. (3.4) nu1 f2 f2bar *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[numuee[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 14 11 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nutauee[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 16 11 -11", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuemumu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 12 13 -13",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nutaumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 16 13 -13", (* tau *)
	(* 2b. (3.4) nu1 f1 f1bar *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueee[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 12 11 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[numumumu[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 22 14 13 -13",
	(* 3a. (3.5) nu1 nu2 nu2bar *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 14 -14",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 16 -16",  (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 12 -12",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 16 -16", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 12 -12",  (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3a[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 14 -14",  (* tau *)
	(* 3b. (3.5) nu1 nu1 nu1bar *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3b[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 12 -12",	
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3b[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 14 -14",	
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[width3b[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 16 -16",  (* tau *)
	(* 4. (3.6) l- P+ *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eD[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 411",  (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eK[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muK[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[ePi[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muPi[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[tauPi[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 211",	(* tau *)
	(* 5. (3.7) nu P0 *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 111",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 111",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 111", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 221",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 221",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 221", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 331",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 331",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 331",  (* tau *)
	(* 6. (3.8) l- V+ *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eKstar[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muKstar[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[erhoc[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 213",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[murhoc[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 213",
	(* 7. (3.9) nu V0 *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 113",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 113",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 113", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 223",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 223",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 223", (* tau *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 333",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 333",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 333"  (* tau *)
};

(* String with all the HNL BRs *)
expbr[n_,mixe_,mixmu_,mixtau_]:=expBrAll[n,mixe,mixmu,mixtau];

(* BR FActor for HNL decay to tau *)
hnl2taubrfactor[n_,mixe_,mixmu_,mixtau_]:=
    enutautau[n,mixe]/totalw[n,mixe,mixmu,mixtau]
  + munutautau[n,mixmu]/totalw[n,mixe,mixmu,mixtau]
  + taunuee[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + taunumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + nutauee[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + nutaumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + width3a[n,mixe]/totalw[n,mixe,mixmu,mixtau]
  + width3a[n,mixmu]/totalw[n,mixe,mixmu,mixtau]
  + width3a[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + width3a[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + width3b[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + eD[n,mixe]/totalw[n,mixe,mixmu,mixtau]
  + tauPi[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + nupi0[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + nueta[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + nuetaprime[n,mixtau]/totalw[n,mixe,mixmu,mixtau]

  + nurho[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + uomega[n,mixtau]/totalw[n,mixe,mixmu,mixtau]
  + nuphi[n,mixtau]/totalw[n,mixe,mixmu,mixtau];

(*  HNL CONFIG  *)
hnlparams = {"2000000001:new",
"2000000001:isResonance = false",
"2000000001:name = HNL",
"2000000001:antiName = HNLbar", (* Esto lo define como Dirac *)
"2000000001:spinType = 2",
"2000000001:chargeType = 0",
"2000000001:colType = 0"}

hnlmass[n_]:={"2000000001:m0 = "<>ToString[AccountingForm[n,16]]}
(* HNL LIFETIME  *)
lifetime[n_,mixe_,mixmu_,mixtau_]:=Quiet[Re[1/totalw[n,mixe,mixmu,mixtau]]*1/tfactor*mmc];
(* PYTHIA8 string *)
explifetime[n_,mixe_,mixmu_,mixtau_]:=Quiet[{
	"2000000001:tau0 = "<>ToString[AccountingForm[Re[1/totalw[n,mixe,mixmu,mixtau]]*1/tfactor*mmc,16]]}]; 

(* Join ALL config for HNL (particleData + BRs) *)
hnldecay[n_,mixe_,mixmu_,mixtau_]:=Join[hnlparams,hnlmass[n],explifetime[n,mixe,mixmu,mixtau],expbr[n,mixe,mixmu,mixtau]];

(******************************************************************************)
(******************************** MESON DECAY *********************************)

(* leptonic decays *)

width2body[meson_,l_,n_,u_]:=If[n<meson[[2]]-l,
	Re[(gf^2*meson[[3]]^2*meson[[2]]*n^2)/(8*Pi)*meson[[4]]^2*u^2
		*(1-n^2/meson[[2]]^2+2l^2/meson[[2]]^2+l^2/n^2(1-l^2/meson[[2]]^2))
		*Sqrt[(1+n^2/meson[[2]]^2-l^2/meson[[2]]^2)^2-4*n^2/meson[[2]]^2]],
	0];

(* Semileptonic Scalar Bondarenko (pp 26, 38) *)
(* 	[411]D+ -> [-311]K0
		[421]D0 -> [-321]K-
		[431]Ds -> [221]eta
		[431]Ds -> [331]eta'	*)

z[meson1_,meson2_,q2_]:=Module[{tp,t0},
tp=(meson1[[2]]+meson2[[2]])^2;
t0=(meson1[[2]]+meson2[[2]])(Sqrt[meson1[[2]]]-Sqrt[meson2[[2]]])^2;
(Sqrt[tp-q2]-Sqrt[tp-t0])/(Sqrt[tp-q2]+Sqrt[tp-t0])];

fp[meson1_,meson2_,q2_]:=Module[
{f0, c, P},
If[(meson1[[5]]==411||meson1[[5]]==421)&& (meson2[[5]]==321||meson2[[5]]==311),
f0=0.7647; c=0.066; P=0.224];
If[(meson1[[5]]==411||meson1[[5]]==421)&& meson2[[5]]==211,
f0=0.6117; c=1.985; P=0.1314];
If[meson1[[5]]!=431,
1/(1-P*q2)*(f0-c(z[meson1,meson2,q2]-z[meson1,meson2,0])(1+(z[meson1,meson2,q2]+z[meson1,meson2,0])/2)),
0.495/((1-q2/2.1121^2)(1-(0.198*q2)/2.1121^2))]
];
f0[meson1_,meson2_,q2_]:=Module[
{f0, c, P},
If[(Abs[meson1[[5]]]==411||Abs[meson1[[5]]]==421)&& (Abs[meson2[[5]]]==321||Abs[meson2[[5]]]==311),
f0=0.7647; c=2.084; P=0];
If[(Abs[meson1[[5]]]==411||Abs[meson1[[5]]]==421)&& Abs[meson2[[5]]]==211,
f0=0.6117; c=1.188; P=0.0342];
If[Abs[meson1[[5]]]!=431,
1/(1-P*q2)*(f0-c(z[meson1,meson2,q2]-z[meson1,meson2,0])(1+(z[meson1,meson2,q2]+z[meson1,meson2,0])/2)),
0.495/((1-q2/2.1121^2)(1-(0.198*q2)/2.1121^2))]
];

Lambda[meson1_,meson2_,l_,n_,xi_]:=Sqrt[lambda[1,(meson2[[2]]/meson1[[2]])^2,xi]]*Sqrt[lambda[xi,(n/meson1[[2]])^2,(l/meson1[[2]])^2]];
gminus[meson1_,meson2_,l_,n_,xi_]:=xi*((n/meson1[[2]])^2+(l/meson1[[2]])^2)-((n/meson1[[2]])^2-(l/meson1[[2]])^2)^2;

ip1[meson1_,meson2_,l_,n_]:=If[n<meson1[[2]]-meson2[[2]]-l,
NIntegrate[1/(3xi^3)*fp[meson1,meson2,xi*meson1[[2]]^2]^2*Lambda[meson1,meson2,l,n,xi]^3,{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
0];
ip2[meson1_,meson2_,l_,n_]:=If[n<meson1[[2]]-meson2[[2]]-l,
NIntegrate[1/(2xi^3)*fp[meson1,meson2,xi*meson1[[2]]^2]^2*Lambda[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi]*lambda[1,(meson2[[2]]/meson1[[2]])^2,xi],{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
0];
ip3[meson1_,meson2_,l_,n_]:=If[n<meson1[[2]]-meson2[[2]]-l,
NIntegrate[1/(2xi^3)*f0[meson1,meson2,meson1[[2]]^2*xi]^2*Lambda[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi]*(1-(meson2[[2]]/meson1[[2]])^2)^2,{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
0];

width3body[meson1_,meson2_,l_,n_,u_]:=Module[{ck},
If[meson2[[5]]==211,ck=1/Sqrt[2],ck=1];
If[n<meson1[[2]]-meson2[[2]]-l,
Re[(gf^2 meson1[[2]]^5)/(64Pi^3)*ck^2*ckm[[2,2]]^2 u^2 (ip1[meson1,meson2,l,n]+ip2[meson1,meson2,l,n]+ip3[meson1,meson2,l,n])],
0]
];


(* Semileptonic Vector Bondarenko (pp 26, 38) *)
(*	[411]D+ -> [-313]Kbar*(892)0
		[421]D0 -> [-323]K*(892)-
		[431]Ds -> [333]phi(1020)	*)

vhh[q2_]:=1.03/((1-q2/2.1121^2)(1-0.27*q2/2.1121^2));
ahh0[q2_]:=0.76/((1-q2/1.96827^2)(1-0.17*q2/2.1121^2));
ahh1[q2_]:=0.66/(1-0.30*q2/2.1121^2-0.20*q2^2/2.1121^4);
ahh2[q2_]:=0.49/(1-0.67*q2/2.1121^2-0.16*q2^2/2.1121^4);
(* 431->333 (2004) Aliev et al *)
ahh1phi[q2_]:=0.54/(1-1.57*q2/d431[[2]]^2+0.16*(q2/d431[[2]]^2)^2);
ahh2phi[q2_]:=0.57/(1-4.7*q2/d431[[2]]^2+5.96*(q2/d431[[2]]^2)^2);
ahh0phi[q2_]:=0.53/(1-0.64*q2/d431[[2]]^2+2.13*(q2/d431[[2]]^2)^2);
vhhphi[q2_]:=0.64/(1-2.81*q2/d431[[2]]^2+1.34*(q2/d431[[2]]^2)^2);

ghh[meson1_,meson2_,q2_]:=If[meson1[[5]]==431&&meson2[[5]]==333,
	vhhphi[q2]/(meson1[[2]]+meson2[[2]]),
	vhh[q2]/(meson1[[2]]+meson2[[2]])];
fhh[meson1_,meson2_,q2_]:=If[meson1[[5]]==431&&meson2[[5]]==333,
	ahh1phi[q2](meson1[[2]]+meson2[[2]]),
	ahh1[q2](meson1[[2]]+meson2[[2]])];
ahhp[meson1_,meson2_,q2_]:=If[meson1[[5]]==431&&meson2[[5]]==333,
	-(ahh2phi[q2]/(meson1[[2]]+meson2[[2]])),
	-(ahh2[q2]/(meson1[[2]]+meson2[[2]]))];
ahhm[meson1_,meson2_,q2_]:=If[meson1[[5]]==431&&meson2[[5]]==333,
	1/q2*(ahh0phi[q2]*2meson2[[2]]-fhh[meson1,meson2,q2]-(meson1[[2]]^2-meson2[[2]]^2)ahhp[meson1,meson2,q2]),
	1/q2*(ahh0[q2]*2meson2[[2]]-fhh[meson1,meson2,q2]-(meson1[[2]]^2-meson2[[2]]^2)ahhp[meson1,meson2,q2])];

F[meson1_,meson2_,l_,n_,xi_]:=(1-xi)^2-2(meson2[[2]]/meson1[[2]])^2 (1+xi)+(meson2[[2]]/meson1[[2]])^4;
gplus[meson1_,meson2_,l_,n_,xi_]:=xi((n/meson1[[2]])^2+(l/meson1[[2]])^2)+((n/meson1[[2]])^2-(l/meson1[[2]])^2)^2

ivg2[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		1/3 meson1[[2]]^2 (meson2[[2]]/meson1[[2]])^2 NIntegrate[1/xi^2 ghh[meson1,meson2,meson1[[2]]^2*xi]^2 Lambda[meson1,meson2,l,n,xi]*F[meson1,meson2,l,n,xi]*(2xi^2-gminus[meson1,meson2,l,n,xi]),{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivf2[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		1/(24meson1[[2]]^2) NIntegrate[1/xi^3 fhh[meson1,meson2,meson1[[2]]^2*xi]^2*Lambda[meson1,meson2,l,n,xi]*(3F[meson1,meson2,l,n,xi](xi^2-((l/meson1[[2]])^2-(n/meson1[[2]])^2)^2)-Lambda[meson1,meson2,l,n,xi]^2+12(meson2[[2]]/meson1[[2]])^2*xi*(2xi^2-gplus[meson1,meson2,l,n,xi])),{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivap2[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		meson1[[2]]^2/24 NIntegrate[1/xi^3 ahhp[meson1,meson2,meson1[[2]]^2*xi]^2*Lambda[meson1,meson2,l,n,xi]*F[meson1,meson2,l,n,xi]*(F[meson1,meson2,l,n,xi](2xi^2-gplus[meson1,meson2,l,n,xi])+3gminus[meson1,meson2,l,n,xi](1-(meson2[[2]]/meson1[[2]])^2)^2),{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivam2[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		meson1[[2]]^2/8 NIntegrate[1/xi ahhm[meson1,meson2,meson1[[2]]^2*xi]^2 Lambda[meson1,meson2,l,n,xi]*F[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi],{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivfap[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		1/12 NIntegrate[1/xi^3 fhh[meson1,meson2,meson1[[2]]^2*xi]*ahhp[meson1,meson2,meson1[[2]]^2*xi]*Lambda[meson1,meson2,l,n,xi]*(3*xi*F[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi]+(1-xi-(meson2[[2]]/meson1[[2]])^2)(3F[meson1,meson2,l,n,xi](xi^2-((l/meson1[[2]])^2-(n/meson1[[2]])^2)^2)-Lambda[meson1,meson2,l,n,xi]^2)),{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivfam[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		1/4 NIntegrate[1/xi^2 fhh[meson1,meson2,meson1[[2]]^2*xi]*ahhm[meson1,meson2,meson1[[2]]^2*xi]*Lambda[meson1,meson2,l,n,xi]*F[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi],{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivapam[meson1_,meson2_,l_,n_]:=Module[{},
	If[n<meson1[[2]]-meson2[[2]]-l,
		meson1[[2]]^2/4 NIntegrate[1/xi^2 ahhp[meson1,meson2,meson1[[2]]^2*xi]*ahhm[meson1,meson2,meson1[[2]]^2*xi]*Lambda[meson1,meson2,l,n,xi]*F[meson1,meson2,l,n,xi]*gminus[meson1,meson2,l,n,xi]*(1-(meson2[[2]]/meson1[[2]])^2),{xi,(l/meson1[[2]]+n/meson1[[2]])^2,(1-meson2[[2]]/meson1[[2]])^2}],
		0]];
ivall[meson1_,meson2_,l_,n_]:=ivg2[meson1,meson2,l,n]+ivf2[meson1,meson2,l,n]+ivap2[meson1,meson2,l,n]+ivam2[meson1,meson2,l,n]+ivfap[meson1,meson2,l,n]+ivfam[meson1,meson2,l,n]+ivapam[meson1,meson2,l,n];

width3bodyv[meson1_,meson2_,l_,n_,u_]:=Module[{ck=1},
If[n<meson1[[2]]-meson2[[2]]-l,
Re[(gf^2 meson1[[2]]^7)/(64Pi^3 meson2[[2]]^2)*ck^2*ckm[[2,2]]^2 u^2 ivall[meson1,meson2,l,n]],
0]
];

d411k311bare[n_,u_]:=width3body[d411, k311, e, n, u];
d411k311barmu[n_,u_]:=width3body[d411, k311, mu, n, u];
d411k313bare[n_,u_]:=width3bodyv[d411, k313, e, n, u];
d411k313barmu[n_,u_]:=width3bodyv[d411, k313, mu, n, u];

d431eta221e[n_,u_]:=width3body[d431, eta221, e, n, u];
d431eta221mu[n_,u_]:=width3body[d431, eta221, mu, n, u];
d431eta331e[n_,u_]:=width3body[d431, eta331, e, n, u];
d431eta331mu[n_,u_]:=width3body[d431, eta331, mu, n, u];
d431phi333e[n_,u_]:=width3bodyv[d431, phi1020, e, n, u];
d431phi333mu[n_,u_]:=width3bodyv[d431, phi1020, mu, n, u];

d421k321bare[n_,u_]:=width3body[d421, k321, e, n, u];
d421k321barmu[n_,u_]:=width3body[d421, k321, mu, n, u];
d421k323bare[n_,u_]:=width3bodyv[d421, k323, e, n, u];
d421k323barmu[n_,u_]:=width3bodyv[d421, k323, mu, n, u];

newwidth2hnl[meson_,n_,mixe_,mixmu_,mixtau_]:=
Switch[meson[[5]],
411,
	Re[
	width2body[meson,e,n,mixe]
	+width2body[meson,mu,n,mixmu]
	+width2body[meson,tau,n,mixtau]
	+d411k311bare[n,mixe]
	+d411k311barmu[n,mixmu]
	+d411k313bare[n,mixe]
	+d411k313barmu[n,mixmu]],
431,
	Re[
	width2body[meson,e,n,mixe]
	+width2body[meson,mu,n,mixmu]
	+width2body[meson,tau,n,mixtau]
	+d431eta221e[n,mixe]
	+d431eta221mu[n,mixmu]
	+d431eta331e[n,mixe]
	+d431eta331mu[n,mixmu]
	+d431phi333e[n,mixe]
	+d431phi333mu[n,mixmu]],
421,
	Re[
	d421k321bare[n,mixe]
	+d421k321barmu[n,mixmu]
	+d421k323bare[n,mixe]
	+d421k323barmu[n,mixmu]
	]
];

newwidth[meson_,n_,mixe_,mixmu_,mixtau_]:=
Re[1/meson[[1]]+newwidth2hnl[meson,n,mixe,mixmu,mixtau]];

(* Meson main Parameters*)
mesonparam[meson_,n_,mixe_,mixmu_,mixtau_]:={
	ToString[meson[[5]]]<>":mayDecay = true",
	ToString[meson[[5]]]<>":tau0 = "<>ToString[AccountingForm[Re[1/newwidth[meson, n, mixe, mixmu, mixtau]]*1/tfactor*mmc,16]],
	ToString[meson[[5]]]<>":m0 = "<>ToString[AccountingForm[meson[[2]],16]]};

(* Meson decays *)	

exp411decay[n_,mixe_,mixmu_,mixtau_]:={
	"411:oneChannel = 1 "<>ToString[AccountingForm[width2body[d411,e,n,mixe]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 0 -11 2000000001",
	"411:addChannel = 1 "<>ToString[AccountingForm[width2body[d411,mu,n,mixmu]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 0 -13 2000000001",
	"411:addChannel = 1 "<>ToString[AccountingForm[width2body[d411,tau,n,mixtau]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 0 -15 2000000001",
	"411:addChannel = 1 "<>ToString[AccountingForm[d411k311bare[n,mixe]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 -311",
	"411:addChannel = 1 "<>ToString[AccountingForm[d411k311barmu[n,mixmu]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 -311",
	"411:addChannel = 1 "<>ToString[AccountingForm[d411k313bare[n,mixe]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 -313",
	"411:addChannel = 1 "<>ToString[AccountingForm[d411k313barmu[n,mixmu]/newwidth2hnl[d411, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 -313"
};

exp431decay[n_,mixe_,mixmu_,mixtau_]:={
	"431:oneChannel = 1 "<>ToString[AccountingForm[width2body[d431,e,n,mixe]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 0 -11 2000000001",
	"431:addChannel = 1 "<>ToString[AccountingForm[width2body[d431,mu,n,mixmu]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 0 -13 2000000001",
	"431:addChannel = 1 "<>ToString[AccountingForm[width2body[d431,tau,n,mixtau]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 0 -15 2000000001",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431eta221e[n,mixe]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 221",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431eta221mu[n,mixmu]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 221",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431eta331e[n,mixe]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 331",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431eta331mu[n,mixmu]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 331",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431phi333e[n,mixe]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 333",
	"431:addChannel = 1 "<>ToString[AccountingForm[d431phi333mu[n,mixmu]/newwidth2hnl[d431, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 333"
};

exp421decay[n_,mixe_,mixmu_,mixtau_]:={
	"421:oneChannel = 1 "<>ToString[AccountingForm[d421k321bare[n,mixe]/newwidth2hnl[d421, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 -321",
	"421:addChannel = 1 "<>ToString[AccountingForm[d421k321barmu[n,mixmu]/newwidth2hnl[d421, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 -321",
	"421:addChannel = 1 "<>ToString[AccountingForm[d421k323bare[n,mixe]/newwidth2hnl[d421, n, mixe, mixmu, mixtau],16]]<>" 22 -11 2000000001 -323",
	"421:addChannel = 1 "<>ToString[AccountingForm[d421k323barmu[n,mixmu]/newwidth2hnl[d421, n, mixe, mixmu, mixtau],16]]<>" 22 -13 2000000001 -323"
};

expmesonbrfactor[meson_,n_,mixe_,mixmu_,mixtau_]:={
	"# "<>ToString[meson[[5]]]<>" BR factor: "<>ToString[AccountingForm[newwidth2hnl[meson, n, mixe, mixmu, mixtau]/newwidth[meson, n, mixe, mixmu, mixtau],16]]}; 

mesonbrfactor[meson_,n_,mixe_,mixmu_,mixtau_]:={
  ToString[n]<>" "<>ToString[AccountingForm[newwidth2hnl[meson, n, mixe, mixmu, mixtau]/newwidth[meson, n, mixe, mixmu, mixtau],16]]};   


(* EXPORTA PYTHIA CONFIG FILES PARA "nfiles" MASAS CON <0,1.9>.  *)

Print["Todos los par\[AAcute]metros y funciones definidos correctamente."]

names={};
brfactors431={};
brfactors411={};
brfactors421={};

Do[
  name=ToString[n]<>".dat";
  Print["Generando "<>name<>" ..."];
  AppendTo[names,name];
  AppendTo[brfactors431,mesonbrfactor[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]]];
  AppendTo[brfactors411,mesonbrfactor[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]]];
  AppendTo[brfactors421,mesonbrfactor[d421,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]]];
  Quiet[Export["./mathbrs/"<>name,
    Join[
      {"#"<>ToString[n]},
      expmesonbrfactor[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		  expmesonbrfactor[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		  expmesonbrfactor[d421,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		  hnldecay[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		  mesonparam[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
  		mesonparam[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
  		mesonparam[d421,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
  		exp411decay[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
  		exp431decay[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
  		exp421decay[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]]
  		]
    ]
  ];
,{n,mhnli,mhnlf,paso}];
AppendTo[names,""];
AppendTo[brfactors431,""];
AppendTo[brfactors411,""];
AppendTo[brfactors421,""];
Export["./mathbrs/confignames.dat",names];
Export["./mathbrs/brfactors431.dat",brfactors431];
Export["./mathbrs/brfactors431bar.dat",brfactors431];
Export["./mathbrs/brfactors411.dat",brfactors411];
Export["./mathbrs/brfactors411bar.dat",brfactors411];
Export["./mathbrs/brfactors421.dat",brfactors421];
Export["./mathbrs/brfactors421bar.dat",brfactors421];
Print["All done"]

