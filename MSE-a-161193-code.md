
    (* Code dump for https://mathematica.stackexchange.com/a/161193/ *)

    (*
     * Setup for `Y` and `U` equations
     *)

    {eqs, ics} = {(1/2) Y'[x]^2 == (1 - Log[Y[x]^2]) Y[x]^2, Y[0] == 1};
    dics = Reduce[{eqs /. x -> 0, ics, Y'[0] > 0}, {Y'[0]}];
    (** RESIDUALS (unused) **)
    res = eqs /. Equal -> Subtract;
    dres = D[eqs, x] /. Equal -> Subtract;
    (** SUBSTITUTION **)
    usub = Y -> Function[x, Y[x]] /. Solve[8 U[x] == (1 - 2 Log[Y[x]]), Y[x], Reals];
    ueqs = {eqs, ics} /. usub // Reduce[#, {U[x], U[0]}, Reals] &;
    dueqs = D[eqs, x] /. usub // Reduce[#, {U[x]}, Reals] & // Simplify[#, U'[x] != 0] &;
    duics = Solve[ueqs /. x -> 0, {U[0], U'[0]}] /. Rule -> Equal;

    (*
     * Code for figure 1 - Y contact manifold
     *)

    singPlus = NDSolveValue[{D[eqs, {x, 1}], Y[0] == Sqrt@E, Y'[0] == 0},
       Y, {x, -15, 15},
       Method -> {"EquationSimplification" -> "Residual"}];
    singMinus = NDSolveValue[{D[eqs, {x, 1}], Y[0] == -Sqrt@E, Y'[0] == 0},
       Y, {x, -15, 15},
       Method -> {"EquationSimplification" -> "Residual"}];
    singZero = NDSolveValue[
       {Simplify[D[eqs /. Equal -> Subtract, x]/Y'[x]] == 0 /.
          Log[_] :> Piecewise[{{Y[x] Log[Y[x]^2], Y[x] != 0}}]/Y[x] //
         Expand,
        Y[0] == 0, Y'[0] == 0},
       Y, {x, -15, 15}];

    Clear[phSol];
    With[{x1 = -15., x2 = 15.},
      phSol[x0_, branch_: 1] :=
       First@NDSolve[{D[eqs, {x, 2}], Y[x0] == branch*Sqrt@E, Y'[x0] == 0,
           Y''[x0] == -2 Sqrt[E], WhenEvent[Y[x] == 0, "StopIntegration"],
           WhenEvent[Y'[x] == 0, "StopIntegration"]},
         Y, {x, x1, x2},
         Method -> {"EquationSimplification" -> "Residual"},
         "ExtrapolationHandler" -> {Indeterminate &,
           "WarningMessage" -> False}];
      phPlus = Show@Table[
         ParametricPlot3D[{t, Y[t], Y'[t]} /. phSol[xt] // Evaluate, {t,
           x1, x2}, PlotStyle -> Gray],
         {xt, x1 - 1, x2 + 1, 1.9}];
      phMinus = Show@Table[
         ParametricPlot3D[{t, -Y[t], -Y'[t]} /. phSol[xt] // Evaluate, {t,
            x1, x2}, PlotStyle -> Gray],
         {xt, x1 - 1, x2 + 1, 1.9}]
      ];

    cplot = Show[
       phPlus, phMinus,
       PlotRange -> All];

    cmanifold = With[{cm = eqs /. {Y'[x] -> p, Y[x] -> y}},
       ContourPlot3D[
        cm,
        {x, -15, 15}, {y, -1.7, 1.7}, {p, -2.7, 2.7},
        AxesLabel -> Automatic, BoxRatios -> {2, 1, 1},
        ContourStyle -> {Opacity[0.35], LightGreen}, Mesh -> None,
        BoundaryStyle -> None,
        PlotRange -> {15 {-1, 1}, 2.7 {-1, 1}, 2.7 {-1, 1}},
        PlotLabel -> cm]
       ];

    ys1 = NDSolve[{D[eqs, x], Y[0] == 1, Y'[0]^2 == 2}, Y, {x, -15, 15},
       AccuracyGoal -> 14];

    With[{pF =
        Function[y,
         Piecewise[{{Sqrt[2] Sqrt[y^2 - y^2 Log[y^2]], y != 0}}]],
       x1 = -15, x2 = 15},
      cpP2Y = Show[
        cmanifold,
        ParametricPlot3D[{{x, Sqrt@E, 0}, {x, -Sqrt@E, 0}}, {x, x1, x2},
         PlotStyle -> {Directive[Red, Thickness[Medium]]}],
        cplot,
        cplot /. {x_Real, y_Real, p_Real} :> {x, y, -2.7},
        ParametricPlot3D[
         {t, Y[t], Y'[t]} /.
           Join[ys1, {{Y -> singMinus}, {Y -> singPlus}, {Y -> singZero}}] // Evaluate,
         {t, -15, 15},
         PlotStyle -> {Directive[Purple, Tube[0.08]],
           Directive[Orange, Tube[0.08]], Directive[Red, Tube[0.05]],
           Directive[Red, Tube[0.05]], Directive[Green, Tube[0.05]]}],
        Lighting -> "Neutral", AxesLabel -> {x, Y, p == Y'}
        ]
      ];

    Show[cpP2Y, ViewPoint -> {-5, -3, 1.5}]

    (*
     * Code for Figure 2 - Y projection
     *)
    qcmanifold =
      ContourPlot3D[
       Evaluate[
        NestList[D[#, x] &, eqs, 1] /. {Y[x] -> y, Y'[x] -> p,
           Y''[x] -> q} // Simplify[#, p != 0] &],
       {y, -1.71, 1.7}, {p, -1.7, 1.701}, {q, -3.5, 3.5},
       ContourStyle -> Opacity[0.5], Mesh -> None, AxesLabel -> Automatic];

    s1 = Quiet@
      NDSolveValue[{D[eqs, x], Y[0] == 1, Y'[0] == Sqrt[2]},
       Y, {x, -50, 50},
       AccuracyGoal -> Infinity,
       Method -> {"Projection", "Invariants" -> {eqs}}]
    (*s2=NDSolveValue[{D[eqs,x],Y[0]==1,Y'[0]==-Sqrt[2]},Y,{x,-15,15}];*)

    qcpP2Y = GraphicsRow[{
        Show[
         qcmanifold,
         ParametricPlot3D[
          {Y[t], Y'[t], Y''[t]} /. {{Y -> s1}} // Evaluate,
          {t, -15, 15}, PlotStyle -> {Red, Tube[0.08]}
          ],
         Lighting -> "Neutral", ViewPoint -> {-1, -4, 2.5}
         ], ListLinePlot[s1, PlotStyle -> Red, PlotRange -> 1.72, InterpolationOrder -> 3]}];

    (*
     * Code for figure 3 - AccuracyGoal
     *)

    ClearAll[stepValuePlot];
    stepValuePlot[s_, scale_: 100, opts___?OptionQ] /; NumericQ[scale] :=
      ListLinePlot[{scale * RealExponent@ Differences@ Flatten@ s["Grid"],
        RealExponent@s["ValuesOnGrid"]},
       opts,
       PlotLabel -> s["Domain"],
       PlotLegends -> {HoldForm[scale*Log10["\[CapitalDelta]x"]],
          HoldForm[Log10[Y]]}];

    ClearAll[agPlot];
    agPlot[ag_, scale_: 100, opts___?OptionQ] /; NumericQ[scale] :=
      With[{s1 = NDSolveValue[{D[eqs, x], dics}, Y, {x, -25, 25},
          AccuracyGoal -> ag,
          Method -> {"Projection", "Invariants" -> {eqs}}]},
       Sow[s1, "Solution"];
       Row[{
         Plot[s1[x], Evaluate@Flatten@{x, s1@"Domain"}, opts,
          PlotStyle -> {AbsoluteThickness[2.5], ColorData[97][2]},
          ImageSize -> {Automatic, 150}, PlotLabel -> s1["Domain"],
          PlotRange -> {{-25., 25.}, {-1.72, 1.72}}],
         stepValuePlot[s1, scale, opts, ImageSize -> {Automatic, 150},
          GridLines -> {None, Cases[{-ag, -162}, _?NumericQ]},
          PlotRange -> {Max@{-40 - ag, -220}, 5}]
         }, Spacer[1]]
       ];

    PrintTemporary@Dynamic@{a, MemoryInUse[], Clock[]};
    movie = Table[
       Labeled[
        MapThread[
          Pane,
          {(*{Round[a,0.1]}~Join~*)
           First@Quiet@
             agPlot[a, If[TrueQ[a >= 100], 100, 10],
              ImageSize -> {Automatic, 140}],
           {(*40,*)200, 350}}
          ] // Row,
        AccuracyGoal -> a, Top],
       {a, 8, 200, 1}];

    Export["foo.gif", movie[[1 ;; ;; 2]], "DisplayDurations" -> 0.2]

    Manipulate[
     movie[[n]],
     {n, 1, Length@movie, 1}]

    (*
     * Code for figure 4 - U contact manifold
     *)

    Clear[uphSol];
    With[{x1 = -21., x2 = 21.},
      uphSol[x0_] :=
       First@NDSolve[{dueqs, U[x0] == 0, U'[x0] == 0,
           WhenEvent[U[x] == 0, "StopIntegration"],
          WhenEvent[U'[x] == 0, "StopIntegration"]},
         U, {x, x1, x2},
         Method -> {"Projection", "Invariants" -> {First@ueqs}},
         Method -> {"EquationSimplification" -> "Residual"},
         "ExtrapolationHandler" -> {Indeterminate &,
           "WarningMessage" -> False}];
      ucplot = Show[Table[
         ParametricPlot3D[{t, U[t], U'[t]} /. uphSol[xt] // Evaluate, {t,
           x1; -15, x2; 15}, PlotStyle -> Gray],
         {xt, x1, x2, 2.}],
        PlotRange -> {{-15, 15}, {-0.5, 10}, {-5, 5}}]
      ];

    Show[
     ContourPlot3D[
      {First@ueqs} /. {U[x] -> u, U'[x] -> p, U''[x] -> q} //
       Evaluate,
      {x, -15, 15}, {u, -0.5, 7.3}, {p, -2.7, 2.6}(*,{q,-25,1}*),
      ContourStyle -> Opacity[0.5], Mesh -> None,
      AxesLabel -> Automatic],
     ParametricPlot3D[{x, 0, 0}, {x, -15, 15},
      PlotStyle -> {Directive[Red, Tube[0.05]]}],
     ucplot /. {x_Real, y_Real, p_Real} :> {x, y, -2.7},
     ucplot,
     ParametricPlot3D[
      {x, U[x], U'[x]} /. usol // Evaluate,
      {x, -15, 15}, PlotStyle -> {Purple, Tube[0.08]}],
     ViewPoint -> {1, -5, 2}, PlotRangeClipping -> True,
     PlotRangePadding -> None
     ]

     (*
      * Code for figure 5 - U projection
      *)

    usol1 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15},  WorkingPrecision -> 16];
    usol2 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15}, WorkingPrecision -> 24];
    usol3 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15}, WorkingPrecision -> 32];
    usol = NDSolve[{dueqs, Last@duics}, U, {x, -15, 15},
       Method -> {"Projection", "Invariants" -> {First@ueqs}}];

    Show[
      surf12 = {First@ueqs, dueqs} /. {U[x] -> u, U'[x] -> p, U''[x] -> q};
      ContourPlot3D[
       surf12 // Evaluate,
       {u, -0.5, 10}, {p, -5, 5}, {q, -22, 2},
       ContourStyle -> Opacity[0.5], Mesh -> None, AxesLabel -> Automatic,
        PlotLegends -> surf12(*{"Order 1","Order 2"}*)],
      ParametricPlot3D[
       {U[x], U'[x], U''[x]} /. Join[usol1, usol2, usol3] // Evaluate,
       {x, -15, 15},
       PlotStyle ->
        Table[Directive[ColorData[97][k + 2], Tube[0.08]], {k, 3}],
       PlotLegends -> HoldForm /@ Thread[WP -> {16, 24, 32}]],
      ViewPoint -> {4, -1, 1}
      ] /. Legended[x_, leg_] :> Legended[x, leg /. Tube[_] :> {}]

    Show[
     ContourPlot3D[
      {First@ueqs, dueqs} /. {U[x] -> u, U'[x] -> p, U''[x] -> q} //
       Evaluate,
      {u, -0.5, 10}, {p, -5, 5}, {q, -25, 1},
      ContourStyle -> Opacity[0.5], Mesh -> None,
      AxesLabel -> Automatic],
     ParametricPlot3D[
      {U[x], U'[x], U''[x]} /. usol // Evaluate,
      {x, -15, 15},
      PlotStyle -> {Blend[{Yellow, Red}, 0.2], Tube[0.08]}],
     ViewPoint -> {4, -3, 2}
     ]

     (*
      * Code for figure 6 - U WorkingPrecision
      *)

    usol1 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15}, WorkingPrecision -> 16];
    usol2 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15}, WorkingPrecision -> 24];
    usol3 = NDSolve[{dueqs, Last@duics} // N[#, 55] &, U, {x, -15, 15}, WorkingPrecision -> 32];
    ListLinePlot[U /. Join[usol1, usol2, usol3], PlotRange -> 7,
     PlotLegends -> HoldForm /@ Thread[WP -> {16, 24, 32}]]

    usol = NDSolve[{dueqs, Last@duics}, U, {x, -15, 15},
       Method -> {"Projection", "Invariants" -> {First@ueqs}}];
    ListLinePlot[U /. usol, PlotRange -> 7]

    ListLinePlot[U /. Join[usol1, usol2, usol3, usol],
     PlotStyle ->
      Join[Table[
        ColorData[97][k + 2], {k, 3}], {Blend[{Yellow, Red}, 0.2]}],
     PlotRange -> 7,
     PlotLegends ->
      Append[HoldForm /@ Thread[WP -> {16, 24, 32}], "Projection"]]

    (*
     * Code for figure 7 - Y from U
     *)

    ysol = Join[
         U["Grid"] /. usol, (* usol is above *)
         NestList[D[#, x] &, Y[x] /. usub, 2] /. First[usol] /. x -> "ValuesOnGrid"
         ] // Transpose // Interpolation;
    Plot[ysol[x], {x, -15, 15}, PlotRange -> All]

     (*
      * Ununsed
      *)

    (* Phase portrait of reduced 2nd-order U equation (cf. fig. 5) *)
    StreamPlot[{p, 1 + 8 p^2 - 8 u}, {u, 0, 5}, {p, -3, 3}]

    (*
     * Utility
     *
     * To plot the surfaces, the following would be convenient to convert
     * the differential equations to equations with ordinary variables,
     * but copy-paste editing meant I never got around to using it. :)
     *)

     toVars[Y_[X_], vars_: {y, p, q, x}] :=
       NestList[D[#, X] &, Y[X], Length@vars - 2] ~Append~ X -> vars // Thread;

     (1/2) Y'[x]^2 == (1 - Log[Y[x]^2]) Y[x]^2 /. toVars[Y[x]]
     toVars[U[x]]
     (*
       p^2/2 == y^2 (1 - Log[y^2])
       {U[x] -> y, U'[x] -> p, U''[x] -> q, x -> x}
     *)
