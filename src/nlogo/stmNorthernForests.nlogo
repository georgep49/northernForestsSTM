__includes [
  "setup.nls"
  "construct-lsp.nls"
  "veg-dynamics.nls"
  "pathogens.nls"
  "fire-routines.nls"
  "climate.nls"
  "list-utils.nls"
  "dist-utils.nls"
  "catena.nls"
]

extensions [
  profiler ;; in case we need to figure out why things are slow
  matrix   ;; matrix maths
  palette  ;; nicer colours
  table    ;; dictionary like tables for storing lots of the input parameters
  csv      ;; easy reading of CSVs
  gis      ;; load FBM Surfaces
  ;;sr       ;; because...
]

;; As per the Ecosystems paper (Perry et al. 2015)
;; classes: 0 = invaded, 1 = manuka, 2 = kanuka, 3 = yng forest (mid-succ), 4 = mature forest (late-succ)

;; New names for STM for northern forests under pathogens, fire. etc.:
;; classes "gr" (0), "d-sh" (1), m-sh (2), "ksh-k" (3) "ksh-nok" (5)  "yf-k" (4)  "yf-nok" (6) "old" (7) "ksh-p" (8) "yf-p" (9) "old-p" (10)

;; set and forget globals...
globals [

  ;; run control aux
  starting-seed
  from-R?
  run-id

  ;; model world structure
  world-size
  patch-grain
  n-states
  record-tag
  fbm-code

  ;; forest composition etc.
  class-init-abund
  class-list
  class-names-list
  forest-classes-list
  forest-idx-list
  forest-gully-cover
  ;;base-changes


  ; succession related
  succ-tpi-wgt
  succ-aspect-wgt
  base-changes-dict
  tr-mtx-bank-yfor
  tr-mtx-bank-ofor

  ;; fire-related
  flammability-dict
  flamm-wind-wgt
  flamm-slope-wgt
  flamm-slope-breaks
  flamm-aspect-wgt
  enso-freq-wgt-list

  fire-front
  fire-size
  fire-year
  rnd-seed
  extinguished?

  ;; pathogen-related
  rust-susc-dict
  phy-susc-dict

  ;; fire history
  ;; two lists to deal with individual fires and their characteristics
  ;; third list keeps track of fire sizes for aggregate stats at end
  fire-stats
  fire-record
  fire-size-list

  ;; climate related
  enso-mtx
  enso-list
  enso-state
  enso-record

  ;; auxillary (Secondary) parameters
  aux-params-dict

  ;; lsp reporters
  abundances
  beyond-flamm-time
  old-growth-abund
  sum-stalled       ;; total stalled list 0 = 2 -> 3, 1 = 3 -> 4
  sum-changes       ;; no. stalled list 0 = 2 -> 3, 1 = 3 -> 4
  n-changes

  ;; these globals control the MRC algorithm
  colour-list        ;; list of spectral colors
  colour-dict
  num-classes       ;; the number of suitability classes
  class-raw-FD      ;; the raw freq distribution of classes requested
  class-FD          ;; the modified/scaled (used!) freq distribution of classes requested
  class-CFD         ;; the modified/scaled (used!) CFD of classes requested
  clusters          ;; a list of patch-sets of patches in each cluster

  ;;  if we need to add gully forest
  forestFront
  forestArea

  ;; topo
  nhb-scalar  ;; this is for the topographic metrics
]


patches-own
[
  fire-history      ;; list of years in which cell burned
  last-change
  next-change
  times-change
  t-colonised

  regenbank-yfor
  regenbank-ofor
  ;ldd-this-tick

  ;; topo and hydrological params
  flow-accum
  flow-to
  aspect
  aspect-arc
  edaphic-grad
  elevation
  slope
  slope-horn
  TPI
  TWI
  hs

  ;; lsp context
  distance-to-coast
  farm-node?
  farm?

  ;; pathogen status
  myrtle-rust?
  myrtle-rust-time

  kauri-mate?
  kauri-mate-time
  kauri-mate-nhb

  ;; fire history
  stalled
  flammability
  burned?

  ;; MRC stuff
  class             ;; used to label initial clusters and to index veg class
  prev-class        ;; class previously occupied

  cluster-leader    ;; patch which 'leads' the cluster (see model lib example)
  in-perc-cluster?  ;; true if patch is in initial percolation cluster
  cluster-id

  ;; topo helpers
  gully-forest?
  valley?
]


to go
  ;while [ticks <= 5]
  ;[
    ;profiler:start         ;; start profiling
    ;if ticks = 0 [print date-and-time]

   ;while [ ticks < max-ticks and old-growth-abund <= max-forest ]
   ;[
   ;; Get the climate conditions for the year
   transition-enso
   let extrinsic random-normal 0 extrinsic-sd

   ;; fire dynamics
   let fire-this-tick false
   if ticks > burn-in-regen
   [
     let enso-wgt 1
     ;if enso-state = "ENL" or enso-state = "EN" [ set enso-wgt enso-freq-wgt ]
     set enso-wgt item (position enso-state enso-list) enso-freq-wgt-list

     if random-float 1 <= (fire-frequency * enso-wgt * (1 + extrinsic))
     [
       let start? ignite-fire extrinsic
       if start? = true
       [
         fire-spread extrinsic
         post-fire
         set fire-this-tick true
       ]
     ]
   ]

   set n-changes 0

   ;; dispersal and regeneration bank dynamics
   dispersal
   regenerate-patch-bank
   if sap-herbivory > 0 [ herbivory-patch-bank ]

   ;; succesional dynamics
   if ticks >= burn-in-regen
   [
     succession

     if n-changes > 0 or fire-this-tick = true
     [
       colour-by-class
       update-abundances
       if abund-flammable < 0.3 and ticks < beyond-flamm-time [ set beyond-flamm-time ticks ]
     ]
   ]

  ;; pathogen dynamics
  if rust-global-inf > 0 [spread-rust]

  if phy-global-inf > 0
  [
    spread-kauri-mate
    kauri-mate-switch
  ]

   ; profiler:stop          ;; stop profiling
   ; print profiler:report  ;; view the results
   ; profiler:reset

   update-abundances

   tick

   ;]

end

;; Various helper functions...
to update-abundances
  let n [class] of patches
  set abundances (map [ i -> occurrences i n ] class-names-list)
  set old-growth-abund count patches with [class = "old-p" or class = "old-f"] / world-size
end


to-report occurrences [x the-list]
  report reduce
    [ [i j] -> ifelse-value (j = x) [i + 1] [i] ] (fput 0 the-list)
end

;; Colours patches by class
to colour-by-class
  ask patches [
    set pcolor table:get colour-dict class
  ]
end

to color-by-stall
  let mx max [stalled] of patches
  ask patches with [ stalled > 0 ] [
    set pcolor scale-color orange stalled 1 mx + 5
  ]
end

to colour-by-lastfire

  if length fire-record > 0  ; check to make sure there have been some fires
  [
    let burned-patches patches with [not empty? fire-history]
    let mx max last [fire-history] of burned-patches

    ask patches
    [
      ifelse not empty? fire-history
      [
        set pcolor scale-color red (last fire-history) mx 1
      ]
      [
        set pcolor grey
      ]
    ]
  ]

end

to color-by-regen3
  ask patches [
    set pcolor scale-color blue (sum regenbank-yfor) 100 0
  ]
end

to color-by-regen4
  ask patches [
    set pcolor scale-color blue (sum regenbank-ofor) 100 0
  ]
end

to color-by-edaphic
ask patches
[
  set pcolor palette:scale-gradient [[236 226 240] [28 144 153]] edaphic-grad 1 0
]
end

;; Abundance flammable
to-report abund-flammable
  report count patches with [flammability > 0.25] / world-size
end

;; Abundance stalled
to-report abund-stalled
  report count patches with [stalled > 0] / world-size
end

;; Mean flammability (lsp)
to-report mean-lsp-flammable
  report mean [flammability] of patches
end

;; b1 upper, b0, lower, b2 steep
to-report decr-limit-fx [b0 b1 b2 x]
  report b1 *  exp(- b2 * x) + b0
end

;; b1 upper, b1 steep
to-report incr-limit-fx [b1 b2 x]
  report b1 * (1 - exp(- b2 * x))
end
@#$#@#$#@
GRAPHICS-WINDOW
241
10
761
531
-1
-1
2.0
1
10
1
1
1
0
0
0
1
0
255
0
255
1
1
1
Year
30.0

SLIDER
30
54
202
87
perc-seed
perc-seed
0
1
0.57
0.01
1
NIL
HORIZONTAL

BUTTON
17
13
80
46
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
84
12
147
45
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
151
12
214
45
step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
806
10
1154
224
Abundances
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"invaded" 1.0 0 -2570826 true "" "plot count patches with [class = \"d-sh\" or class = \"gr\"] / world-size"
"manuka" 1.0 0 -8431303 true "" "plot count patches with [class = \"m-sh\"] / world-size"
"kanuka" 1.0 0 -5509967 true "" "plot count patches with [class = \"ksh-p\" or class = \"ksh-k\" or class = \"ksh-nok\"] / world-size"
"yngfor" 1.0 0 -13840069 true "" "plot count patches with [class = \"yf-p\" or class = \"yf-k\" or class = \"yf-nok\"] / world-size"
"oldfor" 1.0 0 -15575016 true "" "plot count patches with [class = \"old-p\" or class = \"old-f\"] / world-size"

SLIDER
31
93
203
126
fire-frequency
fire-frequency
0
1
0.0
.01
1
NIL
HORIZONTAL

SLIDER
30
130
202
163
fraction-seed-ldd
fraction-seed-ldd
0
1
0.15
.01
1
NIL
HORIZONTAL

SLIDER
31
166
203
199
seed-pred
seed-pred
0
1
0.0
0.01
1
NIL
HORIZONTAL

SLIDER
31
241
203
274
flamm-start
flamm-start
0
1
0.5
0.01
1
NIL
HORIZONTAL

SWITCH
46
388
175
421
track-stalled?
track-stalled?
1
1
-1000

PLOT
806
227
1155
377
Lsp flammability
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"lspflm" 1.0 0 -16777216 true "" "plot mean-lsp-flammable"
"propflm" 1.0 0 -2674135 true "" "plot abund-flammable"
"propstl" 1.0 0 -5825686 true "" "plot abund-stalled"

SLIDER
31
200
203
233
mean-ldd
mean-ldd
0
10
5.0
0.5
1
NIL
HORIZONTAL

TEXTBOX
25
423
225
447
frac-consumed estimated from Dijkstra thesis\nfrac-predated estimated from C & A and W & K
9
0.0
1

SLIDER
32
277
204
310
base-invasion
base-invasion
0
1
0.05
.01
1
NIL
HORIZONTAL

SLIDER
32
312
204
345
fire-slow
fire-slow
0
10
2.0
1
1
NIL
HORIZONTAL

SLIDER
33
348
205
381
fire-invasion
fire-invasion
0
1
0.1
0.01
1
NIL
HORIZONTAL

TEXTBOX
3
254
34
272
FIRE
11
15.0
1

TEXTBOX
5
65
28
83
RGN
11
15.0
1

SWITCH
808
385
921
418
invasion?
invasion?
1
1
-1000

INPUTBOX
1032
388
1307
448
init-composition-file
parameter_files/initial_shrub_composition.dat
1
0
String

TEXTBOX
208
290
240
308
0.05
11
105.0
1

TEXTBOX
207
361
238
379
0.10
11
105.0
1

SLIDER
30
449
202
482
base-seed-prod-of
base-seed-prod-of
0
20
8.0
1
1
NIL
HORIZONTAL

SLIDER
29
486
201
519
base-seed-prod-yf
base-seed-prod-yf
0
20
4.0
1
1
NIL
HORIZONTAL

SLIDER
21
589
193
622
crit-density-yng
crit-density-yng
1
20
12.0
1
1
NIL
HORIZONTAL

SLIDER
21
628
193
661
crit-density-old
crit-density-old
1
20
14.0
1
1
NIL
HORIZONTAL

TEXTBOX
39
550
189
583
Minimum 'density' of juvenile saps to make transition - this is effectively an index of dispersal.
9
0.0
1

BUTTON
928
386
1026
419
show-stalled
color-by-stall
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
928
422
1028
455
show-class
colour-by-class
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
198
628
370
661
burn-in-regen
burn-in-regen
0
50
15.0
5
1
NIL
HORIZONTAL

BUTTON
807
461
927
494
show-r3
color-by-regen3
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
808
498
928
531
show-r4
color-by-regen4
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
808
421
918
454
write-record?
write-record?
1
1
-1000

PLOT
1158
11
1522
221
Fire size history
Year
Fire size
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "(foreach fire-year fire-size-list [ [y f ] -> plotxy y f] )"

TEXTBOX
242
534
648
576
Chocolate - manuka\nBlues - no kauri, Oranges - kauri, Turquoises - pohutakawa\nDark green - old-growth\n
11
0.0
1

BUTTON
936
462
1077
495
highlight-kauri-mate
ask patches with [kauri-mate?] [ set pcolor red ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
936
497
1077
530
highlight-rust
ask patches with [myrtle-rust?] [ set pcolor yellow ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
1311
324
1483
357
extrinsic-sd
extrinsic-sd
0
.5
0.0
.01
1
NIL
HORIZONTAL

SLIDER
1314
449
1486
482
rust-global-inf
rust-global-inf
0
.3
0.0
.005
1
NIL
HORIZONTAL

SLIDER
1316
485
1488
518
phy-global-inf
phy-global-inf
0
.3
0.0
.005
1
NIL
HORIZONTAL

SLIDER
1317
522
1489
555
phy-local-inf
phy-local-inf
0
.3
0.0
.005
1
NIL
HORIZONTAL

SLIDER
1317
558
1489
591
phy-radius-inf
phy-radius-inf
0
3
1.5
0.5
1
NIL
HORIZONTAL

INPUTBOX
1310
256
1518
316
enso-matrix-file
parameter_files/enso_matrix.dat
1
0
String

SLIDER
197
591
369
624
sap-herbivory
sap-herbivory
0
1.0
0.0
.01
1
NIL
HORIZONTAL

SLIDER
1312
360
1484
393
enso-freq-wgt
enso-freq-wgt
0.75
1.25
0.9
.01
1
NIL
HORIZONTAL

SWITCH
954
578
1108
611
farm-edge?
farm-edge?
1
1
-1000

SLIDER
1115
578
1287
611
farm-edge-nodes
farm-edge-nodes
1
40
30.0
1
1
NIL
HORIZONTAL

SLIDER
1116
617
1288
650
mean-farm-depth
mean-farm-depth
1
40
10.0
1
1
NIL
HORIZONTAL

SWITCH
955
615
1109
648
farm-revegetate?
farm-revegetate?
1
1
-1000

BUTTON
811
536
953
569
highlight-fire-mosaic
colour-by-lastfire
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
960
535
1066
568
show edaphic
ask patches [ \n  set pcolor scale-color blue edaphic-grad 0 1]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
538
536
762
596
nlrx-info
NIL
1
0
String

BUTTON
1080
496
1194
529
show elevation
ask patches [ \n  set pcolor scale-color grey elevation 0 150]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
1070
535
1162
568
show lform
ask patches [\nif tpi = 1 [set pcolor blue]\nif tpi = 2 [set pcolor green]\nif tpi = 3 [set pcolor orange]\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
812
626
945
659
forest-gully-prop
forest-gully-prop
0
1
0.5
.01
1
NIL
HORIZONTAL

CHOOSER
809
574
947
619
terrain-type
terrain-type
"flat" "ridge-gully" "dtm"
1

MONITOR
634
601
763
646
seed
starting-seed
0
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -7500403 true true 135 285 195 285 270 90 30 90 105 285
Polygon -7500403 true true 270 90 225 15 180 90
Polygon -7500403 true true 30 90 75 15 120 90
Circle -1 true false 183 138 24
Circle -1 true false 93 138 24

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="fire-invasion-forest-start-OLD" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <postRun>write-fire-record</postRun>
    <timeLimit steps="300"/>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <metric>mean-lsp-flammable</metric>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-seed-ldd">
      <value value="0.15"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fire-frequency" first="0" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="mean-ldd">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perc-seed">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-revegetate?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition-file">
      <value value="&quot;parameter_files/initial_forest_composition.csv&quot;"/>
      <value value="&quot;parameter_files/initial_shrubland_composition.csv&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="baseline-OLD" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>fbm-code</metric>
    <metric>abundances</metric>
    <enumeratedValueSet variable="enso-freq-wgt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-revegetate?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="enso-matrix-file">
      <value value="&quot;parameter_files/enso_matrix.csv&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-farm-depth">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-yf">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-old">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rust-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="burn-in-regen">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-of">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-seed-ldd">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sap-herbivory">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-radius-inf">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-local-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-yng">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-forest">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="extrinsic-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition-file">
      <value value="&quot;parameter_files/initial_forest_composition.csv&quot;"/>
      <value value="&quot;parameter_files/initial_shrubland_composition.csv&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge-nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perc-seed">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="terrain-type">
      <value value="&quot;flat&quot;"/>
      <value value="&quot;ridge-gully&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="baseline-shrub-ridge" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <postRun>write-fire-record</postRun>
    <exitCondition>ticks &gt;= 300</exitCondition>
    <metric>starting-seed</metric>
    <metric>abundances</metric>
    <metric>abund-stalled</metric>
    <metric>mean [flammability] of patches</metric>
    <enumeratedValueSet variable="enso-freq-wgt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-revegetate?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="enso-matrix-file">
      <value value="&quot;parameter_files/enso_matrix.dat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-farm-depth">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-yf">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-old">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rust-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="burn-in-regen">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-of">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-seed-ldd">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sap-herbivory">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-radius-inf">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-local-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-yng">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="extrinsic-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition-file">
      <value value="&quot;parameter_files/initial_forest_composition.dat&quot;"/>
      <value value="&quot;parameter_files/initial_shrub_composition.dat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge-nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perc-seed">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <subExperiment>
      <enumeratedValueSet variable="terrain-type">
        <value value="&quot;flat&quot;"/>
      </enumeratedValueSet>
      <enumeratedValueSet variable="forest-gully-prop">
        <value value="0"/>
      </enumeratedValueSet>
    </subExperiment>
    <subExperiment>
      <enumeratedValueSet variable="terrain-type">
        <value value="&quot;ridge-gully&quot;"/>
      </enumeratedValueSet>
      <enumeratedValueSet variable="forest-gully-prop">
        <value value="0.5"/>
      </enumeratedValueSet>
    </subExperiment>
  </experiment>
  <experiment name="shrub-one" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <postRun>write-fire-record</postRun>
    <exitCondition>ticks &gt;= 300</exitCondition>
    <metric>abundances</metric>
    <metric>abund-stalled</metric>
    <metric>mean [flammability] of patches</metric>
    <enumeratedValueSet variable="enso-freq-wgt">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-revegetate?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="enso-matrix-file">
      <value value="&quot;parameter_files/enso_matrix.dat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-farm-depth">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-yf">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-old">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rust-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="burn-in-regen">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-of">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-seed-ldd">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sap-herbivory">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-radius-inf">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-local-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-yng">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-record?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="extrinsic-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition-file">
      <value value="&quot;parameter_files/initial_shrub_composition.dat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="farm-edge-nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="perc-seed">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phy-global-inf">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="terrain-type">
      <value value="&quot;flat&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="forest-gully-prop">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
