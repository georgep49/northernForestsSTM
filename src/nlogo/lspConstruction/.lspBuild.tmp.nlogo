__includes [
  "../list-utils.nls"
  "../dist-utils.nls"
]

extensions [
  csv
  table
]

globals [

  ;;
  forestFront
  forestArea

  ;; run control aux
  from-R?
  run-id

  ;; model world structure
  world-size
  patch-grain
  n-states
  record-tag
  fbm-code

  ;; forest composition etc.
  class-list
  class-names-list
  forest-classes-list
  ;;base-changes
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
  class-FD          ;; freq distribution of classes requested
  class-CFD         ;; cumulative freq dist of classes requested
  clusters          ;; a list of patch-sets of patches in each cluster

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
  ldd-this-tick

  ;; topo and hydrological params
  flow-accum
  flow-to
  aspect
  edaphic-grad
  elevation
  slope
  TPI
  TWI

  ;; lsp context
  distance-to-coast
  farm-node?
  farm?

  ;; patiogen status
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
  forest?
  f-shell
  valley?
]


to setup
  clear-all
  show "Building parameter dictionaries..."

  ; if from-r? != TRUE [set from-r? FALSE]

  set world-size count patches
  set patch-grain 20
  set n-states 11

  set record-tag ""

  set rnd-seed new-seed

  set colour-dict table:from-list [
    ["gr" 19]
    ["d-sh" 38]
    ["m-sh" 34]
    ["ksh-k" 27]
    ["yf-k" 24]
    ["ksh-nok" 107]
    ["yf-nok" 104]
    ["old-f" 63]
    ["ksh-p" 88]
    ["yf-p" 85]
    ["old-p" 83]
  ]

  set class-list n-values n-states [ i -> i ]
  set class-names-list (list "gr" "d-sh" "m-sh" "ksh-k" "yf-k" "ksh-nok" "yf-nok" "old-f" "ksh-p" "yf-p" "old-p")
  set forest-classes-list (list "yf-k" "yf-nok" "old-f" "yf-p" "old-p")

  ;; load dictionaries (NL tables) from csv files
  set flammability-dict table:from-list csv:from-file "parameter_files/flammability_table.dat"
  set base-changes-dict table:from-list csv:from-file "parameter_files/base_changes_table.dat"
  set rust-susc-dict table:from-list csv:from-file "parameter_files/rust_table.dat"
  set phy-susc-dict table:from-list csv:from-file "parameter_files/kaurimate_table.dat"

  ;; Load other parameter files
  load-aux-parameters

  ;; Set up the transition matrices for seedling/sapling dynamics -> ROW wise
  set tr-mtx-bank-yfor (list [0.8 0] [0.05 0.9])    ;; transitions from 1.5m -> 5 cm dbh individuals in yng forest
  set tr-mtx-bank-ofor (list [0.8 0] [0.05 0.9])    ;; ... and old forest

  ;; Set up fire history
  set fire-record []
  set fire-stats []
  set fire-size-list []
  set fire-year []
  set beyond-flamm-time 99999

  ;; Setup clamate, ENSO, ...
  show "Setting up ENSO dynamics..."
  set enso-record []
  init-enso

  ;; Lsp
  build-topo-landscape
  build-veg-landscape

    reset-ticks
end

;; Colours patches by class
to colour-by-class
  ask patches [
    set pcolor table:get colour-dict class
  ]
end


to build-topo-landscape

  ifelse flat-terrain? = true
  [
    ask patches
    [
      set elevation 0
      set aspect -1
      set TPI -1
      set TWI -1
    ]
  ]
  [
    build-catena
  ]

end


to build-catena

  set nhb-scalar 3 * sqrt 2
  ridge-gully true 5 5 100 0.05 0.1 ;;  [desc valley-width elev-step curr-elevation noise diffuse-strength]

  ask patches [
    set flow-accum 1           ;; Each patch contributes at least itself
    set flow-to nobody
  ]

  calculate-D8-flow-direction
  calculate-slope
  calculate-optim-flow-accum
  calculate-TWI
  calculate-TPI

  ;; edaphic gradient is TWI rescaled from 0-1 with 0 wettest and 1 driest
  let max-TWI max [TWI] of patches
  ask patches [ set edaphic-grad 1 - ( TWI / max-TWI ) ]

end

;; D8 flow direction: Choose the steepest downhill neighbor from 8 directions
to calculate-D8-flow-direction
  ask patches [
    let downhill-neighbors neighbors with [elevation < [elevation] of myself]
    if any? downhill-neighbors [
      set flow-to min-one-of downhill-neighbors [elevation]
    ]
  ]
end

;; Slope calculation based on the steepest descent direction
to calculate-slope
  ask patches [
    if flow-to != nobody [
      let dz elevation - ([elevation] of flow-to)
      let di distance flow-to *   patch-grain
      ;set slope math:atan (dz / di)  ;; Slope in radians
      set slope atan dz di * (pi / 180) ;; Slope in radians
    ]
  ]
end

;; Flow accumulation using D8 routing
to calculate-flow-accum
  let patches-sorted patches with [flow-to != nobody]
  ask patches-sorted [
    let receivers patches with [flow-to = myself]
    ask receivers [
      set flow-accum flow-accum + [flow-accum] of myself  ;; Accumulate flow
    ]
  ]
end

;; Flow accumulation using D8 routing from high to low
to calculate-optim-flow-accum
  let sorted-patches reverse sort-on [elevation] patches

  foreach sorted-patches [ p ->
    ask p [
      if flow-to != nobody [
        ask flow-to [
          set flow-accum flow-accum + [flow-accum] of p
        ]
      ]
    ]
  ]
end

;; Compute Topographic Wetness Index (TWI; 0 dry to ... high wet)
;; https://en.wikipedia.org/wiki/Topographic_wetness_index
to calculate-TWI
  ask patches [
    if slope > 0 [
      set TWI ln (flow-accum / ( tan slope * (180 / pi)))
    ]
    if TWI <= 0 [ set TWI 0 ]  ;; Prevent division errors - this is very steep slopes
;; https://gis.stackexchange.com/questions/104543/getting-negative-topographic-wetness-index-twi-values-in-saga-gis/113535#113535
  ]

  ;; Need to tidy up the zero slope locales (set to high value, 95th% of TWI aribitrarily)
  let sort-twi (sort [twi] of patches with [slope > 0 and twi > 0])
  let n round length sort-twi * 0.95
  let zero-twi item n sort-twi
  ask patches with [slope = 0 ] [ set twi zero-twi ]

end

to calculate-TPI
  ask patches
  [
    let m mean [elevation] of other patches in-radius nhb-scalar
    let rel-elevation elevation - m
    if rel-elevation <= -2 [set TPI 1]   ; gully
    if rel-elevation >= -2 and rel-elevation < 2 [set TPI 2] ; slope
    if rel-elevation >= 2 [set TPI 3] ; ridge
  ]

  ask patches with [valley? = true] [ set TPI 1]
  ;; line above is a hack to stop artificial slopes due to valley bottom z noise
end

to ridge-gully [desc valley-width elev-step curr-elevation noise diffuse-strength]
;  let desc true
;  let valley-width 3
;  let elev-step 5
;  let curr-elevation 100
  let i 0

  ask patches with [pxcor = i] [ set elevation curr-elevation ]

  while [i <= max-pxcor] [
    ask patches with [pxcor = i]
    [
      ifelse desc = true
      [set elevation curr-elevation - elev-step]
      [set elevation curr-elevation + elev-step]
    ]

    ifelse desc = true
    [ set curr-elevation curr-elevation - elev-step ]
    [ set curr-elevation curr-elevation + elev-step ]

    if curr-elevation = 0 [
      set desc false
      set i i + (valley-width - 1) + random-poisson (valley-width / 2)
    ]

    if curr-elevation = 100 [set desc true]
    set i i + 1
  ]

  ask patches with [elevation = 0 ] [
    set valley? true
    set elevation elevation + random-float 0.01
  ]

  ask patches [ set elevation (elevation * (random-normal 1 noise))]
  diffuse elevation diffuse-strength

end


;; drive the MRC landscape building
to build-veg-landscape
  set clusters []

  set sum-stalled (list 0 0)
  set sum-changes (list 0 0)

  ;; set up the MRC landscape
  show "Building starting landscape (MRC)..."
  init-patch-mrc-variables
  setup-mrc true false ;; args:: remove-singleton-patches? sequentially-assign-clusters?

  ask patches [set class item class class-names-list]

  ;; This adds forest patches to gullies (preferentially if requested)
  if forest-cover > 0
  [
    show "Adding forest in gullies..."
    build-gully-forest
  ]

  show "Building starting landscape (ageing via cluster id, etc.)..."
  init-patches
  find-clusters
  id-clusters
  if forest-cover > 0 [ tag-gully-forest ]
  age-by-cluster

  ;; Ensure that there are no cells in the 'invaded' state if invasion is false
  if invasion? = false
  [
    ask patches with [class = "d-sh"] [set class "m-sh"]
  ]

  show "Building starting landscape (coastal & farmland)..."
  starting-coastal-forest
  if farm-edge? [edge-farmland]

  ;update-abundances
  colour-by-class

end


to tag-gully-forest

  let ids remove-duplicates sort [cluster-id] of patches with [forest? = true]

  foreach ids
  [  f ->

    let ft random-float 1

    ask patches with [cluster-id = f] [
      if ft <= 0.3 [ set class "yf-k" ]
      if ft > 0.7 [ set class "yf-nok" ]
    ]
  ]
end

;;;;;;;;;;;
;; Code to build landscapes using the MRC method - adapted from the book.

to init-patch-mrc-variables

  ifelse file-exists? init-composition-file = true
  [
    read-class-file-csv init-composition-file
  ]
  [
    print "Initial class composition file missing!"
  ]

  ask patches [
    set class -1
    set in-perc-cluster? false
    set cluster-leader nobody
  ]

end


;;; Retrieve the list of relative probabilities of
;;; suitability classes from a CSV file with names
to read-class-file-csv [class-file]
  carefully [
    let class-abund csv:from-file class-file
    set num-classes length class-abund
    set class-FD slice unnest class-abund range-from-to-by 1 (num-classes * 2) 2
    let sum-class-p sum class-FD

    set class-CFD cumulative-sum class-FD

    set class-CFD map [i -> i / sum-class-p] class-CFD
    set class-FD map [i -> i / sum-class-p] class-FD

;    show class-abund
;    show num-classes
;    show class-CFD
;    show class-FD

    file-close
  ] [ file-close ]
end

;; Internal setup based on SIMMAP algorithm
;; This is slow -- I believe due to repeated use of long lists
to setup-mrc [remove-singleton-patches? sequentially-assign-clusters?]
  ;; assign patches to percolation cluster
  ask patches [
    if random-float 1 < perc-seed [
      set in-perc-cluster? true
      ; set pcolor white
    ]
  ]
  mark-habitat-patches
  assign-habitat-patches-by-proportion sequentially-assign-clusters?
  set clusters [] ;; throw away the list
  assign-unassigned-to-habitat-patches
  if remove-singleton-patches? [
    let singletons patches with [not any? neighbors4 with [class = [class] of myself]]
    ask singletons [
      set class [class] of one-of neighbors4
    ]
  ]
end


;; identifies and labels sequentially
;; from 0 the connected regions
;; in the percolation cluster
to mark-habitat-patches
  let patch-count 1
  let patches-to-mark sort patches with [in-perc-cluster?]
  while [length patches-to-mark > 0] [
    let this-cluster (patch-set one-of patches-to-mark)
    let current-cluster patch-set nobody
    ;; iteratively grow the patch-set that is this cluster
    while [any? this-cluster] [
      ;; mark by setting habitat-patch value
      ask this-cluster [
        set class patch-count
      ]
      set current-cluster (patch-set current-cluster this-cluster)
      ;; get the next set from the neighbours 4 of the current set
      set this-cluster (patch-set [neighbors4] of this-cluster) with [to-mark?]
    ]
    ;; increment the patch-count and reset the set to mark
    set patch-count patch-count + 1
    set patches-to-mark filter [i -> [class] of i < 0] patches-to-mark
    ;; add the current cluster to the list of all clusters
    set clusters lput current-cluster clusters
  ]
end


to-report to-mark?
  report class < 0 and in-perc-cluster?
end


;; class-FD is a list of the required landscape proportion in each type
;; cells in regions of the percolation cluster should be relabelled
;; as 1, 2, 3... in the required proportions...
to assign-habitat-patches-by-proportion [sequentially-assign-clusters?]
  ;; make a list of the target number of patches in each class
  let num-patches-set count patches with [class > 0]
  let target-counts map [ i -> i * num-patches-set ] class-FD
  ;; iteratively assigning clusters in the percolation cluster
  ;; to a class, always adding to the class which is furthest
  ;; from its target count
  if sequentially-assign-clusters? [
    set clusters reverse sort-by [ [a b ] -> count a > count b] clusters
  ]

  foreach clusters [ i ->
    let biggest-shortfall max target-counts
    let biggest-shortfall-index position biggest-shortfall target-counts
    ask i [
      set class biggest-shortfall-index
    ]
    ;; update the target counts to reflect assignment just made
    set target-counts replace-item biggest-shortfall-index target-counts (biggest-shortfall - count i)
  ]
end


;; patches so far not assigned to a class are
;; now assigned based on the classes of neighbouring assigned patches   let ids remove-duplicates sort [cluster-id] of patches   let ids remove-duplicates sort [cluster-id] of patches   let ids remove-duplicates sort [cluster-id] of patches
to assign-unassigned-to-habitat-patches
  ask patches with [class < 0] [
    let N neighbors with [class >= 0]
    ifelse any? N
    [ set class one-of modes [class] of N ]
    [ let random-draw random-float 1
      let i 0
      let value-set? false
      while [ i < num-classes and not value-set? ] [
        if random-draw <= (item i class-CFD) [
          set class i
          set value-set? true
        ]
        set i i + 1
      ]
    ]
  ]
end


;; routines here label the patches at the end:: modified nl model library
to find-clusters
  loop [
    ;; pick a random patch that isn't in a cluster yet
    let seed one-of patches with [cluster-leader = nobody]

    ;; if we can't find one, then we're done!
    if seed = nobody [ stop ]

    ;; otherwise, make the patch the "leader" of a new cluster
    ;; by assigning itself to its own cluster, then call
    ;; grow-cluster to find the rest of the cluster
    ask seed
    [ set cluster-leader self
      grow-cluster ]
  ]
end

to grow-cluster  ;; patch procedure
  ask neighbors4 with [(cluster-leader = nobody) and (class = [class] of myself)]
  [ set cluster-leader [cluster-leader] of myself
    grow-cluster ]
end

;; once all the clusters have been found, this is called
;; to put numeric labels on them so the user can see
;; that the clusters were identified correctly
to id-clusters
  let counter 0
  loop
  [ ;; pick a random patch we haven't labeled yet
    let pch one-of patches with [cluster-id = 0]
    if pch = nobody [ stop ]

    ;; give all patches in the chosen patch's cluster
    ;; the same label
    ask pch
    [ ask patches with [cluster-leader = [cluster-leader] of myself]
      [ set cluster-id counter ] ]
    set counter counter + 1 ]
end

;; this then gives each cluster an age based on its id
to age-by-cluster
  let ids remove-duplicates sort [cluster-id] of patches

  foreach ids
  [  i ->
     let f-class [class] of one-of patches with [cluster-id = i]
     let half-age ((table:get base-changes-dict f-class ) / 2)
     let nc half-age + random half-age

     ask patches with [cluster-id = i]
     [
       set next-change nc
     ]
  ]
end

;; this is code to switch patch classes around so pohutakawa appears close to coast
;; at start up.
to starting-coastal-forest

  ;let not-coastal-forest patches with
  ;let coastal-forest patches

  let ids remove-duplicates sort [cluster-id] of patches

  ;; loop over patches using their cluster id from MRC
  foreach ids
  [  i ->

    let focal-cl patches with [cluster-id = i]
    let focal-cl-class [class] of one-of focal-cl


    let mean-dist mean [distance-to-coast] of focal-cl

    ;; if closer than critical dist to coast (average) then see if change to *-p
    ;; the one-of lets us use member? and relies on the fact that the entire cluster is one class.

    if (mean-dist <= table:get aux-params-dict "max-poh-coast") and
                   (focal-cl-class = "ksh-k" or focal-cl-class = "yf-k" or focal-cl-class = "ksh-nok" or focal-cl-class = "yf-nok" or focal-cl-class = "old-f")
    [

      let b0 table:get aux-params-dict "poh-coast-b0"
      let b1 table:get aux-params-dict "poh-coast-b1"
      let b2 table:get aux-params-dict "poh-coast-b1"

      if random-float 1 < (decr-limit-fx b0 b1 b2 (mean-dist / table:get aux-params-dict "max-poh-coast"))
      [
        ; print "changing"(mean [edaphic-grad] of focal-cl)
        ask focal-cl
        [
          if class = "ksh-k" [set class "ksh-p" ]
          if class = "yf-k" [set class "yf-p"]
          if class = "ksh-bok" [set class "ksh-p"]
          if class = "yf-nok" [set class "yf-p"]
          if class = "old-f" [set class "old-p"]

          set prev-class class
        ]
      ]
    ]

    ;; and the opposite *-p to X, with X also considering the edaphic conditions

    if (mean-dist > table:get aux-params-dict "max-poh-coast") and
                     (focal-cl-class  = "ksh-p" or focal-cl-class = "yf-p" or focal-cl-class = "old-p")
    [
      let edaphic-flip random-float 1    ;; do it outside ask so whole patch changes
      let mean-ed (mean [edaphic-grad] of focal-cl)
      ask focal-cl
      [
        if class = "ksh-p"
        [
          ifelse edaphic-flip < mean-ed  ;; *-nok or *-k based on edaphics
          [ set class "ksh-k"]
          [ set class "ksh-nok" ]

        ]

        if class = "yf-p"
        [
          ifelse edaphic-flip < mean-ed
          [ set class "yf-k" set prev-class class]
          [ set class "yf-nok" set prev-class class]
        ]

        if class = "old-p" [ set class "old-f" ]

        set prev-class class
      ]
   ]

  ]
end

to edge-farmland

  let edge-patches patches with [pxcor = max-pxcor]

  ask n-of (random-poisson farm-edge-nodes) edge-patches [set farm-node? true]

  ask patches with [farm-node? = true]
  [
    let sz random-poisson mean-farm-depth
    let px pxcor
    let py pycor

    ask patches with [(pycor < py + sz) and (pycor > py - sz) and (pxcor > px - sz)]
    [
      set class "gr"
      set prev-class class
      set farm? true
    ]
  ]
end

;; initialise the patch fire and regen history
to init-patches

  ask patches
  [
    ; set forest? false

    set times-change 0
    set flammability table:get flammability-dict class
    set burned? false
    set regenbank-yfor (list 0 0)
    set regenbank-ofor (list 0 0)

    set fire-history []
    set next-change table:get base-changes-dict class
    set t-colonised (list -1 -1)
    set prev-class class

    set distance-to-coast pxcor * patch-grain
    set farm-node? false

    set myrtle-rust? false
    set kauri-mate? false
    set kauri-mate-nhb other patches in-radius phy-radius-inf
  ]

end

;; b1 upper, b0, lower, b2 steep
to-report decr-limit-fx [b0 b1 b2 x]
  report b1 *  exp(- b2 * x) + b0
end

;; b1 upper, b1 steep
to-report incr-limit-fx [b1 b2 x]
  report b1 * (1 - exp(- b2 * x))
end



to load-aux-parameters

  ;; general parameters for curves etc.
  set aux-params-dict csv:from-file "parameter_files/aux_parameters.dat"
  set aux-params-dict table:from-list aux-params-dict

  ;; fire wind and slope modifiers
  ;; we strip off the title so it is a list of numeric lists
  let X map but-first (csv:from-file "parameter_files/fire_weights.dat")

  set flamm-wind-wgt sublist X 0 3 ;; upper is exclusive
  set flamm-slope-wgt item 3 X
  set flamm-slope-breaks item 4 X
  set flamm-aspect-wgt item 5 X
  set enso-freq-wgt-list item 6 X

end


to init-enso

  set enso-mtx csv:from-file enso-matrix-file   ; load the matrix
  set enso-list first enso-mtx                   ; get the states possible (headers)
  set enso-mtx but-first enso-mtx                ; strip to mtx only

  set enso-state one-of enso-list         ; get the first state

  set enso-record lput enso-state enso-record
  ; show enso-mtx
  ; show enso-state

end

to build-gully-forest
  let idx 0
  while [count patches with [forest? = true] < (forest-cover * max-pxcor * max-pycor)]
  [
    start-forest-grow idx
    set idx idx + 1
  ]

  ask patches with [forest? = true] [set class "old-f"]

end

to start-forest-grow [idx]
  set forestArea 0

  ask one-of patches with [edaphic-grad < 0.4 and tpi = 1] ;; gully
  [
    set forest? true
    set f-shell idx
    set forestFront patch-set self              ;; a new disturb-front patch-set
  ]

  forest-grow idx (forest-cover * max-pxcor * max-pycor)

end


to forest-grow [idx maxArea]
  while [ any? forestFront and forestArea <= maxArea ]                 ;; Stop when we run out of active disturbance front
  [
    let newForestFront patch-set nobody     ;; Empty set of patches for the next 'round' of the disturbance

    ask forestFront
    [
      set f-shell idx
      let N neighbors;;  with [ forest? != true ]

      ask N
      [
        let p 0.125
        if twi < 0.35 [set p 0.4]
        if twi > 0.35 and twi < 0.75 [set p 0.25]

        if (random-float 1) <= p [ set newForestFront ( patch-set newForestFront self) ]

      ]
   ]

   set forestFront newForestFront
   ask newForestFront [set forest? true]
   set forestArea forestArea + ( count newForestFront )
  ]

end



@#$#@#$#@
GRAPHICS-WINDOW
230
10
642
423
-1
-1
4.0
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
100
0
100
1
1
1
ticks
30.0

BUTTON
81
166
220
199
NIL
build-catena
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
757
95
863
128
show ed
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

BUTTON
867
95
938
128
show lf
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

BUTTON
80
203
222
236
NIL
build-veg-landscape
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
761
178
867
211
invasion?
invasion?
1
1
-1000

SWITCH
759
140
878
173
farm-edge?
farm-edge?
0
1
-1000

INPUTBOX
761
226
1035
286
init-composition-file
parameter_files/initial_shrub_composition.dat
1
0
String

SLIDER
1058
93
1230
126
perc-seed
perc-seed
0
1
0.58
0.01
1
NIL
HORIZONTAL

SLIDER
1058
172
1230
205
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
1060
132
1232
165
mean-farm-depth
mean-farm-depth
0
20
10.0
1
1
NIL
HORIZONTAL

SLIDER
1059
235
1231
268
phy-radius-inf
phy-radius-inf
0
3
1.5
.1
1
NIL
HORIZONTAL

BUTTON
113
39
186
72
build-all
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

INPUTBOX
761
289
1033
349
enso-matrix-file
parameter_files/enso_matrix.dat
1
0
String

BUTTON
941
96
1055
129
show veg
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

PLOT
763
366
963
516
EG (sc TWI)
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"set-histogram-num-bars 20" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [edaphic-grad] of patches"

SWITCH
889
139
1010
172
flat-terrain?
flat-terrain?
1
1
-1000

PLOT
976
366
1176
516
TWI
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"set-histogram-num-bars 20" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [twi] of patches"

TEXTBOX
336
480
351
498
*
11
0.0
1

TEXTBOX
469
480
619
498
*
11
0.0
1

SLIDER
1063
301
1235
334
forest-cover
forest-cover
0
1
0.1
.01
1
NIL
HORIZONTAL

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
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

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
