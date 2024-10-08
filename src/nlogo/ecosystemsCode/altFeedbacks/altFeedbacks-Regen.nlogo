__includes["mrc.nls" "regen-mtx.nls"]
extensions [profiler matrix]

;; classes: 0 = invaded, 1 = manuka, 2 = kanuka, 3 = yng forest (mid-succ), 4 = mature forest (late-succ)


globals [ 
  from-R?
  
  run-id
  
  world-size
  record-tag
  
  ;; fire-related
  flammability-list
  fire-front 
  fire-size
  s
  extinguished?
  class-list
  abundances
  
  ;; fire history 
  ;; two lists to deal with individual fires and their characteristics
  ;; third list keeps track of fire sizes for aggregate stats at end
  fire-stats
  fire-record
  fire-size-list
  
  
  beyond-flamm-time
  
  base-changes
  ;;accum-changes
  ;;base-lag
  
  ;; regeneration and succession things
  matrix-3
  matrix-4
  sum-stalled       ;; total stalled list 0 = 2 -> 3, 1 = 3 -> 4
  sum-changes       ;; no. stalled list 0 = 2 -> 3, 1 = 3 -> 4
  n-changes
  
  ;; these globals control the MRC algorithm
  colour-list        ;; list of spectral colors
  num-classes       ;; the number of suitability classes
  class-FD          ;; freq distribution of classes requested
  class-CFD         ;; cumulative freq dist of classes requested  
  clusters          ;; a list of patch-sets of patches in each cluster
] 


patches-own [ 
  fire-history      ;; list of years in which cell burned
  last-change
  next-change
  times-change
  t-colonised
  
  regenbank-3
  regenbank-4
  ldd-this-tick
   
  stalled
  
  flammability
  burned? 
  
  class             ;; used to label initial clusters and to index veg class
  prev-class        ;; class previously occupied
  
  cluster-leader    ;; patch which 'leads' the cluster (see model lib example)
  in-perc-cluster?  ;; true if patch is in initial percolation cluster
  cluster-id
]


  
to setup
  clear-all
  if from-r? != TRUE [set from-r? FALSE]
  
  set world-size count patches
  
  set record-tag ""
  
  set s new-seed
  set colour-list (list 35 43 68 65 62)
  set flammability-list (list 0.65 0.35 0.25 0.1 0.05)
  set base-changes (list 20 30 80 120 -999);; -> this is time spent in each state
  ;;set accum-changes (list 0 20 50 130 250)
  ;; so if start in all manuka [item *1*] and no stalling:: 231 years total
  
  set fire-record []
  set fire-stats []
  set fire-size-list []
  
  ;; set base-lag (list 15 40)   ;; time after min density acheived of transition (legacy from when didn't have stage structure)
  
  set beyond-flamm-time 99999
  
  set matrix-3 matrix:from-column-list [[0.9 0.05][0 0.99]]    ;; transitions from 1.5m -> 5 cm dbh individuals
  set matrix-4 matrix:from-column-list [[0.9 0.05][0 0.99]]    

   
  ;; Build the landscape using the MRC code...
  build-landscape
  
  reset-ticks
end


to build-landscape
  set clusters []
  set class-list n-values 5 [?]
  
  set sum-stalled (list 0 0)
  set sum-changes (list 0 0)
  
  init-patch-mrc-variables
  
  read-class-file init-composition
    
  setup-mrc true false ;; args:: remove-singleton-patches? sequentially-assign-clusters?
  
  ;; Ensure that there are no cells in the 'invaded' state if invasion is false
  if invasion? = false
  [
    ask patches with [class = 0] [set class 1]
  ]  
  
  init-patch-dynamics
  find-clusters
  id-clusters
  age-by-cluster
  
  update-abundances
end  
  

to init-patch-dynamics
  ask patches 
  [
    set times-change 0
    set flammability item class flammability-list 
    set burned? false 
    set regenbank-3 matrix:from-column-list [[0 0]]
    set regenbank-4 matrix:from-column-list [[0 0]]     
        
    set fire-history []
    set next-change item class base-changes
    set t-colonised (list -1 -1)
    set prev-class class
  ]
  
  color-patches
end  


to go
  ;;while [ticks <= 5]
  ;;[
    ;;profiler:start         ;; start profiling
    
    ;if ticks = 0 [print date-and-time]
    ;if ticks = 300 [print date-and-time  stop ]
    
    while [ ticks <= max-ticks and ((item 4 abundances / world-size) <= max-forest)  ]
    [
      
      let fire-this-tick false
      if ticks > burn-in-regen
      [
        if random-float 1 <= fire-frequency
        [
          let start? ignite-fire
          if start? = true
          [
            fire-spread
            post-fire
            set fire-this-tick true
          ]  
        ]
      ]  
      
      set n-changes 0
      
      dispersal
      regenerate-patch-bank
      
      if ticks > burn-in-regen
      [
        succession
      
        if n-changes > 0 or fire-this-tick = true
        [
          color-patches
          update-abundances
          if abund-flammable < 0.3 and ticks < beyond-flamm-time
          [
            set beyond-flamm-time ticks
          ]  
        ]
      ]
;      profiler:stop          ;; stop profiling
;      print profiler:report  ;; view the results
;      profiler:reset    
      
      tick
    ]
  ;write-record  
end

  
to-report ignite-fire
  set fire-size 0
  let succ-ignition? false

  let target nobody
  ifelse any? patches with [ class = 0 or class = 1 ] and random-float 1 <= flamm-start
  [ 
    set target one-of patches with [ class = 0 or class = 1 ]
  ]
  [
    set target one-of patches 
  ]  
  
  ask one-of patch-set target [ ;; with [pcolor = green] [
    
    if random-float 1 <= ([flammability] of target * 2)
    [
      ;; record initial fire conditions
      set fire-stats (list (ticks)  ([class] of target)  (count [neighbors with [class <= 2]] of target ) (mean [flammability] of patches) )
    
      foreach abundances
      [
        set fire-stats lput (? / world-size) fire-stats
      ]  
    
      set burned? true
      set pcolor red
      set fire-front patch-set self               ;; a new fire-front patch-set
      set fire-size 1
      set succ-ignition? true                                 ;; a fire has started    
    ]
  ]
  
  report succ-ignition?  
end

to fire-spread
  while [ any? fire-front ]                 ;; Stop when we run out of flaming fire-front
  [ 
    fire-shell
  ]
end

to fire-shell 
   let new-fire-front patch-set nobody     ;; Empty set of patches for the next 'round' of the fire
   
   ask fire-front [
     
     let N neighbors with [ not burned? ]
     
     ask N
     [
       let edge-wgt count neighbors with [class <= 1] * 0.025
       ;; show edge-wgt
       
       if random-float 1 <= flammability * (1 + edge-wgt) 
       [ 
         set new-fire-front ( patch-set new-fire-front self) 
         set burned? true 
       ]
     ]
   ]
   
   ;; ask fire-front [set pcolor grey]    
   ask new-fire-front [set pcolor red]
   set fire-front new-fire-front
end


;; This sorts out what happens after a fire - this is where pyrophyllic invasion takes place
to post-fire
   
   let new-invasions 0
   ask patches with [burned?]
   [
     set prev-class class
     
     set fire-history lput ticks fire-history
     let last30 length filter [? >= ticks - 30] fire-history
     let local-weeds count neighbors with [class = 0]
     
     set regenbank-3 matrix:from-column-list [[0 0]]
     set regenbank-4 matrix:from-column-list [[0 0]]
     
     ifelse invasion? = true
     [      
       ifelse random-float 1 <= (base-invasion + (fire-invasion * last30)) * (1 + (local-weeds / 8)) 
       [  
         set class 0
         set next-change ticks + (item class base-changes) + (last30 * fire-slow)
         set flammability item class flammability-list 
       ]
       [
         set class 1
         set next-change ticks + (item class base-changes) + (last30 * fire-slow)
         set flammability item class flammability-list 
       ]
     ]
     
     [
         set class 1
         set next-change ticks + (item class base-changes) + (last30 * fire-slow)
         set flammability item class flammability-list 
     ]     
   ]
   
   set new-invasions count patches with [class = 0 and prev-class > 0]
   
   set fire-size count patches with [burned?]
   
   set fire-size-list lput fire-size fire-size-list
   set fire-stats lput fire-size fire-stats


   let types-burned map [count patches with [burned? and prev-class = ?]]  (n-values 5 [?])
   foreach types-burned
   [
     set fire-stats lput (? / world-size) fire-stats
   ]
   
   set fire-stats lput mean [flammability] of patches fire-stats
   set fire-stats lput new-invasions fire-stats
   
   set fire-record lput fire-stats fire-record
   
   ask patches [set burned? false] 
    
   color-patches
end


;; This is the linear 'time elapsed' component of the succession process
to succession
     
  ;; set base-changes (list 20 30 80 120 -999)  -> time spent in class (30 in man, 80 in kan, 120 in yng)
  let change-patches patches with [next-change <= ticks and next-change != -999]
  if any? change-patches
  [
    ask change-patches
    [
      set prev-class class
      ;; invaded -> manuka shrubland
      
      if class = 0
      [
        set class 1
        
        let last30 length filter [? >= ticks - 30] fire-history
        let base item class base-changes    ;; so time to *next* change estimated here
        
        set next-change (ticks + (base + (last30 * 2))) * (0.9 + random-float 0.2)
        ;;next change yr = base change + fire slowing + random noise   
       
        set last-change ticks
        set times-change times-change + 1
        
        set n-changes n-changes + 1
      ]
      
      ;; manuka shrubland -> kanuka  
      if class = 1 and last-change != ticks
      [
        set class 2
        
        let last30 length filter [? >= ticks - 30] fire-history
        let base item class base-changes
        
        set next-change ceiling ((ticks + (base + (last30 * 2)))  * (0.9 + random-float 0.2)) 
        ;;show next-change
        set last-change ticks
        set times-change times-change + 1
        set stalled 0
                
        set n-changes n-changes + 1
      ]
      
      ;; kanuka -> young forest
      if class = 2 and last-change != ticks
      [
        let n-saps max (list matrix:get regenbank-3 1 0 matrix:get regenbank-4 1 0)    ;; n-trees is max of type 3 or 4 to stop spurious stalling
        ifelse n-saps >= crit-density-yng ;; and (ticks - (item 0 t-colonised)) >= (item 0 base-lag)
        [
          
          if track-stalled?
          [
            let n item 0 sum-changes
            let stl item 0 sum-stalled
            
            set sum-changes replace-item 0 sum-changes (n + 1)
            set sum-stalled replace-item 0 sum-stalled  (stl + stalled)    ; replace-item index list value 
          ]
          
          set class 3
          let base item class base-changes * (0.9 + random-float 0.2) 
          
          
          set next-change ticks + base
          set last-change ticks
          set times-change times-change + 1
          set stalled 0
          
          set n-changes n-changes + 1
        ]  
        [
          set stalled stalled + 1
        ]  
      ]
      
      ;; young forest -> old forest
      if class = 3 and last-change != ticks
      [
        let n-saps matrix:get regenbank-4 1 0
        
        ifelse n-saps >= crit-density-old ;; and (ticks - (item 1 t-colonised)) >= (item 1 base-lag)
        [
          if track-stalled?
          [
            let n item 1 sum-changes
            let stl item 1 sum-stalled
            
            set sum-changes replace-item 1 sum-changes (n + 1)
            set sum-stalled replace-item 1 sum-stalled  (stl + stalled)    ; replace-item index list value 
          ]
          
          set class 4
          set stalled 0
          set next-change -999
          set times-change times-change + 1
          set last-change ticks
          
          set n-changes n-changes + 1
        ]
        [
          set stalled stalled + 1
        ]
      ]
      
      set flammability item class flammability-list
      
    ]
  ]
end  

to dispersal
  let target nobody
   
   ;; what is dispersed is 1.5 m high saps (effectively)
     
    ask patches with [ class >= 3 ]
    [
     let seed-rain 0
     if class = 3 [set seed-rain round (random-poisson base-seed-prod-3)]
     if class = 4 [set seed-rain round (random-poisson base-seed-prod-4)]
      
     repeat seed-rain
     [
       
       if random-float 1 >= seed-pred    ;; avoids predation?
       [
         ifelse random-float 1 <= fraction-consumed  ;; consumed and dispersed by pigeon
         [
           let d random-exponential mean-ldd
           set target patch-at-heading-and-distance (random-float 360) d
         ]
         [
           set target one-of (patch-set self neighbors)
         ]
     
         if target != nobody
         [
        ;; assume that seedlings of class 3 can establish under any sort of veg, except invaded shrubland
        ;; but class 4 needs to be under class 2 or 3
        
           if class = 3 and [class] of target >= 1 ;; no establishment of class 3 in invaded shrubland
           [
             let idx class - 3   ; so 3 -> 0 and 4 -> 1
             let n-sdl matrix:get regenbank-3 0 0        ; abundance of sdl of class x in the patch
             
             ask target 
             [ 
               matrix:set regenbank-3 0 0 (n-sdl + 1)       ;; matrix:set matrix row-i col-j new-value
               if (n-sdl + 1) >= crit-density-yng [set t-colonised replace-item idx t-colonised ticks]
             ]
           
           ]
           
           if class = 4 and [class] of target >= 2   ;; class 4 only establishes in kanuka woodland [2] or older
           [
             let idx class - 3
             let n-sdl matrix:get regenbank-4 0 0    ; abundance of sdl of class x in the patch
             
             ask target 
             [ 
               matrix:set regenbank-4 0 0 (n-sdl + 1)       ;; matrix:set matrix row-i col-j new-value
               if (n-sdl + 1) >= crit-density-old [set t-colonised replace-item idx t-colonised ticks]
             ]
           ]    
             
      
           ]
         ]   
     ] 
   ]
     
  
end

;to thin-regenbank
;  
;  if sap-mortality > 0
;  [
;    ask patches 
;    [
;      if sum propagules > 0 
;      [
;        set propagules( map [rbinom ?1 (1 - sap-mortality)] propagules)
;      ]
;    ]     
;  ]
;end

to update-abundances
  let n [class] of patches
   set abundances (map [occurrences ?1 n] class-list)  
end 

to-report occurrences [x the-list]
  report reduce
    [ifelse-value (?2 = x) [?1 + 1] [?1]] (fput 0 the-list)
end

;; Colours patches by class
to color-patches
  ask patches [
    set pcolor item class colour-list
  ]  
end

to color-by-stall
  let mx max [stalled] of patches
  ask patches with [ stalled > 0 ] [
    set pcolor scale-color orange stalled 1 mx + 5
  ]
end


to color-by-regen3
  ask patches [ 
    set pcolor scale-color blue (sum matrix:get-column regenbank-3 0) 100 0
  ]  
end

to color-by-regen4
  ask patches [ 
    set pcolor scale-color blue (sum matrix:get-column regenbank-4 0) 100 0
  ]  
end
  
;; Abundance flammable
to-report abund-flammable
  report count patches with [class <= 2] / world-size
end  

;; Abundance stalled
to-report abund-stalled
  report count patches with [stalled > 0] / world-size
end  


to-report rbinom [n prob]

  let x 0
  let i 0
  while [i < n]
  [
    if random-float 1 < prob [set x x + 1]
    set i i + 1
  ]
  report x
  
end

to write-record
  if write-record? = true
  [
    ;; let exists? file-exists? "fireRecord.txt"
    let tag 0
    if from-r? = FALSE
    [ set run-id behaviorspace-run-number ]
    
    
    let file-name ( word record-tag "fireRecord-" run-id ".txt" )
    file-open file-name
    
    file-print " Run Year start.veg start.hb pre.flamm pre.inv pre.man pre.kan pre.yfor pre.ofor fire.size burn.inv burn.man burn.kan burn.yfor burn.ofor post.flamm new.invasions"  
    
    let n length fire-record
    let idx 0
    
    foreach fire-record
    [
      let data ?
      file-write run-id
      
      foreach data
      [
        file-write precision ? 3 
      ]
      file-newline
      set idx idx + 1
    ]  
    
    file-close
  ]  
end    

to file-newline
  file-print ""
end

to-report random-binomial [sz pr]
  report length filter [?] n-values sz [random-float 1 < pr]
end  

;to-report get-binomial [sz pr]
;  r:put "p" pr
;  r:put "s" sz
;  
;  report r:get "rbinom(1,s,p)"
;end  

to calibrate-bank
 file-open "calibrate-bank.txt"
 file-print "class last.change n5.d3 n5.d4 n10.d3 n10.d4 dist.d3 last.d3 last.d4 dist.d4  r3.sdl r3.sap r4.sdl r4.sap"
 
 ask patches
 [
   let d3 patches in-radius 10 with [ class = 3 ]
   let d4 patches in-radius 10 with [ class = 4 ]
   
             
   file-write class
   file-write last-change
   
   file-write count d3 in-radius 5 with [ class = 3 ]
   file-write count d4 in-radius 5 with [ class = 4 ]
      
   file-write count d3
   file-write count d4
   
   
   ifelse any? d3
     [ 
       let nn-d3 min-one-of d3 [distance myself] 
       file-write precision distance nn-d3 3
       file-write precision (mean [last-change] of d3) 2
     ]     
     [ repeat 2 [ file-write -1 ] ]
 
   ifelse any? d4
     [ 
       let nn-d4 min-one-of d4 [distance myself] 
       file-write precision distance nn-d4 3
       file-write precision (mean [last-change] of d4) 2
     ]     
     [ repeat 2 [ file-write -1 ] ]
   
   file-write matrix:to-column-list regenbank-3
   file-write matrix:to-column-list regenbank-4
   
   file-print ""
 ]
 file-close
end 
   
;; ret <- data.frame(x, NLReport("max-ticks"), NLReport("burn-in-regen"), NLReport("fire-frequency"), NLReport("item 3 abundances")) 
to-report report-evaluation-variables
    ;; LENGTH = 11 + 1 for X
    
    
    let dump lput beyond-flamm-time abundances
    set dump lput ticks dump
    ifelse length fire-size-list > 0
    [
      set dump lput (length fire-size-list) dump
      set dump lput (max fire-size-list) dump
      set dump lput (sum fire-size-list) dump
      set dump lput (mean fire-size-list) dump
    ]
    [
      set dump lput 0 dump
      set dump lput 0 dump
      set dump lput 0 dump
      set dump lput 0 dump
    ]
    
    
    report dump
end  

to-report mle-exponent [size-list xmin]
  carefully
  [
    let b 1
    let mle-est map [(log ? 2) / xmin] size-list
    set b 1 + ((sum mle-est) ^ -1) 
    report b
  ]
  [
    report -999
  ]  
end
   
   
  
@#$#@#$#@
GRAPHICS-WINDOW
241
10
701
491
-1
-1
3.0
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
149
0
149
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
p
p
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
706
10
1054
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
"invaded" 1.0 0 -6459832 true "" "plot item 0 abundances / world-size"
"manuka" 1.0 0 -7171555 true "" "plot item 1 abundances / world-size"
"kanuka" 1.0 0 -5509967 true "" "plot item 2 abundances / world-size"
"yngfor" 1.0 0 -13840069 true "" "plot item 3 abundances / world-size"
"oldfor" 1.0 0 -15575016 true "" "plot item 4 abundances / world-size"

SLIDER
31
93
203
126
fire-frequency
fire-frequency
0
1
0.161069144494832
.01
1
NIL
HORIZONTAL

SLIDER
30
130
202
163
fraction-consumed
fraction-consumed
0
1
0.266056878282689
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
0.88986198581988
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
0.161238820757717
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
706
227
1055
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
"lspflm" 1.0 0 -16777216 true "" "plot mean [flammability] of patches"
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
4.90936847426929
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

MONITOR
1057
10
1114
55
BFT
beyond-flamm-time
0
1
11

SLIDER
32
277
204
310
base-invasion
base-invasion
0
1
0.102435674885055
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
1
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
0.0897684312309139
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
710
386
816
419
invasion?
invasion?
0
1
-1000

INPUTBOX
937
382
1056
442
init-composition
classes.txt
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
base-seed-prod-4
base-seed-prod-4
0
20
8
1
1
NIL
HORIZONTAL

SLIDER
29
486
201
519
base-seed-prod-3
base-seed-prod-3
0
20
4
1
1
NIL
HORIZONTAL

SLIDER
239
495
411
528
crit-density-yng
crit-density-yng
1
20
12
1
1
NIL
HORIZONTAL

SLIDER
241
532
413
565
crit-density-old
crit-density-old
1
20
14
1
1
NIL
HORIZONTAL

TEXTBOX
47
526
197
559
Minimum 'density' of juvenile saps to make transition - this is effectively an index of dispersal.
9
0.0
1

BUTTON
818
386
919
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
828
422
928
455
show-class
color-patches
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
711
462
883
495
max-ticks
max-ticks
1
4000
250
10
1
NIL
HORIZONTAL

SLIDER
710
498
882
531
burn-in-regen
burn-in-regen
0
50
10
5
1
NIL
HORIZONTAL

BUTTON
887
463
1007
496
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
888
500
1008
533
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

SLIDER
710
534
882
567
max-forest
max-forest
0
1.01
1.01
.01
1
NIL
HORIZONTAL

SWITCH
711
422
821
455
write-record?
write-record?
0
1
-1000

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
NetLogo 5.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="predation" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <final>write-record</final>
    <exitCondition>count patches with [class = 4] &gt;= (0.95 * world-size)</exitCondition>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed-pred" first="0" step="0.05" last="1"/>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="ldd-by-pred" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>count patches with [class = 4] &gt;= (0.95 * world-size)</exitCondition>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fraction-consumed" first="0" step="0.1" last="1"/>
    <steppedValueSet variable="seed-pred" first="0" step="0.1" last="0.9"/>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fire-noinvasion" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="500"/>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fire-frequency" first="0" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fire-invasion" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="500"/>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fire-frequency" first="0" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p">
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
  </experiment>
  <experiment name="consumed" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>count patches with [class = 4] &gt;= world-size * 0.95</exitCondition>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fraction-consumed" first="0.05" step="0.05" last="1"/>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fire-invasion-regen" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="500"/>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.2"/>
      <value value="0.4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="fire-frequency" first="0" step="0.01" last="0.2"/>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0"/>
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="test" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <exitCondition>count patches with [class = 4] &gt;= (0.95 * world-size)</exitCondition>
    <metric>abundances</metric>
    <metric>beyond-flamm-time</metric>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="burn-in-regen" first="0" step="10" last="50"/>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-3">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition">
      <value value="&quot;classes.txt&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-old">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-yng">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-4">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p">
      <value value="0.56"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="start-test" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>write-record</final>
    <timeLimit steps="35"/>
    <metric>abundances</metric>
    <enumeratedValueSet variable="flamm-start">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-4">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="invasion?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-consumed">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-yng">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-invasion">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-composition">
      <value value="&quot;classes.txt&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="track-stalled?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-slow">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-ldd">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-forest">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-pred">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="crit-density-old">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="burn-in-regen">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fire-frequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="35"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-seed-prod-3">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="base-invasion">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p">
      <value value="0.57"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="write-record?">
      <value value="true"/>
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
