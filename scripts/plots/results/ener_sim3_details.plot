
reset

set term qt size 300,260 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [0:6]
set yrange [32:36]
set ytics 33,1,35 right scale .5
set grid xtics ytics

NOXTICS = "set xtics ('' 0, '' 2, '' 4, '' 6) scale .5; \
          unset xlabel"
XTICS = "set xtics ('0' 0, '2' 2, '4' 4, '6' 6) scale .5 offset 0,.3;\
          set xlabel 'Time (sec)'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'Power (W)'"

TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.58"
BMARGIN = "set tmargin at screen 0.56; set bmargin at screen 0.18"
LMARGIN = "set lmargin at screen 0.16; set rmargin at screen 0.49"
RMARGIN = "set lmargin at screen 0.56; set rmargin at screen 0.89"


set multiplot layout 2,2 rowsfirst

@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
LABEL='I'
set ylabel offset 0,-3.4;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.818,0.86
plot '../../../data/simulation3/raw5/updated/energy_simulation3D.csv' using 1:2 w l ls 1 lw .6, \
     '' using 1:3 w l ls 1 lw .6 lc rgb '#FF0000' 

@TMARGIN; @RMARGIN
@NOXTICS;

set yrange [35:39]
set ytics 36,1,38 right scale .5
unset ylabel

LABEL='II'
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.79,0.86
set key left
set key spacing .5
plot '../../../data/simulation3/raw1/updated/energy_simulation3.csv' using 1:2 w l title 'data' ls 1 lw .6, \
     '' using 1:3 w l ls 1 lw .6 lc rgb '#FF0000' title 'KF' 

@BMARGIN; @LMARGIN
@XTICS; @YTICS
LABEL='III'
unset key
unset ylabel
set xtics ('0' 0, '2' 2, '4' 4, '6' 6) scale .5 offset 0,.3
set xlabel offset 12.4,.9
set yrange [30.6:39]
set ytics 32,2,38 scale .5
set obj 1 rect at graph 0.85,0.86 size char strlen(LABEL)+.7, char 1
set label 1 LABEL at graph 0.8,0.86
plot '../../../data/simulation3/raw3/updated/energy_simulation3B.csv' using 1:2 w l ls 1 lw .6,\
     '' using 1:3 w l ls 1 lw .6 lc rgb '#FF0000'

@BMARGIN; @RMARGIN
@XTICS; 
LABEL='IV'
set xtics ('0' 0, '2' 2, '4' 4, '6' 6) scale .5 offset 0,.3
unset xlabel
set yrange [46.6:55]
set ytics 48,2,54 scale .5
set obj 1 rect at graph 0.858,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.8,0.86 font ',8'
plot '../../../data/simulation3/raw4/updated/energy_simulation3C.csv' using 1:2 w l ls 1 lw .6,\
     '' using 1:3 w l ls 1 lw .6 lc rgb '#FF0000'

unset label 1
unset obj 1
unset margin
unset ytics
unset xtics
set grid x2tics,ytics

set yrange [0:45]
set ytics ('37.7' 37.7) scale .5
set x2range [0:20000]
set x2tics ('65' 6538) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.225, 0.56
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../data/simulation3/raw5/updated/perioddata_simulation3D.csv' using 1:2 axis x2y1 w l ls 1 lw .6 

set yrange [0:46]
set ytics ('38.4' 38.36) scale .5
set x2range [0:20000]
set x2tics ('69' 6923) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.625, 0.56
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../data/simulation3/raw1/updated/perioddata_simulation3.csv' using 1:2 axis x2y1 w l ls 1 lw .6 

set yrange [0:55]
set ytics ('44.6' 44.56) scale .5
set x2range [0:20000]
set x2tics ('83' 8264) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.225, 0.16
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../data/simulation3/raw3/updated/perioddata_simulation3B.csv' using 1:2 axis x2y1 w l ls 1 lw .6 

set yrange [0:82]
set ytics ('68.4' 68.43) scale .5
set x2range [0:20000]
set x2tics ('116' 11578) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.625, 0.16
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../data/simulation3/raw4/updated/perioddata_simulation3C.csv' using 1:2 axis x2y1 w l ls 1 lw .6 


unset multiplot

pause -1
