
reset

set term qt size 260,240 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [0:6]
set yrange [22.6:31.1]
set ytics 24,3,30 right scale .5
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
set key left
set key spacing .5
set ylabel offset 0,-3.4;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.818,0.86
plot '../../../../data/simulation3/raw5/new_physics/static/energy_simulation3Ds.csv' using 1:2 w l title 'data' ls 1 lw .6, \
     '' using 1:3 w l title 'KF' ls 1 lw .6 lc rgb '#FF0000' 


@BMARGIN; @LMARGIN
@XTICS; @YTICS
LABEL='II'
unset key
unset ylabel
set xtics ('0' 0, '2' 2, '4' 4, '6' 6) scale .5 offset 0,.3
set xlabel offset .4,.9
set yrange [26:31.5]
set ytics 27,2,31 right scale .5
unset ylabel
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.79,0.86
plot '../../../../data/simulation3/raw1/new_physics/static/energy_simulation3s.csv' using 1:2 w l ls 1 lw .6,\
     '' using 1:3 w l ls 1 lw .6 lc rgb '#FF0000'

unset xlabel
unset label 1
unset obj 1
unset margin
unset ytics
unset xtics
set grid x2tics,ytics

set yrange [0:58]
set ytics ('46' 46.32) scale .5
set x2range [0:20000]
set x2tics ('83' 8258) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.225, 0.56
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../../data/simulation3/raw5/new_physics/static/perioddata_simulation3Ds.csv' using 1:2 axis x2y1 w l ls 1 lw .6 

set yrange [0:60]
set ytics ('48' 48.02) scale .5
set x2range [0:20000]
set x2tics ('85' 8507) offset 0,-.5 scale .5
set size 0.3, 0.3
set origin 0.225, 0.16
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
plot '../../../../data/simulation3/raw1/new_physics/static/perioddata_simulation3s.csv' using 1:2 axis x2y1 w l ls 1 lw .6 


unset multiplot

pause -1
