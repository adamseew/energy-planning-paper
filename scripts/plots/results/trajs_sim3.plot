
reset

set term qt size 260,260 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [-150:250]
set yrange [-100:300]
set ytics 100 right scale .5
set grid xtics ytics

NOXTICS = "set xtics ('' -150, '' -50, '' 50, '' 150, '' 250) scale .5; \
          unset xlabel"
XTICS = "set xtics -150,100 scale .5 offset 0,.3;\
          set xlabel 'x (m)'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'y (m)'"

TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.58"
BMARGIN = "set tmargin at screen 0.56; set bmargin at screen 0.18"
LMARGIN = "set lmargin at screen 0.16; set rmargin at screen 0.54"
RMARGIN = "set lmargin at screen 0.56; set rmargin at screen 0.94"

array point[1]

set multiplot layout 2,2 rowsfirst

@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
LABEL='I'
set ylabel offset 0,-3.4;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.828,0.86
plot '../../../data/simulation3/raw5/position_simulation3D_resized.csv' using 1:2 w l ls 1 lw .6, point us (-59.824):(12.337) pt 1 lc rgb '#FF0000'

@TMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
LABEL='II'
set ytics -100,100,200 scale .5
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.795,0.86
set key left
plot '../../../data/simulation3/raw1/position_simulation3_resized.csv' using 1:2 w l title 'path' ls 1 lw .6, point us (-65.032):(104.03) pt 1 lc rgb '#FF0000' notitle

@BMARGIN; @LMARGIN
@XTICS; @YTICS
LABEL='III'
unset key
unset ylabel
set xtics -150,100,150 scale .5 offset 0,.3
set xlabel offset 10.2,.9
set obj 1 rect at graph 0.88,0.86 size char strlen(LABEL)+.7, char 1
set label 1 LABEL at graph 0.83,0.86
plot '../../../data/simulation3/raw3/position_simulation3B_resized.csv' using 1:2 w l ls 1 lw .6, point us (35.108):(-53.611) pt 1 lc rgb '#FF0000'

@BMARGIN; @RMARGIN
@XTICS; @NOYTICS
LABEL='IV'
set xtics -150,100,250 scale .5 offset 0,.3
unset xlabel
set obj 1 rect at graph 0.858,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.81,0.86 font ',8'
plot '../../../data/simulation3/raw4/position_simulation3C_resized.csv' using 1:2 w l ls 1 lw .6, point us (-32.82):(228.46) pt 1 lc rgb '#FF0000'

unset multiplot

pause -1
