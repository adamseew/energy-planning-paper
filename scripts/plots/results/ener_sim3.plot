
reset

set term qt size 340,260 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [0:200]
set yrange [29.3:41]
set ytics 30,4,38 right scale .5
set grid xtics ytics

NOXTICS = "set xtics ('' 0, '' 60, '' 120, '' 180) scale .5; \
          unset xlabel"
XTICS = "set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3;\
          set xlabel 'Time (min)'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'Power (W)'"

TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.58"
BMARGIN = "set tmargin at screen 0.56; set bmargin at screen 0.18"
LMARGIN = "set lmargin at screen 0.16; set rmargin at screen 0.54"
RMARGIN = "set lmargin at screen 0.56; set rmargin at screen 0.94"


set multiplot layout 2,2 rowsfirst

@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
LABEL='I'
set ylabel offset 0,-3.4;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.828,0.86
plot '../../../data/simulation3/raw5/energy_simulation3D_resized.csv' using 1:2 w l ls 1 lw .6

@TMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
LABEL='II'
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.80,0.86
set key left
plot '../../../data/simulation3/raw1/energy_simulation3_resized.csv' using 1:2 w l title 'energy' ls 1 lw .6

@BMARGIN; @LMARGIN
@XTICS; @YTICS
LABEL='III'
unset key
unset ylabel
set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3
set xlabel offset 12.4,.9
set yrange [28:62]
set ytics 30,8,54 scale .5
set obj 1 rect at graph 0.88,0.86 size char strlen(LABEL)+.7, char 1
set label 1 LABEL at graph 0.83,0.86
plot '../../../data/simulation3/raw3/energy_simulation3B_resized.csv' using 1:2 w l ls 1 lw .6

@BMARGIN; @RMARGIN
@XTICS; @NOYTICS
LABEL='IV'
set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3
unset xlabel
set obj 1 rect at graph 0.858,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.818,0.86 font ',8'
plot '../../../data/simulation3/raw4/energy_simulation3C_resized.csv' using 1:2 w l ls 1 lw .6

unset multiplot

pause -1
