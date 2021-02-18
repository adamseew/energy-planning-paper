
reset

set term qt size 160,200 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [0:200]
set yrange [26:42]
set ytics ('30' 30, '40' 40) right scale .5
set grid xtics ytics

NOXTICS = "set xtics ('' 0, '' 60, '' 120, '' 180) scale .5; \
          unset xlabel"
XTICS = "set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3;\
          set xlabel 'Time (min)'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'"


set multiplot layout 7,1

set tmargin at screen 0.96; set bmargin at screen 0.84
set lmargin at screen 0.16; set rmargin at screen 0.89
@NOXTICS; @YTICS
LABEL='I'
set ylabel 'Value' offset .5,-6
set key at screen 0.86,0.90 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:4 w l ls 1 lw .6 title 'a0';

set tmargin at screen 0.84; set bmargin at screen 0.72
set lmargin at screen 0.16; set rmargin at screen 0.89
@NOXTICS;
set yrange [-28:15]
set ytics ('-10' -10, '10' 10) right scale .5
unset ylabel
set key at screen 0.86,0.78 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:5 w l ls 1 lw .6 title 'a1';

set tmargin at screen 0.72; set bmargin at screen 0.60
set lmargin at screen 0.16; set rmargin at screen 0.89
@YTICS
set key at screen 0.86,0.66 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:6 w l ls 1 lw .6 title 'b1';

set tmargin at screen 0.60; set bmargin at screen 0.48
set lmargin at screen 0.16; set rmargin at screen 0.89
@YTICS
set key at screen 0.86,0.54 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:7 w l ls 1 lw .6 title 'a2';

set tmargin at screen 0.48; set bmargin at screen 0.36
set lmargin at screen 0.16; set rmargin at screen 0.89
@YTICS
set key at screen 0.86,0.42 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:8 w l ls 1 lw .6 title 'b2';

set tmargin at screen 0.36; set bmargin at screen 0.24
set lmargin at screen 0.16; set rmargin at screen 0.89
@YTICS
set yrange [-40:25]
set ytics ('-20' -20, '20' 20) right scale .5
set key at screen 0.86,0.30 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:9 w l ls 1 lw .6 title 'a3';


set tmargin at screen 0.24; set bmargin at screen 0.12
set lmargin at screen 0.16; set rmargin at screen 0.89
@XTICS; @YTICS
set xlabel offset .7,.9
set key at screen 0.86,0.18 samplen .5
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:10 w l ls 1 lw .6 title 'b3';

pause -1

unset multiplot

