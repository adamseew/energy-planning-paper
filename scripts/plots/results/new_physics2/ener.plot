
reset

set term qt size 280,340 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set style line 2 lc rgb '#00FF00' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [0:200]
set yrange [31:44]
set ytics 32,4,42 right scale .5
set grid xtics ytics

set multiplot layout 2,2 rowsfirst

set tmargin at screen 0.96; set bmargin at screen 0.72; 
set lmargin at screen 0.20; set rmargin at screen 0.51
set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3
set format y '%.0f'
set ylabel 'Power (W)'
LABEL='I'
set key left
set key spacing .5
set ylabel offset 0,-6.8;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.828,0.86
plot '../../../../data/simulation3/raw5/new_physics/static2/energy_simulation3Ds_resized.csv' using 1:2 title 'data' w l ls 1 lw .6,\
     '../../../../data/simulation3/raw5/new_physics/static2/evol_simulation3Ds_resized.csv' using 1:2 title 'mod' w l ls 1 lw .6 lc rgb '#FF0000'

set tmargin at screen 0.96; set bmargin at screen 0.72; 
set lmargin at screen 0.59; set rmargin at screen 0.90
unset xlabel 
unset ylabel
set xrange [0:200]
set yrange [26:39]
set ytics 27,4,37 right scale .5
LABEL='II'
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.792,0.86
unset key
plot '../../../../data/simulation3/raw1/new_physics/static2/energy_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6,\
     '../../../../data/simulation3/raw1/new_physics/static2/evol_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6 lc rgb '#FF0000

set xrange [0:403.1]
set tmargin at screen 0.66; set bmargin at screen 0.41; 
set lmargin at screen 0.20; set rmargin at screen 0.90
set format y '%.0f'
unset obj 1
unset label 1
unset ylabel
set xtics ('0' 0, '1' 60, '2' 120, '3' 180, '4' 240, '5' 300, '6' 360) scale .5 offset 0,.3
set yrange [25.9:41]
set ytics 27,4,41 scale .5
set key left bottom
set key spacing .5
LABEL='I'
set key left
set obj 1 rect at graph 0.923,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.914,0.86
plot '../../../../data/simulation3/raw5/new_physics/dynamic/energy_simulation3D_resized.csv' using 1:2 w l notitle ls 1 lw .6,\
     '../../../../data/simulation3/raw5/new_physics/dynamic/bat_simulation3D.csv' using 1:3 w l title 'bat' ls 2 lw .6

set xrange [0:661.99]
set tmargin at screen 0.35; set bmargin at screen 0.10; 
set lmargin at screen 0.20; set rmargin at screen 0.90
set xlabel 'Time (min)' 
set format y '%.0f'
unset obj 1
unset label 1
unset ylabel
set xtics ('0' 0, '2' 120, '4' 240, '6' 360, '8' 480, '10' 600) scale .5 offset 0,.3
set xlabel offset 0.4,.9
set yrange [25.9:41]
set ytics 27,4,41 scale .5
LABEL='II'
set obj 1 rect at graph 0.908,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.895,0.86
plot '../../../../data/simulation3/raw1/new_physics/dynamic2/energy_simulation3_resized.csv' using 1:2 w l ls 1 lw .6,\
     '../../../../data/simulation3/raw1/new_physics/dynamic2/bat_simulation3.csv' using 1:3 w l title 'bat' ls 2 lw .6

pause -1 

unset multiplot
