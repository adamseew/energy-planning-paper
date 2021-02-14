
reset

set term qt size 280,260 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set style line 2 lc rgb '#000000' pt 5 ps 0.2 lt 1 dt 2 lw .6
set border lw .6

unset key

set macros

set xrange [200:400]
set yrange [21:34]
set ytics 22,3,31 right scale .5
set grid xtics ytics

set multiplot layout 2,2 rowsfirst

set tmargin at screen 0.96; set bmargin at screen 0.68; 
set lmargin at screen 0.20; set rmargin at screen 0.54
set xtics ('4' 240, '5' 300, '6' 360) scale .5 offset 0,.3
set format y '%.0f'; set ylabel 'Power (W)'
LABEL='I'
set key left
set key spacing .5
set ylabel offset 0,-4.4;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.828,0.86
plot '../../../../data/simulation3/raw5/new_physics/static/energy_simulation3Ds_resized.csv' using 1:2 title 'data' w l ls 1 lw .6,\
     '../../../../data/simulation3/raw5/new_physics/static/evol_simulation3Ds_resized.csv' using 1:2 title 'mod' w l ls 1 lw .6 lc rgb '#FF0000'

set tmargin at screen 0.96; set bmargin at screen 0.68; 
set lmargin at screen 0.56; set rmargin at screen 0.90
set xtics ('4' 240, '5' 300, '6' 360) scale .5 offset 0,.3
unset xlabel 
set format y ''; unset ylabel
LABEL='II'
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.80,0.86
unset key
plot '../../../../data/simulation3/raw1/new_physics/static/energy_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6,\
     '../../../../data/simulation3/raw1/new_physics/static/evol_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6 lc rgb '#FF0000

set xrange [0:403.1]
set tmargin at screen 0.60; set bmargin at screen 0.18; 
set lmargin at screen 0.20; set rmargin at screen 0.90
set xtics ('0' 0, '1' 60, '2' 120, '3' 180) scale .5 offset 0,.3
set xlabel 'Time (min)' 
set format y '%.0f'; set ylabel 'Power (W)'
unset obj 1
unset label 1
unset ylabel
set xtics ('0' 0, '1' 60, '2' 120, '3' 180, '4' 240, '5' 300, '6' 360) scale .5 offset 0,.3
set xlabel offset 0.4,.9
set yrange [26.5:41]
set ytics 27,4,41 scale .5

plot '../../../../data/simulation3/raw5/new_physics/dynamic/energy_simulation3D_resized.csv' using 1:2 w l ls 1 lw .6

unset xlabel
unset label 1
unset obj 1
unset margin

set yrange [28:40]
set xrange [90:170]
set size 0.42, 0.54
set origin 0.61, 0.18
set xtics ('' 120, '' 180)
set ytics ('' 31, '' 35, '' 39)
set style rect fc lt -1 fs solid 0.15 noborder
set obj 1 rect from graph 0.0,0.0 to graph 1,1 behind fc rgb "#FFFFFF" 
set key right samplen .5
set key spacing .5
plot '../../../../data/simulation3/raw5/new_physics/dynamic/energy_simulation3D.csv' using 1:2 w l title 'data' ls 1 lw .6,\
     '../../../../data/simulation3/raw5/new_physics/dynamic/bat_simulation3D.csv' using 1:3 w l title 'bat' ls 2 lw .6

pause -1 

unset multiplot
