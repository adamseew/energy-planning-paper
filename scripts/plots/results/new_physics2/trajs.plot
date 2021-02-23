
reset

set term qt size 280,340 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

array point[1]

set multiplot layout 3,2 rowsfirst

set xrange [-150:250]
set yrange [-100:300]
set grid xtics ytics
set tmargin at screen 0.96; set bmargin at screen 0.70 
set lmargin at screen 0.16; set rmargin at screen 0.51
set xtics ('' -150, '' -50, '' 50, '' 150, '' 250) scale .5
LABEL='I'
set ytics ('-100' -100, '0' 0, '100' 100, '200' 200, '300' 300) right scale .5
set ylabel 'y (m)' offset 1.2,-7;
set obj 1 rect at graph 0.853,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.83,0.86
set key left
set key spacing .5
plot '../../../../data/simulation3/raw5/new_physics/static2/position_simulation3Ds_resized.csv' using 1:2 w l title 'path' ls 1 lw .6 

set tmargin at screen 0.96; set bmargin at screen 0.70
set lmargin at screen 0.53; set rmargin at screen 0.88
LABEL='II'
unset ylabel
set ytics ('' -100, '' 0, '' 100, '' 200, '' 300) right scale .5
set obj 1 rect at graph 0.838,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.795,0.86
plot '../../../../data/simulation3/raw1/new_physics/static2/position_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6



set tmargin at screen 0.68; set bmargin at screen 0.38; 
set lmargin at screen 0.16; set rmargin at screen 0.64
set xtics ('' -150, '' -50, '' 50, '' 150, '' 250) scale .5
unset key
unset ylabel
LABEL='I'
set obj 1 rect at graph 0.883,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.868,0.86
unset xlabel 
set yrange [-100:250]
set ytics -100,100,200 right scale .5
plot '../../../../data/simulation3/raw5/new_physics/dynamic_revised/position_simulation3D_resized.csv' using 1:2 w l ls 1 lw .6

unset label
unset obj
set tmargin at screen 0.68; set bmargin at screen 0.54; 
set lmargin at screen 0.66; set rmargin at screen 0.88
set yrange [1:11]
set ytics ('' 2, '' 6, '' 10) mirror scale .2
set format y ''
unset ylabel
set xrange [0:403.01] 
set xtics ('' 0, '' 120, '' 240, '' 360) scale .2
unset xlabel
set y2range [1:11]
set y2tics 2,4,10 offset -.6,0 scale .2
set y2label 'cn' offset -2,0
plot '../../../../data/simulation3/raw5/new_physics/dynamic_revised/ctl_simulation3D_resized.csv' using 1:3 w l axis x1y2 ls 1 lw .6 

set tmargin at screen 0.54; set bmargin at screen 0.40; 
set lmargin at screen 0.66; set rmargin at screen 0.88
set yrange [-1100:100]
set ytics ('' 0, '' -500, '' -1000) offset .1,0 mirror scale .2
set xrange [0:403.01] 
set xtics ('0' 0, '2' 120, '4' 240, '6' 360) offset 0,.3 scale .2
unset xlabel
unset ylabel
set y2range [-1100:100]
set y2tics ('0' 0, '-.5' -500, '-1' -1000) offset .1,0 scale .2
set y2label 'cn' offset -2,0
plot '../../../../data/simulation3/raw5/new_physics/dynamic_revised/ctl_simulation3D_resized.csv' using 1:2 w l axis x1y2 ls 1 lw .6



unset y2tics
unset y2range
unset y2label
set xrange [-150:250]
set tmargin at screen 0.36; set bmargin at screen 0.06; 
set lmargin at screen 0.16; set rmargin at screen 0.64
unset key
unset ylabel
LABEL='II'
set obj 1 rect at graph 0.888,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.865,0.86
set xtics -150,100,150 scale .5
set xlabel 'x (m)' offset .2,.9;
set yrange [-100:250]
set ytics ('-100' -100, '0' 0, '100' 100, '200' 200) right scale .5
plot '../../../../data/simulation3/raw1/new_physics/dynamic_revised/position_simulation3_resized.csv' using 1:2 w l ls 1 lw .6

set tmargin at screen 0.34; set bmargin at screen 0.20; 
set lmargin at screen 0.66; set rmargin at screen 0.88
set yrange [1:11]
set ytics ('' 2, '' 6, '' 10) mirror scale .2
set format y ''
unset label
unset obj
unset ylabel
set xrange [0:661.99] 
set xtics ('' 0, '' 180, '' 360, '' 540) offset 0,.3 scale .2
unset xlabel
set y2range [1:11]
set y2tics 2,4,10 offset -.6,0 scale .2
set y2label 'cn' offset -2,0
plot '../../../../data/simulation3/raw1/new_physics/dynamic_revised/ctl_simulation3_resized.csv' using 1:3 w l axis x1y2 ls 1 lw .6 

set tmargin at screen 0.20; set bmargin at screen 0.06; 
set lmargin at screen 0.66; set rmargin at screen 0.88
set ytics mirror scale .2
set yrange [-1100:100]
set ytics ('' 0, '' -500, '' -1000) offset .1,0 scale .2
set xrange [0:661.99] 
set xtics ('0' 0, '3' 180, '6' 360, '9' 540) offset 0,.3 scale .2
set xlabel 'Time (min)' offset 0,.9
unset ylabel
set y2range [-1100:100]
set y2tics ('0' 0, '-.5' -500, '-1' -1000) offset .1,0 scale .2
set y2label 'cn' offset -2,0
plot '../../../../data/simulation3/raw1/new_physics/dynamic_revised/ctl_simulation3_resized.csv' using 1:2 w l axis x1y2 ls 1 lw .6



pause -1

unset multiplot

