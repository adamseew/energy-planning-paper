
reset

set term qt size 240,240 font 'Times,8'
set datafile separator comma

set style line 1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw .6
set border lw .6

unset key

set macros

set xrange [-150:250]
set yrange [-100:300]
set grid xtics ytics

NOXTICS = "set xtics ('' -150, '' -50, '' 50, '' 150, '' 250) scale .5; \
          unset xlabel"
XTICS = "set xtics -150,100 scale .5 offset 0,.3;\
          set xlabel 'x (m)'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.0f'; set ylabel 'y (m)'"

array point[1]

set multiplot layout 2,2 rowsfirst

set tmargin at screen 0.96; set bmargin at screen 0.58 
set lmargin at screen 0.16; set rmargin at screen 0.54
@NOXTICS; @YTICS
LABEL='I'
set ytics 100 right scale .5
set ylabel 'y (m)' offset 1.2,0;
set obj 1 rect at graph 0.843,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.828,0.86
set key left
set key spacing .5
plot '../../../../data/simulation3/raw5/new_physics/static/position_simulation3Ds_resized.csv' using 1:2 w l title 'path' ls 1 lw .6 

set tmargin at screen 0.96; set bmargin at screen 0.58
set lmargin at screen 0.56; set rmargin at screen 0.94
@NOXTICS; @NOYTICS
LABEL='II'
unset ylabel
set xlabel 'x (m)' offset 1.2,1.4;
set xtics -150,100,250 scale .5 offset 0,.5 scale .5
set obj 1 rect at graph 0.828,0.86 size char strlen(LABEL)+1.2, char 1
set label 1 LABEL at graph 0.795,0.86
plot '../../../../data/simulation3/raw1/new_physics/static/position_simulation3s_resized.csv' using 1:2 w l ls 1 lw .6

set tmargin at screen 0.56; set bmargin at screen 0.10; 
set lmargin at screen 0.16; set rmargin at screen 0.62
@XTICS; @YTICS
LABEL='III'
unset key
unset ylabel
unset label
unset obj
unset xlabel
set xtics -150,100,250
set ytics 100 right scale .5
set zrange [20:30]
set ztics 20,10,30 offset 1.2,0 scale .5
set zlabel 'h (m)' offset 1.1,.1 rotate by 90
set view 10,10,1
splot '../../../../data/simulation3/raw5/new_physics/dynamic/position_simulation3D_resized.csv' w l ls 1 lw .6

set tmargin at screen 0.25; set bmargin at screen 0.09; 
set lmargin at screen 0.76; set rmargin at screen 0.94
@XTICS; @NOYTICS
set ytics mirror scale .2
set xrange [0:403.01] 
set xtics ('0' 0, '2' 120, '4' 240, '6' 360) offset 0,.5 scale .2
set xlabel 'Time (min)' offset 0,1.4
unset ylabel
set y2range [-1000:0]
set y2tics ('0' 0, '-500' -500, '-1000' -1000) offset .1,0 scale .2
set key right
set key spacing .2 
set key samplen .1
plot '../../../../data/simulation3/raw5/new_physics/dynamic/ctl_simulation3D_resized.csv' using 1:2 w l axis x1y2 title 'c1' ls 1 lw .6

set tmargin at screen 0.45; set bmargin at screen 0.29; 
set lmargin at screen 0.76; set rmargin at screen 0.94
@XTICS; @NOYTICS
set ytics mirror scale .2
set xrange [0:403.01] 
set xtics ('' 0, '' 120, '' 240, '' 360) scale .2
unset xlabel
set y2range [2:10]
set y2tics 2,4,10 offset -.6,0 scale .2
plot '../../../../data/simulation3/raw5/new_physics/dynamic/ctl_simulation3D_resized.csv' using 1:3 w l axis x1y2 title 'c2' ls 1 lw .6 



pause -1

unset multiplot

