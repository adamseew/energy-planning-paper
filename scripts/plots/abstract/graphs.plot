
reset

set term qt size 180,120 font 'Times,8'

set datafile separator ','
set tics nomirror scale 0.5    
unset title

set border lw .5

set xrange [0:320]

set xtics 100 center offset -.5,.5

set yrange [19:31]
set ytics 5 center offset -.9,0

set grid xtics ytics

set xlabel 'Time (sec)' rotate parallel center offset -.7,1.1
set ylabel 'Power (W)' rotate parallel center offset 2.5,-.1

plot '../periodic_energy.csv' using 1:2 w l notitle lc rgb 'black' lw .5

pause -1

set term qt size 160,120 font 'Times,8'

unset ylabel
unset ytics

set xtics ('-2' -0.02, '0' 0, '2' 0.02) center offset -.1,.5
set y2tics ('0' 0, '2' 200, '4' 400) center offset .3,0

set grid xtics y2tics

set xlabel 'Frequency (cHz)' rotate parallel center offset .1,1.1
set y2label 'Spectrum (hdB)' rotate parallel center offset -2,-.1
set xrange [-0.03:0.03]
set y2range [0:750]

plot '../spectrum2.csv' using 1:2 w l notitle lc rgb 'black' lw .5 axes x1y2

pause -1

set y2tics ('4000' 400000) center offset 1.8,0
set y2range [750:440000]

plot '../spectrum2.csv' using 1:2 w l notitle lc rgb 'black' lw .5 axes x1y2

pause -1
