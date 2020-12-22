#!/bin/bash

gnuplot -e "
    set datafile separator ',';
    set term qt size 180,120;
    set tics nomirror out scale 0.75;    
    unset title;
    set border lw .5;
    set xrange [0:320];
    set xtics 100 center offset -.0,.5 font 'Times,8';
    set yrange [19:31];
    set ytics 5 center offset -.3,0 font 'Times,8';
    set xlabel 'Time (sec)' rotate parallel center offset -.7,1.1 font 'Times,8';
    set ylabel 'Power (W)' rotate parallel center offset 2.5,-.1 font 'Times,8';
    set key font 'Times,8';
    plot 'periodic_energy.csv' using 1:2 w l notitle lc rgb 'black' lw .5;
    pause -1;
    set term qt size 160,120;
    unset ylabel;
    unset ytics;
    set xtics ('-2' -0.02, '0' 0, '2' 0.02) center offset -.1,.5 font 'Times,8';
    set y2tics ('0' 0, '2' 200, '4' 400) center offset -.3,0 font 'Times,8';
    set xlabel 'Frequency (cHz)' rotate parallel center offset .1,1.1 font 'Times,8';
    set y2label 'Spectrum (hdB)' rotate parallel center offset -2,-.1 font 'Times,8';
    set xrange [-0.03:0.03];
    set y2range [0:500];
    plot 'spectrum.csv' using 1:2 w l notitle lc rgb 'black' lw .5 axes x1y2;
    pause -1;
"
