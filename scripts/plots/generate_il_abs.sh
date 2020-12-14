#!/bin/bash

gnuplot -e "
    set datafile separator ',';
    set term qt size 180,120;
    set tics nomirror out scale 0.75;    
    unset title;
    set border lw .5;
    set xrange [0:320];
    set xtics 100 center offset .6,.5 font 'Times,8';
    set yrange [19:41];
    set ytics 5 center offset -.3,0 font 'Times,8';
    set xlabel 'Time (sec)' rotate parallel center offset -.7,1.1 font 'Times,8';
    set ylabel 'Power (W)' rotate parallel center offset 2.5,-.1 font 'Times,8';
    set key font 'Times,8';
    set multiplot;
    plot 'periodic_energy.csv' using 1:2 w l notitle lc rgb 'black' lw .5;
    set size 0.5, 0.5;
    set origin 0.5, 0.5;
    unset xtics;
    unset ytics;
    unset xlabel;
    unset ylabel;
    set xrange [-0.03:0.03];
    set yrange [0:500];
    plot 'spectrum.csv' using 1:2 w l notitle lc rgb 'black' lw .5;
    unset multiplot;
    pause -1;
"
