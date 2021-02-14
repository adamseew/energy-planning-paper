#!/bin/bash

gnuplot -e "
    set datafile separator ',';
    set term qt size 280,280;
    set border linewidth 1.5;
    set style line 11 lc rgb '#808080' lt 1;
    set border 3 back ls 11;
    set tics nomirror out scale 0.75;
    set style line 12 lc rgb'#808080' lt 0 lw 1;
    set grid back ls 12;
    unset title;
    set xrange [-133:95];
    set xtics 20*3 center offset .6,.5 font 'Times,8';
    set yrange [-76:260];
    set ytics 60 center offset -.3,0 font 'Times,8';
    set grid ytics,xtics;
    set xlabel 'x (m)' rotate parallel center offset -1.4,.9 font 'Times,11';
    set ylabel 'y (m)' rotate parallel center offset 1.8,-1.1 font 'Times,11';    set key font 'Times,8';
    plot  '../../data/simulation2/path_sim.csv' using 1:2 w l title 'Trajectory' lc rgb 'black';
    pause -1;
"

#for i in generated/*.svg; do svg2tikz --figonly $i > $i.tikz; done

#sed -i -- 's/y\ (m)/\\rotatebox{90}{y\ (m)}/g' generated/*.tikz
