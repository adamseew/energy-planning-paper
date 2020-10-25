#!/bin/bash
rm generated/*

gnuplot -e "
    set datafile separator ',';
    set term svg size 320,320;
    set border linewidth 1.5;
    set style line 11 lc rgb '#808080' lt 1;
    set border 3 back ls 11;
    set tics nomirror out scale 0.75;
    set style line 12 lc rgb'#808080' lt 0 lw 1;
    set grid back ls 12;
    unset title;
    set xrange [-133:95];
    set xtics 20*3 center offset -.4,.5;
    set yrange [-76:231];
    set ytics 60 center offset -1.6,0;
    set grid ytics,xtics;
    set xlabel 'x (m)' rotate parallel center offset -1.4,.9;
    set ylabel 'y (m)' rotate parallel center offset 1.8,-1.1;
    set output 'generated/trajectory.svg';
    plot  '../../data/simulation2/path_sim.csv' using 1:2 w l notitle lc rgb 'black';
"


for i in generated/*.svg; do svg2tikz --figonly $i > $i.tikz; done

sed -i -- 's/y\ (m)/\\rotatebox{90}{y\ (m)}/g' generated/*.tikz
