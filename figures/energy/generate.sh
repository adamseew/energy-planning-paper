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
    set xrange [0:320];
    set xtics 60 center offset -.4,.5;
    set yrange [25:45];
    set ytics 5 center offset -1.6,0;
    set grid ytics,xtics;
    set xlabel 'Time (sec)' rotate parallel center offset -1.4,.9;
    set ylabel 'Power (W)' rotate parallel center offset 1.8,-1.1;
    set output 'generated/max_qos_tees.svg';
    plot '../../data/simulation2/max_qos_tees.csv' using 1:2 w l title 'estimated' lc rgb 'black', '' using 1:3 w l title 'data' lc rgb 'red';
"


for i in generated/*.svg; do svg2tikz --figonly $i > $i.tikz; done

sed -i -- 's/y\ (m)/\\rotatebox{90}{y\ (m)}/g' generated/*.tikz
