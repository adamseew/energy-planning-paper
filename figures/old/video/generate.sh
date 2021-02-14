#!/bin/bash
rm generated/*

gnuplot -e "
    set datafile separator ',';
    set term svg size 200,200;
    unset title;
    set xrange [-133:95];
    set xtics 20*3 center offset -.4,.5;
    set yrange [-76:231];
    set ytics 60 center offset -1.6,0;
    set grid ytics,xtics;
    set xlabel 'x (m)' rotate parallel center offset -1.4,.9;
    set ylabel 'y (m)' rotate parallel center offset 1.8,-1.1;
    stats '../../data/simulation2/path_sim.csv' nooutput;
    print floor(STATS_records);
    do for [i=1:4] {
        outfile = sprintf('generated/path%d.svg',i);
        set output outfile;
        set object 1 ellipse center -45,146  size 140,140 front fs empty border lc rgb 'green';
        plot  '../../data/simulation2/path_sim.csv' using 1:2 every ::1::i w l notitle lc rgb 'red';
    };
"

gnuplot -e "
    set datafile separator ',';
    set term svg size 200,200;
    unset title;
    set xrange [-133:95];
    set xtics 20*3 center offset -.4,.5;
    set yrange [-76:231];
    set ytics 60 center offset -1.6,0;
    set grid ytics,xtics;
    set xlabel 'x (m)' rotate parallel center offset -1.4,.9;
    set ylabel 'y (m)' rotate parallel center offset 1.8,-1.1;
    stats '../../data/simulation2/path_sim.csv' nooutput;
    print floor(STATS_records);
    do for [i=5:11] {
        outfile = sprintf('generated/path%d.svg',i);
        set output outfile;
        set arrow from -115, graph 0 to -115, graph 1 nohead lc rgb 'green';
        plot  '../../data/simulation2/path_sim.csv' using 1:2 every ::1::i w l notitle lc rgb 'red';
    };
"

for i in generated/*.svg; do svg2tikz --figonly $i > $i.tikz; done

sed -i -- 's/y\ (m)/\\rotatebox{90}{y\ (m)}/g' generated/*.tikz

#plot '../../data/simulation2/varphi0_gdn.csv' using 1:2:3:4 with vectors head filled lt 2 notitle;

