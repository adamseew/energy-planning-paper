#!/bin/bash

gnuplot -e "
    set datafile separator ',';
    set term qt size 320,220;
    set border linewidth 1.5;
    set style line 11 lc rgb '#808080' lt 1;
    set border 3 back ls 11;
    set tics nomirror out scale 0.75;
    set style line 12 lc rgb'#808080' lt 0 lw 1;
    set grid back ls 12;
    unset title;
    set xrange [0:320];
    set xtics 60 center offset .6,.5 font 'Times,11';
    set yrange [29:44];
    set ytics 3 center offset -.3,0 font 'Times,11';
    set grid ytics, xtics;
    set xlabel 'Time (sec)' rotate parallel center offset -1.4,.9 font 'Times,11';
    set ylabel 'Power (W)' rotate parallel center offset .5,-.1 font 'Times,11';
    set key font 'Times,8';
    plot '../../data/simulation2/max_qos_tees.csv' using 1:2 w l title 'Data' lc rgb 'black', '' using 1:3 w l title 'Estimated' lc rgb 'red' lw .5;
    pause -1;
"

gnuplot -e "
    set datafile separator ',';
    set term qt size 320,220;
    set border linewidth 1.5;
    set style line 11 lc rgb '#808080' lt 1;
    set border 3 back ls 11;
    set tics nomirror out scale 0.75;
    set style line 12 lc rgb'#808080' lt 0 lw 1;
    set grid back ls 12;
    unset title;
    set xrange [0:320];
    set xtics 60 center offset .6,.5 font 'Times,11';
    set yrange [18000:25000];
    set ytics ('18' 18000, '20' 20000, '22' 22000, '24' 24000) center offset -.3,0 font 'Times,11';
    set grid ytics, xtics;
    set xlabel 'Time (sec)' rotate parallel center offset -1.4,.9 font 'Times,11';
    set ylabel 'Energy (kJ)' rotate parallel center offset .5,-.1 font 'Times,11';
    set key font 'Times,8';
    set style fill solid;
    set boxwidth 1;
    plot '../../data/simulation2/est_vs_joules.csv' using 1:2 w impulses title 'Estimates' lc rgb 'black', '' using 1:3 w l title 'True value' lc rgb 'red' lw .5;
    pause -1;
"

#2 graph on top each other and a graph at begin
#in deep anaysis in section 4

#just skype the comutational energy sans says how you can vary 

#- graph as it is
#- graph with what happpens to figure 1 if from one point on th efigure is not anymore periodic.
#- graph actually adapting

#plot under the actual plot of energy joules rather than the power to show what happend with the predicted data rather than the actual
