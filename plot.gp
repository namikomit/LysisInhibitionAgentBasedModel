set datafile separator ","
set title "Population Dynamics"
set xlabel "Time (minutes)"
set ylabel "Population (/ml)"
set logscale y
set format y "10^{%L}"

plot 'Culture_Figures_Paper/Btimeseries_filtered_LIC(true)_LICT(100).csv' using 1:2 with lines title 'Bacteria','Culture_Figures_Paper/Itimeseries_filtered_LIC(true)_LICT(100).csv' using 1:2 with lines title 'Infected Bacteria'