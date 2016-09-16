#!/usr/bin/gnuplot

# set output
set output 'Populations_All_16_High_Pt.eps'


set terminal postscript enhanced color 'Helvetica' 20

set xlabel 'Time (fs)'
set ylabel 'Change in occupation'


plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):16  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):16 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):16 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):16 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):16 w l lw 2 title '121 331 hi'


set output 'Populations_All_15_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):15  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):15 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):15 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):15 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):15 w l lw 2 title '121 331 hi'


set output 'Populations_All_14_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):14  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):14 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):14 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):14 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):14 w l lw 2 title '121 331 hi'

set output 'Populations_All_13_High_Pt.eps'
plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):13  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):13 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):13 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):13 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):13 w l lw 2 title '121 331 hi'

set output 'Populations_All_12_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):12  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):12 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):12 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):12 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):12 w l lw 2 title '121 331 hi'


set output 'Populations_All_11_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):11  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):11 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):11 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):11 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):11 w l lw 2 title '121 331 hi'

set output 'Populations_All_10_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):10  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):10 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):10 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):10 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):10 w l lw 2 title '121 331 hi'

set output 'Populations_All_9_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):9  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):9 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):9 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):9 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):9 w l lw 2 title '121 331 hi'

set output 'Populations_All_8_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):8  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):8 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):8 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):8 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):8 w l lw 2 title '121 331 hi'

set output 'Populations_All_7_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):7  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):7 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):7 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):7 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):7 w l lw 2 title '121 331 hi'


set output 'Populations_All_6_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):6  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):6 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):6 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):6 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):6 w l lw 2 title '121 331 hi'

set output 'Populations_All_5_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):5  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):5 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):5 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):5 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):5 w l lw 2 title '121 331 hi'


set output 'Populations_All_4_High_Pt.eps'

plot 'Population_Trace_118.1_0_322_Pt.txt' u ($1*0.00228789):4  w l lw 2 title '118 322 hi',\
'Population_Trace_119.1_0_324_Pt.txt' u ($1*0.00228789):4 w l  lw 2 title '119 324 hi',\
'Population_Trace_120.1_0_326_Pt.txt' u ($1*0.00228789):4 w l lw 2 title '120 326 hi',\
'Population_Trace_120.1_0_329_Pt.txt' u ($1*0.00228789):4 w l lw 2 title '120 329 hi',\
'Population_Trace_121.1_0_331_Pt.txt' u ($1*0.00228789):4 w l lw 2 title '121 331 hi'

