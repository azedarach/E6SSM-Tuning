set terminal pngcairo enhanced dashed

filename="test_rge_run.txt"

set output 'essm_soft_higgsmass_rge_flow.png'

set title "ESSM RG Flow of Soft Higgs Masses"
set ylabel "m^2 (GeV^2)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(5):NaN) t "m_{Hu}^2" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(5):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(6):NaN) t "m_{Hd}^2" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(6):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) == 0 ? column(7):NaN) t "m_s^2" with lines lt 1 lw 1.7 lc rgb 'green', \
filename using (log10(column(1))):(column(12) != 0 ? column(7):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'green'

reset

set terminal pngcairo enhanced dashed

set output 'essm_soft_squarkmass_rge_flow.png'

set title "ESSM RG Flow of Soft Squark Masses"
set ylabel "m^2 (GeV^2)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(8):NaN) t "m_{qL3}^2" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(8):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(9):NaN) t "m_{tR}^2" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(9):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue'

reset 

set terminal pngcairo enhanced dashed

set output 'essm_stopmass_rge_flow.png'

set title "ESSM RG Flow of Stop Masses"
set ylabel "m (GeV)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(10):NaN) t "| {t_1} |" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(10):NaN) notitle with lines lt 2 lc rgb 'red' lw 1.7, \
filename using (log10(column(1))):(column(12) == 0 ? column(11):NaN) t "| {t_2} |" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(11):NaN) notitle with lines lt 2 lc rgb 'blue' lw 1.7 

reset

set terminal pngcairo enhanced dashed

set output 'essm_soft_gauginomass_rge_flow.png'

set title "ESSM RG Flow of Soft Gaugino Masses"
set ylabel "M (GeV)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(13):NaN) t "M_1" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(13):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(14):NaN) t "M_2" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(14):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) == 0 ? column(15):NaN) t "M_3" with lines lt 1 lw 1.7 lc rgb 'green', \
filename using (log10(column(1))):(column(12) != 0 ? column(15):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'green', \
filename using (log10(column(1))):(column(12) == 0 ? column(16):NaN) t "M_1'" with lines lt 1 lw 1.7 lc rgb 'yellow', \
filename using (log10(column(1))):(column(12) != 0 ? column(16):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'yellow'

reset 

set terminal pngcairo enhanced dashed

set output 'essm_soft_Aterm_rge_flow.png'

set title "ESSM RG Flow of Soft A Terms"
set ylabel "A (GeV)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(3):NaN) t "A_{/Symbol l}" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(3):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(4):NaN) t "A_t" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(4):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue'

reset

set terminal pngcairo enhanced dashed

set output 'essm_ewsbconditions_rge_flow.png'

set title "ESSM RG Flow of EWSB Conditions"
set ylabel "f (GeV^2)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(17):NaN) t "f_1" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(17):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(18):NaN) t "f_2" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(18):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) == 0 ? column(19):NaN) t "f_3" with lines lt 1 lw 1.7 lc rgb 'green', \
filename using (log10(column(1))):(column(12) != 0 ? column(19):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'green', \

reset 

set terminal pngcairo enhanced dashed

set output 'essm_msusycondition_rge_flow.png'

set title "ESSM RG Flow of M_{SUSY} Condition"
set ylabel "f_4"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(20):NaN) notitle with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(20):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red'

reset

set terminal pngcairo enhanced dashed

set output 'essm_ewsbconditiondiffs_rge_flow.png'

set title "ESSM RG Flow of EWSB Conditions"
set ylabel "{/Symbol D}f (GeV^2)"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? (column(17)-column(18)):NaN) t "f_1-f_2" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? (column(17)-column(18)):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? (column(17)-column(19)):NaN) t "f_1-f_3" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? (column(17)-column(19)):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) == 0 ? (column(18)-column(19)):NaN) t "f_2-f_3" with lines lt 1 lw 1.7 lc rgb 'green', \
filename using (log10(column(1))):(column(12) != 0 ? (column(18)-column(19)):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'green'

reset 

set terminal pngcairo enhanced dashed

set output 'essm_singletyukawas_rge_flow.png'

set title "ESSM RG Flow of Singlet Yukawas"
set xlabel "log_{10}(Q/GeV)"
set grid

plot filename using (log10(column(1))):(column(12) == 0 ? column(2):NaN) t "{/Symbol l}" with lines lt 1 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) != 0 ? column(2):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'red', \
filename using (log10(column(1))):(column(12) == 0 ? column(21):NaN) t "{/Symbol k}" with lines lt 1 lw 1.7 lc rgb 'blue', \
filename using (log10(column(1))):(column(12) != 0 ? column(21):NaN) notitle with lines lt 2 lw 1.7 lc rgb 'blue'