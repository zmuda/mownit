set terminal svg
set output"output.svg"
plot 'outEuler.tmp' using 1:2 with lines lw 2, 'outGSL.tmp' using 1:2 with lines lw 1
