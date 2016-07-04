#!/bin/bash
declare -a nt=("0" "1" "2" "3" "4" "5")
for i in {0,1,2,3,4,5}
do
pdfcrop deformingAdvection_onOrthogW_c1_480x240_${nt[i]}_T.pdf
rm deformingAdvection_onOrthogW_c1_480x240_${nt[i]}_T.pdf
done
for i in {0,1,2,3,4,5}
do
mv deformingAdvection_onOrthogW_c1_480x240_${nt[i]}_T-crop.pdf deformingAdvection_onOrthogW_c1_480x240_${nt[i]}_T.pdf
done