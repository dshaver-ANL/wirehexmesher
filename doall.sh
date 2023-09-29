MULTI=1  #set to 1 for multiple spans

[ -e converter ] && rm converter

W2N=compile_wire2nek_single
[[ $MULTI -gt 0 ]] && W2N=compile_wire2nek_multi

cd matlab
matlab -nodisplay -r wire_mesher
mv *.out ../
cd ..
python bintoascii.py
cd wire2nek
./$W2N
cp wire2nek ../converter
cp base.rea ../
cd ..
./converter
reatore2 << EOF
wire_out
wirehex
EOF
# clean up working dir
#rm wire_out.rea
#rm wirehex.rea
#rm base.rea
#rm *.out
#rm converter
