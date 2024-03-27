[ -e converter ] && rm converter

rm *.out 2>/dev/null #remove old files to eliminate false success
rm wire_out 2>/dev/null

cd matlab
rm *.mat 2>/dev/null #remove old files to eliminate false success
matlab -nodisplay -r wire_mesher </dev/null
mv *.out ../
cd ..
python bintoascii.py
cd wire2nek
./compile_wire2nek
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
