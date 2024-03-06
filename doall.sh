[ -e converter ] && rm converter

cd matlab
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
