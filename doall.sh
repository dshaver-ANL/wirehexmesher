[ -e converter ] && rm converter

#cd matlab
#matlab -nodisplay -r wire_mesher
#mv *.out ../
#cd ..
#python bintoascii.py
cd wire2nek
./compile_wire2nek
cp wire2nek ../converter
cp base.rea ../
cd ..
./converter
