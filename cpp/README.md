
To compile and run:
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./assignment3
```
Note that the paths to the input and output file are hardcoded, and the expected input file is `data/WellnessCenterSama_simplified.obj`

To replicate the conversion, use `IfcConvert`, available at https://docs.ifcopenshell.org/ifcconvert/installation.html. Unzip the archive. 
```
cd path/to/unzipped-folder
./IfcConvert --use-element-names --include entities IfcWall IfcRoof IfcFooting IfcWindow IfcDoor IfcStairFlight IfcSlab IfcMember IfcPlate --exclude attribute Name "Floor:Generic 150mm - Filled:348347" attribute Name "Floor:Street:348398" -v path/to/2022020320211122Wellness center Sama.ifc path/to/WellnessCenterSama_simplified.obj"
```





