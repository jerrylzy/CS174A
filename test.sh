#!/bin/bash

cd 'Ray Tracer'
make
make copy
cd ..

if [ $# -eq 1 ]
then
    ./raytrace "assignment3-tests-and-results/test$1.txt"
    mv "test$1.ppm" test_output/
    exit
fi

testcases=( "Sky" "Ambient" "Background" "Behind" "Diffuse" "Illum" "ImgPlane"\
	   "Intersection" "Parsing" "Reflection" "Sample" "Shadow" "Specular" )

for i in "${testcases[@]}"
do
    ./raytrace "assignment3-tests-and-results/test$i.txt"
done

./raytrace "assignment3-tests-and-results/tcb.txt"

mkdir test_output
mv *.ppm test_output/