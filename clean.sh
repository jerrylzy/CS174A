#!/bin/bash

rm -f raytrace  *.ppm
find . -type f -name "*~" -delete

cd test_output
rm *.ppm