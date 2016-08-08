#!/bin/sh

cd clapp
./app.Rscript 0.0.0.0 15000 &

cd ../pcapp
./app.Rscript 0.0.0.0 20000 &

cd ..
exit 0
