#!/bin/sh

#start apps
cd clapp
./app.Rscript 0.0.0.0 15123 &

cd ../pcapp
./app.Rscript 0.0.0.0 20321 &

cd ..

exit 0