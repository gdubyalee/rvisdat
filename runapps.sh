#!/bin/sh
cd /srv/rvisdat

#start apps
cd clapp
./app.Rscript 0.0.0.0 15123 &

cd ../pcapp
./app.Rscript 0.0.0.0 20321 &

cd ../getinit
./app.Rscript 0.0.0.0 25456 &

cd ../chk
./checkup.py &

cd ..

exit 0
