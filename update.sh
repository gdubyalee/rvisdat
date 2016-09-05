#!/bin/sh
cd /srv/rvisdat
git pull
echo 'New version will be run on next startup.'
echo 'Alternatively kill the running scripts and run /srv/rvisdat/runapps.sh'
exit 0
