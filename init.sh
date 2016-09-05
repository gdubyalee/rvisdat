#!/bin/sh

#modify firewall
ufw allow 15123
ufw allow 20321

#Create startup task (on osx)
cp osxlaunchd.xml /System/Library/LaunchDaemons

su george -c ./runapps.sh
