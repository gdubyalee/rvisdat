#!/bin/sh

#modify firewall
ufw allow 15123
ufw allow 20321

#Create startup task (on osx)
cp rvisdatstartup.plist /Library/LaunchAgents/rvisdatstartup.plist

su george -c ./runapps.sh
