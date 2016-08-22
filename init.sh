#!/bin/sh

#modify firewall
ufw allow 15123
ufw allow 20321

su george -c ./runapps.sh
