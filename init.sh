#!/bin/sh

#modify firewall
ufw allow 15000
ufw allow 20000

su george -c ./runapps.sh
