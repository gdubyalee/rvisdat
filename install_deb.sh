#!/bin/sh
#Install dependencies, something like
apt-get -y install gcc r-base-dev git
cd /
mkdir /srv
cd srv
#Pull down program...
if [ ! -f rvisdat ]; then
  git clone https://github.com/gdubyalee/rvisdat.git
fi
cd rvisdat
./installdeps.Rscript


#install DriftR, InferCryptDrift
cd misc/CryptDrift1D
R CMD INSTALL DriftR
R CMD INSTALL InferCryptDrift
cd ../..

##modify firewall
#ufw allow 15123
#ufw allow 20321

#Edit /etc/rc.local and add startup task
sed -i 's/exit 0/\n\/srv\/rvisdat\/runapps.sh\nexit 0/' /etc/rc.local
chmod +x /etc/rc.local

exit 0
