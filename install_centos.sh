#!/bin/sh
#Install dependencies, something like
yum groupinstall -y 'development tools'
yum install -y wget which libcurl-devel openssl-devel git
rpm -Uvh https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum install -y R

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

#modify firewall
firewall-cmd --zone=public --add-port=15123/tcp --permanent
firewall-cmd --zone=public --add-port=20321/tcp --permanent
firewall-cmd --reload

#Edit /etc/rc.local and add startup task
sed -i -z 's/exit 0/\n\/srv\/rvisdat\/runapps.sh\nexit 0/'

