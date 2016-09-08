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
iptables -A INPUT -m state --state NEW -p tcp --dport 15123 -j ACCEPT
iptables -A INPUT -m state --state NEW -p tcp --dport 20321 -j ACCEPT
/etc/init.d/iptables restart

#Edit /etc/rc.local and add startup task
echo '\n\n\/srv\/rvisdat\/runapps.sh\n\nexit 0\n' >>/etc/rc.d/rc.local
chmod +x /etc/rc.d/rc.local

