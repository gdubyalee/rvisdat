#!/usr/bin/python
import time as t
import urllib2 as ul
import os

#Poll every ten minutes forever
pollInt=600
t.sleep(pollInt)
corrRes=ul.urlopen(ul.Request('http://localhost:20321')).read()
while 1:
  t.sleep(pollInt)
  #Check that the response is what we expect
  res=ul.urlopen(ul.Request('http://localhost:20321')).read()
  if res!=corrRes:
    with open('broken.log','w') as f:
      f.write('\n\nCrash at time '+t.strftime('%Y-%m-%d %H:%M:%S :',t.gmtime()))
      os.system('ps aux --sort rss>>broken.log')
      f.write('\n'+res)
    os.system('reboot')

