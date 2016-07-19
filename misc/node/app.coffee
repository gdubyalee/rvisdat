#!/usr/local/bin/coffee
require 'shelljs/global'
express=require 'express'
app=express()

app.get '/',(req,res)->
  res.redirect 'main.htm'

app.get '/r',(req,res)->
  res.send ls('R')

app.use(express.static 'www')

PORT=6666
app.listen PORT, ()->
  console.log 'Listening on port '+PORT
