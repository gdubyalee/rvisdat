#!/usr/local/bin/coffee
express=require 'express'
app=express()

app.get '/',(req,res)->
  res.send 'Hello world!'

app.use(express.static 'www')

PORT=6666
app.listen PORT, ()->
  console.log 'Listening on port '+PORT
