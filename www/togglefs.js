
$(document).ready(function(){
  $('#mcmc').click(function(){
    window.open($('#mcmc img')[0].src())
  })
  //$('#mcmc').css('resize','both')
  //$('#mcmc').css('border','1px solid black')
  //$('#mcmc img').css('max-width','100%')
  //$('#mcmc img').css('height','auto')
  //$('#mcmc img').css('display','block')
})
/*
setTimeout(function(){
  var viewFS=document.getElementById('mcmc')
  alert(viewFS)
  var isFS=0

  if(viewFS){
    viewFS.addEventListener('click',function(){
      if(!isFS){
        var docEl=viewFS//document.documentElement
        if(docEl.requestFullscreen)docEl.requestFullscreen()
        else if(docEl.msRequestFullscreen)docEl.msRequestFullscreen()
        else if(docEl.mozRequestFullScreen)docEl.mozRequestFullScreen()
        else if(docEl.webkitRequestFullScreen)docEl.webkitRequestFullScreen()
      }else{
        if(document.exitFullscreen) document.exitFullscreen()
        else if(document.msExitFullscreen) document.msExitFullscreen()
        else if(document.mozCancelFullScreen) document.mozCancelFullScreen()
        else if(document.webkitCancelFullScreen) document.webkitCancelFullScreen()
      }
      isFS=!isFS
    })
  }
},
2000
)
*/
