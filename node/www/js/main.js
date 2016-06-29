var myFun=(x)=>{
  return x*x
}
var myCube=(x)=>{
  return x*x*x
}

var xVals=[0,1,2,3,4,5]
var f=[]
xVals.forEach((x)=>{f.push({x:x,y:myFun(x)})})

var cv1=$('#chart1'),cv2=$('#chart2')
var chart1,chart2
//graphF =[{x:...,y:...},...]
drawChart=(canvas,graphF,label)=>{
  return(new Chart(
    canvas,
    {
      type:'line',
      data:{
        datasets:[{
          label:label,
          data:graphF 
        }]
      },
      options:{
        scales:{
          xAxes:[{
            type:'linear',
            position:'bottom'
          }]
        }
      }
    }
  ))
}

updateChartVals=(chart,newVals,datasetIndex)=>{
  //Coerce dataset index to number, 0 by default
  datasetIndex|=0
  if(chart.data.datasets[datasetIndex])chart.data.datasets[datasetIndex].data=newVals
  else chart.data.datasets[datasetIndex]={data:newVals}
  chart.update()
}

var params={}
updateParam=(param)=>{
  params[param]=parseFloat($('#'+param).val())
}


//Evaluate analysic solution
analyticNeutralDriftSolution=(N,lambda,tau,n,t)=>{
  var ret=0
  var angle
  if(t<tau)return(0);
  else t-=tau;
  if(n===0){
    for(i=1;i<N;i++){
      angle=.5*Math.PI*i/N
      ret+=Math.pow(Math.cos(angle),2)*Math.exp(
        -4*lambda*Math.pow(.5*Math.sin(angle),2)*t
      )
    }
  }else if(n===N){
    for(i=1;i<N;i++){
      angle=.5*Math.PI*i/N
      ret+=Math.pow(Math.cos(angle),2)*(2*(i%2)-1)*Math.exp(
        -4*lambda*Math.pow(.5*Math.sin(angle),2)*t
      )
    }
  }else{
    for(i=1;i<N;i++){
      angle=.5*Math.PI*i/N
      ret+=Math.sin(n*angle)*Math.sin(angle)*Math.exp(
        -4*lambda*Math.pow(.5*Math.sin(.5*angle),2)*t
      )
    }
  }
  return((2/N)*ret)
}

analyticSolutionGraphPoints=(N,lambda,tau,n,t,min,max,h)=>{
}

//Testing stuff

var pulseChaseChart=drawChart(cv1,f,'A nice graph')

f[2]={x:2,y:1}
setTimeout(()=>{
  updateChartVals(pulseChaseChart,f)
},5000)

var f2=[]
xVals.forEach((x)=>{f2.push({x:x,y:myCube(x)})})

setTimeout(()=>{updateChartVals(pulseChaseChart,f2,1)},6000)
