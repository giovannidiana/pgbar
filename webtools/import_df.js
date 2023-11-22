const inputElement = document.getElementById("file");
inputElement.addEventListener("change", handleFiles, false);

function handleFiles(event) {
    var reader = new FileReader();
    reader.onload = onReaderLoad;
    reader.readAsText(event.target.files[0]);
}

function onReaderLoad(event){
    var obj = JSON.parse(event.target.result);
    update_plot(obj);
}

d3.select("body")
    .on("keypress", function(event){
        if(event.keyCode==102){
            console.log("F pressed.")
            document.getElementById("file").click()
        }
    })
            

function update_plot(data){

    var traces = [];

    var trace_f = {
        name:'data',
        type:'scatter',
		x:data.time,
		y:data.Y
	};

    let ymin = d3.min(data.Y),
        ymax = d3.max(data.Y),
        tmin = d3.min(data.time),
        tmax = d3.max(data.time);

    traces.push(trace_f);
    console.log(trace_f);


    var n = data.time.length;

    var ribbon_x = Array(2*n);
    var ribbon_y = Array(2*n);

    for(let i=0; i<n;i++){
        ribbon_x[i] = data.time[i];
        ribbon_x[n+i] = data.time[n-i-1]
        ribbon_y[i] = data.filter_mean[i]-data.filter_sd[i];
        ribbon_y[n+i] = data.filter_mean[n-i-1]+data.filter_sd[n-i-1];
    }

    traces.push({
        name:"denoised calcium+baseline",
        type:'scatter',
        mode:'lines',
        line:{color:"#39e75f"},
        x:data.time,
        y:data.filter_mean
    })

    var ribbon = {
        name:"stat",
        type:'scatter',
        fill:'tozerox',
        fillcolor: "rgba(0,100,80,0.2)",
        line: {color: "transparent"},
        x:ribbon_x,
        y:ribbon_y
    }

    traces.push(ribbon);

    if(data.spiketimes != undefined){
        if ( data.spiketimes instanceof Array){
            for(let i=0;i<data.spiketimes.length;i++){
                traces.push({
                    type:'scatter',
                    mode:'lines',
                    showlegend: false,
                    line: {color:"rgb(100,100,0)"},
                    x:[data.spiketimes[i],data.spiketimes[i]],
                    y:[ymin-(ymax-ymin)/10,ymax]
                })
            }
        } else {
            traces.push({
                type:'scatter',
                mode:'lines',
                showlegend: false,
                line: {color:"rgb(100,100,0)"},
                x:[data.spiketimes,data.spiketimes],
                y:[ymin,ymax]
            })
        }

    }

    traces.push({
        name:'baseline',
        type:'scatter',
        mode:'lines',
        line: {color:"orange", dash: 'dot'},
        x:data.time,
        y:data.B,
    })
    
    if(data.gtbase != undefined){
        traces.push({
            name:'gt baseline',
            type:'scatter',
            mode:'lines',
            line: {color:"orange"},
            x:data.time,
            y:data.gtbase
        })
    }


    // Now add second plot
    //
    //var intprob=data.prob.map((d,i,a)=>d3.sum(a.slice(Math.max(i-1,0),Math.min(i+2,n))));

    traces.push({
        name:'spike probability',
        type:'scatter',
        mode:'lines',
        line: {color:"rgb(255,0,0)"},
        x:data.time,
        //y:intprob,
        y:data.prob,
        xaxis:'x',
        yaxis:'y2'
    });


    if(data.stimtimes != undefined){
        for(let i=0;i<data.stimtimes.length;i++){
            traces.push({
                type:'scatter',
                mode:'lines',
                showlegend: false,
                line: {color:"rgb(100,100,0)"},
                x:[data.stimtimes[i],data.stimtimes[i]],
                y:[d3.min(data.prob),d3.max(data.prob)],
                xaxis:'x',
                yaxis:'y2'
            })
        }
    }

    // Now add third plot
    //

    if(data.burst != undefined){
        traces.push({
            name:'firing state',
            type:'scatter',
            mode:'lines',
            line: {color:"purple"},
            x:data.time,
            y:data.burst,
            xaxis:'x',
            yaxis:'y3'});
    }


    if(data.spiketimes_init != undefined){
        traces.push({
            name:'initialization',
            type:'scatter',
            mode:'lines',
            line: {color:"rgb(0,255,0)"},
            x:data.time,
            y:data.spiketimes_init,
            xaxis:'x',
            yaxis:'y3'
        });
    }

    var layout={
        height:800,
        xaxis:{
            title:{text:"time (s)"},
            range: [tmin,tmax]
        },

        yaxis:{
            range: [ymin-(ymax-ymin)/10,ymax]
        },

        grid: {
            rows: 3,
            columns: 1
        },

        margin: {
            l:60,
            r:5,
            b:20,
            t:60,
            pad:0
        }
    }

	Plotly.newPlot("plot",traces,layout,{'responsive':true});

}
