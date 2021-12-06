// xtitle
// ytitle
// read at from predictms
// fix Edge Issue

// Empty divs & clear timer

$("#GraphAreaStacked").empty()
$("#RadioButtons").empty()
if (typeof intervalListener != "undefined") clearInterval(intervalListener);

// Margins and heights of graphs 
margin = {top: 10, right: 20, bottom: 50, left: 50}
	, svgwidth = parseInt(d3.select('#GraphAreaStacked').style('width'), 10)*0.9
	, svgheight = 0.8*svgwidth
	, graphwidth = svgwidth - margin.left - margin.right
	, graphheight = svgheight - margin.top - margin.bottom;

// Initial values
CurrentAt = 0
CurrentTime = maxt
CurrentIndex = mspred.timevar.length	
StateColour = d3.scaleOrdinal(d3.schemeCategory10);

// X and Y scales
xtime = d3.scaleLinear()
			.domain([0,maxt])
			.range([ 0, graphwidth]);

yprob = d3.scaleLinear()
			.domain([0,1])
			.range([graphheight, 0 ]);


PlotArea = d3.area()
	.x(function(d,i) {return xtime(mspred.timevar[i]);})
	.y0(function(d) {return yprob(d[0]);})
	.y1(function(d) {return yprob(d[1]);})
/*	
	.interpolate("linear")
	.defined(function(d,i) { if(isFinite(d[0])) {
								return 1}
							 else {
								return 0}})	
*/			
// Set up SVG and plotarea		  
SVGStacked = d3.select("#GraphAreaStacked")
	.append("svg")
	.attr("width", svgwidth)
	.attr("height", svgheight)
	.attr('class', 'svg')

PlotAreaStacked = SVGStacked.append('g')
	.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')
	.attr('width', graphwidth)
	.attr('height', graphheight)
	.attr('class', 'plot')	

DrawGraphs()

/////////////////////////////////////////////////////////////////////////////////
function DrawGraphs() {
	DrawAxes()
	DrawSliderStacked()
	InitialShadedAreasStacked()
	AddRadioButtons()
	//	AutoRunTime()
}
	
function DrawAxes() {
	// Add and draw the X & Y Axes
	PlotAreaStacked.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + graphheight + ")")
		.call(d3.axisBottom(xtime));
	PlotAreaStacked.append("g")
		.attr("class", "y axis")
		.call(d3.axisLeft(yprob));
	PlotAreaStacked.append("text")
		.attr("x",xtime(maxt/2))
		.attr("y",yprob(-0.08))
		.attr("class", "x axis")
		.text("Time")
	PlotAreaStacked.append("text")
		.attr("x",xtime(1))
		.attr("y",yprob(0.5))
		.text("Probability")
		
}

function InitialShadedAreasStacked() {
	var areaupper = Array.apply(0, Array(CurrentIndex)).map(function (x, y) { return 0; })
	var ones = Array.apply(0, Array(CurrentIndex)).map(function (x, y) { return 1; })	
	for(var i=1;i<=msboxes.Nstates;i++) {
		var StateIndex = i - 1
		if(i<msboxes.Nstates) {
			var arealower = []
			for(j=0;j<areaupper.length;j++) {
				arealower[j] = areaupper[j] + mspred["P"+i][CurrentAt][j]
			}
		}
		else var arealower = ones	
		var AreaData = d3.zip(areaupper,arealower)
		
		PlotAreaStacked.append("path")
			.attr("class","area")
			.attr("id","AreaState"+i)
			.datum(AreaData)
			.attr("d",PlotArea)
			.style("opacity",0.7)
			.style("stroke",StateColour(StateIndex))
			.style("fill",StateColour(StateIndex))

		var areaupper = arealower
	}
}

function RemoveShadedAreasStacked() {
	d3.selectAll(".area").remove()
}

function AddRadioButtons() {
	var buttons = d3.select("#RadioButtons").append("form")
	buttons.selectAll("input")
		.data(mspred.atlist)
		.enter()	
		.append("div")
		.attr("class","radio")
		.text(function(d) {return d+" ";})
		.insert("input")
		.attr("type","radio")
		.attr("name","atradiobuttons")
		.attr("class", "atradiobuttons")
		.attr("value", function(d, i) {return i;})
		.property("checked", function(d, i) {return i===CurrentAt;})
		.on("change", ChangeAt);
}

function ChangeAt() {
	CurrentAt = this.value
	TransitionShadedArea()
}

function TransitionShadedArea() {
	var areaupper = Array.apply(0, Array(CurrentIndex)).map(function (x, y) { return 0; })
	var ones = Array.apply(0, Array(CurrentIndex)).map(function (x, y) { return 1; })	
	for(var i=1;i<=msboxes.Nstates;i++) {
		var StateIndex = i - 1
		if(i<msboxes.Nstates) {
			var arealower = []
			for(j=0;j<areaupper.length;j++) {
				arealower[j] = areaupper[j] + mspred["P"+i][CurrentAt][j]
			}
		}
		else var arealower = ones	
		var AreaData = d3.zip(areaupper,arealower)
		
		d3.select("path#AreaState"+i)
			.transition().duration(1500).ease(d3.easeCubic)
			.attr("d",PlotArea(AreaData))
		
		var areaupper = arealower
	}	
}
	




		

// Draw Slider
function DrawSliderStacked() {
	bisecttime = d3.bisector(function(d, x) { return d - x; }).left;
// Drag function
	var SliderDrag = d3.drag()
		.on("start", function() {
			d3.event.sourceEvent.stopPropagation();
		})
		.on("drag", function(){
				var newx = d3.event.x;
				if(inrange(newx,xtime(0),xtime(maxt))) {
					d3.select(this)
						.attr("cx", newx)
						CurrentIndex = bisecttime(mspred.timevar,xtime.invert(newx))
						CurrentTime = mspred.timevar[CurrentIndex]
						RemoveShadedAreasStacked()
						InitialShadedAreasStacked()

				}
		})
		.on("end", function(){
				d3.select(this)
					.attr("opacity",0.8);

		});		

// SliderBlob		
	PlotAreaStacked.append("circle")
		.attr("class","SliderBlob")
		.attr("id","SliderBlob")
		.attr("cx",xtime(CurrentTime))
		.attr("cy",yprob(0))
		.attr("r",SliderBlobWidth)
		.call(SliderDrag)

	function inrange(x, min, max) {
		return x >= min && x <= max;
	}	
		
}
	