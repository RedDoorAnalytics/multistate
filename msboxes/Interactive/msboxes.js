// To do
/////////
// Arrow colours
// move arrows
// default colours for at options// 
// make percentages optional
// hover over transitions to see hazard functions


// Empty divs & clear timer
$("#GraphAreaBoxes").empty()
$("#TimeSliderBoxes").empty();
$("#TimeSliderButtonBoxes").empty();


for(i=1;i<=8;i++) {
	$("#HazardBoxes"+i).empty()	
}

if (typeof intervalListener != "undefined") clearInterval(intervalListener);

// Margins and heights of graphs 
margin = {top: 50, right: 20, bottom: 40, left: 50}
	, svgwidth = parseInt(d3.select('#GraphAreaBoxes').style('width'), 10)*0.9
	, svgheight = 0.7*svgwidth
	, graphwidth = svgwidth - margin.left - margin.right
	, graphheight = svgheight - margin.top - margin.bottom;

SliderMargin = {top: 10, right: 50, bottom: 100, left: 10}
SliderSvgWidth = parseInt(d3.select('#TimeSliderBoxes').style('width'), 10)*0.9
SliderHeight =10 
SliderWidth = SliderSvgWidth - SliderMargin.left - SliderMargin.right
SliderSvgHeight = SliderHeight + SliderMargin.top + SliderMargin.bottom;	

HazardCol = ["HazardLHS","HazardRHS"]

// Add Hazard Divs
for(var i=1;i<=msboxes.Ntransitions;i++) {
	ColSelect = 1 - i%2
	d3.select("#"+HazardCol[ColSelect])
		.append("div")
		.attr("id","HazardBoxes"+i)
}
// add Key div
d3.select("#HazardLHS")
	.append("div")
	.attr("id","BoxesKey")	

HazardMargin = {top: 5, right: 5, bottom: 5, left: 25}
HazardSVGHeight = svgheight/Math.ceil(msboxes.Ntransitions/2)
HazardSVGWidth = parseInt(d3.select('#HazardBoxes1').style('width'), 10)
HazardWidth = HazardSVGWidth - HazardMargin.left - HazardMargin.right
HazardHeight = HazardSVGHeight - HazardMargin.top - HazardMargin.bottom

// Initial values	
CurrentTime = 0	
CurrentIndex = 0
duration = 8000
NTimePoints = mspred.timevar.length
TimeMoving = 0
atcolours = d3.scaleOrdinal(d3.schemeCategory10);
legendRectSize = 30

// X and Y scales
xboxes = d3.scaleLinear()
		.domain([msboxes.xmin, msboxes.xmax])
		.range([ 0, graphwidth]);

yboxes = d3.scaleLinear()
		  .domain([msboxes.ymin, msboxes.ymax])
		  .range([graphheight, 0 ]);

// Set up SVG and plotarea		  
SVG = d3.select("#GraphAreaBoxes")
	.append("svg")
	.attr("width", svgwidth)
	.attr("height", svgheight)
	.attr('class', 'svg')

PlotArea = SVG.append('g')
	.attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')
	.attr('width', graphwidth)
	.attr('height', graphheight)
	.attr('class', 'plot')	

defs = SVG.append("defs")


defs.append("marker")
	.attr("id","arrow")
	.attr("viewBox","0 -5 10 10")
	.attr("refX",10)
	.attr("refY",0)
	.attr("markerWidth",10)
	.attr("markerHeight",10)
	.attr("orient","auto")
	.append("path")
		.attr("d", "M0,-5L10,0L0,5")
		.attr("class","arrowHead");

		
		
DrawGraphs()

/////////////////////////////////////////////////////////////////////////////////
function DrawGraphs() {
	DrawBoxes()
	DrawTransitions()
	DrawSlider()
	DrawPlayButton()
	ShadeArea()
	UpdateBoxPercentage()
	StateNames()
	HazardPlot()
	BoxesAddKey()
}
	 
// Add boxes
function DrawBoxes() {
	for (var state=1;state<=msboxes.Nstates;state++) {
		var sindex = state - 1
		PlotArea.append("rect")
			.attr("class","boxes")
			.attr("height", yboxes(msboxes.ymax - msboxes.boxheight))
			.attr("width", xboxes(msboxes.boxwidth))
			.attr("x", xboxes(msboxes.xvalues[sindex]))
			.attr("y", yboxes(msboxes.yvalues[sindex]))
			.attr("rx", 10)
			.attr("ry", 10)	
		}
}

// Draw Transitions
function DrawTransitions() {
	for(var trans=1;trans<=msboxes.Ntransitions;trans++) {
		var tindex = trans - 1
		PlotArea.append("line")	
			.attr("x1", xboxes(msboxes.arrows.x1[tindex]))
			.attr("y1", yboxes(msboxes.arrows.y1[tindex]))
			.attr("x2", xboxes(msboxes.arrows.x2[tindex]))
			.attr("y2", yboxes(msboxes.arrows.y2[tindex]))
			.attr("marker-end","url(#arrow)") 			
			.style(	"stroke", "black")	
			
		PlotArea.append("circle")
			.attr("class","TextCircle")
			.attr("r",xboxes(0.01))
			.attr("cx", xboxes((msboxes.arrowstext.x[tindex])))
			.attr("cy", yboxes(msboxes.arrowstext.y[tindex]+0.02))
			
		PlotArea.append("text")
			.attr("class","TextCircleText")
			.attr("x", xboxes((msboxes.arrowstext.x[tindex])))
			.attr("y", yboxes(msboxes.arrowstext.y[tindex]+0.015	))
			.attr("text-anchor","middle")
			.text(trans)
	}
}

// Initial Shaded areas
function ShadeArea() {
	for (var state=1;state<=msboxes.Nstates;state++) {
		var sindex = state - 1
		for (var at=1;at<=mspred.Nats;at++) {
			atindex = at - 1 
			probarea = mspred["P"+state][atindex][CurrentIndex]
			PlotArea.append("rect")
				.attr("id","fill"+state+at)
				.attr("class","fillboxes")
				.attr("fill",atcolours(atindex))
				.attr("height", yboxes(msboxes.ymax - msboxes.boxheight)*probarea)
				.attr("width", xboxes(msboxes.boxwidth)/mspred.Nats)
				.attr("y", yboxes(msboxes.yvalues[sindex]))
				.attr("x", xboxes(msboxes.xvalues[sindex] + (at - 1)*msboxes.boxwidth/mspred.Nats))
				.attr("rx", 10)
				.attr("ry", 10)	
		}
	}
}

// Update Shaded areas
function UpdateShadedAreas() {
	for (var state=1;state<=msboxes.Nstates;state++) {
		var sindex = state - 1
		for(var at=1;at<=mspred.Nats;at++) {
			atindex = at - 1 
			probarea = mspred["P"+state][atindex][CurrentIndex]
			d3.select("rect#fill"+state+at)
				.attr("y", yboxes(msboxes.yvalues[sindex] - msboxes.boxheight*(1-probarea)))
				.attr("height", yboxes(msboxes.ymax - msboxes.boxheight*probarea))
		}
	}
}

// State Names
function StateNames() {
	for (var state=1;state<=msboxes.Nstates;state++) {
		var sindex = state - 1
		PlotArea.append("text")
			.attr("class","Statenames")
			.text(msboxes.statenames[sindex])
			.attr("x", xboxes(msboxes.xvalues[sindex] + msboxes.boxwidth/2))
			.attr("y", yboxes(msboxes.yvalues[sindex] - msboxes.boxheight/2));
	}
}

// Percentage in Boxes	
function UpdateBoxPercentage() {
	d3.selectAll(".boxpercentages").remove()
	for (var state=1;state<=msboxes.Nstates;state++) {
		var sindex = state - 1
		for (var at=1;at<=mspred.Nats;at++) {
			atindex = at - 1 
			probarea = d3.format(".0%")(mspred["P"+state][atindex][CurrentIndex])		
			PlotArea.append("text")
				.attr("class","boxpercentages")
				.attr("x",xboxes(msboxes.xvalues[sindex] + (at - 1)*msboxes.boxwidth/mspred.Nats + msboxes.boxwidth/(2*mspred.Nats)))
				.attr("y",yboxes(msboxes.yvalues[sindex] - 0.2*msboxes.boxheight))
				.attr("id","boxpercentage"+state)
				.text(probarea)
		}
	}
}

// Draw Slider
function DrawSlider() {
	maxt = d3.max(mspred.timevar)
	SliderSVG = d3.select("#TimeSliderBoxes")
		.append("svg")
		.attr("width", SliderSvgWidth)
		.attr("height", SliderSvgHeight)

	SliderPlotArea = SliderSVG.append('g')
		.attr('transform', 'translate(' + SliderMargin.left + ',' + SliderMargin.top + ')')
		.attr('width', SliderWidth)
	.attr('height', SliderHeight)
	.attr('class', 'SliderBox')	

	xSlider = d3.scaleLinear().range([0, SliderWidth]).domain([0,maxt]);
	SliderBlobWidth = SliderHeight*1.1
	
// BoxArea
	SliderPlotArea.append("rect")
		.attr("class","SliderRect")
		.attr("width",SliderWidth)
		.attr("height", SliderHeight)
		.attr("x", 0)
 		.attr("y", 0)
		
// Axis
	SliderPlotArea.append("g")
		.attr("class", "SliderAxis")
		.attr("transform", "translate(0," + SliderHeight + ")")
		.call(d3.axisBottom(xSlider).ticks(10))
		
bisecttime = d3.bisector(function(d, x) { return d - x; }).left;
// Drag function
var SliderDrag = d3.drag()
		.on("start", function() {
			d3.event.sourceEvent.stopPropagation();
		})
		.on("drag", function(){
				var newx = d3.event.x;
				if(inrange(newx,xSlider(0),xSlider(maxt))) {
					d3.select(this)
						.attr("cx", newx)
						CurrentIndex = bisecttime(mspred.timevar,xSlider.invert(newx))
						CurrentTime = mspred.timevar[CurrentIndex]
						UpdateShadedAreas()
						UpdateBoxPercentage()
						MoveHazardCircles()

				}
		})
		.on("end", function(){
				d3.select(this)
					.attr("opacity",0.8);

		});		

// SliderBlob		
	SliderPlotArea.append("circle")
		.attr("class","SliderBlob")
		.attr("id","SliderBlob")
		.attr("cx",xSlider(CurrentTime))
		.attr("cy",SliderHeight/2)
		.attr("r",SliderBlobWidth)
		.call(SliderDrag)
	
	function inrange(x, min, max) {
		return x >= min && x <= max;
	}	
		
}


// PlayButton
function DrawPlayButton() {
	SliderButton = d3.select("#TimeSliderButtonBoxes")
		.append("button")
		.attr("id","play-button")
		.text("Play")

	SliderButton.on("click", function() {
			var button = d3.select(this);
				if (button.text() == "Pause") {
					clearInterval(intervalListener);
					TimeMoving = false;
					button.text("Play");
				} 
				else {
					if(CurrentIndex ==(NTimePoints-1)) CurrentIndex=0
					console.log(duration)
					TimeMoving = true;
					intervalListener = setInterval(TimeUpdate, duration/NTimePoints);
					button.text("Pause");
				}
		});
}

function DurationChange() {
	d3.selectAll(("input[name='speed']"))
		.on("change", function(){
			duration = this.value;
			if(TimeMoving) {
				clearInterval(intervalListener);
				intervalListener = setInterval(TimeUpdate, duration/NTimePoints);
			}
		});
}

function TimeUpdate() {
	CurrentIndex = CurrentIndex + 1

	if(CurrentIndex==(NTimePoints-1)) {
		clearInterval(intervalListener);
		SliderButton.text("Play")
		TimeMoving = false;
	}
	else {
		CurrentTime = mspred.timevar[CurrentIndex]
		newx = xSlider(CurrentTime)
		d3.select("circle#SliderBlob")
			.attr("cx", newx)
		UpdateShadedAreas()
		UpdateBoxPercentage()
		MoveHazardCircles()
	}
}

// Hazard Plots
function HazardPlot() {
	maxhazard = 0
	for(var i=1;i<=msboxes.Ntransitions;i++) {
		for(var j=0;j<mspred.Nats;j++) {
			maxhazard = Math.max(d3.max(mspred["h"+i][j]),maxhazard)
		}
	}

	xHazard = d3.scaleLinear().range([0, HazardWidth]).domain([0,maxt]);
	yHazard = d3.scaleLinear()
			.domain([0,maxhazard])
			.range([HazardHeight, 0 ]);
	
	HazardLine = d3.line()
		.defined(function(d) {return d; })
		.x(function(d,i) { return xHazard(mspred.timevar[i]);})
		.y(function(d) {return yHazard(d);});
		
	HazardSVG = []
	HazardPlotArea = []

	for(i=1;i<=msboxes.Ntransitions;i++) {
		trans_index = i - 1
		HazardSVG[trans_index] = d3.select("#HazardBoxes"+i)
			.append("svg")
			.attr("width", HazardSVGWidth)
			.attr("height", HazardSVGHeight)

		HazardPlotArea[trans_index] = HazardSVG[trans_index].append('g')
			.attr('transform', 'translate(' + HazardMargin.left + ',' + HazardMargin.top + ')')
			.attr('width', HazardWidth)
			.attr('height', HazardHeight)
			.attr('class', 'HazardPlot')
			
		HazardPlotArea[trans_index].append("g")
			.attr("class", "HazardAxis")
			.attr("transform", "translate(0," + HazardHeight + ")")
			.call(d3.axisBottom(xHazard).ticks([]))
			
		HazardPlotArea[trans_index].append("g")
			.attr("class", "HazardAxis")
			.call(d3.axisLeft(yHazard).ticks([]));	
			
		HazardPlotArea[trans_index].selectAll("#Hazardpath"+i)
			.attr("class","#HazardPath"+i)
			.attr("id","Hazard"+i)
			.data(mspred["h"+i])
			.enter()
			.append("path")
			.attr("d",HazardLine)
			.style("opacity",0.7)
			.style("stroke",function(d,i) {return(atcolours(i));})
			.style("fill","none")

		for(var j=0;j<mspred.Nats;j++) {
			HazardPlotArea[trans_index].append("circle")
				.attr("class","HazardCircles")
				.attr("id","Hazard"+i+j)
				.attr("cx",xHazard(CurrentTime))
				.attr("cy",yHazard(mspred["h"+i][j][CurrentIndex]))
				.attr("r",3)
				.attr("opacity",0.5)
				.attr("fill",atcolours(j))
		}
		
		HazardPlotArea[trans_index].append("circle")
			.attr("class","TextCircle")
			.attr("r", 10)
			.attr("cx", xHazard(-maxt/20))
			.attr("cy", yHazard(maxhazard/2))
			
		HazardPlotArea[trans_index].append("text")
			.attr("class","TextCircleText")
			.attr("x", xHazard(-maxt/20))
			.attr("y", yHazard(maxhazard/2))
			.attr("dominant-baseline", "central")			
			.text(i)		
		
	}
}

function MoveHazardCircles() {
	for(var i=1;i<=msboxes.Ntransitions;i++) {
		var trans_index = i - 1
		for(var j=0;j<mspred.Nats;j++) {
			HazardPlotArea[trans_index].select("#Hazard"+i+j)
				.attr("cx",xHazard(CurrentTime))
				.attr("cy",yHazard(mspred["h"+i][j][CurrentIndex]))
		}	
	}		
}

// Add Legend
function BoxesAddKey() {
	var legendColumns = Math.ceil(mspred.Nats/3)
	var legend = d3.select("#BoxesKey")
			.append("svg")
			.attr("width", HazardSVGWidth)
			.attr("height", HazardSVGHeight)

  tooltip = legend.append("g")
    .attr("class", "tooltip")
    .style("display", "none");
      
  tooltip.append("rect")
    .attr("width", 60)
    .attr("height", 60)
    .attr("fill", "blue")
    .style("opacity", 0.5);

  tooltip.append("text")
    .style("text-anchor", "middle")
    .attr("font-size", "12px")
    .attr("font-weight", "bold")
	.attr("id","tooltiptext")
	.text("HELLO");

		
	legend.selectAll("rect")
		.data(mspred.atlist)
		.enter()
		.append("rect")
		.attr("id",function(d,i) {return(i)})
		.attr("width", legendRectSize)
		.attr("height", legendRectSize)
		.attr("y",function(d,i) {return((i%legendColumns)*(legendRectSize*1.1));})
		.attr("x",function(d,i) {return(Math.floor(i/legendColumns))*100;})
		.style("fill", function(d,i) {return(atcolours(i));})
		.style("stroke", function(d,i) {return(atcolours(i))})
		.style("opacity",0.5)
		.on("mouseover",BoxesLegendMouseOver)

	legend.selectAll("text")
		.data(mspred.atlist)
		.enter()
		.append("text")
		.text(function(d,i) {return("at"+(i+1));})
		.attr("y",function(d,i) {return((i%legendColumns)*(legendRectSize*1.1) + legendRectSize/2);})
		.attr("x",legendRectSize*1.2)
		.attr("x",function(d,i) {return(Math.floor(i/legendColumns)*100 + legendRectSize*1.2) ;})
}

function BoxesLegendMouseOver() {
	console.log(-99)
	var xPosition = parseFloat(d3.select(this).attr("x")) 
	var yPosition = parseFloat(d3.select(this).attr("y")) 
	console.log(xPosition,yPosition)
	var atvalue = this.id
	console.log(mspred.atlist[atvalue])
	d3.select("#tooltiptext")
		.attr("x", xPosition)
		.attr("y", yPosition)
		.text(mspred.atlist[atvalue])	
		
	tooltip.style("display", null);		
}
