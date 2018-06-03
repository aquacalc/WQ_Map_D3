// initially, code to re-create pH isopleths on WQ map

// # nb: pH(NBS) = pH(Free) -
// #               log10(ahFreeToSwsFactor(S, T(Kelvin), 0)) -
// #               log10(ahSwsToNbsFactor())
let phNbsToPhFree = (ph, sal, temp, p) => {
  return(ph +
           Math.log10(ahSwsToNbsFactor(sal, temp, p)) +
           Math.log10(ahFreeToSwsFactor(sal, temp, p))
         );
}

let ahSwsToNbsFactor = (sal, temp, p) => {
  return(calcProtonActivityCoeffZg(temp, sal, 0.0) / ahMolalToMolinforSalinity(sal));
}

let calcProtonActivityCoeffZg = (temp, sal, p) => {

  let rootGamma = Math.sqrt(calcIonicStrength(sal));

  let myHActivityCoeff = 1820000.0 * Math.pow((79 * temp), -1.5);

  myHActivityCoeff = myHActivityCoeff * ((rootGamma / (1 + rootGamma)) - 0.2 * calcIonicStrength(sal));

  myHActivityCoeff = Math.pow(10, -myHActivityCoeff);

  return(myHActivityCoeff);
}

let calcIonicStrength = (sal) => {
  let myIS = 19.924 * sal / (1000.0 - 1.005 * sal); // mole/kg-H2O (molal)
  return(myIS);
}

let ahMolalToMolinforSalinity = (sal) => {
  return(1.0 - 0.001005 * sal);
}


let ahFreeToSwsFactor = (sal, temp, p) => {
  return(1 + (calcTS(sal) / calcKsDickson(temp, sal, 0.0)) +
           (calcTF(sal) / calcKfDickson(temp, sal, 0.0)));
}

let calcTF = (sal) => {
  return(0.0000019522 * sal);
}

let calcTS = (sal) => {
  return(0.0008067267 * sal);
}

// # pH scale: FREE
// # concentration scale: mole/kg-H2O ?? -> molin...?
let calcKfDickson = (temp, sal, p) => {

  let fluorFactor1 = 1590.2 / temp;
  let fluorFactor2 = -12.641;
  let fluorFactor3 = 1.525 * Math.sqrt(calcIonicStrength(sal));

  // # ** nb: molalToMolin factor in sqrt()...
  // #    double fluorFactor3 = 1.525 * sqrt([self calcIonicStrength:s] * (1.0 - 0.001005 * s));

  let molal2molin = Math.log(1.0 - 0.001005 * sal);

  let KF = fluorFactor1 + fluorFactor2 + fluorFactor3 + molal2molin;

  return(Math.exp(KF))
}

// # bisulfate dissociation
// # Dickson (1990); DOE (1994), ch. 5 p. 13; Z & W-G (2001) p. 260
// # pH scale: Free
// # concentration scale: mol/kg-H2O, CONVERTED TO AND RETURNED AS mol/kg-soln (molin)
// # **** nb **** if called from within pH conversion, T already C -> K *******
// # **** nb **** else if called from getKS(), must add 273.15 in call ****************
let calcKsDickson = (temp,sal,p) => {

  // # ** NB: ionic strength calc now returns molal, so
  // # **     change to molin (mol/kg-soln) here for this calc
  // #    double myIS = [self calcIonicStrength:s] * (1.0 - 0.001005 * s);
  let myIS = calcIonicStrength(sal);

  let sulfFactor1 = 141.328 - 4276.1 / temp;
  let sulfFactor2 = -23.093 * Math.log(temp);
  let sulfFactor3 = ( 324.57 - 47.986 * Math.log(temp) - 13856 / temp) * Math.sqrt(myIS)
  let sulfFactor4 = (-771.54 + 114.723 * Math.log(temp) + 35474 / temp) * myIS
  let sulfFactor5 = -2698.0 * Math.pow(myIS, 1.5) / temp
  let sulfFactor6 =  1776.0 * myIS * myIS / temp

  molal2molin = Math.log(1.0 - 0.001005 * sal)

 // # ** in mol/kg-soln (MOLIN)
  let KS = sulfFactor1 + sulfFactor2 + sulfFactor3 + sulfFactor4 + sulfFactor5 + sulfFactor6 + molal2molin;

 // # ** in mol/kg-H2O (MOLAL)
 // #	double KS = sulfFactor1 + sulfFactor2 + sulfFactor3 + sulfFactor4 + sulfFactor5 + sulfFactor6;

  return(Math.exp(KS))
}


let phLineIntercept = (temp, sal, ph) => {
  // console.log(`calcHydroxide = ${calcHydroxide(ph, temp, sal)}`);
  // console.log(`calcHydronium = ${calcHydronium(ph)}`);
  // console.log(`   calcBorate = ${calcBorate(ph, temp, sal)}`);
  return(calcHydroxide(ph, temp, sal) - calcHydronium(ph) + calcBorate(ph, temp, sal));
}

let calcHydroxide = (ph, t, sal) => {
  // if(sal === null)
  //   return();
  let kWToTheTen = Math.log10(calcKWMehrbach(t, sal));

  return(Math.pow(10, (kWToTheTen + ph)));
}

let calcKWMehrbach = (temp, sal) => {
  // if(sal === null)
  //   return()

  let expSum = 148.9652;     // # nb: "148.965 02" in Zeebe & Wolf-Gladrow code
  expSum = expSum - 13847.26 / temp;
  expSum = expSum - 23.6521 * Math.log(temp);
  expSum = expSum + (-5.977 + (118.67 / temp) + 1.0495 * Math.log(temp)) * Math.sqrt(sal);
  expSum = expSum - 0.01615 * sal;

  let KW = Math.exp(expSum);   // # still on TOTAL scale
  let pKW = -Math.log10(KW);   // # still on TOTAL scale

  // # ** nb: convert to FREE pH scale, as per AquaEnv
  pKW = pKW + Math.log10(ahFreeToTotFactor(sal, temp, 0.0));

  return(Math.pow(10, -pKW));
}

let ahFreeToTotFactor = (sal, temp, p) => {
  return(1 + (calcTS(sal) / calcKsDickson(temp, sal, 0.0)));
}

let calcHydronium = (ph) => {
  return(Math.pow(10, -ph));
}

let phLineSlope = (temp, sal, ph) => {
  return(alphaOne(temp, sal, ph) + 2 * alphaTwo(temp, sal, ph));
}

// # Borate ----
// # ** FREE scale?
let calcBorate = (ph, t, sal) => {

  let concB = calcBorateConcOfSalinity(sal);

  let myKB = calcBorateFactor(t, sal);

  let borate = myKB * concB / (myKB + Math.pow(10, -ph));

  return(borate);
}


let calcBorateConcOfSalinity = (sal) => {

  let concB = 0.000232 * sal / (10.811 * 1.80655);

  return(concB);
}

let calcBorateFactor = (t, sal) => {
  // if(sal === null)
  //   return();

  let A = 148.0248 + 137.1942 * Math.sqrt(sal) + 1.62142 * sal;
  let B = -8966.90 - 2890.53 * Math.sqrt(sal) - 77.942 * sal + 1.728 * Math.pow(sal, 1.5) - 0.0996 * sal * sal;
  let C = -24.4344 - 25.085 * Math.sqrt(sal) - 0.2474 * sal;
  let D = 0.053105 * Math.sqrt(sal);

  let K_BOH3 = Math.exp(A + B/t + C * Math.log(t) + D * t);

  let ans = Math.pow(10, -(-Math.log10(K_BOH3) + Math.log10(ahFreeToTotFactor(sal, t, 0))));

  return(ans);
}


let alphaOne = (temp, sal, ph) => {

  // # nb: define p LOCALLY until incorporate in calcs
  let p = 0.0;

  let h = calcHydronium(ph);

  let k1 = getK1(temp, sal, 0);

  let k2 = getK2(temp, sal, 0);

  let numerator = h * k1;

  return (numerator / calcAlphaDenom(h, k1, k2));
}

alphaTwo = function(temp, sal, ph) {

  // # nb: define p LOCALLY until incorporate in calcs
  let p = 0.0;

  let h = calcHydronium(ph);

  let k1 = getK1(temp, sal, 0);

  let k2 = getK2(temp, sal, 0);

  let numerator = k1 * k2;

  return (numerator / calcAlphaDenom(h, k1, k2));
}

let calcAlphaDenom = (h, k1, k2) => {
  return(h * h + k1 * h + k1 * k2);
}

// # T- and S-dependent K1 from Millero et al. (2006)
// # ** pH scale: SWS for calculation, FREE returned
let getK1 = (temp, sal, p) => {
  // if(is.null(sal))
  //   return()

  sqrtS = Math.sqrt(sal);
  lnT   = Math.log(temp);

  let pK1z = -126.34048 + (6320.813 / temp) + 19.568224 * lnT;  // # was "19.56822.."
  let A =   13.4191 * sqrtS + 0.0331 * sal - 0.0000533 * sal * sal;
  let B = -530.123 * sqrtS - 6.103 * sal;
  let C =   -2.06950 * sqrtS;
  let pK1 = pK1z + A + (B / temp) + C * lnT;

  // # ** nb: THIS conversion puts it on the FREE pH scale *from* the SWS scale, as per AquaEnv
  // # **     factor for sws2free = 1.0 / ahFreeToSwsFactorForSalinity:temp:pressure:
  // # **     when dealing with -log10, -log10(1 / ah) = -(-log10(ah)) = +log10(ahFreeToSWSFactor...)
  pK1 = pK1 + Math.log10(ahFreeToSwsFactor(sal, temp, 0.0));

  // # ?? nb: no concentration scale conversion needed, as both on molinity [sic]

  return(Math.pow(10, -pK1));  // # K1
}


// # T- and S-dependent K2 from Millero et al. (2006)
// # ** pH scale: SWS for calculation
// # ** Return: K2 -- not pK2 -- on FREE pH scale

// # ----> nb: Millero (2010) slightly changes some coefficients <----

let getK2 = (temp, sal, p) => {

  // if(is.null(sal))
  //   return()

  let sqrtS = Math.sqrt(sal);
  let lnT   = Math.log(temp);

  let pK2z = -90.18333 + (5143.692 / temp) + 14.613358 * lnT;

  let A = 21.0894 * sqrtS + 0.1248 * sal - 0.0003687 * sal * sal; // # was "21.08945"
  let B = -772.483 * sqrtS - 20.051 * sal;
  let C = -3.3336 * sqrtS; //  # was "-3.3336" or "3.32254"

  // # ** nb: THIS pK2 is on the SWS pH scale
  let pK2 = pK2z + A + (B / temp) + C * lnT;

  pK2 = pK2 + Math.log10(ahFreeToSwsFactor(sal, temp, 0.0));

  return(Math.pow(10, -pK2));
}


let calcDicOfAlk = (alk, ph, temp, sal) => {

  let m = phLineSlope(temp, sal, ph);

  let dic = (alk - calcHydroxide(ph, temp, sal) -
                  calcBorate(ph, temp, sal) +
                  calcHydronium(ph)) / m;
  return(dic);
}



// ************************************************ //
// ************************************************ //
// ************************************************ //
// ************************************************ //

// import CarbCalc from 'CarbCalc';

// $(document).ready(() => {

// d3.json('tz.2012.1.json').then(data => {
//   console.log('RW-I');
//   const cleanedData = data.filter(d => d.temp != null);
//   console.log(cleanedData);
// }).catch(err => console.log(`Blimey! ${err}`));


// programmatically click 'calc' button on 'ready'
  setTimeout(() => {
    $('#btn-calc').click();
  },1);


// ------- COLLECT POINTS FOR pH ISOPLETHS ---------

// temp & sal data
  let temp    = +$('#temperature').val();
  temp += 273.15;
  let sal     = +$('#salinity').val();

  // nb: dicMax to 'help' not drawing beyond chart boundaries
  let dicMax = 5;
  let alkMax = 5;

// see: https://stackoverflow.com/questions/3746725/create-a-javascript-array-containing-1-n
  const phByQuarters = Array.from({length: 30}, (v, k) => 3.75 + (k + 5) / 4);

  let calcPhIsopleths = (temp, sal) => {

    return phByQuarters.map((ph, idx) => {
    // let pHIsopleths = phByQuarters.map((ph, idx) => {

      let ph_free       = phNbsToPhFree(ph, sal, temp, 0);
      let ph_intercept  = 1000 * phLineIntercept(temp, sal, ph_free);
      let ph_slope      = phLineSlope(temp, sal, ph_free);
      let alk_at_dicmax = ph_intercept + ph_slope * dicMax;

      // if(idx <= 3) {
      //   // console.log(`${idx} (${ph}): slope = ${ph_slope}, int = ${ph_intercept}`);
      //   console.log(`${idx} (${ph}): alk_at_dicmax = ${alk_at_dicmax}`);
      // }

      // max DIC point to plot ALK UNLESS resulting ALK > alkMax
      let dicMaxPlot = dicMax;

      // when alk_at_max > alkMax, d3 draws beyond chart group, into margin
      // nb: ...AND ADJUST VALUE OF dicMax KEY...
      if(alk_at_dicmax > alkMax) {
        // find dic at alkMax
        dicMaxPlot = 1000 * calcDicOfAlk(alkMax / 1000, ph_free, temp, sal);
        // re-calc alk_at_max
        alk_at_dicmax = ph_intercept + ph_slope * dicMaxPlot;
      }

      // console.log(`${idx} (${ph}): alk_at_dicmax = ${alk_at_dicmax}`);

      return [ { dic: 0.01, alk: ph_intercept }, { dic: dicMaxPlot, alk: alk_at_dicmax } ];

    });
  }

  let pHIsopleths = calcPhIsopleths(temp, sal);


  // ------------------------------------------------

  // JQuery button event handler
  $('#btn-calc').on('click', () => {
    let temp    = +$('#temperature').val();
    temp += 273.15;
    let sal     = +$('#salinity').val();
    // let ph      = +$('#ph').val(); // NBS
    // let ph_free = phNbsToPhFree(ph, sal, temp, 0);
    // let dic     = +$('#alk_at_dic').val();

    // console.log(`Temp    = ${temp} K`);
    // console.log(`Sal     = ${sal} ppt`);
    // console.log(`pH      = ${ph} NBS`);
    // console.log(`pH_free = ${ph_free} FREE`);
    // console.log(`dic     = ${dic} mg/kg`);

    // ph_intercept = 1000 * phLineIntercept(temp, sal, ph_free);
    // ph_slope     = phLineSlope(temp, sal, ph_free);
    // alk_at_dic   = ph_intercept + ph_slope * dic;
    //
    // $('#answer1').html(ph_intercept);
    // $('#answer2').html(ph_slope);
    // $('#answer3').html(alk_at_dic);

    pHIsopleths = calcPhIsopleths(temp, sal);
    update(pHIsopleths);

  });

// ** EVENT LISTENERS **

// Event listeners
// $("#temperature").on("change", update(pHIsopleths));

// $("#temperature").on("change", update(pHIsopleths));
// $("#salinity").on("change", update)

// ** TEMPERATURE INPUT
  // $('#temperature').on('change', () => {
  //   let temp    = +$('#temperature').val();
  //   temp += 273.15;
  //   let sal     = +$('#salinity').val();
  //
  //   pHIsopleths = calcPhIsopleths(temp, sal);
  //   update(pHIsopleths);
  // });

  $('#temperature').on('keyup', () => {
    let temp    = +$('#temperature').val();
    temp += 273.15;
    let sal     = +$('#salinity').val();

    pHIsopleths = calcPhIsopleths(temp, sal);
    update(pHIsopleths);
  });

// ** SALINITY INPUT
  $('#salinity').on('change', () => {
    let temp    = +$('#temperature').val();
    temp += 273.15;
    let sal     = +$('#salinity').val();

    pHIsopleths = calcPhIsopleths(temp, sal);
    update(pHIsopleths);
  });

  $('#salinity').on('keyup', () => {
    let temp    = +$('#temperature').val();
    temp += 273.15;
    let sal     = +$('#salinity').val();

    pHIsopleths = calcPhIsopleths(temp, sal);
    update(pHIsopleths);
  });




// D3 CODE ******************************************************

// ------------ DIMENSIONS --------------------- //
	const margins = { top: 5, right: 75, bottom: 60, left: 75 };

	let widthTotal  = 700;
	let heightTotal = 550;

	let widthPlot   = widthTotal  - (margins.left + margins.right);
	let heightPlot  = heightTotal - (margins.top  + margins.bottom);

// -------------------------- SELECTION  ------------------------------- //
	let g = d3.select('#chart-area')
							.append('svg')
								.attr('width', widthTotal)
								.attr('height', heightTotal)
							.append('g')
								.attr('transform', 'translate(' + margins.left + ', ' + margins.top + ')');


// // see: https://bl.ocks.org/mbostock/1087001
//       let div = d3.select("#chart-area").append("div")
//         // .attr("class", "tooltip")
//         .attr('position', 'absolute')
//         .attr('text-align', 'center')
//         .attr('width', 60)
//         .attr('height', 12)
//         .attr('padding', 8)
//         .attr('margin-top', -20)
//         .attr('font', '10px sans-serif')
//         .attr('background', '#ddd')
//         .attr('pointer-events', 'none')
//         .style("display", "none");
//
//       // see: https://bl.ocks.org/mbostock/1087001
//       let mouseover = () => {
//         // div.style("display", null);
//         div.style("display", "inline");
//       };
//
//       let mousemove = () => {
//         // console.log(x.invert(d3.event.pageX) + ", " + y.invert(d3.event.pageY));
//         div.text(d3.event.pageX + ", " + d3.event.pageY)
//         .style("left", (d3.event.pageX - 34) + "px")
//         .style("top", (d3.event.pageY - 12) + "px");
//       };
//
//       let mouseout = () => {
//         div.style("display", "none");
//       };
      //
      // g.append("rect")
      //   .attr("width", widthPlot)
      //   .attr("height", heightPlot)
      //   .attr('opacity', 0)
      //   .on("mouseover", mouseover)
      //   .on("mousemove", mousemove)
      //   .on("mouseout", mouseout);


// ------------ INITIALIZE TOOL-TIP --------------------- //
// let tip = d3.tip().attr('class', 'd3-tip')
//     .html(d => {
//         var text = "<strong>DIC:</strong> <span style='color:yellow'>" + d.dic + "</span><br>";
//
//         text += "<strong>[Alk]:</strong> <span style='color:yellow;text-transform:capitalize'>" + d.alk + "</span><br>";
//         return text;
//     });
//
// 	g.call(tip);


// --------------------------------- SCALES ---------------------------- //

// ** X-AXIS SCALE **
	let x = d3.scaleLinear()
		// .base(2)
		.domain([ 0, 5 ])
		.range([ 0, widthPlot ]);

	// ** Y-AXIS-LEFT SCALE **
	let y = d3.scaleLinear()
		.domain([ 0, 5 ])
		.range([ heightPlot, 0 ]);

	// ** Y-AXIS-RIGHTSCALE **
  // 100.086 mg/mmol CaCO3, 50.043 mg/meq
	let y_right = d3.scaleLinear()
		.domain([ 0, 250.215 ])
		// .domain([ 0, 5 ])
		.range([ heightPlot, 0 ]);


// -------------------------------- AXES ------------------------------- //
	// x-AXIS GENERATOR
	let xAxisGroup = d3.axisBottom(x);
    // .tickValues([400, 4000, 40000]);
    // .tickFormat(d3.format("$"));
	// append group and call(x-AXIX GEMERATOR)
	g.append('g')
		.attr('class', 'x axis')
		.attr('transform', 'translate(0, ' + heightPlot + ')')
		.call(xAxisGroup);

	let yAxisGroup = d3.axisLeft(y);
				// .ticks(8);
	g.append('g')
		.attr('class', 'y axis')
		.call(yAxisGroup);

	let yAxisGroupR = d3.axisRight(y_right)
    .ticks(6);
	g.append('g')
		.attr('class', 'y axis')
		.attr('transform', 'translate(' + widthPlot + ', 0)')
		.call(yAxisGroupR);


// --------------------------------- LABELS ---------------------------- //
	g.append('text')
			.attr('class', 'x class-label')
			.attr('x', widthPlot / 2)
			.attr('y', heightPlot + 50)
			.attr('font-size', '16')
			.attr('font-weight', 'bold')
			.attr('text-anchor', 'middle')
			.text('DIC (mmol/kg)');

	let yAxisLabel = g.append('text')
			.attr('class', 'y class-label')
			.attr('x', -heightPlot / 2)
			.attr('y', -53)
			.attr('font-size', '16')
			.attr('font-weight', 'bold')
			.attr('text-anchor', 'middle')
			.attr('transform', 'rotate(-90)')
			.text('Total Alkalinity (meq/kg)');


// ------------------------- pH ISOPLETHS ------------------------------- //
// see: https://stackoverflow.com/questions/8689498/drawing-multiple-lines-in-d3-js

// ----------------- based on BOSTOCK EXAMPLE ------------------------------//

  let update = data => {

    let t = () => d3.transition().duration(1000);

    let line = d3.line()  // lineGenerator
        .x(d => x(d.dic))
        .y(d => y(d.alk));
        // .y(function(d, i) {
        //   console.log(`${i}: ${d.alk}`);
        //   return y(d.alk);
        // }); // set the y values for the line generator


    let lines = g.selectAll(".line")
        .data(data, d => d.alk);

        // console.log('BIND...');
        // console.log(lines);
        // console.log('-----');

        lines.exit()
          .attr('class', 'exit')
          // .transition(t())
          .transition().duration(500)
          .attr("stroke", 'white')
          .style('opacity', 0)
          .remove();

          // console.log('EXIT...');
          // console.log(lines);
          // console.log('-----');

        // data.map(d => console.log(d));

        lines.enter().append("path")
    			// .on('mouseover', tip.show)
    			// .on('mouseout', tip.hide)
// nb-nb-nb: MUST MERGE, or FIRST DATUM LOST ON RE-RENDER
          .merge(lines)
            .attr("class", "line")
            .attr('fill', 'none')
            .transition(t())
            // .transition().duration(5000)
              .attr("stroke", (d, i) => (i % 4 === 0) ? 'black' : 'lightgrey')
              .attr("opacity", (d, i) => (i % 4 === 0) ? '1' : '0.9')
              .attr("d", line);
            // .attr("d", line(data));

          // console.log('ENTER...');
          // console.log(lines);
          // console.log('-----');
  }
// }); // END document ready...

$("#temperature").on("change", update(pHIsopleths));


// see: https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API/Using_Fetch
  // let myHeaders = new Headers();
  // let myInit = { method: 'GET',
  //                headers: myHeaders,
  //                mode: 'cors',
  //                cache: 'default' };
  // let myRequest = new Request('pHIsopleths', myInit);

// see: https://stackoverflow.com/questions/17214293/importing-local-json-file-using-d3-json-does-not-work
  // d3.json(myRequest).then(data => {
  //   data.map(d => console.log(d));
  // }).catch(err => console.log(`Blimey! ${err}`));

// ------------------------------------------------
