var counter = 0;
var limit = 100;
var epsilon = [];
var mu = [];
var length = [];
var crystal;
var oldOverlayLength;

function quickDraw() {
    addNumLayers('dynamicInput', 3);
    fillSampleValues();
    buildCrystal();
    printFieldsChart('linechart2_material', 900, 500);
    getIncomingMode();
    printDispersionChart('dispersionChart', 900, 500);

}



function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}

function addStruct(htmlClass, x, color, width, height) {
    //change addhere to id, customize offset based on which div is id'd

    if(htmlClass=="addHere")
        var overlayDiv = "<div class='overlay'><svg width='" + (width+1) + "' height='" + height + "'><defs><pattern id='pattern-stripe' width='4' height='4' patternUnits='userSpaceOnUse' patternTransform='rotate(45)'><rect width='2' height='4' transform='translate(0,0)' fill='white'></rect></pattern><mask id='mask-stripe'><rect x='0' y='0' width='100%' height='100%' fill='url(#pattern-stripe)' /></mask></defs><rect class='struct rect" + x + "' width='" + (width + 1) + "' height='" + height + "' fill='" + color + "'></svg></div>"
    if(htmlClass=="addHere2")
        var overlayDiv = "<div class='overlay2'><svg width='" + (width+1) + "' height='" + height + "'><defs><pattern id='pattern-stripe' width='4' height='4' patternUnits='userSpaceOnUse' patternTransform='rotate(45)'><rect width='2' height='4' transform='translate(0,0)' fill='white'></rect></pattern><mask id='mask-stripe'><rect x='0' y='0' width='100%' height='100%' fill='url(#pattern-stripe)' /></mask></defs><rect class='struct rect" + x + "' width='" + (width + 1) + "' height='" + height + "' fill='" + color + "'></svg></div>"
    if(htmlClass=="addHere3")
        var overlayDiv = "<div class='overlay3'><svg width='" + (width+1) + "' height='" + height + "'><defs><pattern id='pattern-stripe' width='4' height='4' patternUnits='userSpaceOnUse' patternTransform='rotate(45)'><rect width='2' height='4' transform='translate(0,0)' fill='white'></rect></pattern><mask id='mask-stripe'><rect x='0' y='0' width='100%' height='100%' fill='url(#pattern-stripe)' /></mask></defs><rect class='struct rect" + x + "' width='" + (width + 1) + "' height='" + height + "' fill='" + color + "'></svg></div>"

    $("."+htmlClass).append(overlayDiv);
}

function addStruct1() {
    addStruct(3, 'red', 100, 300);

}


function addNumLayers2(divName, numInputs) {
    numberOfLayers = parseInt(numInputs) + 2;
    var layersToAdd = numberOfLayers;
    if (counter < layersToAdd) {
        if (counter == 0) {

            // document.getElementsByClassName("h"+counter).remove();
            $("#h" + counter).remove();
            var newdiv = document.createElement('div');
            newdiv.innerHTML = "<div id='h" + counter + "'>Ambient Left " + " <br><b>Epsilon:</b> <input type='text' size='6' id='e" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Mu:</b> <input type= 'text' size='6'  id = 'm" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Length:</b> <input type = 'text' size='6' id = 'l" + counter + "'></div>";
            document.getElementById(divName).appendChild(newdiv);
            counter++;
            addNumLayers2(divName, numInputs);
        } else if ((counter + 1) < layersToAdd) {
            $("#h" + counter).remove();
            // document.getElementsByClassName("h"+counter).remove();
            var newdiv = document.createElement('div');
            newdiv.innerHTML = "<div id='h" + counter + "'>Layer" + (counter) + " <br><b>Epsilon:</b> <input type='text' size='6' id='e" + counter + "' class='e" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Mu:</b> <input type= 'text' size='6'  id = 'm" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Length:</b> <input type = 'text' size='6' id = 'l" + counter + "'> </div>";
            document.getElementById(divName).appendChild(newdiv);
            counter++;
            addNumLayers2(divName, numInputs);
        } else if ((counter + 1) == layersToAdd) {
            // document.getElementsByClassName("h"+counter).remove();
            $("#h" + counter).remove();
            var newdiv = document.createElement('div');
            newdiv.innerHTML = "<div id='h" + counter + "'>Ambient Right" + " <br><b>Epsilon:</b> <input type='text' size='6' id='e" + counter + "' class='e" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Mu:</b> <input type= 'text' size='6'  id = 'm" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Length:</b> <input type = 'text' size='6' id = 'l" + counter + "'> </div>";
            document.getElementById(divName).appendChild(newdiv);
            counter = 0;
        }
    }
}

function addNumLayers(divName, numInputs) {
    numberOfLayers = parseInt(numInputs) + 2;
    if (counter < numberOfLayers) {

        var newdiv = document.createElement('div');
        newdiv.innerHTML = "Layer " + (counter + 1) + " <br><b>Epsilon:</b> <input type='text' size='6' id='e" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Mu:</b> <input type= 'text' size='6'  id = 'm" + counter + "'> &nbsp;&nbsp;&nbsp; <b>Length:</b> <input type = 'text' size='6' id = 'l" + counter + "'>";
        document.getElementById(divName).appendChild(newdiv);
        counter++;
        if (counter < numberOfLayers) {
            addNumLayers(divName, numInputs);
        }


    }
    var myElements = document.querySelectorAll(".hiddenButtons");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }
}


function fillSampleValues(numOfLayers) {
    eps = [1, 5, 10, 5, 1];
    muu = [1, 10, 1, 10, 1];
    leng = [6, 3, 3, 3, 6]; //error if ambient size is below six or if layers are size 1
    om = 1;
    kk1 = .2;
    kk2 = .4;
    jj1 = 1;
    jj2 = 0;
    jj3 = 0;
    jj4 = 0;
    numberOfLayers = parseInt(numOfLayers) + 2;
    for (var i = 0; i < numberOfLayers; i++) {
        document.getElementById("e" + i).value = eps[i];
        document.getElementById("m" + i).value = muu[i];
        document.getElementById("l" + i).value = leng[i];
    }
    document.getElementById("omega").value = om;
    document.getElementById("k1").value = kk1;
    document.getElementById("k2").value = kk2;
    document.getElementById("j1").value = jj1;
    document.getElementById("j2").value = jj2;
    document.getElementById("j3").value = jj3;
    document.getElementById("j4").value = jj4;
}

function printStructureChart(htmlClass, divName, width, height) {

      numOfLayers = parseInt(document.getElementById("numLayers").value);
    layers = numOfLayers;
    epsilon = [];
    mu = [];
    length = [];
    var totalLength = 0;
    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
        totalLength = length[i]+totalLength;
    }
    var o = parseFloat(document.getElementById("omega").value);
    var k1 = parseFloat(document.getElementById("k1").value);
    var k2 = parseFloat(document.getElementById("k2").value);
    var j1 = parseFloat(document.getElementById("j1").value);
    var j2 = parseFloat(document.getElementById("j2").value);
    var j3 = parseFloat(document.getElementById("j3").value);
    var j4 = parseFloat(document.getElementById("j4").value);

    console.log("o:" + o + " k1:" + k1 + " k2:" + k2 + " j1:" + j1 + " j2:" + j2 + " j3:" + j3 + " j4:" + j4);
    crystal1 = new emScattering.PhotonicStructure1D(epsilon, mu, length);

    fields = crystal1.determineField(o, k1, k2, [j1, j2, j3, j4]);
    var interfaces = crystal1.materialInterfaces();
    var interfaceLength = interfaces[interfaces.length-1];

    var data = google.visualization.arrayToDataTable([
          [{
              f: 'Date',
              type: 'number' // wont work whithout this
          }, {
              f: 'Line',
              type: 'number' // wont work whithout this
          }], ]);

      var options = {

        chart: {
            title: 'dispersion relationship'
        },
        chartArea: {
            left: 40,
            top: 40
        },
        width: width,
        height: height,
        hAxis: {
            gridlines: {count: interfaceLength },
            viewWindow: {
                min: 0,
                max: totalLength//modify max by adding layer lengths, then apply to buildstructure button, then set up value placement of layers, then custom overlays
            }
        },
        vAxis: {
            gridlines: {count: 10},
            viewWindow: {
                min: -1,
                max: 1
            }
        }
    };

      if (data.getNumberOfRows() == 0) { // if you have no data, add a data point and make the series transparent
          data.addRow([0, 0])
          options.series = {
              0: {
                  color: 'transparent'
              }
          }
      }
     function printInterfaces(dataTable) {
        var cli = this.getChartLayoutInterface();
        var chartArea = cli.getChartAreaBoundingBox();
        var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
        var w = cli.getXLocation(interfaces[1]) - cli.getXLocation(interfaces[0]);
        var y = cli.getChartAreaBoundingBox().height;
        console.log(interfaces);
        console.log("w:" + w);
        console.log("y:" + y);
        console.log("gety:" + Math.floor(cli.getYLocation(1)));
        console.log("bounding:" + cli.getChartAreaBoundingBox().top);
        var yBound = cli.getChartAreaBoundingBox().top;
        //       Element.prototype.remove = function() {
        // this.parentElement.removeChild(this);
        // }
        // //allows removal of elements without parents
        // NodeList.prototype.remove = HTMLCollection.prototype.remove = function() {
        // for(var i = this.length - 1; i >= 0; i--) {
        //     if(this[i] && this[i].parentElement) {
        //         this[i].parentElement.removeChild(this[i]);
        //     }
        //  }
        //  }
        if(htmlClass=="addHere")
            $(".overlay").remove();
        if(htmlClass=="addHere2")
            $(".overlay2").remove();
        if(htmlClass=="addHere3")
            $(".overlay3").remove();
        for (var i = 0; i < interfaces.length - 1; i++) {
            var w = cli.getXLocation(interfaces[i + 1]) - cli.getXLocation(interfaces[i]);

            // document.getElementsByClassName('overlay' + i).remove();
            addStruct(htmlClass, i, cols[i % 5], w, y);

            // document.querySelector('.overlay').style.position = 'absolute';
            // document.querySelector('.overlay').style.opacity = '.5';
            // document.querySelector('.overlay').style.top = Math.floor(cli.getChartAreaBoundingBox().top) + "px";
            // document.querySelector('.overlay').style.left = Math.floor(cli.getXLocation(interfaces[i])) + 15 + "px";
        };

        if(htmlClass=="addHere"){
            var overlays = document.querySelectorAll('.overlay');
             for (j = 0; j < overlays.length; j++) {
            overlays[j].style.position = 'absolute';
            overlays[j].style.opacity = '.5';
            overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 141 + "px";
            overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 17 + "px";
            }
        }else if(htmlClass=="addHere2") {
            var overlays = document.querySelectorAll('.overlay2');
            for (j = 0; j < overlays.length; j++) {
            overlays[j].style.position = 'absolute';
            overlays[j].style.opacity = '.5';
            overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 15 + "px";
            overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 15 + "px";
            }


        }
        if(htmlClass=="addHere3") {
            var overlays = document.querySelectorAll('.overlay3');
            for (j = 0; j < overlays.length; j++) {
            overlays[j].style.position = 'absolute';
            overlays[j].style.opacity = '.5';
            overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 87 + "px";
            overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 17 + "px";
            }
        }




    }

    var chart = new google.visualization.LineChart(document.getElementById(divName));
    google.visualization.events.addListener(chart, 'ready', printInterfaces.bind(chart, data));

    chart.draw(data, options);

}

function printTransmissionChart(divName, width, height) {
    numOfLayers = parseInt(document.getElementById("numLayers").value);
    layers = numOfLayers;
    epsilon = [];
    mu = [];
    length = [];
    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
    }

    var crystal = new emScattering.PhotonicStructure1D(epsilon, mu, length);
    var k1 = parseFloat(document.getElementById("k1").value);
    var k2 = parseFloat(document.getElementById("k2").value);
    var omegLow = 0; //document.getElementById("oLow").value; //hardcode
    var omegHigh = 5; //document.getElementById("oHi").value; //hardcode
    var kZsList = [1, 2, 3];
    // var kZs = document.getElementById("kzList").value;
    // for(var i = 0; i < kZs.length; i++) {
    //  console.log(kZs.charAt(i));
    //  if(isNumeric(kZs.charAt(i))) {
    //    var str = ""+kZs.charAt(i);
    //    kZsList.push(parseInt(str));
    //  }
    // }
    console.log(kZsList);
    var transmissionGraph = crystal.transmission(kZsList, k1, k2, omegLow, omegHigh, 100);
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'omega');
    for (var i = 0; i < transmissionGraph.kzList.length; i++) {
        data.addColumn('number', 'kz' + transmissionGraph.kzList[i]);
    }
    var dataArray = new Array(transmissionGraph.omegaRange.length);
    for (var i = 0; i < dataArray.length; i++) {
        dataArray[i] = new Array(transmissionGraph.kzList.length);
    }

    for (var i = 0; i < transmissionGraph.omegaRange.length; i++) {
        dataArray[i][0] = transmissionGraph.omegaRange[i];
        for (var j = 1; j <= transmissionGraph.kzList.length; j++) {
            dataArray[i][j] = transmissionGraph.transmissionCoeffArrays[j - 1][i];
        }
    }


    for (var i = 0; i < dataArray.length; i++) {
        data.addRows([
            dataArray[i]
        ]);
    }
    var options = {
        chart: {
            title: 'dispersion relationship'
        },
        chartArea: {
            left: 40,
            top: 5
        },
        width: width,
        height: height
    };

    var chart = new google.visualization.LineChart(document.getElementById(divName));

    chart.draw(data, options);
    console.log(dataArray[0][1]);
    //emScattering.printDispersion(dispersion);

    var myElements = document.querySelectorAll(".hiddenChart1");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }



}

function printDispersionChart(divName, width, height) {
    numOfLayers = parseInt(document.getElementById("numLayers").value);
    layers = numOfLayers;
    epsilon = [];
    mu = [];
    length = [];
    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
    }
    var range = 4; //document.getElementById("range").value; //hardcode
    var crystal = new emScattering.PhotonicStructure1D(epsilon, mu, length);
    var k1 = parseFloat(document.getElementById("k1").value);
    var k2 = parseFloat(document.getElementById("k2").value);
    dispersion = crystal.dispersionRelationship(k1, k2, range, 100);
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'kz');
    for (var i = 0; i < dispersion.layersDispersions.length; i++) {
        data.addColumn('number', 'Layer ' + i);
    }
    var dataArray = new Array(dispersion.kz.length);
    for (var i = 0; i < dataArray.length; i++) {
        dataArray[i] = new Array(dispersion.layersDispersions.length + 1);
    }

    for (var i = 0; i < dispersion.kz.length; i++) {
        dataArray[i][0] = dispersion.kz[i];
        for (var j = 1; j < dispersion.layersDispersions.length + 1; j++) {
            dataArray[i][j] = dispersion.layersDispersions[j - 1][i];
            //this is where it fails
            //j not including first dispersion
        }
    }
    for (var i = 0; i < dataArray.length; i++) {
        data.addRows([
            dataArray[i]
        ]);
    }
    var options = {
        chart: {
            title: 'dispersion relationship'
        },
        chartArea: {
            left: 40,
            top: 5
        },
        width: width,
        height: height,
        vAxis: {
            viewWindow: {
                // min: -50,
                // max: 50
            }
        }
    };

    var chart = new google.visualization.LineChart(document.getElementById(divName));

    chart.draw(data, options);
    console.log(dataArray[0][1]);
    //emScattering.printDispersion(dispersion);

    var myElements = document.querySelectorAll(".hiddenChart1");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }



}

function printFieldsChart(htmlClass, divName, width, height) {
    numOfLayers = parseInt(document.getElementById("numLayers").value);
    layers = numOfLayers;
    epsilon = [];
    mu = [];
    length = [];
    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
    }
    var o = parseFloat(document.getElementById("omega").value);
    var k1 = parseFloat(document.getElementById("k1").value);
    var k2 = parseFloat(document.getElementById("k2").value);
    var j1 = parseFloat(document.getElementById("j1").value);
    var j2 = parseFloat(document.getElementById("j2").value);
    var j3 = parseFloat(document.getElementById("j3").value);
    var j4 = parseFloat(document.getElementById("j4").value);

    console.log("o:" + o + " k1:" + k1 + " k2:" + k2 + " j1:" + j1 + " j2:" + j2 + " j3:" + j3 + " j4:" + j4);
    crystal1 = new emScattering.PhotonicStructure1D(epsilon, mu, length);

    fields = crystal1.determineField(o, k1, k2, [j1, j2, j3, j4]);
    var interfaces = crystal1.materialInterfaces();



    var data = new google.visualization.DataTable();
    data.addColumn('number', 'z');
    data.addColumn('number', document.getElementById("shownVal").value);
    //Iterate through fields values

    for (var i = 0, N = fields.z.length; i < N; i++) {
        if (document.getElementById("shownVal").value == 'Ex') {
            data.addRows([
                [fields.z[i], fields.Ex[i]]
            ]);
        }
        if (document.getElementById("shownVal").value == 'Ey') {
            data.addRows([
                [fields.z[i], fields.Ey[i]]
            ]);
        }
        if (document.getElementById("shownVal").value == 'Hx') {
            data.addRows([
                [fields.z[i], fields.Hx[i]]
            ]);
        }
        if (document.getElementById("shownVal").value == 'Hy') {
            data.addRows([
                [fields.z[i], fields.Hy[i]]
            ]);
        }
        // console.log(fields.z[i] + " " + fields.Ex[i] + " " + fields.Ey[i] + " " + fields.Hx[i] + " " + fields.Hy[i]);
    }

    var options = {
        chart: {
            title: document.getElementById("shownVal").value + ' Values in Relation to Z'
        },
        width: width,
        height: height,
        chartArea: {
            left: 40,
            top: 40
        }

    };

    function printInterfaces(dataTable) {
        var cli = this.getChartLayoutInterface();
        var chartArea = cli.getChartAreaBoundingBox();
        var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
        var w = cli.getXLocation(interfaces[1]) - cli.getXLocation(interfaces[0]);
        var y = cli.getChartAreaBoundingBox().height;
        console.log(interfaces);
        console.log("w:" + w);
        console.log("y:" + y);
        console.log("gety:" + Math.floor(cli.getYLocation(1)));
        console.log("bounding:" + cli.getChartAreaBoundingBox().top);
        var yBound = cli.getChartAreaBoundingBox().top;
        //       Element.prototype.remove = function() {
        // this.parentElement.removeChild(this);
        // }
        // //allows removal of elements without parents
        // NodeList.prototype.remove = HTMLCollection.prototype.remove = function() {
        // for(var i = this.length - 1; i >= 0; i--) {
        //     if(this[i] && this[i].parentElement) {
        //         this[i].parentElement.removeChild(this[i]);
        //     }
        //  }
        //  }
       if(htmlClass=="addHere2")
            $(".overlay2").remove();
        else if(htmlClass=="addHere")
            $(".overlay").remove();
        else if(htmlClass=="addHere3")
            $(".overlay3").remove();
        for (var i = 0; i < interfaces.length - 1; i++) {
            var w = cli.getXLocation(interfaces[i + 1]) - cli.getXLocation(interfaces[i]);

            // document.getElementsByClassName('overlay' + i).remove();
            addStruct(htmlClass, i, cols[i % 5], w, y);

            // document.querySelector('.overlay').style.position = 'absolute';
            // document.querySelector('.overlay').style.opacity = '.5';
            // document.querySelector('.overlay').style.top = Math.floor(cli.getChartAreaBoundingBox().top) + "px";
            // document.querySelector('.overlay').style.left = Math.floor(cli.getXLocation(interfaces[i])) + 15 + "px";
        };

        if(htmlClass=="addHere"){
            var overlays = document.querySelectorAll('.overlay');
             for (j = 0; j < overlays.length; j++) {
            overlays[j].style.position = 'absolute';
            overlays[j].style.opacity = '.5';
            overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 141 + "px";
            overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 16 + "px";

            }
        }else if(htmlClass=="addHere2") {
            var overlays = document.querySelectorAll('.overlay2');
            for (j = 0; j < overlays.length; j++) {
            overlays[j].style.position = 'absolute';
            overlays[j].style.opacity = '.5';
            overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 15 + "px";
            overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 30 + "px";
            }
        }



    }

    var chart = new google.visualization.LineChart(document.getElementById(divName));
    google.visualization.events.addListener(chart, 'ready', printInterfaces.bind(chart, data));
    //var chart = new google.charts.Line(document.getElementById('linechart_material'));
    chart.draw(data, options);

    var myElements = document.querySelectorAll(".hiddenChart");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }

}



function getIncomingMode() {
    numOfLayers = parseInt(document.getElementById("numLayers").value);
    layers = numOfLayers;
    epsilon = [];
    mu = [];
    length = [];
    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
    }
    omega = document.getElementById("omega").value;
    k1 = document.getElementById("k1").value;
    k2 = document.getElementById("k2").value;
    crystal = new emScattering.PhotonicStructure1D(epsilon, mu, length);
    incoming = crystal.incomingModes(omega, k1, k2);
    mBack1x = parseFloat(incoming[0].eigenvalue.x.toFixed(8));
    mBack1y = parseFloat(incoming[0].eigenvalue.y.toFixed(8));
    mBack2x = parseFloat(incoming[1].eigenvalue.x.toFixed(8));
    mBack2y = parseFloat(incoming[1].eigenvalue.y.toFixed(8));
    mFor1x = parseFloat(incoming[2].eigenvalue.x.toFixed(8));
    mFor1y = parseFloat(incoming[2].eigenvalue.y.toFixed(8));
    mFor2x = parseFloat(incoming[3].eigenvalue.x.toFixed(8));
    mFor2y = parseFloat(incoming[3].eigenvalue.y.toFixed(8));

    var newdiv = document.createElement('div');
    var newdiv2 = document.createElement('div');
    $("#mode1").remove();
    $("#mode2").remove();
    newdiv.innerHTML = "<div  id='mode1'><i class='fa fa-arrow-right' aria-hidden='true'></i><u>" + mFor1x.toString() + " + " + mFor1y.toString() + "i<br><i class='fa fa-arrow-right' aria-hidden='true'></i>" + mFor2x + " + " + mFor2y + "i</u></div>";
    newdiv2.innerHTML = "<div  id='mode2'><u>" + mBack1x + " + " + mBack1y + "i<i class='fa fa-arrow-left' aria-hidden='true'></i><br> " + mBack2x + " + " + mBack2y + "i</u><i class='fa fa-arrow-left' aria-hidden='true'></i></div>";
    document.getElementById("forwardModes").appendChild(newdiv);
    document.getElementById("backwardModes").appendChild(newdiv2);
    //emScattering.printModes(incoming);

    var myElements = document.querySelectorAll("#hiddenModes");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }

}


function buildCrystal() {
    layers = parseInt(numOfLayers);

    for (var i = 0; i < layers; i++) {
        epsilon[i] = parseInt(document.getElementById("e" + i).value);
        mu[i] = parseInt(document.getElementById("m" + i).value);
        length[i] = parseInt(document.getElementById("l" + i).value);
    }
    var myElements = document.querySelectorAll(".hiddentab");
    for (var i = 0; i < myElements.length; i++) {
        myElements[i].style.opacity = 1;
    }
    console.log(epsilon);
    console.log(mu);
    console.log(length);
    crystal = new emScattering.PhotonicStructure1D(epsilon, mu, length);
    //alert("Variables assigned. Crystal created; layers including ambients = "+layers);
}
