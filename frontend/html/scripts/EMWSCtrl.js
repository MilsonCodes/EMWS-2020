// code style: https://github.com/johnpapa/angular-styleguide 
angular.module('myApp', []).controller('EMWSCtrl', function($scope) {






        google.charts.load('current', { packages: ['corechart'] });             //Loads Google Charts packages
        $scope.context;
        $scope.NumLayers = 3;                                                   //Number of Layers (Default)
        $scope.o = 0.398;                                                           //Omega (Default)
        $scope.k1 = 0.5;                                                          //k1 (Default)
        $scope.k2 = 0.22;                                                          //k2 (Default)
        $scope.wLeft = -5;                                                      //Left bound for transmissions graph (Default)
        $scope.wRight = 5;                                                      //Right bound for transmissions graph (Default)
        $scope.wPoints = 10;                                                   //Number of points for transmissions graph (Default)
        $scope.zPoint = 0;                                                      //Point for each array to pull for transmissions graph
        $scope.incoming = [1, 0, 0, 0];                                         //Incoming coefficients (Defaults)
        $scope.eArray = [];                                                     //Epsilon Array for Layers
        $scope.muArray = [];                                                    //Mu Array for Layers
        $scope.lArray = [];                                                     //Length Array for Layers
        $scope.totalLength = 0;                                                 //Total Length of All Layers
        $scope.outputModes;
        $scope.modesBack = new Array();                                         //Array of the backward modes for the GUI
        $scope.modesForward = new Array();                                      //Array of the forward modes for the GUI

        /*          OLD VARIABLES FOR MODES - ARRAYS USED NOW
        $scope.mBack1x;                                                         //Real value used in modes (for all mBackix and mForix values, i being a value of 1-4)
        $scope.mBack1y;                                                         //Imaginary value used in modes (for all mBackiy and mForiy values, i being a value of 1-4)
        $scope.mBack2x;
        $scope.mBack2y;
        $scope.mBack3x; 
        $scope.mBack3y;
        $scope.mBack4x;
        $scope.mBack4y;
        $scope.mFor1x;
        $scope.mFor1y;
        $scope.mFor2x;
        $scope.mFor2y;
        $scope.mFor3x;
        $scope.mFor3y;
        $scope.mFor4x;
        $scope.mFor4y;
        $scope.mFor1;                                                           //Complex mode (for all mBacki and mFori, i being a value of 1-4)
        $scope.mFor2;
        $scope.mFor3;
        $scope.mFor4;
        $scope.mBack1; 
        $scope.mBack2;
        $scope.mBack3; 
        $scope.mBack4;
        */

        $scope.crystal;                                                         //The Photonic Crystal created in emScattering3.js
        $scope.field;                                                           //The field determined using the Photonic Crystal
        $scope.dispersion;                                                      
        $scope.EX = 'Eₓ';                                                       //Label for Ex
        $scope.EY = 'Eᵧ';                                                       //Label for Ey
        $scope.HX = 'Hₓ';                                                       //Label for Hx
        $scope.HY = 'Hᵧ';                                                       //Label for Hy


        //Defined default values of the layers
        $scope.Layers = [{
                "layerName": "Ambient Left",
                "epsilon": [1, 2, 3, 4],
                "epsilonA": [[1.5, 0, 0], [0, 8, 0], [0, 0, 1]],
                "mu": [1, 2, 3, 4],
                "muA": [[4,0,0],[0,1,0], [0, 0, 1]],
                "length": 10
            },
            {
                "layerName": "Layer 1",
                "epsilon": [2, 2, 3, 6],
                "epsilonA": [[8, 0, 0], [0, 1.5, 0], [0, 0, 1]],
                "mu": [2, 2, 3, 6],
                "muA": [[1,0,0],[0,4,0], [0, 0, 1]],
                "length": 7
            },
            {
                "layerName": "Ambient Right",
                "epsilon": [1, 2, 3, 4],
                "epsilonA": [[1.5, 0, 0], [0, 8, 0], [0, 0, 1]],
                "mu": [1, 2, 3, 4],
                "muA": [[4,0,0],[0,1,0], [0, 0, 1]],
                "length": 10
            },
        ];

        /**Initializes the website application. */
        $scope.init = function() {     
            getArrays();
            updateAll();
            google.charts.setOnLoadCallback(createStructureChart);
            google.charts.setOnLoadCallback(createFieldChart);
            google.charts.setOnLoadCallback(createTransmissionChart)
            //google.charts.setOnLoadCallback(createDispersionChart);
            //google.charts.setOnLoadCallback(createTransmissionChart);
            createAnim();
            //printDispersionChart("dispView", "100%", "100%");
            
            //$(".p2").css("visibility", "hidden");



            //$(".p1").addClass("ng-hide");

            setTimeout(hideIfInDevMode, 250);
        }

        function getQueryVariable(variable)
        {
            var query = window.location.search.substring(1);
            var vars = query.split("&");

            for (var i = 0;i < vars.length; i++) {
                var pair = vars[i].split("=");

                if(pair[0] == variable) return pair[1];
            }
            
            return false;
        }

        function hideIfInDevMode() {
            var devMode = true;        //Change if making developmental changes!!

            if(getQueryVariable("devMode") === "true") devMode = true;

            $scope.runMathBoxField = devMode;

            if(!devMode) {
                document.getElementById("mathbox-field-box").classList.add("ng-hide");
                document.getElementById("nav-tab3").classList.add("ng-hide");
                document.getElementById("nav-tab4").classList.add("ng-hide");

                document.getElementsByClassName("dropdown-menu")[0].getElementsByClassName("col-xs-6")[2].classList.add("ng-hide");
                document.getElementsByClassName("dropdown-menu")[0].getElementsByClassName("col-xs-6")[3].classList.add("ng-hide");
            }
        }

        function addHide() {
            // $(".p2").addClass("ng-hide");
            // $(".p3").addClass("ng-hide");
            // $('.p2').css('visibility', 'visible');
            // $('.p3').css('visibility', 'visible');
        }

        /** Adds a layer to the current set of layers. */
        $scope.addLayer = function() {
            // getIncomingMode();
            $scope.Layers[$scope.NumLayers - 1].layerName = ("Layer " + ($scope.NumLayers - 1));    //Renames last ambient right layer 
            $scope.NumLayers++;                                                                     //Adds to the number of layers
            $scope.Layers.push({
                "layerName": "Ambient Right",
                "epsilon": [1, 2, 3, 4],
                "epsilonA": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "mu": [1, 2, 3, 4],
                "muA": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                "length": 6
            });                                                                                     //Pushes layer to end
        };

        /** Removes a layer from the current set of layers. */
        $scope.removeLayer = function() {
            $scope.Layers.pop();                                                                    //Removes last layer
            $scope.NumLayers--;                                                                     //Subtracts to the number of layers
            $scope.Layers[$scope.NumLayers - 1].layerName = ("Ambient Right");                      //Renames last layer to Ambient Right
        }

        //Not called
        $scope.printLengths = function() {
            angular.forEach($scope.Layers, function(value, index) {

                console.log("length" + value.length);
            })
        }

        /** Updates current set of layers. */
        $scope.updateLayers = function() {
            if ($scope.NumLayers == null || $scope.NumLayers == 0) { console.log("scope is null or 0"); } else if ($scope.NumLayers > $scope.Layers.length) {
                if ($scope.Layers[$scope.NumLayers - 1]) {
                    $scope.Layers[$scope.NumLayers - 1].layerName = 
                            ("Layer " + ($scope.NumLayers - 1)); } 
                else 
                    $scope.Layers[$scope.Layers.length - 1].layerName = 
                        ("Layer " + ($scope.Layers.length - 1));
                if ($scope.NumLayers - 2 == 0) {
                    $scope.Layers[$scope.NumLayers - 2].layerName = "Ambient Left"; }
                $scope.Layers.push({
                    "layerName": "Ambient Right",
                    "epsilon": [1, 2, 3, 4],
                    "epsilonA": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "mu": [1, 2, 3, 4],
                    "muA": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                    "length": 6
                })
                if ($scope.NumLayers > $scope.Layers.length)
                    $scope.updateLayers();
            } else if ($scope.NumLayers < $scope.Layers.length) {
                $scope.Layers.pop();
                $scope.Layers[$scope.NumLayers - 1].layerName = ("Ambient Right");
                if ($scope.NumLayers < $scope.Layers.length)
                    $scope.updateLayers();
            }
        }
        
        /** Builds chart for the Structure tab. */
        $scope.buildStruct = function() {
            getArrays();
            updateCrystal();
            createStructureChart();
        }

        /** Calculates the modes. */
        $scope.buildModes = function() {
            getArrays();
            updateAll();
        }

        /** Checks the state of the check boxes in the Field tabs and sets the order of the selected modes.
         * 
         * @param fieldTab - Boolean if tab being checked is fieldTab (ONLY USE TRUE OR FALSE)
         */
        $scope.checkBoxes = function(fieldTab) {
            var backChecked = 0;                //Create variable for amount of leftward mode boxes checked
            var forChecked = 0;                 //Create variable for amount of rightward mode boxes checked

            var backChkStr = "backModeChk";     //Create variable for base of element id names
            var forChkStr = "forModeChk";       //Create variable for base of element id names
            var incStr = "incoming";            //Create variable for base of element id names

            if(!fieldTab) {
                backChkStr = backChkStr + "T";          //If tab is the transmission tab being checked, add T
                forChkStr = forChkStr + "T";            //If tab is the transmission tab being checked, add T
                incStr = incStr + "T";                  //If tab is the transmission tab being checked, add T
            }

            for(let i = 1; i <= 4; i++){
                if(document.getElementById(backChkStr + i).checked == true) backChecked++;   //If box is checked, add one to amount of checked boxes
                if(document.getElementById(forChkStr + i).checked == true) forChecked++;     //If box is checked, add one to amount of checked boxes
            }

            //Check if exactly two boxes are checked.
            if(backChecked == 2){
                let j = 0;                                              //Counter for setting text in coefficients section (should be only 0 and then 1)

                var backArr = $scope.crystal.Struct.eigenvalues[0];     //Gets eigenvalues of first layer
                var backVecArr = $scope.crystal.Struct.eigenvectors[0]._data;

                //Loops through check boxes
                for(let i = 1; i <= 4; i++){
                    if(document.getElementById(backChkStr + i).checked == false) { document.getElementById(backChkStr + i).disabled = true; }         //If check box is not checked, disable it
                    else {                                                                                                                                  //Runs this code if the box is checked
                        //arr.splice(new_index, 0, arr.splice(old_index, 1)[0]);

                        backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);            //Swap array around based on order of checked boxes
                        backVecArr.splice(j, 0, backVecArr.splice(i-1, 1)[0]);      //Swap eigenvectors around based on order of checked boxes

                        //Depending on i value, change text of incomingj element to that mode
                        document.getElementById(incStr + j).innerHTML = $scope.modesBack[i-1].toString();

                        /*if(i == 1) {
                            backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mBack1.toString();
                        } else if(i == 2) {
                            backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mBack2.toString();
                        } else if(i == 3) {
                            backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mBack3.toString();
                        } else if(i == 4) {
                            backArr.splice(j, 0, backArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mBack4.toString();
                        }*/

                        j++;                                                        //Increase j counter for next loop
                    }
                }

                //console.log(backArr);    //Logs first layer eigenvalues (commented out, uncomment when needed)
            }else if(backChecked < 2){
                //Runs loop through check boxes
                for(let i = 1; i <= 4; i++){
                    if(document.getElementById(backChkStr + i).disabled == true) { document.getElementById(backChkStr + i).disabled = false; }            //If check box is disabled, enable it.
                }

                document.getElementById(incStr + "0").innerHTML = $scope.modesBack[0].toString();                          //Set incoming0 element to default mode
                document.getElementById(incStr + "1").innerHTML = $scope.modesBack[1].toString();                          //Set incoming1 element to default mode
            }

            //Check if exactly two boxes are checked
            if(forChecked == 2){
                let j = 2;                                                          //Counter for setting text in coefficients section (should be only 2 and then 3)

                var lastLayer = $scope.crystal.Struct.eigenvalues.length - 1;  //Gets value of last set of eigenvalues
                var forArr = $scope.crystal.Struct.eigenvalues[lastLayer];     //Gets eigenvalues of last layer
                var forVecArr = $scope.crystal.Struct.eigenvectors[lastLayer]._data;     //Gets eigenvalues of last layer

                //Loop through the check boxes
                for(let i = 1; i <= 4; i++){
                    if(document.getElementById(forChkStr + i).checked == false) { document.getElementById(forChkStr + i).disabled = true;  }      //If box is not checked, disable it
                    else {
                        forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);              //Swap array around based on order of checked boxes
                        forVecArr.splice(j, 0, forVecArr.splice(i-1, 1)[0]);        //Swap eigenvectors around based on order of checked boxes

                        //Depending on i value, change text of incomingj element to that mode
                        document.getElementById(incStr + j).innerHTML = $scope.modesForward[i-1].toString();

                        /*
                        if(i == 1){
                            forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mFor1.toString();
                        } else if(i == 2) {
                            forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mFor2.toString();
                        } else if(i == 3) {
                            forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mFor3.toString();
                        } else if(i == 4) {
                            forArr.splice(j, 0, forArr.splice(i-1, 1)[0]);
                            document.getElementById("incoming" + j).innerHTML = $scope.mFor4.toString();
                        }*/

                        j++;                                            //Increments j counter for next loop
                    }
                }

                //console.log(forArr);      //Logs last layer eigenvalues (commented out, uncomment when needed)
            }else if(forChecked < 2){
                //Loops through check boxes
                for(let i = 1; i <= 4; i++){
                    if(document.getElementById(forChkStr + i).disabled == true) { document.getElementById(forChkStr + i).disabled = false; }          //If check box is disabled, enable it
                }

                document.getElementById(incStr + "2").innerHTML = $scope.modesForward[2].toString();               //Set incoming2 element text to default mode
                document.getElementById(incStr + "3").innerHTML = $scope.modesForward[3].toString();               //Set incoming3 element text to default mode
            }
        }
        
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################


        /** Runs the experiment in the Field tab. */
        $scope.runExp = function() {
            $("canvas").remove();
            getArrays();
            updateAll();

            createFieldChart();
            createAnim();
            //$scope.buildFieldsWithAnim();
        };
        
        /** Updates the arrays based on the layers. */
        function getArrays() {
            $scope.totalLength = 0;
            for (var layer in $scope.Layers){
                $scope.eArray[layer] = ($scope.Layers[layer].epsilonA);
                $scope.muArray[layer] = ($scope.Layers[layer].muA);
                $scope.lArray[layer] = parseFloat($scope.Layers[layer].length);
                $scope.totalLength += parseFloat($scope.Layers[layer].length);
            }
        }

        /** Updates the modes using the structure's Eigensystems and eigenvalues. */
        function updateModes() {
            var lastEigensystem = $scope.crystal.Struct.Eigensystems.length - 1;
            var eigensystems = $scope.crystal.Struct.Eigensystems;

            for(let i = 0; i < 4; i++){
                $scope.modesBack[i] = math.complex(
                    parseFloat(math.re(eigensystems[lastEigensystem][i].eigenvalue)).toFixed(4),
                    parseFloat(math.im(eigensystems[lastEigensystem][i].eigenvalue)).toFixed(4)
                );
                  
                $scope.modesForward[i] = math.complex(
                    parseFloat(math.re(eigensystems[0][i].eigenvalue)).toFixed(4),
                    parseFloat(math.im(eigensystems[0][i].eigenvalue)).toFixed(4)
                );
            }
        }
        /** Updates the Photonic Crystal. */
        async function updateCrystal(){
            var k = [Number($scope.k1),Number($scope.k2),Number($scope.o)];
            $scope.crystal = emScattering3.Driver($scope.eArray,$scope.muArray,$scope.lArray,$scope.NumLayers,k,$scope.incoming);

            $scope.structure = backendAPI.createStructureObject($scope.o, $scope.k1, $scope.k2, $scope.Layers)

            await $scope.structure.buildStructure()

            console.log($scope.structure)
        };
        
        /** Updates the field using the Photonic Crystal. */
        function updateFields(){
            $scope.crystal.Struct.updateScattering();
            $scope.field = $scope.crystal.determineField();
        }
        
        /** Updates the crystal, modes, and fields, and checks the check boxes on the Field tab. */
        function updateAll(){
            updateCrystal();
            updateModes();
            $scope.checkBoxes(true);
            $scope.checkBoxes(false);
            updateFields();
        }
        
        /** Creates the line chart on the Field tab. */
        function createFieldChart() {
            console.log($scope.crystal);
            var fields = $scope.field;                                  //Takes field and puts it to a variable
            console.log(fields);
            var interfaces = $scope.crystal.Struct.materialInterfaces();       //Takes interfaces and puts it to a variable

            var data = new google.visualization.DataTable();            //Creates a data table
            data.addColumn('number', 'z');                              //Adds z column to data table
            //data.addColumn('number', document.getElementById("shownVal").value);
            data.addColumn('number', $scope.EX);                        //Adds Ex column to data table
            data.addColumn('number', $scope.EY);                        //Adds Ey column to data table
            data.addColumn('number', $scope.HX);                        //Adds Hx column to data table
            data.addColumn('number', $scope.HY);                        //Adds Hy column to data table
            
            //Iterate through fields values
            for (var i = 0, N = fields.z.length; i < N; i++) {
                data.addRows([
                    [fields.z[i], fields.Ex[i], fields.Ey[i], fields.Hx[i], fields.Hy[i]]
                ]);
            }

            function printInterfaces(dataTable) { //prints the colored squares on the top of the chart
                var cli = this.getChartLayoutInterface();
                var chartArea = cli.getChartAreaBoundingBox();
                var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
                var oddColors = '#FFFF00';
                var evenColors = '#888888';
                var w = cli.getXLocation(interfaces[1]) - cli.getXLocation(interfaces[0]);
                var y = cli.getChartAreaBoundingBox().height;
                // console.log(interfaces);
                // console.log("w:" + w);
                // console.log("y:" + y);
                // console.log("gety:" + Math.floor(cli.getYLocation(1)));
                // console.log("bounding:" + cli.getChartAreaBoundingBox().top);
                var yBound = cli.getChartAreaBoundingBox().top;
              
                var htmlClass = "addHere"; //confusing set of addHere variables. Essentially, they were set up originally to allow for me to add the same chart to two different places. 
                                            //Will work to remove in the next couple weeks.
                if (htmlClass == "addHere2")
                    $(".overlay2").remove();
                else if (htmlClass == "addHere")
                    $(".overlay").remove();
                else if (htmlClass == "addHere3")
                    $(".overlay3").remove();
                for (var i = 0; i < interfaces.length - 1; i++) {
                    var w = cli.getXLocation(interfaces[i + 1]) - cli.getXLocation(interfaces[i]);

                    // document.getElementsByClassName('overlay' + i).remove();
                    addStruct(htmlClass, i, (i % 2 === 0 ? evenColors : oddColors), w, y);
                    
                    // document.querySelector('.overlay').style.position = 'absolute';
                    // document.querySelector('.overlay').style.opacity = '.5';
                    // document.querySelector('.overlay').style.top = Math.floor(cli.getChartAreaBoundingBox().top) + "px";
                    // document.querySelector('.overlay').style.left = Math.floor(cli.getXLocation(interfaces[i])) + 15 + "px";

                    //axisTicks[j] = (cli.getXLocation(interfaces[j]) / 1511);
                };

                if (htmlClass == "addHere") {
                    var overlays = document.querySelectorAll('.overlay');
                    for (var j = 0; j < overlays.length; j++) {
                        overlays[j].style.position = 'absolute';
                        overlays[j].style.pointerEvents = 'none';
                        overlays[j].style.opacity = '.5';
                        overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 65 + "px";
                        overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 16 + "px";
                    }
                } else if (htmlClass == "addHere2") {
                    var overlays = document.querySelectorAll('.overlay2');
                    for (var j = 0; j < overlays.length; j++) {
                        overlays[j].style.position = 'absolute';
                        overlays[j].style.pointerEvents = 'none';
                        overlays[j].style.opacity = '.5';
                        overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 15 + "px";
                        overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 30 + "px";
                    }
                }
            }
            
            var chart = new google.visualization.LineChart(document.getElementById('structView'));      //Create line chart in structView section
            google.visualization.events.addListener(chart, 'ready', printInterfaces.bind(chart, data)); //all from Google's tutorial https://developers.google.com/chart/interactive/docs/overlays
            //var chart = new google.charts.Line(document.getElementById('linechart_material'));

            //Chart Options
            var options = {
                chart: {
                    title: 'Ex, Ey, Hx, and Hy Values in Relation to Z'
                },
                hAxis: {
                    gridlines: {
                        color: 'transparent'  
                    },
                    ticks: interfaces
                },
                vAxis:{
                    gridlines: {
                        color: 'transparent'  
                    },
                    viewWindowMode: 'maximized'
                },
                width: "100%",
                height: "100%",
                chartArea: {
                    left: 40,
                    top: 40
                },
                backgroundColor: 'transparent',
                colors: [ "red", "orange", "blue", "green" ]
            };

            chart.draw(data, options);          //Draw the chart

            //These lines make the elements with .hiddenChart class visible
            var myElements = document.querySelectorAll(".hiddenChart");
            for (var i = 0; i < myElements.length; i++) {
                myElements[i].style.opacity = 1;
            }

            document.getElementById('structView').children[0].style.zIndex = 1;         //Takes the chart div and puts it above the overlays

            var chartView = new google.visualization.DataView(data);                    //Creates a new data view variable based on the data

            function removeElem(arr, elem) {
                for(var i = 0; i < arr.length; i++) {
                    if(arr[i] === elem) arr.splice(i, 1);
                }
            }

            /** Shows and hides columns based on id.
             * 
             * @param id - The column id to show/hide.
             */
            $scope.toggleLine = function(id) {
                var columns = chartView.getViewColumns();               //Gets the array of visible columns
                var fullColorsArr = [ "red", "orange", "blue", "green" ];
                var isHidden = true;                                    //Creates flag if column specified is hidden

                //Loop through visible columns
                for(var i = 0; i < columns.length; i++){
                    if(columns[i] == id) isHidden = false;              //If column is visible, set flag to false to signify that it is visible
                }

                if(isHidden == true){
                    columns.splice(id, 0, id);                          //If column is not visible, add it back to the visible column array.
                    options.colors.splice(id - 1, 0, fullColorsArr[id - 1]);    //Add color back to visible arr
                    chartView.setColumns(columns);                      //Set visible columns to column array
                }else{
                    chartView.hideColumns([id]);                        //If column is visible, hide the column
                    removeElem(options.colors, fullColorsArr[id - 1]);
                }

                chart.draw(chartView, options);                         //Draw the updated chart
            }

            /** Resets the columns of the field chart. */
            $scope.resetFieldChart = function() {
                chartView.setColumns([0, 1, 2, 3, 4]);                  //Reset to show all visible columns
                options.colors = [ "red", "orange", "blue", "green" ];  //Reset colors
                chart.draw(chartView, options);                         //Draw the updated chart
            }
        }

        /** Creates the MathBox animation on the Field tab. */
        function createAnim() {
            $scope.runMathBoxField = true;
            var E1 = $scope.crystal.mathboxSetupEf();
            var H1 = $scope.crystal.mathboxSetupHf();
            // var hMaxCalc = $scope.crystal.mathboxHf($scope.crystal.Struct.lengths,1,1,H1.HxR,H1.HxPhi,H1.HyR,H1.HyPhi);
            // var eMaxCalc = $scope.crystal.mathboxEf($scope.crystal.Struct.lengths,1,1,E1.ExR,E1.ExPhi,E1.EyR,E1.EyPhi);
            // console.log(hMaxCalc); console.log(eMaxCalc);
            var eXmax = 1;
            var hXmax = 1;
            var endRange = $scope.totalLength;
            var canvasElement = "testcanvas";
            var interfaces = $scope.crystal.Struct.materialInterfacesStartZero();
            var elem = document.getElementById(canvasElement);
            var jelem = $("#" + canvasElement);
            var rgbColor = jelem.parent().css("background-color");
            var WIDTH = elem.offsetWidth;
            var HEIGHT = elem.offsetHeight;
            var w1 = WIDTH;
            var h1 = HEIGHT;
            var three = THREE.Bootstrap({
                plugins: ['core', 'controls'],
                controls: {
                    klass: THREE.OrbitControls
                },
                size: {
                    width: w1,
                    height: h1,
                },
            });            
            var renderer = three.renderer;
            var scene = three.scene;
            var camera = three.camera;
            // Insert into document
            elem.appendChild(renderer.domElement);
            // MathBox $scope.context
            $scope.context = new MathBox.Context(renderer, scene, camera).init();
            var mathbox = $scope.context.api;
            // Set size
            renderer.setSize(WIDTH, HEIGHT);
            $scope.context.resize({ viewWidth: WIDTH, viewHeight: HEIGHT });
            // Place camera and set background
            camera.position.set(-1, 2, 4);
            renderer.setClearColor(new THREE.Color(rgbColor), 1.0);
            var view = mathbox
                .set({
                    focus: 3,
                })
                .cartesian({
                    range: [
                        [0, endRange],
                        [-eXmax, eXmax],
                        [-hXmax, hXmax]
                    ],
                    scale: [2, 2, 2],
                });

            view.axis({
                detail: 30,
            });

            view.axis({
                axis: 2,
            });

            view.scale({
                divide: 10,
            })
            view.ticks({
                classes: ['foo', 'bar'],
                width: 2
            });

            view.scale({ //adds "X-Axis" to the graph
                divide: 1,
                origin: [0, 1, 0, 0],
                axis: "x",
            }).text({
                live: false,
                data: ["Electric Field"]
            }).label({
                color: 0x0074D9
            })

            view.scale({ //adds "Y-Axis" to the graph
                divide: 1,
                origin: [25, 0, 0, 0],
                axis: "y"
            }).text({
                live: false,
                data: ["Magnetic Field"]
            }).label({
                color: 0xFF4136
            })

            view.grid({
                    axes: "xy",
                    divideX: endRange,
                    divideY: 10
                }) //makes two axes in space
                .grid({
                    axes: "xz",
                    divideX: endRange,
                    divideY: 10
                })

            // var colorCoords = []; //possibly remove, replace with just applying interfaces to arrays
            var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
            for (var i = 0; i < interfaces.length - 1; i++) {
                var array1 = [ //only x changes on all shapes coordinates
                    //seperate 1st and last arrays. also, make each array part of an array[][]                    
                    0 + interfaces[0 + i], -2, 0, 0 + interfaces[1 + i], -2, 0, 0 + interfaces[1 + i], 2, 0, 0 + interfaces[0 + i], 2, 0, //first column
                ];
                var array2 = [0 + interfaces[0 + i], 0, 2, 0 + interfaces[0 + i], 0, -2, 0 + interfaces[1 + i], 0, -2, 0 + interfaces[1 + i], 0, 2, ];
                // console.log("Array1: " + array1 + "\nArray2: " + array2 + "\n")
                view.voxel({
                    data: array1,
                    width: 4,
                    height: 2,
                    depth: 1,
                    items: 4,
                    channels: 3,
                });

                view
                    .transform({
                        pass: 'eye',
                        position: [0, 0, 0],
                        scale: [1, 1, 1],
                    })
                    .face({
                        color: cols[i % 5],
                        width: 3,
                        line: false,
                        shaded: true,
                        opacity: .3,
                        zOrder: 9,
                    })
                    .face({
                        color: cols[i % 5],
                        width: 3,
                        fill: false,
                        opacity: .3,
                        zOrder: 9,
                    })
                    .array({
                        data: [0, 0, 0],
                        channels: 3,
                    })
                view.voxel({
                    data: array2,
                    width: 4,
                    height: 2,
                    depth: 1,
                    items: 4,
                    channels: 3,
                });

                view
                    .transform({
                        pass: 'eye',
                        position: [0, 0, 0],
                        scale: [1, 1, 1],
                    })
                    .face({
                        color: cols[i % 5],
                        width: 3,
                        line: false,
                        shaded: true,
                        opacity: .3,
                        zOrder: 9,
                    })
                    .face({
                        color: cols[i % 5],
                        width: 3,
                        fill: false,
                        line: true,
                        opacity: .3,
                        zOrder: 9,
                    })
                    .array({
                        data: [0, 0, 0],
                        channels: 3,
                    })
            }


            var runsetup = true;
            var E1,E2,H1,H2;
            view.interval({
                id: 'ElectricFieldPlot',
                width: endRange, //fields.Ex.length,
                expr: function(emit, z, i, t) {
                    if(z<=endRange){
                    if (runsetup) {
                        E1 = $scope.crystal.mathboxSetupEf();
                        runsetup = false;
                    };                    
                    E2 = $scope.crystal.mathboxEf($scope.crystal.Struct.lengths,t,z,E1.ExR,E1.ExPhi,E1.EyR,E1.EyPhi);
                    // console.log("z: " + z + "\nEx: " + E2.Ex + "\nEy: " + E2.Ey + "\n")
                    emit(z, E2.Ex, E2.Ey);
                    emit(z, 0, 0);
                    }
                },
                items: 2,
                channels: 3,
                live: true,
            });
            view.vector({
                points: '#ElectricFieldPlot',
                color: 0x0074D9,
                width: 1,
                start: true,
            });

            view.interval({
                id: 'MagneticFieldPlot',
                width: endRange, //fields.Ex.length,
                expr: function(emit, z, i, t) {
                    if(z <= endRange){
                    if (runsetup) {
                        H1 = $scope.crystal.mathboxSetupHf();
                        runsetup = false;
                    };
                    H2 = $scope.crystal.mathboxHf($scope.crystal.Struct.lengths,t,z,H1.HxR,H1.HxPhi,H1.HyR,H1.HyPhi);
                    // console.log("z: " + z + "\nHx: " + H2.Hx + "\nHy: " + H2.Hy + "\n")
                    emit(z, H2.Hx, H2.Hy);
                    emit(z, 0, 0);
                    }
                },
                items: 2,
                channels: 3,
                live: true,
            });

            view.vector({
                points: '#MagneticFieldPlot',
                color: 0xFF4136,
                width: 1,
                start: true,
            });
            var visible = false;
            var madeVisible = false;
            var frame = function() {
                var parentVisibility = jelem.parent().css("visibility");
                    //console.log("doings things: ");
                if(parentVisibility != 'hidden'){
                    requestAnimationFrame(frame);
                    visible = true;
                    
                }
                else if(parentVisibility == 'hidden' || !$scope.runMathBoxField){
                    visible = false;

                    requestAnimationFrame(frame);
                    renderer.domElement.style.visibility = parentVisibility;
                    return;
                }
                
                if(!madeVisible) {
                    renderer.domElement.style.visibility = parentVisibility; //wheretohide
                    madeVisible = true;
                }
                
//?
                if(!visible)
                renderer.setSize(WIDTH, HEIGHT);
                $scope.context.frame();
                renderer.render(scene, camera);
                rgbColor = jelem.parent().css("background-color"); //0xFF851B;
                renderer.setClearColor(new THREE.Color(rgbColor), 1.0);
            };
            requestAnimationFrame(frame);
        }

        
 
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################
        //#########################################################################################
        
        /** Builds the information for the Transmissions tab. WIP */
        $scope.buildTransmission = function() {
            getArrays();
            updateAll();
        }

        /** Runs the experiment in the Transmissions tab. WIP */
        $scope.runTransmissionExp = function() {
            getArrays();
            updateAll();

            createTransmissionChart();
        }

        /** Creates the chart in the Transmissions tab. WIP */
        function createTransmissionChart() {
            //var kZsList = [1, 2, 3];
            // var kZs = document.getElementById("kzList").value;
            // for(var i = 0; i < kZs.length; i++) {
            //  console.log(kZs.charAt(i));
            //  if(isNumeric(kZs.charAt(i))) {
            //    var str = ""+kZs.charAt(i);
            //    kZsList.push(parseInt(str));
            //  }
            // }
            var divName = "transmissionView";
            //console.log(kZsList);
            console.time("Transmission Graph");
            var transmission = emScattering3.createTransmissionArrays($scope.eArray, $scope.muArray, $scope.lArray, $scope.NumLayers, $scope.k1, $scope.k2, $scope.incoming ,$scope.wLeft, $scope.wRight, $scope.wPoints, $scope.zPoint);        //Method needs to be created in emScattering3!
            console.timeEnd("Transmission Graph");
            console.log(transmission);
            var data = new google.visualization.DataTable();
            data.addColumn('number', 'omega');                          //Adds Omega column to data table
            data.addColumn('number', 'c');                              //Adds const transmission to data table
            /*
            data.addColumn('number', $scope.EX);                        //Adds Ex column to data table
            data.addColumn('number', $scope.EY);                        //Adds Ey column to data table
            data.addColumn('number', $scope.HX);                        //Adds Hx column to data table
            data.addColumn('number', $scope.HY);                        //Adds Hy column to data table
            */

            /*
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
            */


            for (var i = 0; i < $scope.wPoints; i++) {
                data.addRows([
                    [transmission.omegas[i], transmission.cArr[i] /*transmission.Ex[i], transmission.Ey[i], transmission.Hx[i], transmission.Hy[i]*/]
                ]);
            }
            var options = {
                chart: {
                    title: 'transmission graph'
                },
                chartArea: {
                    left: 40,
                    top: 5
                },
                hAxis: {
                    gridlines: {
                        color: 'transparent',
                        count: 10
                   }
                },
                vAxis: {
                    viewWindowMode: 'maximized'
                },
                width: '100%',
                height: '100%',
                colors: [ "blue", "green", "red", "orange" ]
            };

            var chart = new google.visualization.LineChart(document.getElementById(divName));

            chart.draw(data, options);
            //console.log(dataArray[0][1]);

            var myElements = document.querySelectorAll(".hiddenChart1");
            for (var i = 0; i < myElements.length; i++) {
                myElements[i].style.opacity = 1;
            }
        }

        // function _createAnim_DEPRECATED() {
        //     var eXmax = 1;
        //     var hXmax = 1;
        //     var endRange = $scope.totalLength;
        //     var canvasElement = "testcanvas";
        //     var crystal1 = $scope.crystal;
        //     var fields = $scope.field;
        //     var epsilon = $scope.epsilon2D;
        //     var mu = $scope.mu2D;
        //     var interfaces = $scope.crystal.Struct.materialInterfaces();
        //     var elem = document.getElementById(canvasElement);
        //     var jelem = $("#" + canvasElement);
        //     console.log("canvas element jquer:" + jelem.parent().css("background-color"));
        //     var rgbColor = jelem.parent().css("background-color");
        //     var WIDTH = elem.offsetWidth;
        //     var HEIGHT = elem.offsetHeight;
        //     var w1 = WIDTH;
        //     var h1 = HEIGHT;
        //     var three = THREE.Bootstrap({
        //         plugins: ['core', 'controls'],
        //         controls: {
        //             klass: THREE.OrbitControls
        //         },
        //         size: {
        //             width: w1,
        //             height: h1,
        //             // left: "400",
        //         },
        //     });
        //     // expect(three.renderer.domElement.style.visibility).toBe('hidden');

        //     var renderer = three.renderer;
        //     var scene = three.scene;
        //     var camera = three.camera;
        //     // renderer.domElement.style.visibility = 'hidden';

        //     console.log("styling canvas elem");
        //     // Insert into document

        //     elem.appendChild(renderer.domElement);

        //     // MathBox $scope.context
        //     $scope.context = new MathBox.Context(renderer, scene, camera).init();
        //     var mathbox = $scope.context.api;
        //     // mathbox2 = mathBox({
        //     //   plugins: ['core', 'controls', 'cursor'],
        //     //   controls: {
        //     //     klass: THREE.OrbitControls
        //     //   },    });
        //     // mathbox.install('controls');
        //     // Set size
        //     renderer.setSize(WIDTH, HEIGHT);
        //     $scope.context.resize({ viewWidth: WIDTH, viewHeight: HEIGHT });

        //     // Place camera and set background
        //     camera.position.set(-1, 2, 4);
        //     renderer.setClearColor(new THREE.Color(rgbColor), 1.0);
        //     // renderer.domElement.style.visibility = 'hidden';

        //     console.log("made it here too");
        //     var view = mathbox
        //         .set({
        //             focus: 3,
        //         })
        //         .cartesian({
        //             range: [
        //                 [0, endRange],
        //                 [-eXmax, eXmax],
        //                 [-hXmax, hXmax]
        //             ],
        //             scale: [2, 2, 2],
        //         });

        //     view.axis({
        //         detail: 30,
        //     });

        //     view.axis({
        //         axis: 2,
        //     });

        //     view.scale({
        //         divide: 10,
        //     })
        //     view.ticks({
        //         classes: ['foo', 'bar'],
        //         width: 2
        //     });

        //     // view.grid({
        //     //   divideX: 30,
        //     //   width: 1,
        //     //   opacity: 0.5,
        //     //   zBias: -5,
        //     // });
        //     view.scale({ //adds "X-Axis" to the graph
        //         divide: 1,
        //         origin: [0, -2, 0, 0],
        //         axis: "x",
        //     }).text({
        //         live: false,
        //         data: ["Electric Field"]
        //     }).label({
        //         color: 0x0074D9,
        //         offset: [0, 0]
        //     })
            
        //     view.scale({ //adds "Y-Axis" to the graph
        //         divide: 1,
        //         origin: [0, 0, 2, 0],
        //         axis: "y",
        //     }).text({
        //         live: false,
        //         data: ["Magnetic Field"],
        //         zIndex: 1
        //     }).label({
        //         color: 0xFF4136,
        //         offset: [0, 0]
        //     })

        //     view.grid({
        //             axes: "xy",
        //             divideX: endRange,
        //             divideY: 10
        //         }) //makes two axes in space
        //         .grid({
        //             axes: "xz",
        //             divideX: endRange,
        //             divideY: 10,
        //         })

        //         // three.canvas.style.visibility = "hidden";

        //     var colorCoords = []; //possibly remove, replace with just applying interfaces to arrays
        //     var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
        //     for (var i = 0; i < interfaces.length - 1; i++) {
        //         var array1 = [ //only x changes on all shapes coordinates

        //             //seperate 1st and last arrays. also, make each array part of an array[][]
        //             //
        //             0 + interfaces[0 + i], -2, 0, 0 + interfaces[1 + i], -2, 0, 0 + interfaces[1 + i], 2, 0, 0 + interfaces[0 + i], 2, 0, //first column
        //             //0+interfaces[0+i], 0, 1, 0+interfaces[0+i], 0, -1, 0+interfaces[1+i], 0, -1, 0+interfaces[1+i], 0, 1,  //first row
        //             // 0+interfaces[interfaces.length-2], -1, 0, 0+interfaces[interfaces.length-1], -1, 0, 0+interfaces[interfaces.length-1], 1, 0, 0+interfaces[interfaces.length-2], 1, 0, //last column
        //             // 0+interfaces[interfaces.length-2], 0, 1, 0+interfaces[interfaces.length-2], 0, -1, 0+interfaces[interfaces.length-1], 0, -1, 0+interfaces[interfaces.length-1], 0, 1,  //last row


        //         ];
        //         var array2 = [0 + interfaces[0 + i], 0, 2, 0 + interfaces[0 + i], 0, -2, 0 + interfaces[1 + i], 0, -2, 0 + interfaces[1 + i], 0, 2, ]

        //         // colorCoords.push(array1); //remove push, place the voxel code here
        //         // colorCoords.push(array2);?
        //         console.log("test array1:" + array1)
        //         view.voxel({
        //             data: array1,
        //             width: 4,
        //             height: 2,
        //             depth: 1,
        //             items: 4,
        //             channels: 3,
        //         });

        //         view
        //             .transform({
        //                 pass: 'eye',
        //                 position: [0, 0, 0],
        //                 scale: [1, 1, 1],
        //             })
        //             .face({
        //                 color: cols[i % 5],
        //                 width: 3,
        //                 line: false,
        //                 shaded: true,
        //                 opacity: .3,
        //                 zOrder: 9,
        //             })
        //             .face({
        //                 color: cols[i % 5],
        //                 width: 3,
        //                 fill: false,
        //                 opacity: .3,
        //                 zOrder: 9,
        //             })
        //             .array({
        //                 data: [0, 0, 0],
        //                 channels: 3,
        //             })
        //         view.voxel({
        //             data: array2,
        //             width: 4,
        //             height: 2,
        //             depth: 1,
        //             items: 4,
        //             channels: 3,
        //         });

        //         view
        //             .transform({
        //                 pass: 'eye',
        //                 position: [0, 0, 0],
        //                 scale: [1, 1, 1],
        //             })
        //             .face({
        //                 color: cols[i % 5],
        //                 width: 3,
        //                 line: false,
        //                 shaded: true,
        //                 opacity: .3,
        //                 zOrder: 9,
        //             })
        //             .face({
        //                 color: cols[i % 5],
        //                 width: 3,
        //                 fill: false,
        //                 line: true,
        //                 opacity: .3,
        //                 zOrder: 9,
        //             })
        //             .array({
        //                 data: [0, 0, 0],
        //                 channels: 3,
        //             })
        //     }

        //     var ourFunc;
        //     var runsetup = true;
        //     var EigVals = new Array();
        //     var EigVecs = new Array();
        //     var ExrArray = new Array();
        //     var ExthetaArray = new Array();
        //     var EyrArray = new Array();
        //     var EythetaArray = new Array();

        //     // var HxrArray = new Array();
        //     // var HxthetaArray = new Array();
        //     // var HyrArray = new Array();
        //     // var HythetaArray = new Array();

        //     view.interval({
        //         id: 'ElectricFieldPlot',
        //         width: endRange * 10, //fields.Ex.length,
        //         expr: function(emit, z, i, t) {
        //             if (runsetup) {
        //                 var U = crystal1.scattering($scope.o, $scope.k1, $scope.k2, $scope.incoming);
        //                 //Electric field
        //                 var cxE = [$scope.incoming[0], $scope.incoming[1], U.x[0], U.x[1]];
        //                 var cyE = [0, 0, U.y[0], U.y[1]];
        //                 var cE = new numeric.T(cxE, cyE);
                        
        //                 //Magnetic field
        //                 // var cxH = [j1, j2, U.x[2], U.x[3]];
        //                 // var cyH = [0, 0, U.y[2], U.y[3]];
        //                 // var cH = new numeric.T(cxH, cyH);

        //                 for (var q = 0; q < epsilon.length; q++) {
        //                     EigVals[q] = emScattering.eigenvaluesIsotropic(epsilon[q], mu[q], $scope.k1, $scope.k2);
        //                     EigVecs[q] = emScattering.eigenvectorsIsotropic(epsilon[q], mu[q], $scope.k1, $scope.k2);
                            
        //                     // console.log("EigVecs[q]: "+EigVecs[q].y[0][0]);//EigVecs[layernum].y[wavetype][vectornum]
        //                     // .x or .y for real and imaginary parts
        //                     // wavetype is 0 for ex, 1 for ey, 2 for hx, 3 for hy
        //                     // Electric field
        //                     var ExR = new Array();
        //                     var ExTheta = new Array();
        //                     var EyR = new Array();
        //                     var EyTheta = new Array();
        //                     // Magnetic field
        //                     // var HxR = new Array();
        //                     // var HxTheta = new Array();
        //                     // var HyR = new Array();
        //                     // var HyTheta = new Array();
        //                     for (var j = 0; j < 4; j++) {
        //                         //Electric field
        //                         var ExA = cE.x[j] * EigVecs[q].x[0][j];
        //                         var EyA = cE.x[j] * EigVecs[q].x[1][j];
        //                         var ExB = cE.y[j] * EigVecs[q].y[0][j];
        //                         var EyB = cE.y[j] * EigVecs[q].y[1][j];
        //                         ExR[j] = Math.sqrt(Math.pow(ExA, 2) + Math.pow(ExB, 2));
        //                         EyR[j] = Math.sqrt(Math.pow(EyA, 2) + Math.pow(EyB, 2));
        //                         ExTheta[j] = Math.atan2(ExB, ExA);
        //                         EyTheta[j] = Math.atan2(EyB, EyA);
        //                         //Magnetic field
        //                         // var HxA = cH.x[j]*EigVecs[q].x[2][j];
        //                         // var HyA = cH.x[j]*EigVecs[q].x[3][j];
        //                         // var HxB = cH.y[j]*EigVecs[q].y[2][j];
        //                         // var HyB = cH.y[j]*EigVecs[q].y[3][j];
        //                         // HxR[j] = Math.sqrt(Math.pow(HxA,2)+Math.pow(HxB,2));
        //                         // HyR[j] = Math.sqrt(Math.pow(HyA,2)+Math.pow(HyB,2));
        //                         // HxTheta[j] = Math.atan2(HxB,HxA);
        //                         // HyTheta[j] = Math.atan2(HyB,HyA);
        //                     }
        //                     //electric field
        //                     ExrArray[q] = ExR;
        //                     EyrArray[q] = EyR;
        //                     ExthetaArray[q] = ExTheta;
        //                     EythetaArray[q] = EyTheta;
        //                     //magnetic field
        //                     // HxrArray[q] = HxR;
        //                     // HyrArray[q] = HyR;
        //                     // HxthetaArray[q] = HxTheta;
        //                     // HythetaArray[q] = HyTheta;

        //                     if (q + 1 < epsilon.length) {

        //                 // console.log("here tooman22233124");
        //                     try {
        //                         cE = $scope.crystal.crystal.transferMatrices[q].dot(cE); // generates an exception
        //                     }
        //                     catch (e) {
        //                     // statements to handle any exceptions
        //                     console.log("errored out");
        //                     logMyErrors(e); // pass exception object to error handler
        //                     }
                                
        //                         // console.log("here tooman222332");
        //                         // cH = crystal1.crystal.transferMatrices[q].dot(cH);
        //                         // console.log("c: "+c.x);
        //                     }

        //                 };
        //             };
                    
        //             // runsetup = false;
        //             var layerNumber;
        //             var lowerL = 0;
        //             var upperL = $scope.length2D[0];
        //             for (var ii = 0; ii < $scope.length2D.length; ii++) {
        //                 if (lowerL <= z && z <= upperL) {
        //                     layerNumber = ii;
        //                     break
        //                 } else {
        //                     lowerL = upperL;
        //                     upperL = upperL + $scope.length2D[ii + 1];
                        
        //                 }
        //             }
                    
        //             var Ex = ExrArray[layerNumber][0] * Math.cos(ExthetaArray[layerNumber][0] + EigVals[layerNumber][0].y * z - $scope.o * t) +
        //                 ExrArray[layerNumber][1] * Math.cos(ExthetaArray[layerNumber][1] + EigVals[layerNumber][1].y * z - $scope.o * t) +
        //                 ExrArray[layerNumber][2] * Math.cos(ExthetaArray[layerNumber][2] + EigVals[layerNumber][2].y * z - $scope.o * t) +
        //                 ExrArray[layerNumber][3] * Math.cos(ExthetaArray[layerNumber][3] + EigVals[layerNumber][3].y * z - $scope.o * t);

        //            var  Ey = EyrArray[layerNumber][0] * Math.cos(EythetaArray[layerNumber][0] + EigVals[layerNumber][0].y * z - $scope.o * t) +
        //                 EyrArray[layerNumber][1] * Math.cos(EythetaArray[layerNumber][1] + EigVals[layerNumber][1].y * z - $scope.o * t) +
        //                 EyrArray[layerNumber][2] * Math.cos(EythetaArray[layerNumber][2] + EigVals[layerNumber][2].y * z - $scope.o * t) +
        //                 EyrArray[layerNumber][3] * Math.cos(EythetaArray[layerNumber][3] + EigVals[layerNumber][3].y * z - $scope.o * t);

        //             // Hx = ExrArray[layerNumber][0]*Math.cos(HxthetaArray[layerNumber][0] + EigVals[layerNumber][0].y*z - o*t) +
        //             //     ExrArray[layerNumber][1]*Math.cos(HxthetaArray[layerNumber][1] + EigVals[layerNumber][1].y*z - o*t) +
        //             //     ExrArray[layerNumber][2]*Math.cos(HxthetaArray[layerNumber][2] + EigVals[layerNumber][2].y*z - o*t) +
        //             //     ExrArray[layerNumber][3]*Math.cos(HxthetaArray[layerNumber][3] + EigVals[layerNumber][3].y*z - o*t);

        //             // Hy = EyrArray[layerNumber][0]*Math.cos(HythetaArray[layerNumber][0] + EigVals[layerNumber][0].y*z - o*t) +
        //             //     EyrArray[layerNumber][1]*Math.cos(HythetaArray[layerNumber][1] + EigVals[layerNumber][1].y*z - o*t) +
        //             //     EyrArray[layerNumber][2]*Math.cos(HythetaArray[layerNumber][2] + EigVals[layerNumber][2].y*z - o*t) +
        //             //     EyrArray[layerNumber][3]*Math.cos(HythetaArray[layerNumber][3] + EigVals[layerNumber][3].y*z - o*t);
        //             // x = Math.cos(-5*z);
        //             // y = Math.cos(2*x);

        //             if (Ex > eXmax) {
        //                 // console.log("Ex" + Ex);
        //                 //eXmax++;
        //                // addAnim("testcanvas", endRange);
        //             }
        //             emit(z, Ex, Ey);
        //             emit(z, 0, 0);
        //         },
        //         items: 2,
        //         channels: 3,
        //         live: true,
        //     });

        //     // view.line({
        //     //   points: '#ElectricFieldPlot',
        //     //   color: 0x0074D9,
        //     //   width: 2,
        //     // });
        //     view.vector({
        //         points: '#ElectricFieldPlot',
        //         color: 0x0074D9,
        //         width: 1,
        //         start: true,
        //     });

        //     var HxrArray = new Array();
        //     var HxthetaArray = new Array();
        //     var HyrArray = new Array();
        //     var HythetaArray = new Array();

        //     view.interval({
        //         id: 'MagneticFieldPlot',
        //         width: endRange * 10, //fields.Ex.length,
        //         expr: function(emit, z, i, t) {
        //             if (runsetup) {
        //                 var U = crystal1.scattering($scope.o, $scope.k1, $scope.k2, $scope.incoming);
        //                 //Magnetic field
        //                 var cxH = [$scope.incoming[0], $scope.incoming[1], U.x[2], U.x[3]];
        //                 var cyH = [0, 0, U.y[2], U.y[3]];
        //                 var cH = new numeric.T(cxH, cyH);

        //                 for (i = 0; i < epsilon.length; i++) {
        //                     EigVals[i] = emScattering.eigenvaluesIsotropic(epsilon[i], mu[i], $scope.k1, $scope.k2);
        //                     EigVecs[i] = emScattering.eigenvectorsIsotropic(epsilon[i], mu[i], $scope.k1, $scope.k2);
        //                     // console.log("EigVecs[i]: "+EigVecs[i].y[0][0]);//EigVecs[layernum].y[wavetype][vectornum]
        //                     // .x or .y for real and imaginary parts
        //                     // wavetype is 0 for ex, 1 for ey, 2 for hx, 3 for hy
        //                     var HxR = new Array();
        //                     var HxTheta = new Array();
        //                     var HyR = new Array();
        //                     var HyTheta = new Array();
        //                     for (var j = 0; j < 4; j++) {
        //                         var HxA = cH.x[j] * EigVecs[i].x[2][j];
        //                         var HyA = cH.x[j] * EigVecs[i].x[3][j];
        //                         var HxB = cH.y[j] * EigVecs[i].y[2][j];
        //                         var HyB = cH.y[j] * EigVecs[i].y[3][j];
        //                         HxR[j] = Math.sqrt(Math.pow(HxA, 2) + Math.pow(HxB, 2));
        //                         HyR[j] = Math.sqrt(Math.pow(HyA, 2) + Math.pow(HyB, 2));
        //                         HxTheta[j] = Math.atan2(HxB, HxA);
        //                         HyTheta[j] = Math.atan2(HyB, HyA);
        //                     }
        //                     HxrArray[i] = HxR;
        //                     HyrArray[i] = HyR;
        //                     HxthetaArray[i] = HxTheta;
        //                     HythetaArray[i] = HyTheta;

        //                     if (i + 1 < epsilon.length) {
        //                         cH = crystal1.crystal.transferMatrices[i].dot(cH);
        //                     }
        //                 };
        //                 runsetup = false;
        //             };
        //             var layerNumber;
        //             var lowerL = 0;
        //             var upperL = $scope.length2D[0];
        //             for (var ii = 0; ii < $scope.length2D.length; ii++) {
        //                 if (lowerL <= z && z <= upperL) {
        //                     layerNumber = ii;
        //                     break
        //                 } else {
        //                     lowerL = upperL;
        //                     upperL = upperL + $scope.length2D[ii + 1];
        //                 }
        //             }

        //            var Hx = HxrArray[layerNumber][0] * Math.cos(HxthetaArray[layerNumber][0] + EigVals[layerNumber][0].y * z - $scope.o * t) +
        //                 HxrArray[layerNumber][1] * Math.cos(HxthetaArray[layerNumber][1] + EigVals[layerNumber][1].y * z - $scope.o * t) +
        //                 HxrArray[layerNumber][2] * Math.cos(HxthetaArray[layerNumber][2] + EigVals[layerNumber][2].y * z - $scope.o * t) +
        //                 HxrArray[layerNumber][3] * Math.cos(HxthetaArray[layerNumber][3] + EigVals[layerNumber][3].y * z - $scope.o * t);

        //             var Hy = HyrArray[layerNumber][0] * Math.cos(HythetaArray[layerNumber][0] + EigVals[layerNumber][0].y * z - $scope.o * t) +
        //                 HyrArray[layerNumber][1] * Math.cos(HythetaArray[layerNumber][1] + EigVals[layerNumber][1].y * z - $scope.o * t) +
        //                 HyrArray[layerNumber][2] * Math.cos(HythetaArray[layerNumber][2] + EigVals[layerNumber][2].y * z - $scope.o * t) +
        //                 HyrArray[layerNumber][3] * Math.cos(HythetaArray[layerNumber][3] + EigVals[layerNumber][3].y * z - $scope.o * t);
        //             // x = Math.cos(-5*z);
        //             // y = Math.cos(2*x);
        //             // if(Hy>hXmax){ console.log("Hx"+Hy);
        //             // hXmax++;
        //             // addAnim("testcanvas",endRange);}

        //             emit(z, Hx, Hy);
        //             emit(z, 0, 0);
        //             // emit(z, Hy, Hx);
        //         },
        //         items: 2,
        //         channels: 3,
        //         live: true,
        //     });

        //     // view.line({
        //     //   points: '#MagneticFieldPlot',
        //     //   color: 0xFF4136,
        //     //   width: 2,
        //     // });

        //     view.vector({
        //         points: '#MagneticFieldPlot',
        //         color: 0xFF4136,
        //         width: 1,
        //         start: true,
        //     });
        //     var visible = false;
        //     var madeVisible = false;
        //     var frame = function() {
        //         var parentVisibility = jelem.parent().css("visibility");
        //             //console.log("doings things: ");
        //         if(parentVisibility != 'hidden'){
        //             requestAnimationFrame(frame);
        //             visible = true;
                    
        //         }
        //         else if(parentVisibility == 'hidden'){
        //             visible = false;

        //             requestAnimationFrame(frame);
        //             renderer.domElement.style.visibility = parentVisibility;
        //             return;
        //         }
                
        //         if(!madeVisible) {
        //             renderer.domElement.style.visibility = parentVisibility; //wheretohide
        //             madeVisible = true;
        //         }
                

        //         if(!visible)
        //             console.log("looping while visible is false");
        //         //WIDTH = elem.offsetWidth - 10;
        //         //HEIGHT = elem.offsetHeight - 5;
        //         renderer.setSize(WIDTH, HEIGHT);

        //         //$scope.context.resize({ viewWidth: WIDTH, viewHeight: HEIGHT });
        //         $scope.context.frame();
        //         renderer.render(scene, camera);
        //         // console.log("this is happening");
        //         rgbColor = jelem.parent().css("background-color"); //0xFF851B;
        //         renderer.setClearColor(new THREE.Color(rgbColor), 1.0);

                
        //     };
        //     requestAnimationFrame(frame);
        //     //var off = $(elem).offset();
        //     //off.left += 5;
        //     //$("canvas").offset(off);
        //     // $("canvas").offset()

        // }

        /** Creates the chart for the Dispersion tab. WIP */
        function createDispersionChart() {
            var dispersion = $scope.dispersion;
            var divName = "dispView";
            var width = "100%";
            var height = "100%";
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
        }

        function _createFieldChart_DEPRECATED() {
            var fields = $scope.field;
            var interfaces = $scope.crystal.Struct.materialInterfaces();



            var data = new google.visualization.DataTable();
            data.addColumn('number', 'z');
            data.addColumn('number', document.getElementById("shownVal").value);
            //Iterate through fields values

            for (var i = 0, N = fields.z.length; i < N; i++) {
                if (document.getElementById("shownVal").value == $scope.EX) {
                    data.addRows([
                        [fields.z[i], fields.Ex[i]]
                    ]);
                }
                if (document.getElementById("shownVal").value == $scope.EY) {
                    data.addRows([
                        [fields.z[i], fields.Ey[i]]
                    ]);
                }
                if (document.getElementById("shownVal").value == $scope.HX) {
                    data.addRows([
                        [fields.z[i], fields.Hx[i]]
                    ]);
                }
                if (document.getElementById("shownVal").value == $scope.HY) {
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
                width: "100%",
                height: "100%",
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
                var htmlClass = "addHere";
                if (htmlClass == "addHere2")
                    $(".overlay2").remove();
                else if (htmlClass == "addHere")
                    $(".overlay").remove();
                else if (htmlClass == "addHere3")
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

                if (htmlClass == "addHere") {
                    var overlays = document.querySelectorAll('.overlay');
                    for (var j = 0; j < overlays.length; j++) {
                        overlays[j].style.position = 'absolute';
                        overlays[j].style.opacity = '.5';
                        overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 86 + "px";
                        overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 16 + "px";

                    }
                } else if (htmlClass == "addHere2") {
                    var overlays = document.querySelectorAll('.overlay2');
                    for (var j = 0; j < overlays.length; j++) {
                        overlays[j].style.position = 'absolute';
                        overlays[j].style.opacity = '.5';
                        overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 15 + "px";
                        overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 30 + "px";
                    }
                }



            }

            var chart = new google.visualization.LineChart(document.getElementById('structView'));
            google.visualization.events.addListener(chart, 'ready', printInterfaces.bind(chart, data));
            //var chart = new google.charts.Line(document.getElementById('linechart_material'));
            chart.draw(data, options);

            var myElements = document.querySelectorAll(".hiddenChart");
            for (var i = 0; i < myElements.length; i++) {
                myElements[i].style.opacity = 1;
            }
        }

        /** Creates the chart for the Structures tab. */
        function createStructureChart() {
            var interfaces = $scope.crystal.Struct.materialInterfaces();               //Gets the interfaces
            var interfaceLength = interfaces[interfaces.length - 1];            //Gets highest interface

            //Creates a data table
            var data = google.visualization.arrayToDataTable([
                [{
                    f: 'Date',
                    type: 'number' // wont work whithout this
                }, {
                    f: 'Line',
                    type: 'number' // wont work whithout this
                }],
            ]);
            var jelem = $("#structureView");        //Gets structureView element using jQuery
            console.log("canvas element jquer:" + jelem.parent().css("background-color"));
            var rgbColor = jelem.parent().css("background-color");
            
            //Chart Options
            var options = {

                chart: {
                    title: 'Dispersion Relationship'
                },
                chartArea: {
                    left: 40,
                    top: 40
                },
                width: '100%',
                height: '100%',
                hAxis: {
                    gridlines: { 
                        count: (interfaceLength / 2),
                        color: 'transparent'
                    },
                    viewWindow: {
                        min: interfaces[0],
                        max: interfaces[interfaces.length - 1] //modify max by adding layer lengths, then apply to buildstructure button, then set up value placement of layers, then custom overlays
                    },
                    ticks: interfaces
                },
                vAxis: {
                    //title: 'z',
                    textPosition: 'none',
                    gridlines: {
                         count: 10,
                         color: 'transparent'  
                    },
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

            function printInterfaces(dataTable) { //prints the colored squares on top  of the charts
                var cli = this.getChartLayoutInterface();
                var chartArea = cli.getChartAreaBoundingBox();
                var cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple'];
                var oddColors = '#FFFF00';
                var evenColors = '#888888';
                var w = cli.getXLocation(interfaces[1]) - cli.getXLocation(interfaces[0]);
                var y = cli.getChartAreaBoundingBox().height;
                // console.log(interfaces);
                // console.log("w:" + w);
                // console.log("y:" + y);
                // console.log("gety:" + Math.floor(cli.getYLocation(1)));
                // console.log("bounding:" + cli.getChartAreaBoundingBox().top);
                var yBound = cli.getChartAreaBoundingBox().top;
                var htmlClass = "addHere3";
                if (htmlClass == "addHere3")
                    $(".overlay3").remove();
                for (var i = 0; i < interfaces.length - 1; i++) {
                    var w = cli.getXLocation(interfaces[i + 1]) - cli.getXLocation(interfaces[i]);

                    addStruct(htmlClass, i, (i % 2 === 0 ? evenColors : oddColors), w, y);
                };
                if (htmlClass == "addHere3") {
                    var overlays = document.querySelectorAll('.overlay3');
                    for (var j = 0; j < overlays.length; j++) {
                        overlays[j].style.position = 'absolute';
                        overlays[j].style.opacity = '.5';
                        overlays[j].style.top = Math.floor(cli.getChartAreaBoundingBox().top) + 64 + "px";
                        overlays[j].style.left = Math.floor(cli.getXLocation(interfaces[j])) + 16 + "px";
                    }
                }
            }

            var chart = new google.visualization.LineChart(document.getElementById('structureView'));       //Create a line chart in the structureView element
            google.visualization.events.addListener(chart, 'ready', printInterfaces.bind(chart, data)); //all from Google's tutorial https://developers.google.com/chart/interactive/docs/overlays

            chart.draw(data, options);                      //Draws the chart
        }


    
});
