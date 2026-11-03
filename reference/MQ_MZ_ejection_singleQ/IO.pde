/**
 ***********************************************************************************************************************
 *
 *                                                 PROGRAM OPTIONS
 *
 ***********************************************************************************************************************
 */  

// Function that handles key press events and updates flags or actions accordingly
void keyPressed(){
   
  // Check if "0" key is pressed
  if(keyCode==48 ) { // "0" key pressed
    flag = 0;  // Reset flag
  }
   
  // Check if "7" key is pressed
  if(keyCode==55 ) { // "7" key pressed
    flag = 7;  // Set flag to 7
  }
   
  // Inject a singly-charged ion when "i" key is pressed
  if(keyCode==73) { // "i" key pressed
    // Create a new ion object with specific properties (m, z, Ion energy in eV)
    Ion ion = new Ion(500, 1, ION_ENERGY); // m, z, Vin[eV], Vrf[V], RF[kHz], phase[degrees], P[mTorr]
    CalculateIonTrajSingleIon(ion);  // Calculate the ion trajectory for this ion
    flag = 6;  // Set flag to 6 (indicating ion injection)
  }
  
  // Inject a multiply-charged ion when "m" key is pressed
  if(keyCode==77) { // "m" key pressed
    // Create a new ion object with different mass, charge state, and energy
    Ion ion = new Ion(2500, 5, ION_ENERGY); // m, z, Vin[eV], Vrf[V], RF[kHz], phase[degrees], P[mTorr]
    CalculateIonTrajSingleIon(ion);  // Calculate the ion trajectory for this ion
    flag = 6;  // Set flag to 6
  }
  
  // Increase resonance frequency when "f" key is pressed
  if(keyCode==70) { // "f" key pressed
    RESONANCE_FR = RESONANCE_FR+1;  // Increment resonance frequency
  }
  
   
  // ====================================== FLAG OPTIONS (Change Display Options) ===================================
  // Set different flag values based on key presses (for controlling display options)
  if(keyCode==KeyEvent.VK_1 ) {
    flag = 1;  // Set flag to 1
  }
    
  if(keyCode==KeyEvent.VK_2 ) {
    flag = 2;  // Set flag to 2
  }
    
  if(keyCode==KeyEvent.VK_3 ) {
    flag = 3;  // Set flag to 3
  }
  
  if(keyCode==KeyEvent.VK_4 ) {
    flag = 4;  // Set flag to 4
  }
  
  if(keyCode==KeyEvent.VK_5 ) {
    flag = 5;  // Set flag to 5
  }

  if(keyCode==KeyEvent.VK_6 ) {
    flag = 6;  // Set flag to 6
  }
  
  // Show options related to flag 1 when F1 or other specific keys are pressed
  if(keyCode==97 || keyCode==112) { // "F1" key or "p" key pressed
    showOptions_1();  // Display options for flag 1
  }
  
  // Show options related to flag 2 when F2 is pressed
  if(keyCode==98 ) { // "F2" key pressed
    showOptions_2();  // Display options for flag 2
  }
  
  // Optionally, a place for further key handling, such as "F3", can be added here
  // if(keyCode==99 ) { // "F3" key pressed
  //   showOptions_3();  // Display options for flag 3 (currently commented out)
  // }

  
  // ============================ MOVE MULTI-QUAD =================================
  // These key presses rotate the MultiQ electrodes in various directions by a fixed amount (PI/16)

  if(keyCode==KeyEvent.VK_UP ){
    a1 += PI/16;  // Rotate along axis 1 (upward)
  }
  
  if(keyCode==KeyEvent.VK_DOWN){
    a1 += -PI/16;  // Rotate along axis 1 (downward)
  }
  
  if(keyCode==KeyEvent.VK_LEFT){
    a2 += PI/16;  // Rotate along axis 2 (leftward)
  }
  
  if(keyCode==KeyEvent.VK_RIGHT){
    a2 += -PI/16;  // Rotate along axis 2 (rightward)
  }
  
  if(keyCode==KeyEvent.VK_PAGE_UP){
    a3 += PI/16;  // Rotate along axis 3 (upward)
  }
  
  if(keyCode==KeyEvent.VK_PAGE_DOWN){
    a3 += -PI/16;  // Rotate along axis 3 (downward)
  }
  
  // Navigate through z-axis positions when "," or "." are pressed (adjust zcs variable)
  if(keyCode==44) { // "," key pressed
    if (zcs>=2 && zcs<=Nz-1) { 
      zcs -=1;  // Decrease zcs (move up in z-axis)
      redraw();  // Redraw the view with updated zcs
    }
  }
  
  if(keyCode==46) { // "." key pressed
    if (zcs>=1 && zcs<=Nz-2) { 
      zcs +=1;  // Increase zcs (move down in z-axis)
      redraw();  // Redraw the view with updated zcs
    }
  }

} // End of keyPressed function

/**
 *************************************************************************************************************************
 *
 *                                               1. Geometry & Potentials
 *
 *************************************************************************************************************************
 */

// This method provides a graphical interface for the user to select actions related to geometry and potentials.
public void showOptions_1() {
  
  SwingUtilities.invokeLater(new Runnable() { // run as a separate thread      
    public void run() { 

      // Define the actions the user can select from.
      String[] plays = new String[] { 
                                      "Create Geometry",     // Create electrode geometry
                                      "Calculate Potential",  // Calculate potential distribution
                                      "Save Potential",      // Save calculated potential to file
                                      "Read Potential"       // Read potential from a file                                    
                                    };
      
      // Show a dialog where the user selects an action
      String input = (String) JOptionPane.showInputDialog(
        new JFrame(),
        "Please select the action",
        "Geometry & Potentials", JOptionPane.INFORMATION_MESSAGE,
        new ImageIcon("java2sLogo.GIF"), plays, "Demiurge Actions");
        
        // Perform the corresponding action based on user input
        if(input == "Create Geometry") CreateElectrodeGeometry();
        if(input == "Calculate Potential") CalculatePotentialDistribution();  
        if(input == "Save Potential") SavePotential();
        if(input == "Read Potential") ReadPotential();
    }
  });
}

//-------------------------------------- CREATE ELECTRODE GEOMETRY ----------------------------------------

// This method calculates the electrode geometry and updates the flag to display the appropriate result.
public void CreateElectrodeGeometry() {
  println("Will calculate RF newgeometry");
  Electrodes = calculateElectrodeGeometry();  // Calculate the geometry and store in Electrodes
  flag = 1;  // Set the flag to 1 to update the display
}

//------------------------------------ CALCULATE POTENTIAL DISTRIBUTION -----------------------------------

// This method calculates the potential distribution for the given electrode geometry.
public void CalculatePotentialDistribution() {
  println("Calculate potential for a given electrode geometry");
  POT = calculatePotential(0.001, POT );  // Call calculatePotential with a precision of 0.001
}

//---------------------------------------------- I/O  --------------------------------------------------------

// This method prompts the user to select a file for saving the potential data.
public void SavePotential() {
  selectOutput("Select a file to write to:", "savePotentialIntoFile");
}

// SAVE POTENTIAL

// This method saves the calculated potential distribution into a selected file.
void savePotentialIntoFile(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
  } else {
    println("User selected " + selection.getAbsolutePath());
    
    // Create a writer to write to the selected file
    PrintWriter output = createWriter(selection.getAbsolutePath());
    
    int e = 0;
    // Loop through the potential grid and write each point to the file
    for (int s = 0; s < Nz; s++) {
      if (s == 0) {
        // Write the first row in the file with device dimensions and scale
        output.println("X: " + xDevice + " [mm], Y: " + yDevice + " [mm], Z: " + zDevice + " [mm], SCALE: " + scl);
      } else {
        // Loop through all x, y, z coordinates and write potential data
        for (int i = 0; i < Nx; i++) {
          for (int j = 0; j < Ny; j++) {
            if(POT[i][j][s].electrode) e = 1; else e = 0;
            output.println(i + "," + j + "," + s + "," + POT[i][j][s].potential + "," + e);
          }
        }
      }
    }
        
    output.flush();  // Write the remaining data to the file
    output.close();  // Close the file after writing
 
  }
}

//------------------------------------ SELECTING FILE TO OPEN -------------------------------------------

// This method prompts the user to select a file to read potential data from.
void ReadPotential() {
  selectInput("Select a file to process:", "openPotentialFile");  
}

// This method opens and processes the selected potential file.
void openPotentialFile(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
  } else {
    Electrodes = openPotential(selection.getAbsolutePath(), POT);  // Open and parse the potential file
  }  
}

// This method opens the potential file and reads the data into an array of Potential objects.
ArrayList<Potential> openPotential(String fileName, Potential[][][] pot) {

  startTimer();  // Start a timer to track the file processing time
  
  File chosenFile =  new File(fileName);
  
  try {
     Scanner fileScanner = new Scanner(chosenFile);  // Initialize a scanner to read the file
     String firstLine = fileScanner.nextLine();  // Read the first line (device dimensions and scale)
     String[] fileData = split(firstLine, ",");
     
     // Parse the first line to extract dimensions and scaling factor
     float dimx = Float.parseFloat(split(fileData[0], " ")[1]);
     float dimy = Float.parseFloat(split(fileData[1].trim(), " ")[1]);
     float dimz = Float.parseFloat(split(fileData[2].trim(), " ")[1]);
     float sc   = Float.parseFloat(split(fileData[3].trim(), " ")[1]); 
     int dim = floor((dimx * sc) * (dimy * sc) * (dimz * sc));  // Calculate the total number of points 
         
     int k = 0;
     int progress = 0;
     
     ArrayList<Potential> potout = new ArrayList<Potential>();  // List to store potential data points
  
     // Loop through the rest of the file and parse each line
     while (fileScanner.hasNext()) {   
       int job = floor(float(k)/dim * 100);  // Calculate progress percentage
       if (job == progress) {
         print("|");
         if (progress % 10 == 0) { print(progress); }
         if (progress == 99) { println("100%"); println(""); }
         progress++;
       }  
    
       String line = fileScanner.nextLine();  // Read the next line
       String[] linedata = split(line, ",");
       
       // Parse the data from the current line
       float nx = Float.parseFloat(linedata[0]);
       float ny = Float.parseFloat(linedata[1]);
       float nz = Float.parseFloat(linedata[2]);
       float potential = Float.parseFloat(linedata[3]);
       boolean el = (Float.parseFloat(linedata[4]) > 0);  // Check if this is an electrode potential
       
       // Store the parsed data into the potential grid
       pot[int(nx)][int(ny)][int(nz)] = new Potential(nx * sc, ny * sc, nz * sc, potential, el);
      
       // If it's an electrode, add it to the output list
       if (el) {
         potout.add(new Potential(nx, ny, nz, potential, el)); 
       } 
      
       k++;  
     }
     
     // Print the size of the output list
     println("The size of the array is: " + potout.size());
     Electrodes = potout;  // Update the Electrodes list with the parsed data
     fileScanner.close();  // Close the scanner
     stopTimer();  // Stop the timer
     
     flag = 1;  // Update the flag to trigger the display update
     return potout;  // Return the list of potential data points

  } catch (IOException e) {
    System.err.println("Caught IOException: " + e.getMessage());  // Handle any file reading errors
    return null;  // Return null in case of error
  }
} // close openPotentialFile



/**
 *************************************************************************************************************************
 *
 *                                                2. Ion Motion
 *
 *************************************************************************************************************************
 */
public void showOptions_2() {
    // Method to show the options for ion motion-related actions
    SwingUtilities.invokeLater(new Runnable() { // run as a separate thread      
        public void run() { 
            String[] plays = new String[] { 
                "Load Precalculated Potentials",
                "Calculate ion trajectory for 1+",
                "Calculate ion trajectory for 5+",
                "Find resonance frequency",
                "Spectrum (find transmission of testM at given RESONANCE_FR)"
            };
            // Show a dialog to select the desired action
            String input = (String) JOptionPane.showInputDialog(
                new JFrame(),
                "Display Options",
                "Ion Motion", JOptionPane.INFORMATION_MESSAGE,
                new ImageIcon("java2sLogo.GIF"), plays, "DisplayOptions");
            
            // Perform the selected action
            if(input == "Load Precalculated Potentials") {
                println( "...loading RF potential from file: " + DEFAULT_RFPOT_FILE);
                openPotential(DEFAULT_RFPOT_FILE, RFPOT);
                println( "...loading DC potential from file: " + DEFAULT_DCPOT_FILE);
                openPotential(DEFAULT_DCPOT_FILE, DCPOT);
                println( "...loading AC potential from file: " + DEFAULT_ACPOT_FILE);
                openPotential(DEFAULT_ACPOT_FILE, ACPOT);
            }   
            if(input == "Calculate ion trajectory for 1+") CalculateIonTrajectoryForOnePlus();
            if(input == "Calculate ion trajectory for 5+") CalculateIonTrajectoryForFivePlus();
            if(input == "Find resonance frequency") findResonanceFrequency();
            if(input == "Spectrum (find transmission of testM at given RESONANCE_FR)") findTransmission();
        }
    });
}

/**
 ********************************************************************************************
 *
 *              Opens the existing file defining the geometry of the potential 
 *
 ********************************************************************************************
 */


/** -------------------------- CALCULATE ION TRAJECTORIES FOR 1 ION  ------------------------------------*/
public void CalculateIonTrajSingleIon(Ion ion) {
    println("Calculate ion trajectory");
    // Calls a method to calculate the trajectory of a single ion
    ionTrajectory = calculateSingleIonTrajectory(0.0001, ion); // spacialResolution 0.1 mm
    flag = 6;
}

/** -------------------------- CALCULATE ION TRAJECTORIES FOR 1 ION (Single) ------------------------------------*/
public void CalculateIonTrajectoryForOnePlus() {
    println("Calculate ion trajectory for 1+ ion");
    // Create an Ion object for a +1 charge state with a specific energy
    Ion ion = new Ion(500, 1, ION_ENERGY);
    ionTrajectory = calculateSingleIonTrajectory(0.0001, ion); // spacialResolution 0.1 mm
    flag = 6;
}

/** -------------------------- CALCULATE ION TRAJECTORIES FOR 5 IONS ------------------------------------*/
public void CalculateIonTrajectoryForFivePlus() {
    println("Calculate ion trajectory for 5+ ion");
    // Create an Ion object for a +5 charge state with a specific energy
    Ion ion = new Ion(2500, 5, ION_ENERGY);
    ionTrajectory = calculateSingleIonTrajectory(0.0001, ion); // spacialResolution 0.1 mm
    flag = 6;
}

/** -------------------------- CALCULATE ION TRAJECTORIES FOR MANY IONS --------------------------------*/
int[] CalculateIonTrajForManyIons(int Noi, int z, float m) {
    ArrayList<Trajectory> summary = new ArrayList<Trajectory>();  // Stores all ion trajectories
    int[] result = new int[3];  // Stores the number of ions that exited, returned, or were crushed

    Ion ion = new Ion(m, z, ION_ENERGY);
    
    // Initialize global counters
    Nreturned = 0;
    Nexited = 0;
    Ncrushed = 0;
  
    // Loop through the number of ions to calculate their trajectories
    for(int i = 1; i <= Noi; i++) {    
        summary.addAll(calculateSingleIonTrajectory(0.0004, ion)); // spacial resolution 1 mm
        println( "Calculated trajectory of " + i + " ion");
        ion = new Ion(m, z, ION_ENERGY); 
        ionCount++;
    }
    ionTrajectory = summary;

    // Store the results of the trajectory calculation
    result[0] = Nexited; 
    result[1] = Nreturned; 
    result[2] = Ncrushed; 

    flag = 6;
    return result;
}

/*
 ***************************************************************************************************************************************
 * Method to find the resonance frequency of ions based on their trajectories
 ***************************************************************************************************************************************
 */
public void findResonanceFrequency() {
    println("lets find resonance frequency for " + testM);

    int maxExited = 0;
    int maxReturned = 0;
    int maxCrushed = 0;
    float RES_FR = 0;
  
    // Scan for the resonance frequency within the defined frequency range
    for(float f = FREQ_SCAN_LOW; f <= FREQ_SCAN_HIGH; f = f + FREQ_SCAN_STEP){
        RESONANCE_FR = f;
        int result[] = CalculateIonTrajForManyIons(1000, testCharge, testM); 
        excitationSummary.add(result);
        int Nexit = result[0];
        
        // Update maximum results when necessary
        if(Nexit >= maxExited) {
            maxExited = result[0]; 
            maxReturned = result[1];
            maxCrushed = result[2]; 
            RES_FR = RESONANCE_FR ;
        }
    }
    
    println("Resonance frequency is " + RES_FR + ", number of ions exited: " + maxExited + " returned: " + maxReturned + " crushed: " + maxCrushed);
}

/*
 ***************************************************************************************************************************************
 * Method to find transmission based on m/z ratio of ions within a defined range
 ***************************************************************************************************************************************
 */
public void findTransmission() {
    int maxExited = 0;
    int maxReturned = 0;
    int maxCrushed = 0;

    float deltaMoZ = MOZ_DELTA;
    float rangeMoZ = MOZ_RANGE;
    float centerMoZ = testM/testCharge;
    float startMoZ = centerMoZ - rangeMoZ;
    float endMoZ = centerMoZ + rangeMoZ;  

    float RES_M = 0;

    // Scan across the m/z ratio range
    for(float moz = startMoZ; moz < endMoZ; moz = moz + deltaMoZ){
        int result[] = CalculateIonTrajForManyIons(1000, testCharge, moz * testCharge); 
        excitationSummary.add(result);
        int Nexit = result[0];
        
        // Update the best results for the transmission
        if(Nexit >= maxExited) {
            maxExited = result[0]; 
            maxReturned = result[1];
            maxCrushed = result[2]; 
            RES_M = moz ;
        }
    }

    println("Mass is " + RES_M + ", number of ions exited: " + maxExited + " returned: " + maxReturned + " crushed: " + maxCrushed);
}

/** ---------------------- CALCULATE ION TRAJECTORIES FOR DIFFERENT CHARGES --------------------------------*/
public void CalculateIonTrajForDifCharges(int N) {
    ArrayList<Trajectory> summary = new ArrayList<Trajectory>();
    Ion ion = null;

    // Loop to calculate ion trajectories for different charge states
    for(int i = 1; i <= N; i++) {
        if (i % 2 == 0) ion = new Ion(2500, 5, ION_ENERGY);  // For even ions, use +5 charge state
        else { ion = new Ion(500, 1, ION_ENERGY); }  // For odd ions, use +1 charge state
        
        summary.addAll(calculateSingleIonTrajectory(0.0004, ion)); // spacial resolution 1 mm
        ionCount++;
        println("Calculated trajectory of " + i + " ion");
    }
    ionTrajectory = summary;

    flag = 6;
}

// -------------------------- CALCULATE ION TRAJECTORIES FOR DIFFERENT IONS --------------------------------
public void CalculateIonTrajForDifIons(int N) {
    ArrayList<Trajectory> summary = new ArrayList<Trajectory>();
    int sc = second();  // Current second
    int mc = minute();  // Current minute
    int hc = hour();    // Current hour
    println("Time: " + hc + ":" + mc + ":" + sc);

    Ion ion = null;

    // Loop to calculate ion trajectories for different ions
    for(int i = 1; i <= N; i++) {
        if (i % 2 == 0) ion = new Ion(576, 2, ION_ENERGY);  // Create an ion with specific mass and charge
        else { ion = new Ion(testM, 1, ION_ENERGY); }  // Create a different ion for odd iterations
        
        summary.addAll(calculateSingleIonTrajectory(0.0004, ion)); // spacial resolution 1 mm
        ionCount++;
    }
    ionTrajectory = summary;

    int es = second();  // Final second
    int em = minute();  // Final minute
    int eh = hour();    // Final hour
    int tof = eh * 3600 + em * 60 + es - (hc * 3600 + mc * 60 + sc); // Calculate time of flight
    println("End Time: " + eh + ":" + em + ":" + es + ", Total time of flight: " + tof + " seconds");
}
