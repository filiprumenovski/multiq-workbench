/**
 ***********************************************************************************************************************
 *
 *                                                 PROGRAM OPTIONS
 *
 ***********************************************************************************************************************
 */  

// Handle key presses for various program modes and actions
void keyPressed(){
   
  // Example: Pressed "0" to set flag to 0
  if(keyCode == 48) { // "0" key
    flag = 0;
  }
   
  // Example: Pressed "7" to set flag to 7
  if(keyCode == 55) { // "7" key
    flag = 7;
  }
   
  // INJECT one singly-charged ion when "i" is pressed
  if(keyCode == 73) { // "i" key
    // Ion parameters: m (mass), z (charge), Vin (initial energy), Vrf (RF voltage), RF (frequency), phase, P (pressure)
    Ion ion = new Ion(500, 1, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0); 
    CalculateIonTrajSingleIon(ion);  // Calculate the trajectory of the injected ion
    flag = 6;  // Set flag to 6 (ion trajectory mode)
  }
  
  // INJECT one multiply-charged ion when "m" is pressed
  if(keyCode == 77) { // "m" key
    // Ion parameters: m (mass), z (charge), Vin (initial energy), Vrf (RF voltage), RF (frequency), phase, P (pressure)
    Ion ion = new Ion(1500, 3, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0); 
    CalculateIonTrajSingleIon(ion);  // Calculate the trajectory of the multiply-charged ion
    flag = 6;  // Set flag to 6 (ion trajectory mode)
  }

  //==================================== FLAG OPTIONS (for changing display modes) ================================
  if(keyCode == KeyEvent.VK_1) {  // Pressed "1"
    flag = 1;
  }
    
  if(keyCode == KeyEvent.VK_2) {  // Pressed "2"
    flag = 2;
  }
    
  if(keyCode == KeyEvent.VK_3) {  // Pressed "3"
    flag = 3;
  }
  
  if(keyCode == KeyEvent.VK_4) {  // Pressed "4"
    flag = 4;
  }
  
  if(keyCode == KeyEvent.VK_5) {  // Pressed "5"
    flag = 5;
  }
  
  if(keyCode == KeyEvent.VK_6) {  // Pressed "6"
    flag = 6;
  }

  //============ MODIFY GEOMETRY (add/remove sides) =====================
  if(keyCode == 45) { // Pressed "-" (subtract side from geometry)
    println("keyCode is " + keyCode);
    cut = 0;  // Disable geometry modification
  }
  
  if(keyCode == 61) { // Pressed "+" (add side to geometry)
    println("keyCode is " + keyCode);
    cut = 1;  // Enable geometry modification
  }
  
  //============= FLAG OPTIONS (for displaying additional settings) =====
  if(keyCode == 97 || keyCode == 112) { // Pressed F1
    showOptions_1();  // Show options for geometry and potentials
  }
  
  if(keyCode == 98) { // Pressed F2
    showOptions_2();  // Show additional options
  }

  //============================= MOVE MultiQ using arrow keys =====================
  if(keyCode == KeyEvent.VK_UP) {
    a1 += PI / 16;  // Move MultiQ in the positive direction (UP)
  }
  
  if(keyCode == KeyEvent.VK_DOWN) {
    a1 -= PI / 16;  // Move MultiQ in the negative direction (DOWN)   
  }
  
  if(keyCode == KeyEvent.VK_LEFT) {
    a2 += PI / 16;  // Move MultiQ in the positive direction (LEFT)  
  }
  
  if(keyCode == KeyEvent.VK_RIGHT) {
    a2 -= PI / 16;  // Move MultiQ in the negative direction (RIGHT)
  }
  
  if(keyCode == KeyEvent.VK_PAGE_UP) {
    a3 += PI / 16;  // Move MultiQ in the positive direction (PAGE UP)
  }
  
  if(keyCode == KeyEvent.VK_PAGE_DOWN) {
    a3 -= PI / 16;  // Move MultiQ in the negative direction (PAGE DOWN)
  }
  
  //========================= ADJUST Z-COORDINATE using "," or "." =======================
  if(keyCode == 44) { // Pressed "," (decrease z-coordinate)
    if(zcs >= 2 && zcs <= Nz - 1) { 
      zcs -= 1;  // Move in the negative z-direction
      redraw();   // Redraw the screen
    }
  }
  
  if(keyCode == 46) { // Pressed "." (increase z-coordinate)
    if(zcs >= 1 && zcs <= Nz - 2) { 
      zcs += 1;  // Move in the positive z-direction
      redraw();   // Redraw the screen
    }
  }

} // End of keyPressed()


/**
 *************************************************************************************************************************
 *
 *                                               1. Geometry & Potentials
 *
 *************************************************************************************************************************
 */

public void showOptions_1() {
  SwingUtilities.invokeLater(new Runnable() { // Run as a separate thread      
    public void run() { 
      String[] options = new String[] { 
        "Create Geometry",
        "Calculate Potential",
        "Save Potential",
        "Read Potential",
        "Calculate Stability Diagrams" 
      };
      
      String input = (String) JOptionPane.showInputDialog(
        new JFrame(),
        "Please select the action",
        "Geometry & Potentials", JOptionPane.INFORMATION_MESSAGE,
        new ImageIcon("java2sLogo.GIF"), options, "Demiurge Actions");

      // Handle user selection and call corresponding methods
      if(input == "Create Geometry") {
        CreateElectrodeGeometry();
      }
      if(input == "Calculate Potential") {
        CalculatePotentialDistribution();  
      }
      if(input == "Save Potential") {
        SavePotential();
      }
      if(input == "Read Potential") {
        ReadPotential();
      }
      if(input == "Calculate Stability Diagrams") {
        // Calculate stability for different RF amplitudes
        CalculateApex(RF_AMPLITUDE/2); 
        CalculateApex(RF_AMPLITUDE);
        CalculateApex(RF_AMPLITUDE*2);
        flag = 8;
      }
    }
  });
}

//-------------------------------------- CREATE ELECTRODE GEOMETRY ----------------------------------------
public void CreateElectrodeGeometry() {
  println("Will calculate RF new geometry");
  Electrodes = calculateElectrodeGeometry();  // Calculate geometry of electrodes
  flag = 1;  // Set flag for geometry display mode
}

//------------------------------------ CALCULATE POTENTIAL DISTRIBUTION -----------------------------------
public void CalculatePotentialDistribution() {
  println("Calculating potential for a given electrode geometry");
  POT = calculatePotential(0.001, POT);  // Calculate potential distribution with specified precision
}

//---------------------------------------------- I/O OPERATIONS --------------------------------------------------------

/**
 * Selecting file for saving the potential data.
 */
public void SavePotential() {
  selectOutput("Select a file to write to:", "savePotentialIntoFile");
}

/**
 * Save potential data to the selected file.
 */
void savePotentialIntoFile(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
  } else {
    println("User selected " + selection.getAbsolutePath());

    PrintWriter output = createWriter(selection.getAbsolutePath());

    // Write potential data to file
    int e = 0;
    for (int s = 0; s < Nz; s++) {
      if (s == 0) {
        // Write the first row in the file with device settings
        output.println("X: " + xDevice + " [mm], Y: " + yDevice + " [mm], Z: " + zDevice + " [mm], SCALE: " + scl);
      } else {
        for (int i = 0; i < Nx; i++) {
          for (int j = 0; j < Ny; j++) {
            // Determine if the electrode exists for this point
            e = POT[i][j][s].electrode ? 1 : 0;
            output.println(i + "," + j + "," + s + "," + POT[i][j][s].potential + "," + e);
          }
        }
      }
    }
    
    output.flush();  // Write remaining data
    output.close();  // Close the file
  }
}

/**
 * Selecting file to read potential data from.
 */
void ReadPotential() {
  selectInput("Select a file to process:", "openPotentialFile");  
}

/**
 * Open and read the potential data from the selected file.
 */
void openPotentialFile(File selection) {
  if (selection == null) {
    println("Window was closed or the user hit cancel.");
  } else {
    Electrodes = openPotential(selection.getAbsolutePath(), POT);
    // println("User selected " + selection.getAbsolutePath());
  }  
}

/**
 * Opens a potential file and parses its contents.
 */
ArrayList<Potential> openPotential(String fileName, Potential[][][] pot) {
  startTimer();

  File chosenFile = new File(fileName);
  
  try {
    Scanner fileScanner = new Scanner(chosenFile);
    String firstLine = fileScanner.nextLine();
    String[] fileData = split(firstLine, ",");
    
    // Parse first line for device dimensions and scale factor
    float dimx = Float.parseFloat(split(fileData[0], " ")[1]);
    float dimy = Float.parseFloat(split(fileData[1].trim(), " ")[1]);
    float dimz = Float.parseFloat(split(fileData[2].trim(), " ")[1]);
    float sc = Float.parseFloat(split(fileData[3].trim(), " ")[1]);
    int dim = floor((dimx * sc) * (dimy * sc) * (dimz * sc));  // Calculate total number of points
    
    int k = 0;
    int progress = 0;
    ArrayList<Potential> potout = new ArrayList<Potential>();  // To store potential data
    
    // Read the file line by line
    while (fileScanner.hasNext()) {
      int job = floor(float(k) / dim * 100);
      if (job == progress) {
        print("|");
        if (progress % 10 == 0) { print(progress); }
        if (progress == 99) { println("100%"); println(""); }
        progress++;
      }
      
      String line = fileScanner.nextLine();
      String[] linedata = split(line, ",");
      
      // Parse the data for each point
      float nx = Float.parseFloat(linedata[0]);
      float ny = Float.parseFloat(linedata[1]);
      float nz = Float.parseFloat(linedata[2]);
      float potential = Float.parseFloat(linedata[3]);
      boolean el = Float.parseFloat(linedata[4]) > 0;
      
      // Store the potential data for each point
      pot[int(nx)][int(ny)][int(nz)] = new Potential(nx * sc, ny * sc, nz * sc, potential, el);
      
      // Count electrode points (where potential > 0)
      if (el) {
        potout.add(new Potential(nx, ny, nz, potential, el));
      }
      
      k++;  // Increment counter
    }

    println("The size of the array is: " + potout.size());
    Electrodes = potout;  // Store the electrodes data
    fileScanner.close();
    stopTimer();

    flag = 1;  // Set flag for geometry mode
    return potout;

  } catch (IOException e) {
    System.err.println("Caught IOException: " + e.getMessage());
    return null;
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
  
 SwingUtilities.invokeLater(new Runnable() { // run as a separate thread      
    public void run() { 
      String[] plays = new String[] { "Load Precalculated Potentials",
                                      "Calculate ion trajectories for Ntarget ions",
                                      "Calculate Cahrge Deplition Effect",
                                      "Save ion trajectories in the file",                                    
                                    };
      String input = (String) JOptionPane.showInputDialog(
        new JFrame(),
        "Display Options",
        "Ion Motion", JOptionPane.INFORMATION_MESSAGE,
        new ImageIcon("java2sLogo.GIF"), plays, "DisplayOptions");
        if(input == "Load Precalculated Potentials") {
          println( "...loading RF potential from file: " + DEFAULT_RFPOT_FILE);
          openPotential(DEFAULT_RFPOT_FILE, RFPOT);
          println( "...loading DC potential from file: " + DEFAULT_DCPOT_FILE);
          openPotential(DEFAULT_DCPOT_FILE, DCPOT);
        }   
        
        if(input == "Calculate ion trajectories for Ntarget ions") CalculateIonTrajForManyIons(Ntarget);
        if(input == "Calculate Cahrge Deplition Effect") CalculateIonTrajForDifCharges(Ntarget);
        if(input == "Save ion trajectories in the file") SaveIonTrajectories();
        // System.out.println("User's input: " + input);
    }
  });

}

/**
 ********************************************************************************************
 *
 *              Opens the existing file defining the geometry of the potental 
 *
 ********************************************************************************************
 */


/** -------------------------- CALCULATE ION TRAJECTORIES FOR 1 ION  ------------------------------------*/
public void CalculateIonTrajSingleIon(Ion ion) {
println("Calculate ion trajectory");
                
 ionTrajectory = calculateSingleIonTrajectory(0.0001, ion); // spacialResolution 0.1 mm
 flag = 6;
}


/** -------------------------- CALCULATE ION TRAJECTORIES FOR N IONS  --------------------------------*/
public void CalculateIonTrajForManyIons(int N) {
  ArrayList<Trajectory> summary = new ArrayList<Trajectory>();
  int sc = second();  // Values from 0 - 59
  int mc = minute();  // Values from 0 - 59
  int hc = hour();
  println( "Time "+hc+":"+mc+":"+sc);

  Ion ion = new Ion(1500, 3, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0);
  for(int i=1; i<=N; i++) {    
    summary.addAll(calculateSingleIonTrajectory(0.0004, ion)); // spacial resolution 1 mm
    println( "Calculated trajectory of " + i + " ion");
    ion = new Ion(1500, 3, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0); 
    ionCount++;
  }
  ionTrajectory = summary;
  
  //println(singlyChargedExited + " ions came out through the designated exit" );
  //println(multiplyChargedExited + " ions came out through the designated exit" );
 
  int es = second();  
  int em = minute();
  int eh = hour();
  int tof = eh*3600+em*60 +es -(hc*3600 + mc*60 +es);
  println( "It took "+tof+" seconds to parse the file");
 
 
  flag = 6;
}


/** ---------------------- CALCULATE ION TRAJECTORIES FOR DIFFERENT IONS  --------------------------------*/
public void CalculateIonTrajForDifCharges(int N) {
  ArrayList<Trajectory> summary = new ArrayList<Trajectory>();
  int sc = second();  // Values from 0 - 59
  int mc = minute();  // Values from 0 - 59
  int hc = hour();
  println( "Time "+hc+":"+mc+":"+sc);

  Ion ion = null;

  for(int i=1; i<=N; i++) {
    
    if (i%2==0) ion = new Ion(2500, 5, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0);
    else { ion = new Ion(500, 1, ION_ENERGY, RF_AMPLITUDE, RF_FREQUENCY, 0); }
    
    summary.addAll(calculateSingleIonTrajectory(0.0004, ion)); // spacial resolution 1 mm
    ionCount++;
    println( "Calculated trajectory of " + i + " ion");
  }
  ionTrajectory = summary;
  
  //println(singlyChargedExited + " singly charged ions came out through the designated exit" );
  //println(multiplyChargedExited + " multiply charged ions came out through the designated exit" );
 
  int es = second();  
  int em = minute();
  int eh = hour();
  int tof = eh*3600+em*60 +es -(hc*3600 + mc*60 +es);
  println( "It took "+tof+" seconds to parse the file");
 
 
  flag = 6;
}


public void SaveIonTrajectories() {
  println("Save ion trajectories in the file");
} 


/**
 ********************************************************************************************
 *
 *              Opens the existing file defining the geometry of the potental 
 *
 ********************************************************************************************
 */
ArrayList<Potential> readGeometry(String fileName, String swtch) {

  startTimer();
  
  File chosenFile =  new File(fileName);
  
  try {
     Scanner fileScanner = new Scanner(chosenFile);
     String firstLine = fileScanner.nextLine();
     String[] fileData = split(firstLine, ",");
  
     float dimz = Float.parseFloat(split(fileData[0], " ")[1]);
     float dimx = Float.parseFloat(split(fileData[1].trim(), " ")[1]);
     float dimy = Float.parseFloat(split(fileData[2].trim(), " ")[1]);
     float sc   = Float.parseFloat(split(fileData[3].trim(), " ")[1]); 
     int dim = floor((dimz*sc) * (dimx*sc) * (dimy*sc)); 
     
     
     int k =0;
     int progress = 0;
     
     ArrayList<Potential> out = new ArrayList<Potential>(); 
  
     while ( fileScanner.hasNext() ) {
     
       int job = floor(float(k)/dim*100); 
       if (job == progress ) {
         print("|");
         if (progress%10 == 0) { print(progress);}
         if (progress == 99) { println("100%"); println("");}
         progress++;
       }  
    
       String line = fileScanner.nextLine();
       String[] linedata = split(line, ",");
       // println(linedata[k]);
       float nz = Float.parseFloat(linedata[0] );
       float nx = Float.parseFloat( linedata[1] );
       float ny = Float.parseFloat( linedata[2] );   
       float potential = float( linedata[3] );
       boolean el = false;
       float boo = float( linedata[4] );
       if(boo>0) el=true;
       POT[int(nx)][int(ny)][int(nz)] = new Potential(nx*sc, ny*sc, nz*sc, potential, el);
 
       if(el) {
         out.add( new Potential(nx, ny, nz, potential, el) ); 
       } 
      
       k++;  
    }
    println( "the size of the array is: " + out.size());
      
    stopTimer();
    fileScanner.close();
    
    flag = 1;
    //geometryFileExists = true;
 
    return out;
    
  } catch (IOException e) {
    System.err.println("Caught IOException: " + e.getMessage());
    return null;
  }
 
  
} // close v
