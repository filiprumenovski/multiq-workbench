/**
 *************************************************************************************************************************
 *
 *                                                   Ion Trajectories
 *
 *************************************************************************************************************************
 */

/**
 * Calculates the trajectory of a single ion based on its initial conditions and specified resolution.
 * The ion's path is tracked over time until it exits, returns to the MQ, or is crushed.
 * 
 * @param trajectoryResolution The resolution for trajectory sampling
 * @param ion The ion object whose trajectory is to be calculated
 * @return A list of trajectory points for the ion
 */
ArrayList<Trajectory> calculateSingleIonTrajectory(float trajectoryResolution, Ion ion) {
  
  int ionGroup = (int)ion.Nq;  // Group identifier based on ion charge state
  ion.tr = trajectoryResolution;  // Set trajectory resolution
  
  // Initialize counters and variables
  int count = 0;
  samplingCount = 0;
  
  // Store previous ion position for trajectory calculation
  ion.px = ion.x; 
  ion.py = ion.y; 
  ion.pz = ion.z; 
  
  // Arrays to store electric potential components
  float [] Urf = new float[3];
  float [] Udc = new float[3];
  float [] Uac = new float[3];
 
  // List to store the trajectory points
  ArrayList<Trajectory> trajectory = new ArrayList<Trajectory>();
  float kin = 0;
 
  // Clear previous kinetic energy list
  ionEnergies.clear();
 
  double Ex = 0, Ey = 0, Ez = 0;
  float phase = random(0, 2*PI);  // Random phase for AC excitation
  
  //===========================================  CALCULATIONS  ==============================================
  
  // Start calculating the ion trajectory over time
  while (ion.t <= ion.end_time) {
    
    // CHECK FOR DIFFERENT SCENARIOS    
    if (ionExited(ion)) { 
      println("EXITED"); 
      Nexited++; 
      return trajectory; 
    }
    if (ionReturnedToMQ(ion)) { 
      println("RETURNED TO MQ"); 
      Nreturned++; 
      return trajectory; 
    }
    if (ionCrushed(ion)) { 
      println("ION CRUSHED!"); 
      Ncrushed++; 
      return trajectory; 
    } 
    
    // PERFORM CALCULATIONS:
    // Calculate the difference in potentials between adjacent grid points
    Urf = calculatePotentialDifference(ion.x, ion.y, ion.z, RFPOT); 
    Udc = calculatePotentialDifference(ion.x, ion.y, ion.z, DCPOT);
    Uac = calculatePotentialDifference(ion.x, ion.y, ion.z, ACPOT);

    // Calculate the electric field components based on the applied potentials and frequencies
    double RF = A_RF * Math.sin(2 * PI * RF_FREQUENCY * 1000 * ion.t);  // RF field
    double AC = V_AC * Math.sin(2 * PI * RESONANCE_FR * 1000 * ion.t + phase);  // AC field
    double DC = V_DC;  // DC field (constant)
   
    // Electric field components in x, y, and z directions
    Ex = Urf[0] * RF/deltas + Uac[0] * AC/deltas + Udc[0] * DC/deltas;
    Ey = Urf[1] * RF/deltas + Uac[1] * AC/deltas + Udc[1] * DC/deltas;
    Ez = Urf[2] * RF/deltas + Uac[2] * AC/deltas + Udc[2] * DC/deltas;
   
    // Calculate the ion's acceleration based on the electric field and its charge/mass
    ion.ax = (ion.q/ion.m) * Ex; 
    ion.ay = (ion.q/ion.m) * Ey; 
    ion.az = (ion.q/ion.m) * Ez; 
    
    // Update the ion's position based on its velocity and acceleration
    ion.x1 = ion.x + ion.vx * ion.deltat + ion.ax * ion.deltat * ion.deltat/2;
    ion.y1 = ion.y + ion.vy * ion.deltat + ion.ay * ion.deltat * ion.deltat/2;
    ion.z1 = ion.z + ion.vz * ion.deltat + ion.az * ion.deltat * ion.deltat/2;
   
    // Update the ion's velocity
    ion.vx = (ion.x1 - ion.x)/ion.deltat;
    ion.vy = (ion.y1 - ion.y)/ion.deltat; 
    ion.vz = (ion.z1 - ion.z)/ion.deltat; 
       
    // Update the ion's position to the new coordinates
    ion.x = ion.x1;
    ion.y = ion.y1;
    ion.z = ion.z1;
       
    // Update the ion's time
    ion.t = ion.t + ion.deltat;
    count++;  // Increment the time step counter
   
    // Add the ion's position to the trajectory if the movement exceeds the trajectory resolution
    if (Math.sqrt((ion.x1 - ion.px) * (ion.x1 - ion.px) + (ion.y1 - ion.py) * (ion.y1 - ion.py) + (ion.z1 - ion.pz) * (ion.z1 - ion.pz)) >= ion.tr) { 
      trajectory.add(new Trajectory((float)(ion.x1/deltas), (float)(ion.y1/deltas), (float)(ion.z1/deltas), ionGroup)); // Convert to millimeters
      ion.px = ion.x1; 
      ion.py = ion.y1; 
      ion.pz = ion.z1;
   }
   
   //...............................................   COLLISION   ..................................................
   if (count >= ion.steps_bc) { // Check if collision occurred
      count = 0;     
      ion = collision(ion);  // Handle the collision
      
      // Calculate new mean free path (mfp)
      ion.mfp = -log(random(0.0000001, 1.0))/(ion.Np * ion.ccs); 
      ion.deltaTbc = (float) (ion.mfp/Math.sqrt(ion.vx * ion.vx + ion.vy * ion.vy + ion.vz * ion.vz));
      ion.steps_bc = (int) (ion.deltaTbc/ion.deltat);
      if (ion.steps_bc == 0) ion.steps_bc = 1;  // Avoid division by zero
      ion.colcount++;  // Increment the collision counter         
   }
    
  } // close while
  
  // Handle case when the calculation runs out of time
  if (calculationRunOutOfTime(ion)) { 
    recordRunOutOfTime(ion); 
  }
  
  return trajectory;  // Return the calculated trajectory
} // close function

  
 
  
void printDiagnostics (Ion ion) {
   
      println("*********************************************DIAGNOSTICS*********************************************");
      println("Ion final X: " +ion.x+ " [m] final Y: " + ion.y +" [m] final Z: " + ion.z  +" [m] at T = " + ion.t );
      println( "Final Vx: " +ion.vx + "[m/s],  Vy: " +ion.vy + " [m/s],  Vz: " +ion.vz + " [m/s]");
      //println( "Final Kinz: " +Kinz+ "[eV], Kiny " +Kiny + " [eV], Kinz: " +Kinz + " [eV]"); 
      println("*****************************************************************************************************");
      
  }
  
  void printDiagnosticsAfterCollision (Ion ion) {
      //println("VELOCITY AFTER COLLISION: ");
      //println("Vx : " + ion.vx +" m/z Vy: " + ion.vy + " m/s Vz: " +ion.vz  + " m/s");
      // Diagnostic data
      //println("-------------------------------------------------------------------------------------------------------");
      //println("time lapsed: " + ion.t); 
      //println("-------------------------------------------------------------------------------------------------------");
      //println("x : " + ion.x +" y: " + ion.y + " z: " +ion.z  + " [m]  ");
      //println("Vx : " + ion.vx +" Vy: " + ion.vy + " Vz: " +ion.vz  + " [m/s]  ");
      //println("ax : " + ion.ax +" ay: " + ion.ay + " az: " +ion.az  + " m/s^2  ");
      //println("-------------------------------------------------------------------------------------------------------");
      //iii
  }
  
  boolean hasCrushed(Ion ion) {
    return ion.x>xDevice*0.001 || ion.x<0 || ion.y>yDevice*0.001 || ion.y<0 || ion.z>zDevice*0.001 || ion.z<0 ;
  }
  
  float kineticEnergy(Ion ion) {   
    float vnet= (float) Math.sqrt(ion.vx*ion.vx+ion.vy*ion.vy+ion.vz*ion.vz);
    return (float) (ion.m*(vnet*vnet)/(2*1.6e-19));   
  }
  
  float getIonKineticEnergy (Ion ion) {
     return 0.5*ion.m*(float)((ion.vx*ion.vx+ion.vy*ion.vy+ion.vz*ion.vz)/1.6022e-19); // in [eV]
  //println("KINETIC ENERGY IS " + kin);       
  }
  
  boolean hasLeftThroughTheExit(Ion ion, float Xex, float Yex, float Zex) {
   double exrad = Math.sqrt((ion.x-Xex*0.001)*(ion.x-Xex*0.001) + (ion.y-Yex*0.001)*(ion.y-Yex*0.001));
   double exitz = Math.sqrt((ion.z-Zex*0.001)*(ion.z-Zex*0.001));
   return exrad<0.006 && exitz<0.001; // within radious of 6mm and distnace 1 mm
  }
  
  
    
  void recordCrushEvent (Ion ion) {
    println("ION CRUSHED ON THE WALL OF THE  ION TRAP!!!");
    if (ion.Nq>1) { mIonsCrushed++; } else {sIonsCrushed++; }
    printDiagnostics (ion);
   
  }
  
  boolean calculationRunOutOfTime(Ion ion) { 
    return ion.t>=ion.end_time;
  }
  
  void recordRunOutOfTime(Ion ion){
    //println("ION STAYED IN THE ION TRAP UTILL THE END OF TIME!!"); 
    if (ion.Nq>1) { mNtrap++;} else { sNtrap++;}
  }
  
  boolean isTimeToSampleKineticEnergy (Ion ion ) {   
    if(ion.t>=SAMPLING_TIME*samplingCount) {samplingCount++; return true; }
    else return false;
  }
  void sampleKineticEnergy(Ion ion) {   
     ionEnergies.add(new KineticEnergy(getIonKineticEnergy(ion), (int) ion.Nq));
  }
 
   
   void printFinalReport() {
   println("************************************************************************************** ");
   println("                                  FINAL REPORT                                         ");
   println("************************************************************************************** ");
   println("Total number of ions run  " + ionCount);
   println("-------------------------------------------------------------------------------------- ");
   println("Total number of ions styed in the trap:  " + (mNtrap + sNtrap));
   println("Total number of SINGLY CHARGED IONS stayed in the trap:  " + sNtrap);
   println("Total number of MULTIPLY CHARGED IONS stayed in the trap:  " + mNtrap);
    
  //flag = 4; 
   }
   
   // Adds approxiamtion of an electric field crated by a space charge
   double EofSpaceCharge( double d) {
     float x = 1000*(float) d;
     float E = -(xDevice/2 - x)/30; //impirical formula to give ~ 0.1 V/m(~10^8 ions) at the border of ion cloud (see Coulomb_CUBE program)
     float k = 0.1; //0.1; //0.1;
     return (double) k*E;
   }
   
   // Adds a component to an electrical field     
  double addComponent(double E, double Component) {
    return E+Component;
  }
 
   
//------------------------------------------------------------------------------------------------------------------------
//                                              CONDITIONS
//------------------------------------------------------------------------------------------------------------------------

/**
 * Checks if the ion has exited the device region based on its z-coordinate.
 * @param ion The ion object to check
 * @return true if the ion has exited, false otherwise
 */
boolean ionExited(Ion ion) {
     boolean result = ion.z > ((zDevice - 0.5) * 0.001); // Ion exits if z-coordinate is greater than the defined threshold
     return result;
}
    
/**
 * Checks if the ion has returned to the main quadrupole (MQ) region based on its z-coordinate.
 * @param ion The ion object to check
 * @return true if the ion has returned, false otherwise
 */
boolean ionReturnedToMQ(Ion ion) {
     boolean result = ion.z < (Zenter * 0.001); // Ion returns if z-coordinate is less than the Zenter threshold
     return result;
} 

/**
 * Checks if the ion has been crushed (either by hitting an electrode or moving out of bounds).
 * @param ion The ion object to check
 * @return true if the ion is crushed, false otherwise
 */
boolean ionCrushed(Ion ion) {
     // Ion is crushed if it hits an electrode or if it moves outside the defined device bounds
     return ((pointIsElectrode(ion.x, ion.y, ion.z) && Qstart * 0.001 < ion.z && ion.z < (Qstart + Q1) * 0.001) 
             || (ion.x > (xDevice + 1) * 0.001 || ion.y > (yDevice + 1) * 0.001 || ion.x < -0.001 || ion.y < -0.001));
}
