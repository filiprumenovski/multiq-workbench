/**
 * Represents a point in the ion's trajectory within the ion trap.
 * Includes the position in 3D space and group information.
 */
class Trajectory {
    public float x, y, z; // Coordinates of the ion's position [m]
    int group;            // Group identifier for the ion

    /**
     * Constructor for a Trajectory point.
     * @param ix X-coordinate of the ion's position.
     * @param iy Y-coordinate of the ion's position.
     * @param iz Z-coordinate of the ion's position.
     * @param ionGroup Group identifier for the ion.
     */
    Trajectory(float ix, float iy, float iz, int ionGroup) {
        this.x = ix;
        this.y = iy;
        this.z = iz;
        this.group = ionGroup;
    }
}

/**
 * Represents the kinetic energy of an ion at a specific time.
 */
class KineticEnergy {
    public float e; // Kinetic energy in electron volts [eV]
    float t;        // Time at which the energy is sampled [s]

    /**
     * Constructor for a KineticEnergy point.
     * @param ie Kinetic energy in electron volts.
     * @param time Time in seconds when the energy was calculated.
     */
    KineticEnergy(float ie, float time) {
        this.e = ie;
        this.t = time;
    }
}

/**
 * Calculates the trajectory of a single ion in a 3D space.
 *
 * @param trajectoryResolution Minimum distance between sampled trajectory points [m].
 * @param ion The Ion object containing initial properties and dynamic variables.
 * @return A list of sampled trajectory points as Trajectory objects.
 */
 ArrayList<Trajectory> calculateSingleIonTrajectory( float trajectoryResolution, Ion ion) {
  
    int ionGroup = (int) ion.Nq;
    ion.tr = trajectoryResolution;

    // Initialize variables
    int count = 0;
    float kin = 0;
    double Ex=0, Ey=0, Ez=0;
  
    //previus position
    ion.px = ion.x; 
    ion.py = ion.y; 
    ion.pz = ion.z; 
  
    float[] Urf = new float[3];  // Array to store RF potential differences
    float[] Udc = new float[3];  // Array to store DC potential differences
    ArrayList<Trajectory> trajectory = new ArrayList<Trajectory>();  // List to store the trajectory points

    // clear previous kinetic enery list
    //ionEnergies.clear();
 
    samplingCount = 0;
 
   // Main simulation loop: runs until the ion's time exceeds the end time
    while (ion.t <= ion.end_time) {

        // Check if the ion has exited the system through the first or second exit
        if (hasLeftThroughTheExit(ion, Xexit1, Yexit1, Zexit1)) {
            recordExitEvent(ion, 1);  // Record exit event for exit 1
            break;  // Exit the loop if the ion exits through exit 1
        }
        if (hasLeftThroughTheExit(ion, Xexit2, Yexit2, Zexit2)) {
            recordExitEvent(ion, 2);  // Record exit event for exit 2
            break;  // Exit the loop if the ion exits through exit 2
        }

        // Check if the ion has crushed, and handle the event
        if (hasCrushed(ion)) {
            recordCrushEvent(ion);  // Record crush event
            break;  // Exit the loop if the ion is crushed
        }

        // Sample kinetic energy if it's time to do so
        if (isTimeToSampleKineticEnergy(ion)) {
            sampleKineticEnergy(ion);
        }
        
        // Calculate the potential differences at the ion's current position for RF and DC fields
        Urf = calculatePotentialDifference(ion.x, ion.y, ion.z, RFPOT);
        Udc = calculatePotentialDifference(ion.x, ion.y, ion.z, DCPOT);

        // Compute the RF and DC field amplitudes
        double RF = ion.Urf * Math.sin(ion.omega * ion.t);
        double DC = ion.Udc;

        // Calculate the electric field components in each direction (x, y, z)
        Ex = (Urf[0] * RF + Udc[0] * DC) / deltas + EofSpaceCharge(ion.x);
        Ey = (Urf[1] * RF + Udc[1] * DC) / deltas + EofSpaceCharge(ion.y);
        Ez = (Urf[2] * RF + Udc[2] * DC) / deltas + EofSpaceCharge(ion.z);

        // Calculate the ion's acceleration based on the electric field components
        ion.ax = (ion.q / ion.m) * Ex;
        ion.ay = (ion.q / ion.m) * Ey;
        ion.az = (ion.q / ion.m) * Ez;

        // Update the ion's position using the current velocity and acceleration
        ion.x1 = ion.x + ion.vx * ion.deltat + 0.5 * ion.ax * ion.deltat * ion.deltat;
        ion.y1 = ion.y + ion.vy * ion.deltat + 0.5 * ion.ay * ion.deltat * ion.deltat;
        ion.z1 = ion.z + ion.vz * ion.deltat + 0.5 * ion.az * ion.deltat * ion.deltat;
  
        // Update the ion's velocity based on the new position
        ion.vx = (ion.x1 - ion.x) / ion.deltat;
        ion.vy = (ion.y1 - ion.y) / ion.deltat;
        ion.vz = (ion.z1 - ion.z) / ion.deltat;

        // Update the ion's position with the new coordinates
        ion.x = ion.x1;
        ion.y = ion.y1;
        ion.z = ion.z1;

        // Increment time by the time step
        ion.t += ion.deltat;
        count++;
        
        // Sample the ion’s trajectory at regular intervals
        double distance = Math.sqrt(Math.pow(ion.x1 - ion.px, 2) + Math.pow(ion.y1 - ion.py, 2) + Math.pow(ion.z1 - ion.pz, 2));
        if (distance >= ion.tr) {
            // Add a new trajectory point if the ion has moved enough
            trajectory.add(new Trajectory((float) (ion.x1 / deltas), (float) (ion.y1 / deltas), (float) (ion.z1 / deltas), ionGroup));
            ion.px = ion.x1;
            ion.py = ion.y1;
            ion.pz = ion.z1;
        }
  
         // IF COLLISION OCCURED
        if (count >= ion.steps_bc) {
            count = 0;
            ion = collision(ion);  // Handle ion collision

            // Recalculate the mean free path (mfp) and the collision time step
            ion.mfp = -log(random(0.0,1.0)) / (ion.Np * ion.ccs);
            ion.deltaTbc = (float) (ion.mfp / Math.sqrt(ion.vx * ion.vx + ion.vy * ion.vy + ion.vz * ion.vz));
            ion.steps_bc = (int) (ion.deltaTbc / ion.deltat);
            if (ion.steps_bc == 0) ion.steps_bc = 1;
            // Icounter of # of collisions
            ion.colcount++;
        }
    
    } // close while
    
    // Check if the simulation ran out of time before completion
    if (calculationRunOutOfTime(ion)) {
        recordRunOutOfTime(ion);  // Record if the calculation ran out of time
    }

    // Print the final report after the simulation ends
    printFinalReport();

    // Return the list of trajectory points
    return trajectory;

 } 


// Print diagnostics for the ion's final state at the end of the simulation
void printDiagnostics(Ion ion) {
    // Output diagnostic information related to the final position and velocity of the ion
    // Uncomment the lines below to enable detailed diagnostic logging
    //println("*********************************************DIAGNOSTICS*********************************************");
    //println("Ion final X: " + ion.x + " [m], final Y: " + ion.y + " [m], final Z: " + ion.z + " [m] at T = " + ion.t);
    //println("Final Vx: " + ion.vx + " [m/s], Vy: " + ion.vy + " [m/s], Vz: " + ion.vz + " [m/s]");
    //println("Final Kinz: " + Kinz + " [eV], Kiny: " + Kiny + " [eV], Kinz: " + Kinz + " [eV]");
    //println("*****************************************************************************************************");
}

// Print diagnostics for the ion's state after a collision
void printDiagnosticsAfterCollision(Ion ion) {
    // Output diagnostic information related to the ion's state after a collision
    // Uncomment the lines below to enable detailed logging after collision events
    //println("VELOCITY AFTER COLLISION: ");
    //println("Vx: " + ion.vx + " [m/s], Vy: " + ion.vy + " [m/s], Vz: " + ion.vz + " [m/s]");
    // Diagnostic data after collision
    //println("-------------------------------------------------------------------------------------------------------");
    //println("Time lapsed: " + ion.t);  // Output the time elapsed since the start of the simulation
    //println("-------------------------------------------------------------------------------------------------------");
    //println("Position: x = " + ion.x + ", y = " + ion.y + ", z = " + ion.z + " [m]");
    //println("Velocity: Vx = " + ion.vx + " [m/s], Vy = " + ion.vy + " [m/s], Vz = " + ion.vz + " [m/s]");
    //println("Acceleration: ax = " + ion.ax + " [m/s^2], ay = " + ion.ay + " [m/s^2], az = " + ion.az + " [m/s^2]");
    //println("-------------------------------------------------------------------------------------------------------");
}



  
  // Check if the ion has "crushed" (left the device boundaries or entered an electrode)
boolean hasCrushed(Ion ion) {
    // Ion is considered crushed if it moves outside the defined device boundaries
    return ion.x > xDevice * 0.001 || ion.x < 0 || ion.y > yDevice * 0.001 || ion.y < 0 || ion.z > zDevice * 0.001 || ion.z < 0 || inElectrode(ion);
}

// Calculate the ion's kinetic energy (in joules)
float kineticEnergy(Ion ion) {
    // Calculate the net velocity magnitude
    float vnet = (float) Math.sqrt(ion.vx * ion.vx + ion.vy * ion.vy + ion.vz * ion.vz);
    // Return the kinetic energy using the standard formula (1/2 * m * v^2), converted to energy units
    return (float) (ion.m * (vnet * vnet) / (2 * 1.6e-19)); // Energy in [J]
}

// Calculate the ion's kinetic energy in electron volts (eV)
float getIonKineticEnergy(Ion ion) { // Returns kinetic energy in [eV]
    // Use the mass (MW) and velocity components to compute kinetic energy in eV
    return 0.5 * ion.MW * 1.66 * 1e-8 * (float) ((ion.vx * ion.vx + ion.vy * ion.vy + ion.vz * ion.vz) / 1.6022); // [eV]
    // println("KINETIC ENERGY IS " + kin);       
}

// Check if the ion has exited through one of the defined exit points
boolean hasLeftThroughTheExit(Ion ion, float Xex, float Yex, float Zex) {
    // Calculate the radial distance from the exit point (Xex, Yex)
    double exrad = Math.sqrt((ion.x - Xex * 0.001) * (ion.x - Xex * 0.001) + (ion.y - Yex * 0.001) * (ion.y - Yex * 0.001));
    // Calculate the distance along the Z-axis from the exit point (Zex)
    double exitz = Math.sqrt((ion.z - Zex * 0.001) * (ion.z - Zex * 0.001));
    // Return true if the ion is within a radius of 6mm and a distance of 1mm along the Z-axis
    return exrad < 0.006 && exitz < 0.001;
}

// Record an exit event when an ion leaves the trap
void recordExitEvent(Ion ion, int exit) {
    // For multiply charged ions
    if (ion.Nq > 1) { 
        if (exit == 1) { 
            // multiplyChargedExited1++; 
        } else {
            // multiplyChargedExited2++; 
        }
        // Accumulate the time at which the exit event occurred
        averageTimeExit5 = averageTimeExit5 + (float) ion.t;
        // println("5+ ions exit at T = " + ion.t);
    } else { 
        // For singly charged ions
        if (exit == 1) { 
            // singlyChargedExited1++; 
        } else {
            // singlyChargedExited2++; 
        }
        // Accumulate the time at which the exit event occurred
        averageTimeExit1 = averageTimeExit1 + (float) ion.t;
        // println("1+ ions exit at T = " + ion.t);
    }
    // Optionally print diagnostics
    // printDiagnostics(ion);
}

// Record a crash event when an ion is destroyed by the trap walls
void recordCrushEvent(Ion ion) {
    // Print message to indicate ion was crushed
    println("ION CRUSHED ON THE WALL OF THE ION TRAP!!!");
    // Increment the crushed ion counters based on charge state
    if (ion.Nq > 1) { 
        mIonsCrushed++; 
    } else { 
        sIonsCrushed++; 
    }
    // Print diagnostics after the ion crush event
    printDiagnostics(ion);
}

// Check if the simulation has run out of time (i.e., the ion has reached its end time)
boolean calculationRunOutOfTime(Ion ion) { 
    return ion.t >= ion.end_time;
}

// Record an event when the ion has stayed in the trap until the end of the simulation
void recordRunOutOfTime(Ion ion) {
    println("ION STAYED IN THE ION TRAP UNTIL THE END OF TIME!!"); 
    // Increment counters based on the ion's charge state
    if (ion.Nq > 1) { 
        mNtrap++; 
    } else { 
        sNtrap++; 
    }
}

// Check if it’s time to sample the ion's kinetic energy based on the sampling interval
boolean isTimeToSampleKineticEnergy(Ion ion) {   
    // Return true if the current time has reached the next sampling interval
    if (ion.t >= SAMPLING_TIME * samplingCount) { 
        samplingCount++; 
        return true; 
    } else {
        return false; 
    }
}

// Sample and record the ion's kinetic energy
void sampleKineticEnergy(Ion ion) {   
    // Add the current kinetic energy value to the list of sampled energies
    ionEnergies.add(new KineticEnergy(getIonKineticEnergy(ion), (float) ion.t));
}
// Print the final report summarizing the simulation results
void printFinalReport() {
    println("************************************************************************************** ");
    println("                                  FINAL REPORT                                         ");
    println("************************************************************************************** ");
    // Print the total number of ions processed in the simulation
    println("Total number of ions run  " + ionCount);
    println("-------------------------------------------------------------------------------------- ");
    // Print the total number of ions that stayed in the trap, split by charge state
    println("Total number of ions stayed in the trap:  " + (mNtrap + sNtrap));
    println("Total number of SINGLY CHARGED IONS stayed in the trap:  " + sNtrap);
    println("Total number of MULTIPLY CHARGED IONS stayed in the trap:  " + mNtrap);
    println("-------------------------------------------------------------------------------------- ");  
    // Additional statistics on crushed ions, exits, etc., can be added as needed
    // The following lines are currently commented out but can be enabled if needed:
    // println("Total number of ions crushed on the walls of the ion trap: " + (mIonsCrushed + sIonsCrushed));
    // println("Total number of SINGLY CHARGED IONS crushed on the walls of the ion trap: " + sIonsCrushed);
    // println("Total number of MULTIPLY CHARGED IONS crushed on the walls of the ion trap: " + mIonsCrushed);
    // println("-------------------------------------------------------------------------------------- ");
    // println("Total number of IONS exited from the trap: " + (multiplyChargedExited1 + singlyChargedExited1 + multiplyChargedExited2 + singlyChargedExited2));
    // println("Total number of SINGLY CHARGED IONS exited from the trap: " + singlyChargedExited + "(in average " + averageTimeExit1 / (singlyChargedExited1 + singlyChargedExited2) + " s)");
    // println("Total number of MULTIPLY CHARGED IONS exited from the trap: " + multiplyChargedExited + "(in average " + averageTimeExit5 / (multiplyChargedExited1 + multiplyChargedExited2) + " s)");
    // println("-------------------------------------------------------------------------------------- ");
    // println("THROUGH EXIT 1 exited " + singlyChargedExited1 + " singly charged ions and " + multiplyChargedExited1 + " multiply charged ions");
    // println("THROUGH EXIT 2 exited " + singlyChargedExited2 + " singly charged ions and " + multiplyChargedExited2 + " multiply charged ions");
  
    // Flag for additional operations can be added if necessary
    // flag = 4; 
}

// Approximates the electric field created by space charge at a given distance
double EofSpaceCharge(double d) {
    // Convert distance to millimeters
    float x = 1000 * (float) d;
    // Use empirical formula to calculate the electric field strength (approximation)
    float E = -(xDevice / 2 - x) / 30; // Approximation gives ~0.1 V/m for ~10^8 ions at the ion cloud border
    // Factor for scaling the field strength
    float k = 0.1; // Scaling factor for the electric field
    return (double) k * E; // Return the electric field value
}

// Adds a component to an existing electric field value
double addComponent(double E, double Component) {
    // Add the given component (e.g., from another source) to the existing field strength
    return E + Component;
}
     
