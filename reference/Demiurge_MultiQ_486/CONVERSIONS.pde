//-----------------------------------------------------------------------------------------------------------------------
// Convert a coordinate in meters to an array index
int convertMetersToArray(double coordinate) {
   int Npoint = (int)(coordinate * scl * 1000.0);
   return Npoint;
}

// Convert a coordinate in millimeters to an array index
int convertMillimetersToArray(float c) {
   int Npoint = (int)(c * scl);
   return Npoint;
}

//-----------------------------------------------------------------------------------------------------------------------
// Coordinate system translations (in [mm])
float[] translateToCenterCoordinateSystem(float x, float y, float z) { 
   // Translates coordinates to the center of the device coordinate system
   float[] result = new float[3];
   result[0] = x - xDevice/2;
   result[1] = y - yDevice/2;
   result[2] = z - zDevice/2;
   return result; 
}

float[] translateFromCenterCoordinateSystem(float x, float y, float z) { 
   // Translates coordinates from the center back to the device coordinate system
   float[] result = new float[3];
   result[0] = xDevice/2 + x;
   result[1] = yDevice/2 + y;
   result[2] = zDevice/2 + z;
   return result; 
}

//-----------------------------------------------------------------------------------------------------------------------
// Timer functionality
public void startTimer() {
   // Starts the timer and prints the start time
   Sstart = second();
   Mstart = minute();
   Hstart = hour();
   Tstart = millis();
   println("The task started at " + Hstart + ":" + Mstart + ":" + Sstart); 
   println("Dimension: " + Nx + "x" + Ny + "x" + Nz);
}

public void stopTimer() {
   // Stops the timer and prints the elapsed time
   Send = second();
   Mend = minute();
   Hend = hour();
   Tend = millis();
   println("The task ended at " + Hend + ":" + Mend + ":" + Send);
   println("Elapsed time: " + ((Tend - Tstart)/3600000) + "h:" + ((Tend - Tstart)/60000 % 60) + "m:" + ((Tend - Tstart)/1000 % 60) + "s");
}

//-----------------------------------------------------------------------------------------------------------------------
// Apex calculation
void CalculateApex(float rfamplitude) {
   // Calculates optimal m/z and parameters for ion transmission
   float moverz = getMoverZ(0.706, rfamplitude);
   float cutoff = moverz * 0.706/0.905;
   float U = 0.237 * rfamplitude/(2 * 0.706);
   
   println("Optimal m/z for transmission @ RF amplitude: " + rfamplitude + " [V] and q=0.706 is " + moverz + " m/z"); 
   println("Optimal U for selecting this m/z is " + U + " [V]"); 
   println("Cutoff m/z at this RF amplitude is " + cutoff + " m/z");
   println("Getting optimal U for m/z " + moverz + " and RF amplitude " + rfamplitude + "V: " + getU(moverz, rfamplitude));
   println("Check: Vmax for m/z " + moverz + " is " + getV(moverz, 0.706) + " [V]");
   println("Check: Vcutoff for m/z " + moverz + " is " + getV(moverz, 0.905) + " [V]");
   println();
   flag = 8;
}

//-----------------------------------------------------------------------------------------------------------------------
// Helper functions for apex calculation
float getMoverZ(float q, float V) { 
   float k = 1.68e+8/1.661; // 1.68e-19/1.661e-27
   return k * 4 * V/(wwrr(RF_FREQUENCY, De/2) * q);
}

float getU(float moverz, float V) { 
   // Calculates U based on moverz and V
   float U = 0;
   float k = 1.68e+8/1.661; 
   float Vmax = getV(moverz, 0.706);
   float Vcutoff = getV(moverz, 0.905);

   if (V < Vmax) {
      U = 2 * k * 0.47546 * (V * V)/(moverz * wwrr(RF_FREQUENCY, De/2));   
   } else if (V >= Vmax && V <= Vcutoff) {
      U = -1.191 * V/2 + 1.0778 * moverz * wwrr(RF_FREQUENCY, De/2)/(8 * k);
   }
   return U;
}

float getLeftMoverZ(float U, float V) {
   float k = 1.68e+8/1.661;
   return 2 * k * 0.47546 * (V * V)/(U * wwrr(RF_FREQUENCY, De/2)); 
}

float getRightMoverZ(float U, float V) {
   float k = 1.68e+8/1.661;
   return (U + 1.191 * V/2)/(1.0778 * wwrr(RF_FREQUENCY, De/2)/(8 * k));
}

float getV(float moverz, float q) {
   float k = 1.68e+8/1.661;
   return q * moverz * wwrr(RF_FREQUENCY, De/2)/(4 * k);
}

float wwrr(float f, float r0) {
   float w = 2 * PI * f * 1000; // Angular frequency; f is in [kHz]
   float r = 0.001 * r0;       // Radius in [m]; r0 is in [mm]
   return w * w * r * r;
}


 
/*
println("//====================================  ION SUMMARY  ========================================");
      println("ion MW: " + MW +" Da; charge: " + Nq);
      //println("Initial position x : " + x +" y: " + y + " z: " +z  + " m");
      //println("Initial velocity Vx : " + vx +" m/z Vy: " + vy + " m/s Vz: " + vz  + " m/s");  
      println("Pressure of buffer gas is "+ p + " Torr");
      println("The mass of buffer gas is "+ mbuffer);
      println("The initial mean free path is "+ mfp + " [m]");
      println("The RF frequency: "+ freq + " [kHz]" + ", amplitude: " + Urf + " [V]");  
      println("Diameter of electrodes: " + De +  " [mm]");
      println("q: " + 4*q*Urf/(m*(0.001*De/2)*(0.001*De/2)*omega*omega) );  
      println("CCS: " + ccs );
      */
 
 
