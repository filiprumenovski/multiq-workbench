


//---------------------------------------------------------------------------------------------------------------------------------------
int convertMetersToArray(double coordinate, int Nmax){
   int Npoint = (int)(coordinate*scl*1000.0);
   if (Npoint<=2) { Npoint = 2; }
   if (Npoint>=Nmax-2) {Npoint = Nmax-2; }
   return Npoint;
}
 
int convertMillimetersToArray(float c){
   int Npoint = (int)(c*scl);
   return Npoint;
}

 //------------------------------------------------------------------------------------------------------------------------
 // in [mm]
float[] translateToCenterCoordinateSystem(float x, float y, float z) { // all are [mm]!!!!
  float [] result = new float[3];
  result[0] = x - xDevice/2;
  result[1] = y - yDevice/2;
  result[2] = z - zDevice/2;
  return result; 
}

float[] translateFromCenterCoordinateSystem(float x, float y, float z) { // all are [mm]!!!!
  float [] result = new float[3];
  result[0] = (xDevice/2 + x);
  result[1] = (yDevice/2 + y);
  result[2] = (zDevice/2 + z);
  return result; 
}
 
 
 
 
 //-----------------------------------------------------------------------------------------------------------------------
 
 // start timer
 public void startTimer() {
   Sstart = second();  // Values from 0 - 59
   Mstart = minute();  // Values from 0 - 59
   Hstart = hour();
   Tstart = millis();
   println( "The task started at "+Hstart+":"+Mstart+":"+Sstart); 
   println("Dimension: " + Nx + "x" + Ny + "x" + Nz);
 }
 //end timer
 public void stopTimer() {
   Send = second();  // Values from 0 - 59
   Mend = minute();  // Values from 0 - 59
   Hend = hour();
   int se = second();  // Values from 0 - 59
   int me = minute();  // Values from 0 - 59
   int he = hour();
   Tend = millis();
   
   println( "The task ended at "+Hend+":"+Mend+":"+Send);
   println( "It took "+((Tend-Tstart)/3600000)+":"+((Tend-Tstart)/60000)+":"+((Tend-Tstart)/1000));
 }
 
 
 
 
 
 
