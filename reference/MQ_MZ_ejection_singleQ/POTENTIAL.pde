/**
 **************************************************************************************************
 *
 *                                Defines the potential at a given x, y, z point 
 *
 **************************************************************************************************
 */
class Potential {
  
  public float x, y, z;         // Coordinates of the point
  public float potential;       // Electric potential at the point
  public boolean electrode;     // Flag indicating if the point is an electrode
  
  // Constructor to initialize the point with its coordinates, potential, and electrode status
  Potential(float mx, float my, float mz, float mpot, boolean el) {
    x = mx; 
    y = my; 
    z = mz;
    potential = mpot;
    electrode = el;
  }
}

/**
 *************************************************************************************************************************
 *
 *                                                   Calculates the potential across the grid
 *
 *************************************************************************************************************************
 */
Potential [][][] calculatePotential(float precesion, Potential [][][] pot) {
   
  startTimer();  // Start timing the calculation
   
  // Convert check point coordinates to grid indices
  int i_ch = convertMillimetersToArray(x_check);
  int j_ch = convertMillimetersToArray(y_check);
  int k_ch = convertMillimetersToArray(z_check);
  
  checkPointIsElectrode(i_ch, j_ch, k_ch);  // Check if the point is an electrode
  
  float dif = 100;  // Initial difference in potential
  float prevpot = 1;  // Previous potential value
  
  println("Will calculate the potential until the difference between iterations is less than " + precesion + "%");
  
  int Nit = 0;  // Iteration counter
  int Niteration = N_ITERATIONS;  // Maximum number of iterations

  // Iterate to calculate the potential
  while (Nit < Niteration && dif > precesion) {  
     for (int k = 1; k < Nz - 1; k++) {  // Iterate over the z-axis
       for (int i = 1; i < Nx - 1; i++) {  // Iterate over the x-axis
         for (int j = 1; j < Ny - 1; j++) {  // Iterate over the y-axis
           if (!pot[i][j][k].electrode) {  // Skip electrodes as their potential is fixed
             // Average the potential of the neighboring points to update the potential at (i, j, k)
             pot[i][j][k].potential = 
               (pot[i-1][j][k].potential + 
               pot[i+1][j][k].potential + 
               pot[i][j-1][k].potential + 
               pot[i][j+1][k].potential + 
               pot[i][j][k-1].potential + 
               pot[i][j][k+1].potential)/6; 
           }        
         }
       }
     }
     
     // Calculate the difference in potential at the check point
     dif = abs(100 * ((prevpot + 0.0000000001) - pot[i_ch][j_ch][k_ch].potential)/(prevpot + 0.0000000001));  
     prevpot = pot[i_ch][j_ch][k_ch].potential;  // Update the previous potential
     Nit++;  // Increment iteration count
     
     // Print progress every 10 iterations
     if (Nit % 10 == 0) {
       println("After " + Nit + " iterations, the difference is " + dif + "%");
     }
  }  // Close while loop

  // Print the final results after completion
  println("The final number of iterations was " + Nit + " to reach the difference of " + dif + "%");
  
  return pot;  // Return the updated potential array
}

/**
 *************************************************************************************************************************
 *
 *                                      Calculates the potential gradient at a given point (x, y, z)
 *
 *************************************************************************************************************************
 */
float [] calculatePotentialDifference(double x, double y, double z, Potential [][][] P) { // z, x, y 

   float [] VOLT = new float[3];  // Array to store the gradient in x, y, and z directions
   
   // Convert x, y, z to grid indices
   int a = convertMetersToArray(x, Nx);  // Convert x to grid index
   if (a <= 2) { a = 2; }  // Ensure index is within bounds
   if (a >= Nx - 2) { a = Nx - 2; }

   int b = convertMetersToArray(y, Ny);  // Convert y to grid index
   if (b <= 2) { b = 2; }  // Ensure index is within bounds
   if (b >= Ny - 2) { b = Ny - 2; }

   int c = convertMetersToArray(z, Nz);  // Convert z to grid index
   if (c <= 2) { c = 2; }  // Ensure index is within bounds
   if (c >= Nz - 2) { c = Nz - 2; }
   
   // Calculate the finite difference of potential in each direction (x, y, z)
   float dUx = P[a+1][b][c].potential - P[a-1][b][c].potential; 
   float dUy = P[a][b+1][c].potential - P[a][b-1][c].potential;
   float dUz = P[a][b][c+1].potential - P[a][b][c-1].potential;
   
   // Store the gradients in the VOLT array
   VOLT[0] = -dUx/2;  // Gradient in the x-direction
   VOLT[1] = -dUy/2;  // Gradient in the y-direction
   VOLT[2] = -dUz/2;  // Gradient in the z-direction
   
   return VOLT;  // Return the potential gradients
}

/**
 *************************************************************************************************************************
 *
 *                                      Checks if a given point (x, y, z) is an electrode
 *
 *************************************************************************************************************************
 */
boolean pointIsElectrode(double x, double y, double z) {
  
  // Convert the coordinates to grid indices
  int i = convertMetersToArray(x, Nx);
  int j = convertMetersToArray(y, Ny);
  int k = convertMetersToArray(z, Nz);
  
  // Check if the point is an electrode by checking both RF and AC potential grids
  boolean result = RFPOT[i][j][k].electrode || ACPOT[i][j][k].electrode;
  
  return result;  // Return true if the point is an electrode, false otherwise
}

/**
 *************************************************************************************************************************
 *
 *                                      Prints if a given point (i, j, k) is an electrode or not
 *
 *************************************************************************************************************************
 */
void checkPointIsElectrode(int i, int j, int k) {
  
  // Check if the point is an electrode in the POT array
  boolean result = POT[i][j][k].electrode;
  
  // Print the status of the check point
  if (result) {
     println("The check point IS an ELECTRODE!!!");
  } else {
     println("The check point IS NOT an ELECTRODE!!!");
  }
}

/**
 *************************************************************************************************************************
 *
 *                                     Performs an initial rough estimation of potential in the grid
 *
 *************************************************************************************************************************
 */
void feelRough(Potential [][][] pot) {
     for (int k = 1; k < Nz - 1; k = k + 2) {  // Skip every other point in the z-direction
       for (int i = 1; i < Nx - 1; i = i + 2) {  // Skip every other point in the x-direction
         for (int j = 1; j < Ny - 1; j = j + 2) {  // Skip every other point in the y-direction
           if (!pot[i][j][k].electrode) {  // Skip electrodes
             // Perform a rough update of the potential by averaging the neighboring points
             pot[i][j][k].potential = 
             (pot[i-1][j][k].potential + 
              pot[i+1][j][k].potential + 
              pot[i][j-1][k].potential + 
              pot[i][j+1][k].potential + 
              pot[i][j][k-1].potential + 
              pot[i][j][k+1].potential)/6.0;         
           }        
         }
       }
     }
}
