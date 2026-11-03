/**
 **************************************************************************************************
 *                                   
 * Class representing the potential at a given point in 3D space.
 * Includes information about the coordinates, potential value, 
 * and whether the point is part of an electrode.
 *
 **************************************************************************************************
 */
class Potential {
  
  public float x, y, z;       // Coordinates of the point (x, y, z)
  public float potential;     // Potential value at the point
  public boolean electrode;   // Indicates if the point is an electrode
  
  // Constructor for initializing a Potential object
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
 * Iteratively calculates the potential in a 3D space until the difference between 
 * iterations falls below a specified precision.
 *
 * @param precision The stopping criterion based on the percent difference between iterations.
 * @param pot       The 3D array of Potential objects representing the potential field.
 * @return          The updated 3D array of Potential objects.
 *
 *************************************************************************************************************************
 */
Potential[][][] calculatePotential(float precision, Potential[][][] pot) {
  
  startTimer(); // Start timing the calculation
  
  // Convert the given checkpoint coordinates to array indices
  int i_ch = convertMillimetersToArray(x_check);
  int j_ch = convertMillimetersToArray(y_check);
  int k_ch = convertMillimetersToArray(z_check);
  
  // Ensure the checkpoint is not part of an electrode
  checkCheckPointNotElectrode(i_ch, j_ch, k_ch);
  
  float dif = 100;           // Initialize difference percentage
  float prevpot = 1;         // Previous potential at the checkpoint
  
  println("Starting potential calculation until the difference is less than " + precision + "%.");
  
  int Nit = 0;               // Iteration counter
  int Niteration = N_ITERATIONS; // Maximum number of iterations allowed

  // Iterate until precision is achieved or maximum iterations are reached
  while (Nit < Niteration && dif > precision) {  
    for (int k = 1; k < Nz - 1; k++) {
      for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
          if (!pot[i][j][k].electrode) {
            // Update potential using the average of neighboring points
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
    
    // Calculate the percentage difference at the checkpoint
    dif = Math.abs(100 * ((prevpot + 1e-10) - pot[i_ch][j_ch][k_ch].potential)/(prevpot + 1e-10));  
    prevpot = pot[i_ch][j_ch][k_ch].potential;
    Nit++;

    // Log progress every 10 iterations
    if (Nit % 10 == 0) {
      println("After " + Nit + " iterations, the difference is " + dif + "%.");
    }
  }

  println("Final iteration count: " + Nit + " with a difference of " + dif + "%.");
  
  return pot;
}

/**
 *************************************************************************************************************************
 * 
 * Calculates the potential difference at a given point in 3D space.
 *
 * @param x The x-coordinate in meters.
 * @param y The y-coordinate in meters.
 * @param z The z-coordinate in meters.
 * @param P The 3D array of Potential objects.
 * @return  An array containing the potential differences along the x, y, and z axes.
 *
 *************************************************************************************************************************
 */
float[] calculatePotentialDifference(double x, double y, double z, Potential[][][] P) {
  
  float[] VOLT = new float[3]; // Array to store potential differences
  
  // Convert x-coordinate to array index and ensure it's within bounds
  int a = convertMetersToArray(x);
  if (a <= 2) a = 2;
  if (a >= Nx - 2) a = Nx - 2;
  
  // Convert y-coordinate to array index and ensure it's within bounds
  int b = convertMetersToArray(y);
  if (b <= 2) b = 2;
  if (b >= Ny - 2) b = Ny - 2;
  
  // Convert z-coordinate to array index and ensure it's within bounds
  int c = convertMetersToArray(z);
  if (c <= 2) c = 2;
  if (c >= Nz - 2) c = Nz - 2;

  // Calculate potential differences along each axis
  float dUx = P[a+1][b][c].potential - P[a-1][b][c].potential; 
  float dUy = P[a][b+1][c].potential - P[a][b-1][c].potential;
  float dUz = P[a][b][c+1].potential - P[a][b][c-1].potential;
  
  VOLT[0] = -dUx/2; 
  VOLT[1] = -dUy/2; 
  VOLT[2] = -dUz/2; 
  
  return VOLT;
}

/**
 *************************************************************************************************************************
 * 
 * Checks if a given checkpoint in the 3D potential array is part of an electrode.
 * Logs the result to the console.
 *
 * @param i The x-index of the checkpoint.
 * @param j The y-index of the checkpoint.
 * @param k The z-index of the checkpoint.
 *
 *************************************************************************************************************************
 */
void checkCheckPointNotElectrode(int i, int j, int k) { 
  if (POT[i][j][k].electrode) {
    println("The checkpoint IS AN ELECTRODE!");
  } else {
    println("The checkpoint IS NOT AN ELECTRODE.");
  }
}
