
/**
 **************************************************************************************************
 *
 *                                   calculate a GEOMETRY
 *
 **************************************************************************************************
 */

// Method to calculate electrode geometry and return a list of electrode potentials
ArrayList<Potential> calculateElectrodeGeometry() { 
  
  startTimer(); // Start timing for performance measurement
  
  // List to store coordinates of electrodes for fast plotting
  ArrayList<Potential> electrodes = new ArrayList<Potential>();

  // Loop through layers along Z-axis (Nz layers)
  for (int k = 0; k < Nz; k++) {
    float z = k * step_z; // Calculate Z-coordinate for current layer
    if(k % scl == 0) println("calculated layer at " + z + "[mm]"); // Print progress every scl layers
    
    // Loop through X-coordinate (Nx grid points)
    for (int i = 0; i < Nx; i++) {
      float x = i * step_x;
      
      // Loop through Y-coordinate (Ny grid points)
      for (int j = 0; j < Ny; j++) {
        float y = j * step_y;
        
        // Define the potential at current (x, y, z) coordinate
        POT[i][j][k] = define_RF(x, y, z);
        
        // If the point corresponds to an electrode, add it to the electrodes list
        if (POT[i][j][k].electrode) {
          electrodes.add(new Potential(i, j, k, POT[i][j][k].potential, POT[i][j][k].electrode));
        }
      }
    }
  }

  stopTimer(); // Stop timing

  return electrodes; // Return the list of electrode potentials
}


//-------------------------------------------------------------------------------------------
//                               Defines RF potential for mass selective ejection
//                                       at a single quadrupole exit
//-------------------------------------------------------------------------------------------

// Method to define the RF potential at a given (x, y, z) coordinate
Potential define_RF(float x, float y, float z) {
  
  // Determine the slice number for the current coordinates
  int SLICE = mapLayer(x, y, z); 
  
  // Initialize the potential for the given coordinates (default potential = 0)
  Potential f = new Potential(x, y, z, 0, false);

  // Switch statement to define potential based on the slice location
  switch (SLICE) {

    // Inside a MultiQ (no potential)
    case 1:
      f = new Potential(x, y, z, 0, false); // No potential inside MultiQ
      break;
    
    // Quadrupole 1 region
    case 2:
      if (withinSqueezedQuadElectrode(x, y, z)) 
        f = new Potential(x, y, z, SIGN(x, y) * U, true); // Define potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;
      
    // Gap 1 region
    case 3:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, SIGN(x, y) * U, true); // Define potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB 1 region
    case 4:
      if (withinHole(x, y, DH/2))
        f = new Potential(x, y, z, 0, false); // No potential within hole (observation purpose)
      else 
        f = new Potential(x, y, z, U0, true); // Define potential outside hole
      break;

    // Gap 1 region (repeated)
    case 5:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, SIGN(x, y) * U, true); // Define potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // Quadrupole 1 region (repeated)
    case 6:
      if (withinSqueezedQuadElectrode(x, y, z)) 
        f = new Potential(x, y, z, SIGN(x, y) * U, true); // Define potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // Gap 3 region
    case 7:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, SIGN(x, y) * U, true); // Define potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB region (final)
    case 8:
      f = new Potential(x, y, z, U0, true); // Define potential for PCB region
      break;
  } // End switch statement

  return f; // Return the calculated potential
}
//--------------------------------------DC------------------------------------------

// Method to define the DC potential at a given (x, y, z) coordinate
Potential define_DC(float x, float y, float z) {
  
  // Determine the slice based on the (x, y, z) coordinates
  int SLICE = mapLayer(x, y, z); 
  
  // Initialize the potential at the given coordinates (default potential = 0)
  Potential f = new Potential(x, y, z, 0, false);

  // Switch statement to define the DC potential based on the slice location
  switch (SLICE) {

    // Inside a MultiQ (no potential)
    case 1:
      f = new Potential(x, y, z, 0, false); // No potential inside MultiQ
      break;
    
    // Quadrupole 1 region (positive DC potential)
    case 2:
      if (withinSqueezedQuadElectrode(x, y, z)) 
        f = new Potential(x, y, z, U0, true); // Define positive potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;
      
    // Gap 1 region (positive DC potential)
    case 3:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, U0, true); // Define positive potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB 1 region (either observation or defined potential)
    case 4:
      if (withinHole(x, y, DH/2)) 
        f = new Potential(x, y, z, 0, false); // No potential within hole (observation purpose)
      else 
        f = new Potential(x, y, z, U, true); // Define potential outside hole
      break;

    // Gap 2 region (negative DC potential)
    case 5:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, -U, true); // Define negative potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // Quadrupole 1 region (negative DC potential)
    case 6:
      if (withinSqueezedQuadElectrode(x, y, z)) 
        f = new Potential(x, y, z, -U, true); // Define negative potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;
    
    // Gap 3 region (negative DC potential)
    case 7:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, -U, true); // Define negative potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB region (negative potential)
    case 8:
      f = new Potential(x, y, z, -U, true); // Define negative potential for PCB region
      break;
  } // End switch statement

  return f; // Return the calculated DC potential
}


//-----------------------------------AC-----------------------------------

// Method to define the AC potential at a given (x, y, z) coordinate
Potential define_AC(float x, float y, float z) {

  // Determine the slice based on the (x, y, z) coordinates
  int SLICE = mapLayer(x, y, z); 
  
  // Initialize the potential at the given coordinates (default potential = 0)
  Potential f = new Potential(x, y, z, 0, false);

  // Switch statement to define the AC potential based on the slice location
  switch (SLICE) {

    // Inside a MultiQ (no potential)
    case 1:
      f = new Potential(x, y, z, 0, false); // No potential inside MultiQ
      break;
    
    // Inside a Quadrupole 1 region (AC potential with conditions)
    case 2:
      if (withinSqueezedQuadElectrode(x, y, z) && (x >= xDevice/2 && y >= yDevice/2)) 
        f = new Potential(x, y, z, +U, true); // Define positive AC potential for electrode
      else if (withinSqueezedQuadElectrode(x, y, z) && (x >= xDevice/2 && y <= yDevice/2)) 
        f = new Potential(x, y, z, +U, true); // Define positive AC potential for electrode
      else if (withinSqueezedQuadElectrode(x, y, z) && ((x <= xDevice/2 && (y <= yDevice/2) || y >= yDevice/2))) 
        f = new Potential(x, y, z, -U, true); // Define negative AC potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;
      
    // Gap 1 region (AC potential with conditions)
    case 3:
      if (withinQuadElectrode(x, y, 4.6) && (x >= xDevice/2 && y >= yDevice/2)) 
        f = new Potential(x, y, z, +U, true); // Define positive AC potential for electrode
      else if (withinQuadElectrode(x, y, 4.6) && (x >= xDevice/2 && y <= yDevice/2)) 
        f = new Potential(x, y, z, +U, true); // Define positive AC potential for electrode
      else if (withinQuadElectrode(x, y, 4.6) && ((x <= xDevice/2 && (y <= yDevice/2) || y >= yDevice/2))) 
        f = new Potential(x, y, z, -U, true); // Define negative AC potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB 1 region (either observation or defined AC potential)
    case 4:
      if (withinHole(x, y, DH/2)) 
        f = new Potential(x, y, z, 0, false); // No potential within hole (observation purpose)
      else 
        f = new Potential(x, y, z, U0, true); // Define potential outside hole
      break;

    // Gap 1 region (AC potential)
    case 5:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, U0, true); // Define AC potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // Quadrupole 1 region (AC potential)
    case 6:
      if (withinSqueezedQuadElectrode(x, y, z)) 
        f = new Potential(x, y, z, U0, true); // Define AC potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // Gap 3 region (AC potential)
    case 7:
      if (withinQuadElectrode(x, y, 4.6)) 
        f = new Potential(x, y, z, U0, true); // Define AC potential for electrode
      else 
        f = new Potential(x, y, z, 0, false); // No potential outside electrode
      break;

    // PCB region (AC potential)
    case 8:
      f = new Potential(x, y, z, U0, true); // Define AC potential for PCB region
      break;
  } // End switch statement

  return f; // Return the calculated AC potential
}


//***********************************************************************************
//
//***********************************************************************************
int mapLayer( float x, float y, float z) {
  
  int CASE = 0;
  float L1 = 0  + Qstart; // quads begin
  float L2 = L1 + Q1; // inside a quad
  float L3 = L2 + GAP;
  float L4 = L3 + PLATE; 
  float L5 = L4 + GAP;
  float L6 = L5 + Q2;
  float L7 = L6 + GAP;
  float L8 = L7 + PLATE;
  
  // ONION PRINCIPLE:
  // INNER LAYER: containing volume
  if ( 0<z && z<=L1) { return CASE = 1; }
  else if ( L1<z && z<=L2) { return CASE = 2; }
  else if ( L2<z && z<=L3) { return CASE = 3; }
  else if ( L3<z && z<=L4) { return CASE = 4; } 
  else if ( L4<z && z<=L5) { return CASE = 5; } 
  else if ( L5<z && z<=L6) { return CASE = 6; } 
  else if ( L6<z && z<=L7) { return CASE = 7; }
  else if ( L7<z && z<=L8) { return CASE = 8; }
  
  return CASE;
     
}


//***********************************************************************************
//                                 IS POINT NEAR EXIT?
//                                All dimensions in [mm]
//***********************************************************************************
// This function checks if a given point (X, Y, Z) is near the exit, within a specified 
// distance threshold (based on De, spacing, and CubeWidth).
boolean isNearExit(float X, float Y, float Z, float Xex, float Yex, float Zex) {
  // Calculate the radial distance in the XY-plane
  float r = sqrt((X - Xex) * (X - Xex) + (Y - Yex) * (Y - Yex)); 
  // Calculate the Z-axis distance
  float Zdist = sqrt((Z - Zex) * (Z - Zex));

  // Check if the point is within a defined distance range from the exit
  if (r < De + 1.415 * spacing && Zdist < CubeWidth + gap + De) {
    return true;  // Point is near the exit
  } else {
    return false; // Point is not near the exit
  }
}

//***********************************************************************************
// Auxiliary function that evaluates if a point (ion) is between electrodes within a "CUBE"
//***********************************************************************************
// This function checks if the point (X, Y, Z) is outside the defined electrode region, 
// specifically looking at the distance from the exit and whether it is in a cube-like 
// region defined by De and CubeWidth.
boolean isBetweenElectrodes(float X, float Y, float Z) {  
  // Calculate the radial distance to the exit in the XY-plane
  float r = sqrt((X - Xexit) * (X - Xexit) + (Y - Yexit) * (Y - Yexit)); 
  // Calculate the Z-axis distance to the exit
  float Zdist = sqrt((Z - Zexit) * (Z - Zexit));

  // Check if the point is outside the defined electrode region
  if (r < De && Zdist < CubeWidth + gap) {
    return false; // Point is inside the electrodes
  } else {
    return true;  // Point is outside the electrodes
  }
}

//---------------------------------------------------------------------------------------------------------------
// Function to check if a point is inside a quadrupole electrode region
//---------------------------------------------------------------------------------------------------------------
boolean withinQuadElectrode(float x, float y, float d) {
  boolean result = false;
  // Position of quadrupole electrode center along X and Y axis
  float center = spacing/2 + De/2;

  // Coordinates for the four quadrants of the quadrupole electrode
  float Xq1 = xDevice/2 - center;
  float Yq1 = yDevice/2 - center;
  
  float Xq2 = xDevice/2 + center;
  float Yq2 = yDevice/2 - center;
  
  float Xq3 = xDevice/2 - center;
  float Yq3 = yDevice/2 + center;
  
  float Xq4 = xDevice/2 + center;
  float Yq4 = yDevice/2 + center;

  // Check if the point (x, y) lies inside one of the four quadrants of the electrode
  if (x < xDevice/2 && y < yDevice/2) {
    if (sqrt((x - Xq1) * (x - Xq1) + (y - Yq1) * (y - Yq1)) < d/2) result = true;
  } else if (x > xDevice/2 && y < yDevice/2) {
    if (sqrt((x - Xq2) * (x - Xq2) + (y - Yq2) * (y - Yq2)) < d/2) result = true;
  } else if (x < xDevice/2 && y > yDevice/2) {
    if (sqrt((x - Xq3) * (x - Xq3) + (y - Yq3) * (y - Yq3)) < d/2) result = true;
  } else if (x > xDevice/2 && y > yDevice/2) {
    if (sqrt((x - Xq4) * (x - Xq4) + (y - Yq4) * (y - Yq4)) < d/2) result = true;
  }

  return result;
}

//----------------------------------------------------------------------------------
// Placeholder function that currently does not implement any specific checks
// for a screen electrode region.
//----------------------------------------------------------------------------------
boolean withinScreenElectrode(float x, float y, float d) {
  boolean result = false;
  return result;
}

//----------------------------------------------------------------------------------
// Function to return the sign based on position in the device, used for symmetry
// across quadrants.
//----------------------------------------------------------------------------------
int SIGN(float x, float y) {
  // Check in which quadrant the point lies, and return -1 or 1 accordingly
  if ((x <= xDevice/2 && y <= yDevice/2) || (x >= xDevice/2 && y >= yDevice/2)) {
    return -1;
  } else {
    return 1;
  }
}

//----------------------------------------------------------------------------------
// Function that calculates the distance between a point (x, y, z) and the centers 
// of the electrodes along different axes (X, Y, Z).
//----------------------------------------------------------------------------------
float[] distance(float x, float y, float z) { 
  float[] result = new float[3];
  
  // Calculate grid indices for the nearest electrode centers
  int i = ceil(abs(0.1 + x - CubeWidth)/(De + 2 * gap));
  if (i >= Ne) i = Ne;
  int j = ceil(abs(0.1 + y - CubeWidth)/(De + 2 * gap));
  if (j >= Ne) j = Ne;
  int k = ceil(abs(0.1 + z - CubeWidth)/(De + 2 * gap));
  if (k >= Ne) k = Ne;

  // Compute the centers of the quadrupoles in the grid
  float qX = CubeWidth + (gap + De/2) + (De + spacing) * (i - 1);
  float qY = CubeWidth + (gap + De/2) + (De + spacing) * (j - 1);  
  float qZ = CubeWidth + (gap + De/2) + (De + spacing) * (k - 1);
         
  // Calculate distances based on the specific layer along the Z-axis
  if (k == 1 || k == Ne) {
     result[0] = sqrt((x - qX) * (x - qX) + (y - qY) * (y - qY));
     result[1] = sqrt((z - qZ) * (z - qZ));
     result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
  }
      
  // Layer along the Y-axis
  if (1 < k && k < Ne && (j == 1 || j == Ne)) {
     result[0] = sqrt((x - qX) * (x - qX) + (z - qZ) * (z - qZ));
     result[1] = sqrt((y - qY) * (y - qY));
     result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
  }
         
  // Layer along the X-axis
  if (1 < k && k < Ne && 1 < j && j < Ne && (i == 1 || i == Ne)) {
      result[0] = sqrt((y - qY) * (y - qY) + (z - qZ) * (z - qZ));
      result[1] = sqrt((x - qX) * (x - qX));
      result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
   }
    
  return result;
}

//------------------------------------------------------------------------------------------------
// Function to check if a point is inside a hole with radius 'r' at the center of the device.
//------------------------------------------------------------------------------------------------
boolean withinHole(float x, float y, float r) { 
  float X = xDevice/2 - x;
  float Y = yDevice/2 - y;
  return ((X * X + Y * Y) < r * r);  // Point is inside the hole if within the radius
} 

//------------------------------------------------------------------------------------------------
// Function to check if a point is inside a wire region based on its (x, z) coordinates.
//------------------------------------------------------------------------------------------------
boolean insideWire(float x, float z) {
  float d = 0.25;  // Wire radius
  float Z = Qstart + Q1 + GAP + PLATE - 1 - z;  // Z position of the wire
  float X = xDevice/2 - x;
  return ((Z * Z + X * X) < d * d/4);  // Check if point is within the wire region
}

//------------------------------------------------------------------------------------------------
// Function to check if a point is inside a square-shaped region based on the device center.
//------------------------------------------------------------------------------------------------
boolean withinSquare(float x, float y, float d) {
   float X = xDevice/2 - x;
   float Y = yDevice/2 - y;
   float quadCenter = spacing/2 + De/2;
   float b = 2 * quadCenter - 1.41 * d/2;  // Compute square's boundary
   return ((X - b) < Y && (-X - b) < Y && Y < (X + b) && Y < (-X + b));  // Check if point is inside square
}

//------------------------------------------------------------------------------------------------
// Function to check if a point lies within a circular region.
//------------------------------------------------------------------------------------------------
boolean withinCircle(float x, float y, float d) {
   float X = xDevice/2 - x;
   float Y = yDevice/2 - y;
   return ((X * X + Y * Y) < d * d/4);  // Check if point is within the circle
}

//-----------------------------------------------------------------------------------------------
// Function to check if a point lies within extra rods, given specific parameters.
//-----------------------------------------------------------------------------------------------
boolean withinExtraRods(float x, float y, float z, float p, float d, float l) {
  float Xr = xDevice/2 - x;
  float Yr1 = yDevice/2 + p - y;
  float Yr2 = yDevice/2 - p - y; 
  return (((Xr * Xr + Yr1 * Yr1) < d * d/4 || (Xr * Xr + Yr2 * Yr2) < d * d/4)) && z < (Qstart + l) && z > (Qstart); 
}

//-----------------------------------------------------------------------------------------------
// Function to check if a point is inside a squeezed quadrupole electrode, considering Z-axis variation.
//-----------------------------------------------------------------------------------------------
boolean withinSqueezedQuadElectrode(float x, float y, float z) {
  boolean result = false;
  
  // Squeezing factor for quadrupole electrode
  float alfa = 0.02;
  float d = De;
  
  // Position of quadrupole axis
  float center = spacing/2 + De/2;
  
  // Calculate the squeeze factor along Z-axis based on current position
  float sq = alfa * (Qstart + Q1 - z);  // Adjust electrode position based on Z
  
  // Coordinates for the four quadrants of the squeezed quadrupole
  float Xq1 = xDevice/2 - center;
  float Yq1 = yDevice/2 - center + sq;
  
  float Xq2 = xDevice/2 + center;
  float Yq2 = yDevice/2 - center + sq;
  
  float Xq3 = xDevice/2 - center;
  float Yq3 = yDevice/2 + center - sq;
  
  float Xq4 = xDevice/2 + center;
  float Yq4 = yDevice/2 + center - sq;
  
  // Check if the point (x, y) is inside one of the squeezed quadrants
  if (x < xDevice/2 && y < yDevice/2) {
    if (sqrt((x - Xq1) * (x - Xq1) + (y - Yq1) * (y - Yq1)) < d/2) result = true;
  } else if (x > xDevice/2 && y < yDevice/2) {
    if (sqrt((x - Xq2) * (x - Xq2) + (y - Yq2) * (y - Yq2)) < d/2) result = true;
  } else if (x < xDevice/2 && y > yDevice/2) {
    if (sqrt((x - Xq3) * (x - Xq3) + (y - Yq3) * (y - Yq3)) < d/2) result = true;
  } else if (x > xDevice/2 && y > yDevice/2) {
    if (sqrt((x - Xq4) * (x - Xq4) + (y - Yq4) * (y - Yq4)) < d/2) result = true;
  }

  return result;  // Return if the point lies within the squeezed quadrupole
}
