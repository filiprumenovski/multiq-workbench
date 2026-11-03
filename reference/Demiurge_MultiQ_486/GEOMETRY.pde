/**
 **************************************************************************************************
 *
 *                                   Calculate Electrode Geometry
 * to calculate different geometry replace function defining particular geometry, 
 * e.g. definePotential_RF(x, y, z) in POT[i][j][k] = definePotential_RF(x, y, z);
 *
 **************************************************************************************************
 */

ArrayList<Potential> calculateElectrodeGeometry() { // (PrintWriter output_potential) {
  
  startTimer();
  
  // Create an ArrayList to store the coordinates of electrodes for faster plotting
  ArrayList<Potential> electrodes = new ArrayList<Potential>();

  for (int k = 0; k < Nz; k++) {
    float z = k * step_z; // step_z = 1/scl
    if (k % scl == 0) println("Calculated layer at " + z + "[mm]");
    
    for (int i = 0; i < Nx; i++) {
      float x = i * step_x;
      for (int j = 0; j < Ny; j++) {
        float y = j * step_y;
        POT[i][j][k] = definePotential_RF(x, y, z);
        
        // Add electrode to the list if it is an electrode
        if (POT[i][j][k].electrode) {
          electrodes.add(new Potential(i, j, k, POT[i][j][k].potential, POT[i][j][k].electrode)); 
        }
      }
    }
  }

  stopTimer();
  
  return electrodes;
}



/**
 **************************************************************************************************
 **************************  Calculate the Structure of the Potential   **************************
 **************************************************************************************************
 */

// "Main" RF potential on all electrodes
Potential definePotential_RF(float x, float y, float z) {
  // Coordinates x and y start from the top-left corner
    
  int POSITION = mapLayer(x, y, z); 
  Potential f = null;
  float[] restrains = new float[3];
  
  switch (POSITION) {
    
    // Inside layer
    case 1:
      f = new Potential(x, y, z, 0, false);
      break;
    
    // Inside QUAD layer
    case 2:  
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height
      if (restrains[0] <= De/2 && restrains[1] <= De/2) 
        f = new Potential(x, y, z, U * restrains[2], true);  // restrains[2] indicates the sign (+1 or -1)
      else 
        f = new Potential(x, y, z, 0, false);     
      break;
    
    // Inside GAP layer
    case 3:
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height, restrains[2] is the sign
      if (restrains[0] <= stem && restrains[1] <= gap + De/2) 
        f = new Potential(x, y, z, U * restrains[2], true); 
      else 
        f = new Potential(x, y, z, 0, false);   
      break;
    
    // Inside SUPPORTING PLATE layer
    case 4:
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height
      if (restrains[0] <= stem && restrains[1] <= gap + De/2 + CubeWidth)      
        f = new Potential(x, y, z, U * restrains[2], true); 
      else if (restrains[0] <= De/2 && restrains[1] <= CubeWidth + gap + De/2) 
        f = new Potential(x, y, z, 0, false);
      else if (restrains[0] <= De * 2 && restrains[1] <= CubeWidth + gap + De/2) 
        f = new Potential(x, y, z, U0, true);
      else 
        f = new Potential(x, y, z, 0, false);
      break;
  }
  return f; 
}

/**
 **************************************************************************************************
 *******************************************   DC Potential   *************************************
 **************************************************************************************************
 */

// "Main" DC potential on all electrodes
Potential definePotential_DC(float x, float y, float z) {
  // Coordinates x and y start from the top-left corner
    
  int POSITION = mapLayer(x, y, z); 
  Potential f = null;
  float[] restrains = new float[3];
  
  switch (POSITION) {
    
    // Inside layer
    case 1:
      f = new Potential(x, y, z, 0, false);
      break;
    
    // Inside QUAD layer
    case 2:  
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height
      if (restrains[0] <= De/2 && restrains[1] <= De/2) 
        f = new Potential(x, y, z, U0, true);  // restrains[2] indicates the sign (+1 or -1)
      else 
        f = new Potential(x, y, z, 0, false);     
      break;
    
    // Inside GAP layer
    case 3:
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height
      if (restrains[0] <= stem && restrains[1] <= gap + De/2) 
        f = new Potential(x, y, z, U0, true); 
      else 
        f = new Potential(x, y, z, 0, false);   
      break;
    
    // Inside SUPPORTING PLATE layer
    case 4:
      restrains = quadrupoles(x, y, z);           
      // restrains[0] is the radius, restrains[1] is the height
      if (restrains[0] <= stem && restrains[1] <= gap + De/2 + CubeWidth)      
        f = new Potential(x, y, z, U0, true); 
      else if (restrains[0] <= De/2 && restrains[1] <= CubeWidth + gap + De/2) 
        f = new Potential(x, y, z, 0, false);
      else if (restrains[0] <= De * 2 && restrains[1] <= CubeWidth + gap + De/2) 
        f = new Potential(x, y, z, U, true);
      else 
        f = new Potential(x, y, z, 0, false);
      break;
  }
  return f; 
}

//***********************************************************************************
//                        MAPPING OF ELECTRODE LAYERS           
//***********************************************************************************
int mapLayer(float x, float y, float z) {
 
  // 1: "Inside"
  // 2: "Inside QUAD Layer"
  // 3: "Inside GAP Layer"
  // 4: "Inside SUPPORTING PLATE Layer"
  int CASE = 0;
  float L1 = CubeWidth + gap + De; // Width of quad supporting plate (PCB) + gap + height of the electrode (= diameter of the electrode)
  float L2 = CubeWidth + gap;
  float L3 = CubeWidth;
  
  // ONION PRINCIPLE:
  // INNER LAYER: containing volume
  if ((L1 <= x && x <= xDevice - L1) && (L1 <= y && y <= yDevice - L1) && (L1 <= z && z <= zDevice - L1)) { 
    return CASE = 1; 
  }
  
  // QUAD LAYER
  else if ((L2 <= x && x <= xDevice - L2) && (L2 <= y && y <= yDevice - L2) && (L2 <= z && z <= zDevice - L2)) { 
    return CASE = 2; 
  }
  
  // GAP LAYER
  else if ((L3 <= x && x <= xDevice - L3) && (L3 <= y && y <= yDevice - L3) && (L3 <= z && z <= zDevice - L3)) { 
    return CASE = 3; 
  }
  
  // SUPPORTING PLATE LAYER
  else if ((0 <= x && x <= xDevice) && (0 <= y && y <= yDevice) && (0 <= z && z <= zDevice)) { 
    return CASE = 4; 
  } 
  
  return CASE;
}

//***********************************************************************************
//                                 IS POINT NEAR EXIT?
//                                 All dimensions in [mm]
//***********************************************************************************
boolean isNearExit(float X, float Y, float Z, float Xex, float Yex, float Zex) {
  float r = sqrt((X - Xex) * (X - Xex) + (Y - Yex) * (Y - Yex)); 
  float Zdist = sqrt((Z - Zex) * (Z - Zex));
  
  // If the point is within a defined radius and Z distance from the exit
  if (r < De + 1.415 * spacing && Zdist < CubeWidth + gap + De) {
    return true; 
  } else {
    return false;
  }
}

//***********************************************************************************
// Auxiliary function to evaluate if a point (ion) is near the exit within a "CUBE"
//***********************************************************************************
boolean isBetweenElectrodes(float X, float Y, float Z) {  
  float r = sqrt((X - Xexit) * (X - Xexit) + (Y - Yexit) * (Y - Yexit)); 
  float Zdist = sqrt((Z - Zexit) * (Z - Zexit));
  
  // Check if the point is outside the electrode's vicinity
  if (r < De && Zdist < CubeWidth + gap) {
    return false; 
  } else {
    return true;
  }
}

//***********************************************************************************
// Function to sort the point into one of four layers
//***********************************************************************************
int SortPoint(float x, float y, float z) {
  // Layers defined by the following cases:
  // 1 - Inside
  // 2 - Inside QUAD Layer
  // 3 - Inside GAP Layer
  // 4 - Inside SUPPORTING PLATE Layer
  
  int CASE = 10;
  float L1 = CubeWidth + gap + De;
  float L2 = CubeWidth + gap;
  float L3 = CubeWidth;

  // ONION PRINCIPLE:
  // 1. INNER LAYER: containing volume
  if ((L1 <= x && x <= xDevice - L1) && (L1 <= y && y <= yDevice - L1) && (L1 <= z && z <= zDevice - L1)) {
    return CASE = 1;
  }
  
  // 2. QUAD LAYER
  else if ((L2 <= x && x <= xDevice - L2) && (L2 <= y && y <= yDevice - L2) && (L2 <= z && z <= zDevice - L2)) {
    return CASE = 2;
  }
  
  // 3. GAP LAYER
  else if ((L3 <= x && x <= xDevice - L3) && (L3 <= y && y <= yDevice - L3) && (L3 <= z && z <= zDevice - L3)) {
    return CASE = 3;
  }
  
  // 4. SUPPORTING PLATE LAYER
  else if ((0 <= x && x <= xDevice) && (0 <= y && y <= yDevice) && (0 <= z && z <= zDevice)) {
    return CASE = 4; 
  }

  return CASE;
}

//***********************************************************************************
// Function to calculate the quadrupole coordinates and their constraints
//***********************************************************************************
float[] quadrupoles(float x, float y, float z) { // All units in [mm]
  float[] result = new float[3];

  // Calculate the indices based on the coordinates (adjusted for cube width and electrode spacing)
  int i = ceil(abs(0.1 + x - CubeWidth)/(De + 2 * gap));
  if (i >= Ne) i = Ne;
  int j = ceil(abs(0.1 + y - CubeWidth)/(De + 2 * gap));
  if (j >= Ne) j = Ne;
  int k = ceil(abs(0.1 + z - CubeWidth)/(De + 2 * gap));
  if (k >= Ne) k = Ne;

  // Calculate the centers of the quadrupoles
  float qX = CubeWidth + (gap + De/2) + (De + spacing) * (i - 1); 
  float qY = CubeWidth + (gap + De/2) + (De + spacing) * (j - 1);  
  float qZ = CubeWidth + (gap + De/2) + (De + spacing) * (k - 1);

  // Determine the layer along the Z-axis
  if (k == 1 || k == Ne) {
    result[0] = sqrt((x - qX) * (x - qX) + (y - qY) * (y - qY));
    result[1] = sqrt((z - qZ) * (z - qZ));
    result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
  }

  // Determine the layer along the Y-axis
  if (1 < k && k < Ne && (j == 1 || j == Ne)) {
    result[0] = sqrt((x - qX) * (x - qX) + (z - qZ) * (z - qZ));
    result[1] = sqrt((y - qY) * (y - qY));
    result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
  }

  // Determine the layer along the X-axis
  if (1 < k && k < Ne && 1 < j && j < Ne && (i == 1 || i == Ne)) {
    result[0] = sqrt((y - qY) * (y - qY) + (z - qZ) * (z - qZ));
    result[1] = sqrt((x - qX) * (x - qX));
    result[2] = pow(-1, i) * pow(-1, j) * pow(-1, k);
  }

  return result;
}

//------------------------------------------------------------------------------------
// Check if an ion is inside an electrode
//***********************************************************************************
boolean inElectrode(Ion ion) { // In mm
  boolean result = false;
  
  float[] restrains = new float[3];
  float ion_x = (float) ion.x * 1000; 
  float ion_y = (float) ion.y * 1000; 
  float ion_z = (float) ion.z * 1000; 
  
  // Determine the position of the ion
  int POSITION = SortPoint(ion_x, ion_y, ion_z);
  
  switch (POSITION) {
    // "Inside QUAD Layer"
    case 2:  
      restrains = quadrupoles(ion_x, ion_y, ion_z);           
      // Check if the ion is within the bounds of the quadrupole (radius and height)
      if (restrains[0] <= De/2 && restrains[1] <= De/2) {
        result = true;
      } else {
        result = false;
      }
      break;
  }
  return result;
}
