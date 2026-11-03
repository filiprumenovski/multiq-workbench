// =======================================   DRAW  ===============================================
void draw() {
  if (flag == 0) return;

  // ===================================    Draw RF Geometry   ==================================== 
  if (flag == 1) { // Draw 3D geometry of RF field
    background(255);
    
    translate(width/2, height/2, -width); // Move to center of canvas
    rotateX(a1);
    rotateY(a2);
    rotateZ(a3);
    
    scale(DISPLAY_SCALE*2); // Scaling factor for display, defined at the start of the program!
    
    stroke(#dddddd); // Light gray for coordinate axes
    strokeWeight(0.5);
    fill(0);
    line(-100.0, 0, 0, 300, 0, 0); // X axis
    fill(0, 102, 153, 204);
    text("Z", -100, 0, 0); // Label Z-axis
    line(0, -100, 0, 0, 1000, 0); // Y axis
    text("Y", 0, -100, 0); // Label Y-axis
    line(0, 0, -100, 0, 0, 1000); // Z axis
    text("X", 0, 0, -100); // Label X-axis

    // Draw electrodes with different potential values
    if (Electrodes != null) {
      for (Potential el: Electrodes) {
        strokeWeight(4);
        if (el.potential > 0) {
          stroke(#FF0000, 100); // Use a distinguishable red (Tomato) for positive potential
        }
        if (el.potential < 0) {
          stroke(#0000FF, 100); // Use a distinguishable blue (SteelBlue) for negative potential
        }
        if (el.potential == 0) {
          stroke(#BEBEBE, 100); // Gray for zero potential
        }
        point(el.z, el.x, el.y); // Draw electrode point
        strokeWeight(1);
      }
    } else {
      text("No geometry is defined...", -50, 0, 0); // If no electrodes are defined
    }
  }

  // ====================   Draw a cross-section of the potential   ============================= 
  if (flag == 2) { 
    background(255);

    int offsetY = Ny;
    int offsetX = Nx;
    float delta = 0.03;
    float step = 0.2;
    scale(1);
    strokeWeight(2);
    
    println("Drawing equipotentials at the z: " + zcs);
    
    // Draw equipotentials
    for (int i = 1; i < Nx; i++) {
      for (int j = 1; j < Ny; j++) {
        for (float p = -1 + 0.1; p <= 1 - 0.1; p = p + step) { // Loop through equipotentials
          if (p - delta < RFPOT[i][j][zcs].potential && RFPOT[i][j][zcs].potential < p + delta) {
            stroke(#8B4513); // Use a distinguishable brown color for equipotentials
            point((i * 2 + w/2 - offsetX), (2 * j + h/2 - offsetY)); // Draw equipotential points
          }
          if (RFPOT[i][j][zcs].electrode) {
            // Color electrodes based on potential values
            if (RFPOT[i][j][zcs].potential > 0) stroke(#FF6347); // Red for positive potential
            if (RFPOT[i][j][zcs].potential < 0) stroke(#4682B4); // Blue for negative potential
            if (RFPOT[i][j][zcs].potential == 0) stroke(#BEBEBE); // Gray for zero potential
            point((i * 2 + w/2 - offsetX), (j * 2 + h/2 - offsetY)); // Draw electrode points
          }
        }
      }
    }

  } // Close if

  // ==================   Draw 3D Trajectories  =================================================
  if (flag == 5) { 
    background(255);
    
    float TRANSPARENT = 10;
    
    translate(width/2, height/2, -width);
    
    rotateX(a1);
    rotateY(a2);
    rotateZ(a3);
    
    scale(DISPLAY_SCALE*2);
    
    stroke(#dddddd); // Light gray for coordinate axes
    strokeWeight(1);
    fill(0);
    line(-100.0, 0, 0, 300, 0, 0); // X axis
    text("X", -100, 0, 0); // Label X-axis
    line(0, -100, 0, 0, 1000, 0); // Y axis
    text("Y", 0, -100, 0); // Label Y-axis
    line(0, 0, -100, 0, 0, 1000); // Z axis
    text("Z", 0, 0, -100); // Label Z-axis

    if (ionTrajectory != null) {
      for (Trajectory coordinates: ionTrajectory) {
        strokeWeight(0.5);
        stroke(0);  
        point(coordinates.x, coordinates.y, coordinates.z); // Draw trajectory points
        strokeWeight(0.5);
      }
    }

    // Draw electrodes
    if (Electrodes != null) {
      for (Potential el: Electrodes) {
        strokeWeight(1);
        stroke(#BEBEBE, TRANSPARENT); // Light gray for electrodes
        point(el.x, el.y, el.z); // Draw electrode points
        strokeWeight(1);
      }
    } else {
      text("No geometry is defined...", -50, 0, 0); // If no geometry is defined
    }

  }

  // =======================================  Draw 2D (Plane) Trajectories  ==============================
  if (flag == 6) {
    background(255);

    fill(0, 100, 150, 50); // Semi-transparent fill for background elements
    int SCALE = DISPLAY_SCALE;
    float UpperLeftZ = w/2 - zDevice/2 * scl * SCALE;
    float UpperLeftY = h/2 - yDevice/2 * scl * SCALE;

    fill(0);

    // Draw electrodes
    if (Electrodes != null) {
      for (Potential el: Electrodes) {
        strokeWeight(1);       
        stroke(#555555); // Use gray for electrodes in 2D view
        point(UpperLeftZ + el.z * SCALE, UpperLeftY + el.x * SCALE); // Draw electrode points in 2D view
        strokeWeight(1);
      }
    } else {
      text("No geometry is defined...", -50, 0, 0); // If no geometry is defined
    }

    // Draw trajectories (projection on ZY plane)
    if (ionTrajectory != null) {
      for (Trajectory coordinates: ionTrajectory) {
        strokeWeight(1);
        
        // Group 1 ions (blue)
        if (coordinates.group == 1) stroke(#0000FF); // Blue for group 1
        // Group > 1 ions (red)
        if (coordinates.group > 1) stroke(#FF0000); // Red for groups greater than 1
        point(UpperLeftZ + coordinates.z * SCALE, UpperLeftY + coordinates.x * SCALE); // Draw trajectory points
        strokeWeight(1);
        stroke(0);
      }
    }

    // Draw X-Y Trajectories
    float graphOffsetX = 1000;
    float graphOffsetY = 240;
    // Draw electrodes in X-Y plane
    if (Electrodes != null) {
      for (Potential el: Electrodes) {
        if (el.z > 100 && el.z < 102) {
          strokeWeight(1);       
          stroke(#555555); // Gray for electrodes in X-Y plane
          point(graphOffsetX + el.y * SCALE, graphOffsetY + el.x * SCALE); // Draw electrode points
          strokeWeight(1);
        }
      }
    } else {
      text("No geometry is defined...", -50, 0, 0); // If no geometry is defined
    }

    // Draw X-Y trajectories
    if (ionTrajectory != null) {
      for (Trajectory coordinates: ionTrajectory) {
        strokeWeight(1);
        
        // Group 1 ions (blue)
        if (coordinates.group == 1) stroke(#0000FF, 100);// Blue for group 1
        // Group > 1 ions (red)
        if (coordinates.group > 1) stroke(#FF0000, 100); // Red for groups greater than 1
        point(graphOffsetX + coordinates.y * SCALE, graphOffsetY + coordinates.x * SCALE); // Draw trajectory points
        strokeWeight(1);
        stroke(0);
      }
    }

    // Visualization of excitation data
    int c = 0;
    float S = 1;
    float COEF = 0.5;
    float Xstart = 100;
    float Yoffset = 700;
    float l = (FREQ_SCAN_HIGH - FREQ_SCAN_LOW)/FREQ_SCAN_STEP;
    float opa = 180;
    float strokeWeight = 1;
    int binWidth = 4;
    try {
      for (int[] es : excitationSummary) {
        // Draw returned (blue) data
        strokeWeight(strokeWeight);
        stroke(#0000FF);
        //point((Xstart + c) * S, S * (Yoffset - es[1] * COEF)); // Returned data
        for ( int i =1; i<=binWidth; i++) {
          line((Xstart + c+i) * S, S * (Yoffset - es[1] * COEF), (Xstart + c+i) * S, S * (Yoffset - es[0] * COEF)); // Line for returned
        }
        // Draw transmitted (green) data
        strokeWeight(strokeWeight);
        //point((Xstart + c) * S, S * (Yoffset - 0 * COEF)); // Transmitted data
        //point((Xstart + c) * S, S * (Yoffset - 1000 * COEF)); // Transmitted data continuation
        stroke(#00FF00);
        //point((Xstart + c) * S, S * (Yoffset - es[0] * COEF)); // Exited data
        for ( int i =1; i<=binWidth; i++) {
          line((Xstart + c+i) * S, S * (Yoffset - es[0] * COEF), (Xstart + c+i) * S, S * (Yoffset - es[2] * COEF)); // Line for exited data
        }
        // Draw crushed (red) data
        strokeWeight(strokeWeight);
        stroke(#FF0000);
        //point((Xstart + c) * S, S * (Yoffset - es[2] * COEF)); // Crushed data
        for ( int i =1; i<=binWidth; i++) {
          line((Xstart + c+i) * S, S * (Yoffset - es[2] * COEF), (Xstart + c+i) * S, S * (Yoffset - 0 * COEF)); // Line for crushed data
        }
        c = c+binWidth;
      }
      //stroke(#4169E1);
      //line(20 * S, S * (Yoffset + 5), S * (l + 20), S * (Yoffset + 5)); // Scale line
    } catch (Exception e) {}
  }

  // ==================   Draw cross-section of the equipotential lines  ============================= 
  if (flag == 7) { 
    // Draw excitation profile
    if (excitationSummary != null) {
      int c = 0;
      float S = 1;
      float Yoffset = 300;
      float l = (FREQ_SCAN_HIGH - FREQ_SCAN_LOW)/FREQ_SCAN_STEP;
      
      strokeWeight(1);
      stroke(0);
      for (int[] es : excitationSummary) {
        point((20 + c) * S, S * (Yoffset - es[0])); // Exited
        point((20 + c) * S, S * ((Yoffset + 10) + es[2])); // Crushed
        c++;
      }
      stroke(0, 0, 255);
      line(20 * S, S * (Yoffset + 5), S * (l + 20), S * (Yoffset + 5)); // Scale line
    }
  } 
} 
