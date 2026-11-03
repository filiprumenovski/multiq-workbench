//=======================================   DRAW   ===============================================
void draw() {
  if (flag == 0) return; // Exit the function if no flag is set

  //===================================   Draw RF Geometry   ====================================
  if (flag == 1) { // Draw 3D geometry of RF field
    background(255); // Set background to white
    translate(width/2, height/2, -width); // Move origin to center
    rotateX(a1);
    rotateY(a2);
    rotateZ(a3);
    scale(DISPLAY_SCALE); // Scale according to a global constant

    // Draw coordinate axes
    stroke(#666666); // Gray for axes
    strokeWeight(0.5);
    fill(0);
    line(-100.0, 0, 0, 300, 0, 0);
    text("Z", -100, 0, 0);
    line(0, -100, 0, 0, 1000, 0);
    text("Y", 0, -100, 0);
    line(0, 0, -100, 0, 0, 1000);
    text("X", 0, 0, -100);

    // Draw electrodes
    if (Electrodes != null) {
      for (Potential el : Electrodes) {
        strokeWeight(4);
        if (el.potential > 0) stroke(#FF0000, 100); // Positive potential: red
        else if (el.potential < 0) stroke(#0000FF, 100); // Negative potential: blue
        else stroke(#888888, 100); // Zero potential: gray
        point(el.z, el.x, el.y); // Plot point
        strokeWeight(1);
      }
    } else {
      text("No geometry is defined...", -50, 0, 0);
    }
  }

  //=================   Draw a Cross-Section of the Potential   =============================
  if (flag == 2) {
    background(255); // Set background to white
    int offsetY = Ny;
    int offsetX = Nx;
    float delta = 0.03;
    float step = 0.2;
    strokeWeight(2);

    println("Drawing equipotentials at z: " + zcs);
    stroke(#3366CC); // Blue for equipotentials
    for (int i = 1; i < Nx; i++) {
      for (int j = 1; j < Ny; j++) {
        for (float p = -1 + 0.1; p <= 1 - 0.1; p += step) { // Set of equipotentials
          if (p - delta < POT[i][j][zcs].potential && POT[i][j][zcs].potential < p + delta) {
            point((i * 2 + w/2 - offsetX), (2 * j + h/2 - offsetY));
          }
          if (POT[i][j][zcs].electrode) {
            if (POT[i][j][zcs].potential > 0) stroke(#FF0000); // Positive: red
            else if (POT[i][j][zcs].potential < 0) stroke(#0000FF); // Negative: blue
            else stroke(#CCCCCC); // Neutral: light gray
            point((i * 2 + w/2 - offsetX), (j * 2 + h/2 - offsetY));
            stroke(#3366CC); // Reset stroke
          }
        }
      }
    }
  }

  //=================   Draw Crossection of Equipotential Lines   =============================
  if (flag == 4) {
    // Uncomment and implement as needed, ensuring the same color optimization principles
  }

  //=====================================  Draw 3D Trajectories  ==================================
  if (flag == 5) {
    background(255); // Set background to white
    translate(width/2, height/2, -width); // Center the scene
    rotateX(a1);
    rotateY(a2);
    rotateZ(a3);
    scale(DISPLAY_SCALE);

    // Draw coordinate axes
    stroke(#666666); // Gray for axes
    strokeWeight(0.25);
    line(-100.0, 0, 0, 300, 0, 0);
    text("X", -100, 0, 0);
    line(0, -100, 0, 0, 1000, 0);
    text("Y", 0, -100, 0);
    line(0, 0, -100, 0, 0, 1000);
    text("Z", 0, 0, -100);

    // Draw ion trajectories
    if (ionTrajectory != null) {
      for (Trajectory coordinates : ionTrajectory) {
        strokeWeight(0.5);
        stroke(#333333); // Dark gray for trajectories
        point(coordinates.x, coordinates.y, coordinates.z);
      }
    }

    // Draw electrodes
    if (Electrodes != null) {
      for (Potential el : Electrodes) {
        strokeWeight(1);
        stroke(#888888, 50); // Semi-transparent gray
        point(el.x, el.y, el.z);
      }
    } else {
      text("No geometry is defined...", -50, 0, 0);
    }
  }

  //=======================================  Draw 2D (Plane) Trajectories  ==========================
  if (flag == 6) {
    background(255); // Set background to white
    int SCALE = DISPLAY_SCALE;
    float UpperLeftZ = w/2 - zDevice/2 * scl * SCALE;
    float UpperLeftY = h/2 - zDevice/2 * scl * SCALE;

    // Draw the device plane
    fill(255);
    rect(UpperLeftZ + CubeWidth * scl * SCALE, UpperLeftY + CubeWidth * scl * SCALE, 
         (zDevice - 2 * CubeWidth) * scl * SCALE, (xDevice - 2 * CubeWidth) * scl * SCALE);

    // Draw electrodes
    fill(#AAAAAA); // Neutral gray for electrodes
    for (int k = 1; k <= Ne; k++) {
      for (int l = 1; l <= Ne; l++) {
        float Qk = CubeWidth + (gap + De/2) + (De + spacing) * (k - 1);
        float Ql = CubeWidth + (gap + De/2) + (De + spacing) * (l - 1);
        if (k == 1 || k == Ne) {
          rect(UpperLeftZ + (Ql - De/2) * scl * SCALE, UpperLeftY + (Qk - De/2) * scl * SCALE, 
               De * scl * SCALE, De * scl * SCALE);
        } else {
          rect(UpperLeftZ + (CubeWidth + gap) * scl * SCALE, UpperLeftY + (Qk - De/2) * scl * SCALE, 
               De * scl * SCALE, De * scl * SCALE);
          rect(UpperLeftZ + (zDevice - CubeWidth - gap - De) * scl * SCALE, UpperLeftY + (Qk - De/2) * scl * SCALE, 
               De * scl * SCALE, De * scl * SCALE);
        }
      }
    }

    // Draw trajectories
    if (ionTrajectory != null) {
      for (Trajectory coordinates : ionTrajectory) {
        if (coordinates.group == 1) stroke(#1F77B4); // Blue for group 1
        else stroke(#D62728); // Red for other groups
        point(UpperLeftZ + coordinates.z * SCALE, UpperLeftY + coordinates.y * SCALE);
      }
    }

    // Draw kinetic energy visualization
    float START_X = 0;
    float START_Y = 700;
    float tSCALE = E_SCALE;
    float MAX_VAL = 500;
    float ENERGY_SCALE = 10;
    if (ionEnergies != null) {
      try {
        for (KineticEnergy E : ionEnergies) {
          stroke(#333333); // Dark gray for energy bars
          line(UpperLeftZ + START_X + E.t * tSCALE, UpperLeftY + (START_Y - E.e * ENERGY_SCALE), 
               UpperLeftZ + START_X + E.t * tSCALE, UpperLeftY + START_Y);
          line(UpperLeftZ + START_X, UpperLeftY + START_Y, 
               UpperLeftZ + START_X + MAX_VAL, UpperLeftY + START_Y);
        }
      } catch (Exception e) {
        println("ConcurrentModificationException");
      }
    }
  }
}
