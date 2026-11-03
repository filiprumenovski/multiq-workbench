/**
 **************************************************************************************************
 *
 *                               Class representing the trajectory of an ion in the given x, y, z coordinates
 *
 **************************************************************************************************
 */
class Trajectory {
  
  public float x, y, z;  // Ion's position in x, y, z coordinates
  int group;             // Group identifier for the ion
  
  // Constructor to initialize the ion's position and group
  Trajectory(float ix, float iy, float iz, int iongroup) {
    x = ix; 
    y = iy; 
    z = iz;
    group = iongroup;
  }
  
}

/**
 **************************************************************************************************
 *
 *                               Class representing the excitation characteristics of an ion
 *
 **************************************************************************************************
 */
class Excition {
  
  public float m;      // Ion mass (in atomic mass units)
  public int charge;   // Ion charge (e.g., +1, +2, etc.)
  public float f;      // Resonant frequency for the ion
  public int out;      // Number of ions that exited the system
  public int in;       // Number of ions that returned to the system
  public int crushed;  // Number of ions that were crushed
  
  // Constructor to initialize the excitation properties of the ion
  Excition(float mass, int z, float freq, int exited, int returned, int killed) {
    m = mass;
    charge = z;
    f = freq;
    out = exited;
    in = returned;
    crushed = killed;
  }
  
}

/**
 **************************************************************************************************
 *
 *                               Class representing the kinetic energy of an ion
 *
 **************************************************************************************************
 */
class KineticEnergy {
  
  public float e;      // Ion's kinetic energy
  int group;           // Group identifier for the ion
  
  // Constructor to initialize the kinetic energy of the ion and its group
  KineticEnergy(float ie, int iongroup) {
    e = ie;
    group = iongroup;
  }
  
}
