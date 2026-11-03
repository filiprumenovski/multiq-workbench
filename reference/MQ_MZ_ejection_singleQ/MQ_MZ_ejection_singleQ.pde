
/**
 ************************************************************************************************************************
 *
 *                            Simulates the ion ejection process in a single quadrupole exit of MultiQ-IT
 *
 * @author Andrew Krutchinsky
 * @organization: The Rockefeller University
 *
 ************************************************************************************************************************
 */
 
import javax.swing.JOptionPane;
import java.awt.event.KeyEvent;
import javax.swing.SwingUtilities;
import java.util.Scanner;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import java.lang.Math;
import java.util.Random;

// File paths for potential files
String Directory = "C:\\Users\\Andrew\\MultiQ-IT\\PAPER_v2024\\PROGRAMS\\MQ_MZ_ejection_singleQ\\";
String DEFAULT_RFPOT_FILE = Directory + "RF";  // Path to RF potential file
String DEFAULT_DCPOT_FILE = Directory + "DC";  // Path to DC potential file
String DEFAULT_ACPOT_FILE = Directory + "AC";  // Path to AC potential file

// Geometry of the calculation area
float xDevice = 16;       // Width of the calculation area in mm
float yDevice = 16;       // Length of the calculation area in mm
int Ne = 1;               // Dimensionality of the electrode array

// Quadrupole and setup dimensions
float Qstart = 3;         // Position of the quadrupole ejection setup from z = 0 mm  
float Q1 = 12;            // Length of the ejection quadrupole in mm
float GAP = 1;            // Gap between the ejection quadrupole and stopping plate in mm
float PLATE = 1.6;        // Thickness of the stopping plate in mm
float Q2 = 3;             // Length of the RF-only quadrupole after the stopping plate
float zDevice = Qstart + Q1 + GAP + PLATE + GAP + Q2 + GAP + PLATE; // Total length of the calculation volume along the z-axis 
float DH = 7;             // Diameter of the hole in the stopping plate

// Scale factor for rendering the model
float scl = 10;           // Scale: pixels per mm



// Electrode dimensions and configuration
float CubeWidth = 1.5;    // Width of the "CUBE" electrode
float De = 6.35;          // Base spacing of electrodes
float spacing = 0.323 * De; // Gap between electrodes (0.323 is a universal coefficient)
float He = De;            // Height of the electrode
float gap = spacing/2;  // Gap between RF electrode and ground
float stem = 2;           // Radius of the quadrupole screw
// Radius of the field (based on electrode spacing and size)
float Rfield = (De + spacing)/sqrt(2) - De/2; // Field radius

// Electric potential values
float U = 1.0;            // 1 V potential
float U0 = 0;             // Ground potential (0 V)
float E_SPACE_CHARGE = 1000; // Space charge electric field in V/m

// Ion entrance and exit coordinates
float Xenter = xDevice/2; // x-coordinate of the ion entrance [mm]
float Yenter = yDevice/2; // y-coordinate of the ion entrance [mm]
float Zenter = 1;         // z-coordinate of the ion entrance [mm]
float Xexit = xDevice/2;  // x-coordinate of the ion exit [mm]
float Yexit = yDevice/2;  // y-coordinate of the ion exit [mm]
float Zexit = zDevice;    // z-coordinate of the ion exit [mm]
  
// Checkpoint for potential calculation (near the second exit)
float x_check = xDevice/2 + 1;    // x-coordinate of the checkpoint [mm]
float y_check = yDevice/2 + 0.5;  // y-coordinate of the checkpoint [mm]
float z_check = Qstart + Q1 - 1;  // z-coordinate of the checkpoint [mm]
 
// Environmental pressure in mTorr
float PRESSURE = 1; // mTorr
  
// Test ion properties
float testM = 2500;  // Ion mass in atomic mass units (amu)
int testCharge = 5;  // Ion charge state

float deltaM = 10;   // Mass difference for scanning

// RF, DC, and AC potential settings
float A_RF = 100;    // RF amplitude in volts
float V_DC = 5;      // DC potential on the exit plate in volts
float V_AC = 0.5;    // AC potential in volts
 
// RF frequency settings
float RF_FREQUENCY = 500; // RF frequency in kHz
float RESONANCE_FR = 100; // Resonance frequency at 100 V RF
float FREQ_SCAN_LOW = 80; // Lower bound for frequency scanning in kHz
float FREQ_SCAN_HIGH = 120; // Upper bound for frequency scanning in kHz
float FREQ_SCAN_STEP = 1; // Step size for frequency scanning in kHz

// m/z scanning range and step
float MOZ_RANGE = 25;  // m/z scanning range
float MOZ_DELTA = 2;   // m/z scanning step

// List to store excitation summary data
ArrayList<int[]> excitationSummary = new ArrayList<int[]>(3);

// Parameters for ion ejection process
// Ion energy in electron volts (eV)
float ION_ENERGY = 0.1; // Energy of ions in eV
// Angle of ion entry (degrees)
float alfa = 10; // +/- angle at which ions enter the system

// Time duration for calculation (in seconds)
float CALCULATION_TIME = 1; // Time for simulation in seconds

// Number of ions to simulate
int Nions = 10000; // Total number of ions to simulate



// Number of iterations for the simulation
int N_ITERATIONS = 6000; // Number of iterations for the simulation

// Flags for mode switching
int flag = 0; // Flag to switch between different modes

// Rotation angles for the system
float a1 = 0.0; // Rotation angle 1
float a2 = 0.0; // Rotation angle 2
float a3 = 0.0; // Rotation angle 3

// Display scale (controls the resolution of the visualization)
int DISPLAY_SCALE = 2; // Set to 6 for high-resolution screens

// Time interval to sample ion parameters (in seconds)
float SAMPLING_TIME = 0.0001; // Sampling time for ion parameters in seconds

// Counters for different simulation results
int ionCount = 1; // Ion counter

// Counts for trapped ions and ions that were crushed
int sNtrap = 0; // Number of trapped small ions
int mNtrap = 0; // Number of trapped medium ions
int sIonsCrushed = 0; // Small ions crushed
int mIonsCrushed = 0; // Medium ions crushed

// Counts for ion behavior at exits
int Nreturned = 0; // Number of ions that returned to the trap
int Nexited = 0;   // Number of ions that exited the trap
int Ncrushed = 0;  // Number of ions crushed

// Target ions for comparison (1+ vs 5+ ions)
int Ntarget = 2000; // Number of target ions for comparison

// Derived parameters for grid dimensions
int Nx = (int)(xDevice * scl);  // Number of grid points along the x-axis
int Ny = (int)(yDevice * scl);  // Number of grid points along the y-axis
int Nz = (int)(zDevice * scl);  // Number of grid points along the z-axis

// z-coordinate for drawing cross-sections of equipotential lines
int zcs = Nz - 2;

// Step sizes for the spatial grid in each direction (in mm)
float step_z = 1/scl;         // Step size in the z-direction (0.1 mm per grid point)
float step_x = 1/scl;         // Step size in the x-direction (0.1 mm per grid point)
float step_y = 1/scl;         // Step size in the y-direction (0.1 mm per grid point)

// Spatial grid resolution in meters
double deltas = (xDevice/Nx) * 0.001;  // Spatial resolution for the grid in meters

// Arrays for storing potential values and ion trajectory data
ArrayList<Potential> Electrodes; // Array for electrodes
ArrayList<Potential> RFstructure; // Array for RF structure
ArrayList<Potential> DCstructure; // Array for DC structure

// Arrays for output potentials
ArrayList<String> RFoutputPotential = new ArrayList<String>(); // RF output potentials
ArrayList<String> DCoutputPotential = new ArrayList<String>(); // DC output potentials

// Arrays for ion trajectories and kinetic energy data
ArrayList<Trajectory> ionTrajectory; // Array for ion trajectory data
ArrayList<KineticEnergy> ionEnergies = new ArrayList<KineticEnergy>(); // Array for ion energies

// Number of ions that stayed inside the trap
int Ntrap = 0;

// Number of ions lost from the trap
int Nlost = 0;

// Counter for sampling iterations
int samplingCount = 0;

// Arrays for potential data in the grid (Nx, Ny, Nz dimensions)
Potential[][][] POT = new Potential[Nx][Ny][Nz];
Potential[][][] RFPOT = new Potential[Nx][Ny][Nz];
Potential[][][] DCPOT = new Potential[Nx][Ny][Nz];
Potential[][][] ACPOT = new Potential[Nx][Ny][Nz];

// Variables for drawing (used in the "setup" function)
int Sstart, Mstart, Hstart, Send, Mend, Hend, Tstart, Tend;

// Average times for ions exiting through specific exit plates
float averageTimeExit1 = 0; // Average time for ions exiting through Exit 1
float averageTimeExit5 = 0; // Average time for ions exiting through Exit 5

// Constants for physical quantities
float amu = 1.660539e-27; // Atomic mass unit in kilograms
float qe = 1.60217646e-19; // Elementary charge in Coulombs

// Window dimensions for visualization (pixels)
int w = 1300, h = 800;

// Setup function to initialize the simulation environment
void setup() {
    System.setProperty("jogl.disable.openglcore", "false"); // Disable OpenGL core profile
    size(1300, 800, P3D); // Set window size (1300x800) and enable 3D rendering
    background(255); // Set background color to white
    
    keyPressed(); // Call keyPressed function to redraw the visualization
} // End of setup function







 
