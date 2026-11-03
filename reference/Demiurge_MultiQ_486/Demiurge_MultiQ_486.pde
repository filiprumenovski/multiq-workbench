/**
 ************************************************************************************************************************
 *
 * Ad hoc DEMIURGE: Computes ion motion in a 486-quadrupole MultiQ-IT
 *
 * @author Andrew Krutchinsky
 * @organization The Rockefeller University
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

// in windows use \\ to escape \; linux is sane
String DEFAULT_RFPOT_FILE = "C:\\Users\\Andrew\\MultiQ-IT\\PAPER_v2024\\PROGRAMS\\Demiurge_MultiQ_486\\RF";
String DEFAULT_DCPOT_FILE = "C:\\Users\\Andrew\\MultiQ-IT\\PAPER_v2024\\PROGRAMS\\Demiurge_MultiQ_486\\DC";

// Dimensions of the calculation area in mm
float xDevice = 83.0;  // Width
float yDevice = xDevice;  // Length
float zDevice = xDevice;  // Depth
int Ne = 10;  // Number of electrodes in the array

float scl = 2;  // Scale: pixels/mm

// Electrode dimensions
float CubeWidth = 1.5;  // Width of the "CUBE" electrode
float De = ((zDevice - 2 * CubeWidth)/Ne)/(1 + 0.323);  // Spacing (calculated based on dimensions)
float spacing = 0.323 * De;  // Gap between electrodes (0.323 is a universal coefficient)
float He = De;  // Height of the electrode

// Gap and screw dimensions
float gap = spacing/2;  // Gap between RF electrode and ground
float stem = 2;  // Radius of the quadrupole screw

// Potential values
float U = 1.0;  // RF voltage [V]
float U0 = 0;  // Ground potential [V]

float E_SPACE_CHARGE = 1;  // Space charge field [V/m]

// Ion entrance and exit coordinates
float Xenter = xDevice/2;  // X-coordinate of entrance [mm]
float Yenter = yDevice - 5;  // Y-coordinate of entrance [mm]
float Zenter = zDevice/2;  // Z-coordinate of entrance [mm]

float Xexit = xDevice/2;  // Central exit of the SIDE plate [mm]
float Yexit = yDevice/2;
float Zexit = zDevice;

float Xexit1 = xDevice/2;  // Additional exit 1 [mm]
float Yexit1 = 3 * (gap + De + gap) + CubeWidth;
float Zexit1 = zDevice;

float Xexit2 = xDevice/2;  // Additional exit 2 [mm]
float Yexit2 = 7 * (gap + De + gap) + CubeWidth;
float Zexit2 = zDevice;

// Checkpoint for potential calculation near the second exit
float x_check = xDevice/2 + 1;  // [mm]
float y_check = 7 * (gap + De + gap) + CubeWidth + 2;  // [mm]
float z_check = zDevice - 9;  // [mm]

// Flags and modes
int flag = 0;  // Mode switch
int cut = 0;  // Geometry modification switch
float a1 = 0.0;  // Rotation angle (X-axis)
float a2 = 0.0;  // Rotation angle (Y-axis)
float a3 = 0.0;  // Rotation angle (Z-axis)

// Display parameters
int DISPLAY_SCALE = 3;  // Scale for display; set to 6 for high-resolution screens

// Ion trapping and crushing statistics
int sNtrap = 0, mNtrap = 0;
int sIonsCrushed = 0, mIonsCrushed = 0;

// Simulation parameters
float PRESSURE = 0.5;  // Pressure in mTorr
float RF_AMPLITUDE = 50;  // RF amplitude [V]
float RF_FREQUENCY = 500;  // RF frequency [kHz]
float ION_ENERGY = 1;  // Ion energy [eV]
float REPEL_POTENTIAL = 0.1;  // Repelling potential [V]
float CALCULATION_TIME = 1;  // Simulation time [s]
float SAMPLING_TIME = 0.0001;  // Sampling interval [s]
float E_SCALE = 500;  // Scale for energy visualization

int Ntarget = 200;  // Target number of ions
int N_ITERATIONS = 4000;  // Iteration count
int ionCount = 1;  // Current ion count

/**
 *************************************************************************************************************************
 * Parameters derived from main variables
 *************************************************************************************************************************
 */

int Nx = (int) (xDevice * scl);  // X-dimension grid size
int Ny = (int) (yDevice * scl);  // Y-dimension grid size
int Nz = (int) (zDevice * scl);  // Z-dimension grid size

int zcs = Nz - 1;  // Z-value for cross-section plotting

float step_z = 1/scl;  // Step size in Z direction [mm]
float step_x = 1/scl;  // Step size in X direction [mm]
float step_y = 1/scl;  // Step size in Y direction [mm]

double deltas = (xDevice/Nx) * 0.001;  // Spatial grid resolution [m]

// Data structures
ArrayList<Potential> Electrodes;  // Electrode potentials
ArrayList<Potential> RFstructure;  // RF structure
ArrayList<Potential> DCstructure;  // DC structure

ArrayList<String> RFoutputPotential = new ArrayList<>();  // RF potentials
ArrayList<String> DCoutputPotential = new ArrayList<>();  // DC potentials
ArrayList<Trajectory> ionTrajectory;  // Ion trajectories
ArrayList<KineticEnergy> ionEnergies = new ArrayList<>();  // Ion kinetic energies

// Trap statistics
int Ntrap = 0;  // Ions trapped
int Nlost = 0;  // Ions lost

int samplingCount = 0;

// Potential grids
Potential[][][] POT = new Potential[Nx][Ny][Nz];
Potential[][][] RFPOT = new Potential[Nx][Ny][Nz];
Potential[][][] DCPOT = new Potential[Nx][Ny][Nz];
Potential[][][] Q1POT = new Potential[Nx][Ny][Nz];
Potential[][][] Q2POT = new Potential[Nx][Ny][Nz];

// Time tracking
int Sstart, Mstart, Hstart, Send, Mend, Hend, Tstart, Tend;
float averageTimeExit1 = 0;
float averageTimeExit5 = 0;

// Constants
float amu = 1.66e-27;  // Atomic mass unit
float q = 1.60217646e-19;  // Elementary charge [C]
float kB = 1.380649e-23;  // Boltzmann constant

int w = 1000, h = 1000;

//======================================= SETUP ===============================================
void setup() {
    System.setProperty("jogl.disable.openglcore", "false");
    size(1000, 1000, P3D);
    background(255);
    keyPressed();  // Trigger initial redraw
}
