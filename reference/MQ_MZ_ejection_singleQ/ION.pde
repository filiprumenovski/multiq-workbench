/**
 *************************************************************************************************************************
 *
 *                                                         Ion
 *
 *************************************************************************************************************************
 */
class Ion {

    // Physical properties of the ion
    public float m;            // Absolute mass of the ion
    public float MW;           // Molecular weight of the ion
    public float q;            // Absolute charge of the ion (in Coulombs)
    public int Nq;             // Number of charges on the ion
    public double x, y, z;     // Ion coordinates in 3D space
    public double x1, y1, z1;  // Ion coordinates in the center-of-mass system
    public double px, py, pz;  // Previous position coordinates for trajectory calculation
    public double vx, vy, vz;  // Ion velocity components (in x, y, and z directions)
    public float minv;         // Inverse mass-related parameter (used in calculations)
    public double ax, ay, az;  // Ion acceleration components
    public float y0, z0;       // Initial positions along y-axis and z-axis
    public float radiee;       // Ion radius
    public float volts;        // Voltage associated with the ion
    public float vin;          // Initial velocity of the ion
    public float ang, angin;   // Initial angle of ion motion
    public float Kin, Kinr, Kinz;  // Kinetic energy and its components
    public float z_begin, y_begin; // Initial positions (y, z) for ion motion
    public float sees_pressure; // Surrounding pressure (in mTorr)
    public float mbuffer, mbg; // Mass of buffer gas and corresponding mass in atomic units
    public float Np, ccs;      // Number of particles (density) and collision cross-section
    public float mfp;          // Mean free path (ion scattering length)
    public double t;           // Current time
    public double deltat;      // Time step for the simulation
    public float deltaTbc;     // Time step for collisions
    public float end_time;     // End time of the simulation
    public float Urf, freq, omega, phaserf;  // RF field parameters (voltage, frequency, etc.)
    public float Udc;          // DC voltage potential
    public float tr;           // Not used in the provided context
    
    public int steps_bc;       // Steps taken for collision calculations
    public int colcount = 0;   // Collision count

    public float exitX, exitY, exitZ; // Exit point coordinates for the ion
    public float rex;                 // Radius of the exit point

    // Constructor for the Ion class
    Ion(float mass, int Ncharges, float VIN) {
    
        // Initialize physical properties
        m  = mass * amu;  // Convert mass to atomic mass units
        MW = mass;        // Set molecular weight
        q  = Ncharges * 1.60217646e-19;  // Convert charge to Coulombs
        Nq = Ncharges;    // Set number of charges

        // Initialize ion position in 3D space (scaled from given units)
        x = Xenter * 0.001; // [m]
        y = Yenter * 0.001; // [m]
        z = Zenter * 0.001; // [m]

        // Set initial ion velocity based on the given voltage
        float Ein = VIN;  // Initial energy in eV
        Kin = VIN * random(0.75, 1.5);  // Randomized kinetic energy
        vin = sqrt(2 * q * Kin/m);  // Calculate velocity based on kinetic energy
        angin = random(-alfa, +alfa);  // Random angle for ion entry
        ang = angin * 2 * PI/360;    // Convert angle to radians
        vz = vin * Math.cos(ang);      // Calculate velocity along z-axis
        vx = (Math.sqrt(2)/2) * vin * Math.sin(ang); // Velocity along x-axis
        vy = vx;  // Assuming symmetry between x and y components

        // Compute inverse mass-related factor
        minv = sqrt(0.035 * 2 * 1.6e-19/m);

        // Initialize time parameters
        t = 0;
        deltat = 1e-7;  // Time step for the simulation
        end_time = CALCULATION_TIME;  // Set the end time of the simulation

        // Set vacuum conditions
        mbuffer = 28;  // Molecular weight of the buffer gas
        mbg = mbuffer * amu;  // Convert buffer gas mass to atomic units
        sees_pressure = PRESSURE * 0.001 + 1e-10;  // Pressure surrounding the ion (converted to mTorr)
        Np = 2.68 * 1e5 * (sees_pressure/760);  // Calculate number of particles in m^3
        ccs = (2.81e-9) * MW * MW * MW - (3.55e-5) * MW * MW + (2.32e-1) * MW + 41.91; // Cross-section formula
        ccs = ccs * 10 * sqrt(Nq);  // Introduce charge dependence of the cross-section
        mfp = -log(random(0.0000001, 1.0))/(Np * ccs * 10.0);  // Mean free path

        // Calculate time step for collision events
        deltaTbc = mfp/sqrt((float)(vx * vx + vy * vy + vz * vz));
        steps_bc = (int)(deltaTbc/deltat) + 1;  // Number of steps per collision
        colcount = 0;

        // Initialize RF field and DC voltage parameters
        phaserf = random(0, 2 * PI);  // Random phase for RF field
        exitX = Xexit * 0.001;  // Set ion exit point in x-direction
        exitY = Yexit * 0.001;  // Set ion exit point in y-direction
        exitZ = Zexit * 0.001;  // Set ion exit point in z-direction
        rex = 2 * 0.001;  // Set the radius of the exit point
    }

}  // End of Ion class
