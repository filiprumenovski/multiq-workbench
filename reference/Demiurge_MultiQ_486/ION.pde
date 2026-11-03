/**
 *************************************************************************************************************************
 *
 *                                                         Ion Class
 *
 *************************************************************************************************************************
 */
class Ion {
  
    // Physical properties of the ion
    public float m;            // Absolute mass
    public float MW;           // Molecular weight
    public float q;            // Absolute ion charge
    public int Nq;             // Number of charges

    // Ion coordinates and velocity
    public double x, y, z;     // Current ion coordinates
    public double x1, y1, z1;  // Center-of-mass system coordinates
    public double px, py, pz;  // Previous trajectory point
    public double vx, vy, vz;  // Velocity components
    public float minv;         // Minimum velocity magnitude
    public double ax, ay, az;  // Acceleration components

    // Initial position and energy properties
    public float y0, z0;       // Initial position along y- and z-axis
    public float radiee;       // Ion radius (if applicable)
    public float volts;        // Voltage
    public float vin;          // Initial velocity
    public float ang, angin;   // Injection angles
    public float Kin, Kinr, Kinz; // Kinetic energy components
    public float z_begin, y_begin; // Starting positions

    // Buffer gas and collision parameters
    public float p;            // Pressure [mTorr]
    public float mbuffer, mbg; // Mass of buffer gas molecules
    public float Np;           // Number density of particles
    public float ccs;          // Collision cross-section
    public float mfp;          // Mean free path

    // Time parameters
    public double t;           // Current time
    public double deltat;      // Time step
    public float deltaTbc;     // Time to the next boundary collision
    public float end_time;     // Maximum calculation time
    public int steps_bc;       // Steps between collisions
    public int colcount = 0;   // Collision count
    public float tr;           // trajectory sampling dist

    // Electrical fields
    public float Urf, freq, omega, phaserf; // RF field properties
    public float Udc;                       // DC potential
   
    
    // Exit conditions
    public float exitX, exitY, exitZ; // Exit coordinates
    public float rex;                 // Exit radius

    /**
     * Constructor for the Ion class
     * @param mass Mass of the ion
     * @param Ncharges Number of charges
     * @param VIN Initial voltage
     * @param URF RF voltage amplitude
     * @param F RF frequency
     * @param PHASE RF phase
     */
    Ion(float mass, int Ncharges, float VIN, float URF, float F, float PHASE) {

        // Ion properties
        m = mass * amu;
        MW = mass;
        q = Ncharges * 1.60217646e-19;
        Nq = Ncharges;

        // Initial ion coordinates (converted to meters)
        x = Xenter * 0.001;
        y = Yenter * 0.001;
        z = Zenter * 0.001;

        // Initial velocity calculation
        Kin = VIN;
        vin = sqrt(2 * q * Kin/m);
        angin = random(-5, 5);          // Injection angle (random within ±5°)
        ang = angin * 2 * PI/360;     // Convert angle to radians
        vy = -vin * Math.cos(ang);      // Initial y-velocity
        vx = 0.707 * vin * Math.sin(ang); // Initial x-velocity
        vz = vx;                         // Initial z-velocity

        minv = sqrt(0.035 * 2 * 1.6e-19/m); // Minimum velocity

        // Time parameters
        t = 0;
        deltat = 1e-7;
        end_time = CALCULATION_TIME;

        // Buffer gas properties
        mbuffer = 28;                   // Molecular weight of buffer gas
        mbg = mbuffer * amu;            // Mass of buffer gas molecules
        p = PRESSURE * 0.001 + 1e-10;   // Surrounding pressure [mTorr]
        Np = 2.68e5 * (p/760);        // Particle density

        // Collision cross-section calculation
        //ccs = 4*pow(MW, 0.6); // in [A]^2 OLD EMPIRICAL FORMULA (see AK thesis)
        //float CCS = -2.724e-5*MW*MW +2.141e-1*MW+40.80; //FORMULA from D.Clemmer et al " A Database of 660 Peptide..." JASMS 1999,101188-1211  
        float CCS_He = (2.81e-9) * MW * MW * MW - (3.55e-5) * MW * MW + (2.32e-1) * MW + 41.91; // D. Russel et al "A Collision Cross-Section Database.." JASMS2007,18,1232-1238
        float COEF = 10;                // Correction factor for N2 (impirical)
        ccs = COEF * CCS_He;
        mfp = -log(random(0.0, 1.0))/(Np * ccs); // Mean free path

        // Time to boundary collision
        deltaTbc = mfp/sqrt((float)(vx * vx + vy * vy + vz * vz));
        steps_bc = (int)(deltaTbc/deltat) + 1;
        colcount = 0;

        // RF and DC field properties
        Urf = URF;
        freq = F;                      // RF frequency [kHz]
        omega = 2 * PI * freq * 1000;  // Angular frequency
        phaserf = 2 * PI * PHASE/360.0; // Phase in radians
        Udc = REPEL_POTENTIAL;         // DC voltage [V]

        // Exit properties
        exitX = Xexit * 0.001;
        exitY = Yexit * 0.001;
        exitZ = Zexit * 0.001;
        rex = 2 * 0.001;               // Exit radius
        
        println("//----------------------------------- ION SUMMARY  -----------------------------------");
        println("ion MW: " + MW +" Da; charge: " + Nq);
        //println("Initial position x : " + x +" y: " + y + " z: " +z  + " m");
        //println("Initial velocity Vx : " + vx +" m/z Vy: " + vy + " m/s Vz: " + vz  + " m/s");  
        println("Pressure of buffer gas is "+ p + " Torr");
        println("The mass of buffer gas is "+ mbuffer);
        println("The initial mean free path is "+ mfp + " [m]");
        println("The RF frequency: "+ freq + " [kHz]" + ", amplitude: " + Urf + " [V]");  
        println("Diameter of electrodes: " + De +  " [mm]");
        println("q: " + 4*q*Urf/(m*(0.001*De/2)*(0.001*De/2)*omega*omega) );  
        println("CCS: " + ccs );
    }
}
