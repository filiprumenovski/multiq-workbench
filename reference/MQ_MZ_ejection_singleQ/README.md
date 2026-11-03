# **MQ_MZ_Ejection_SingleQ**

This program is an _ad hoc_ derivative of the *Demiurge_MultiQ_486* program developed in **Processing** and **Java**. It simulates ion motion within a single quadrupole of a 486-quadrupole MultiQ-IT ion trap in the presence of a buffer gas. The quadrupole is configured for m/z-dependent ion ejection from the MultiQ-IT. The computer model is based on the work by **A.N. Krutchinsky, I.V. Chernushevich, V.L. Spicer, W. Ens, and K.G. Standing**:  *Collisional Damping Interface for an Electrospray Ionization Time-of-Flight Mass Spectrometer*,  **Journal of the American Society for Mass Spectrometry**, Volume 9, Issue 6, June 1998, Pages 569–579.

---

## **Program Overview**

The program calculates ion trajectories in a single quadrupole. Key features include:

- **Ion Motion:**  
  Ions enter the quadrupole from the MultiQ-IT confinement region. The program tracks their fate: ions are either:
  - Returned to the confinement region of the MultiQ.
  - Resonantly ejected from the quadrupole (and hence from the MultiQ).
  - Collided with electrodes during the ejection process.

- **Geometry:**  
  The quadrupole consists of four rods (6.35 mm in diameter, 12 mm long). Two neighboring rods are squeezed at the entrance by 0.02 radians (*Geometry folder*, line 516) to enhance resonant ion ejection.

- **Operating Modes:**  
  The quadrupole operates in RF-only mode with the addition of a small excitation voltage (see *MQ MZ Ejection SingleQ folder*, line 90).

- **Ejection Mechanism:**  
  - Ions are stopped by a 5V potential difference (*MQ MZ EjectionQ folder*, line 86).  
  - Ions whose oscillation frequency resonates with the excitation frequency experience increased oscillation amplitude.  
  - Excited ions are either:  
    - Pushed towards the plate by the tapered electric field of the squeezed rods, overcoming the potential barrier and exiting the trap.  
    - Returned to the confinement region of the MultiQ.
    - Killed in the collision with the electrodes.

- **Potential Calculations:**  
  Before simulating ion motion, one needs to calculate three potentials:  
  - **RF:** Governs RF-only motion within the quadrupole.  
  - **DC:** Affects motion near the stopping plate at the quadrupole exit.  
  - **AC:** Drives ion excitation between unsqueezed rods.

---

## **Setup Instructions**

### 1. **Download Files**
- Download all `*.pde` files.
- Place them in a folder named `MQ_MZ_ejection_singleQ` on your computer.

### 2. **Install Processing**
- Install **Processing 4.3** (latest version as of December 2024).
- Earlier versions may work but are untested.
- Supported operating systems: **Windows 11**, **Linux**, and **MacOS**.

### 3. **Run the Program**
- Open `MQ_MZ_ejection_singleQ.pde` in Processing.
- Run the program to open the **MQ_MZ_ejection_singleQ** canvas window.

---

## **Usage Instructions**

### **Step 1: Create Geometry**
1. Run the program.
2. Press **F1** (or **Alt+F1**) to open the dialog menu.
3. Select **"Create Geometry"** to define the electrode geometry.
   - Geometry details are in the `"Geometry"` folder (starting from **line 26**).

---

### **Step 2: Calculate Potential**
1. Press **F1** again and select **Calculate Potential**.
2. Save the resulting potential file in the program folder:
   - Default filename for RF: **`RF`** (referenced by the `DEFAULT_RFPOT_FILE` variable).

---

### **Step 3: Load Precalculated Potentials**
1. Place the generated `RF`, `DC`, and `AC` files (each ~120MB) in the `MQ_MZ_ejection_singleQ` folder.
2. Run the program and press **F2** to open the options menu.
3. Select **"Load Precalculated Potentials"** to load these files into memory.

---

### **Step 4: Run Simulations**
1. Use the following key commands to examine geometry, potential, or ion trajectories:
   - **1–6**: Visualize simulation parameters (refer to the `"Draw"` folder).
   - **i**: View the trajectory of a single ion with **m/z = 500/1+**.
   - **m**: View the trajectory of a single ion with **m/z = 1500/3+**.
2. Press **F2** to open the options menu again.
3. Select simulation options such as:
   - Ion trajectory for charges **1+** or **5+**.
   - Resonance frequency calculations.
   - Spectrum generation for ions in the defined **MOZ_RANGE**.

---

## **Additional Notes**
- Modifications or adaptations of the program are encouraged.
- For issues or suggestions, contributions are welcome.

---
