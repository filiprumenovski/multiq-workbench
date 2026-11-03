# **Demiurge_MultiQ-IT_486**

This _ad hoc_ program, developed in **Processing** and **Java**, simulates ion motion within a 486-quadrupole MultiQ-IT ion trap. The simulation is based on the computer model described in the work of **A.N. Krutchinsky, I.V. Chernushevich, V.L. Spicer, W. Ens, and K.G. Standing**:  *Collisional Damping Interface for an Electrospray Ionization Time-of-Flight Mass Spectrometer*,  **Journal of the American Society for Mass Spectrometry**, Volume 9, Issue 6, June 1998, Pages 569–579.

---

## **Setup Instructions**

### 1. **Install Processing**
- Install **Processing 4.3** (latest version as of December 2024).
- Earlier versions may work but are untested.
- Supported operating systems: **Windows 11**, **Linux** and **MacOS**
  
### 2. **Download Files**
- Download all `*.pde` files.
- Place them in a folder named `Demiurge_MultiQ_486` on your computer.

### 3. **Run the Program**
- Open `Demiurge_MultiQ_486.pde` in Processing.
- Run the program.
- A canvas window titled **Demiurge_MultiQ_486** will appear.

---

## **Usage Instructions**

### **Step 1: Create Geometry**
1. Move your mouse cursor over the program window.
2. Press **F1** (or **Alt+F1**) to open the dialog menu.
3. Select **"Create Geometry"** to define the electrode geometry.
   - Geometry details are in the `"Geometry"` folder, starting from **line 26**.
   - Default geometry uses the `definePotential_RF(float x, float y, float z)` function.

---

### **Step 2: Calculate Potential**
1. Press **F1** again and select the **Calculate Potential** option.
2. Once calculations finish, save the potential file in the program folder.
   - Default filename: **`RF`** (referenced by the `DEFAULT_RFPOT_FILE` variable).

---

### **Step 3: Define and Save Additional Geometries**
1. To calculate a different geometry, edit **line 26** in the `"Geometry"` folder.
   - Example: Replace `definePotential_RF(x, y, z)` with `definePotential_DC(x, y, z)` for DC potentials.
2. Repeat **Steps 1 and 2** to save a DC potential file named **`DC`**.

---

### **Step 4: Load Precalculated Potentials**
1. Save the generated `RF` (~100MB) and `DC` (~100MB) files in the `Demiurge_MultiQ_486` folder.
2. Run the program and press **F2** to open the options menu.
3. Select **"Load Precalculated Potentials"** to load the potential arrays into memory.

---

### **Step 5: Run Simulations**
1. Use the following key commands to examine geometry, potential, or ion trajectories:
   - **1–6**: Set flags to visualize various simulation parameters (refer to the `"Draw"` folder).
   - **i**: Examine the trajectory of a single ion with **m/z = 500/1+**.
   - **m**: Examine the trajectory of a single ion with **m/z = 1500/3+**.
2. Press **F2** to open the options menu again.
3. Select simulation options such as:
   - Simulating ion motion for multiple ions.
   - Computing ion depletion effects.

---

## **Additional Notes**
- You are welcome to modify or adapt the program as needed.
- If you encounter issues or have suggestions, contributions are encouraged.

