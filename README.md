# 🫁 Zaphyr - Frugal Arduino / ESP32-based Spirometer

![Python](https://img.shields.io/badge/Python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![C++](https://img.shields.io/badge/C++-00599C?style=for-the-badge&logo=c%2B%2B&logoColor=white)
![ESP32](https://img.shields.io/badge/ESP32-000000?style=for-the-badge&logo=espressif&logoColor=white)
![Healthcare](https://img.shields.io/badge/HealthTech-red?style=for-the-badge)

**Zaphyr** is a startup project aiming to develop a low-cost (frugal) spirometer for low-resource environments. It combines a differential pressure sensor (MXP5010DP), a microcontroller (ESP32), and an advanced Python script to acquire, process, and analyze patients' respiratory data to assist in the diagnosis of pulmonary diseases.

## ✨ Key Features

* **High-Frequency Acquisition:** Precise sampling at 1000 Hz for 10 seconds via the ESP32.
* **Advanced Signal Processing:** * Dynamic cut-off frequency calculation using Fast Fourier Transform (FFT).
  * 4th-order low-pass Butterworth filtering with zero-phase distortion.
* **Fluid Mechanics:** Conversion of pressure to volumetric flow rate (Bernoulli's Principle) with dynamic compensation for friction-induced pressure drops (Solving the Colebrook-White equation).
* **Automated Clinical Analysis:** * Extraction of standard parameters: **FVC** (Forced Vital Capacity), **FEV1** (Forced Expiratory Volume in 1 second), **PEF** (Peak Expiratory Flow).
  * Calculation of the Tiffeneau-Pinelli index (FEV1/FVC ratio).
  * Generation of a **pre-diagnosis** (normal, obstructive, or restrictive profile) by comparing results to theoretical values based on the patient's age, gender, and height.
* **Medical Visualization:** Automatic generation of standard spirometry curves (Volume-Time and Flow-Volume graphs).
* **Simulation Mode:** Built-in data generator simulating different clinical profiles to test the software without physical hardware.

## 🛠️ Hardware Requirements

* **Microcontroller:** ESP32 board (or compatible Arduino).
* **Sensor:** **MXP5010DP** differential pressure sensor.
* **Mechanics:** 3D-printed spirometry tube (Venturi effect with R1 = 13.9mm and R2 = 4.9mm).

## 📁 Project Architecture

* `codeArduino/codeArduino.ino`: Embedded script for the ESP32. Handles the 12-bit ADC, strict timing (1 ms), and serial communication at 115200 bauds.
* `zaphire.py`: Main Python script running on the computer. Manages serial communication, mathematical processing (Numpy/Scipy), console user interface, and graph plotting (Matplotlib).
* `A_frugal_arduino-based_spirometer_for_low-resource.pdf`: Documentation and research associated with the project.

## 🚀 Installation & Usage

### 1. Hardware Setup (ESP32)
1. Open `codeArduino.ino` in the Arduino IDE.
2. Connect the `Vout` pin of your sensor to pin `34` of the ESP32 (can be modified in the code).
3. Upload the code to the board.

### 2. Software Setup (Python)
Ensure Python 3 is installed. Install the required dependencies:
```bash
pip install pyserial numpy scipy matplotlib
```

### 3. Running the Analysis
1. Connect the ESP32 to your computer via USB.
2. Identify your board's port:
   * **Windows:** `COM3`, `COM4`, etc. (Check Device Manager).
   * **Mac/OSX:** `/dev/tty.usbserial-...` or `/dev/cu.usbserial-...`
3. If necessary, modify the `port` parameter in the `acquerir_donnees_serie()` function at the bottom of the `zaphire.py` file.
4. Run the script:
```bash
python zaphire.py
```
5. Follow the on-screen instructions (enter the patient's data, then blow into the tube when prompted!).

> 💡 **Tip:** To test the code without the ESP32, comment out the serial acquisition line in the main block (`__main__`) and uncomment the `generer_donnees_simulees(profil="obstructif")` function.

## ✍️ Author
* **Echchorfi Charaf** - *Computer Science and Information Systems Engineering student - UIASS*
