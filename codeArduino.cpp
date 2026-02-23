// --- CONFIGURATION ---
// Sur un ESP32, les broches ADC courantes sont 32, 33, 34, 35, 36, 39.
// À adapter selon le branchement de la broche Vout de votre capteur MXP5010DP.
const int SENSOR_PIN = 34; 

const unsigned long INTERVALLE_US = 1000; // 1000 microsecondes = 1 milliseconde = 1000 Hz
unsigned long previousMicros = 0;

bool isRunning = false;
int sampleCount = 0;
const int MAX_SAMPLES = 10000; // 10 secondes d'enregistrement

void setup() {
  // Initialisation de la communication série à 115200 bauds (doit correspondre au script Python)
  Serial.begin(115200);
  
  // Configuration de la résolution de l'ADC pour l'ESP32 (12 bits par défaut = 0-4095)
  analogReadResolution(12);
  
  while (!Serial) {
    ; // Attendre que le port série soit prêt
  }
}

void loop() {
  // 1. Écoute des commandes venant de l'ordinateur
  if (Serial.available() > 0) {
    String command = Serial.readStringUntil('\n');
    command.trim(); // Nettoyer les caractères invisibles (\r)

    if (command == "START") {
      isRunning = true;
      sampleCount = 0;
      previousMicros = micros(); // Initialisation du chronomètre de précision
    } 
    else if (command == "STOP") {
      isRunning = false;
    }
  }

  // 2. Boucle d'acquisition cadencée
  if (isRunning) {
    unsigned long currentMicros = micros();
    
    // Si exactement 1 milliseconde (ou plus) s'est écoulée
    if (currentMicros - previousMicros >= INTERVALLE_US) {
      // On ajoute l'intervalle strict pour éviter toute dérive temporelle cumulée
      previousMicros += INTERVALLE_US; 
      
      // Lecture de la tension du capteur et envoi direct
      int sensorValue = analogRead(SENSOR_PIN);
      Serial.println(sensorValue);
      
      sampleCount++;
      
      // Arrêt automatique de sécurité une fois les 10 000 points atteints
      if (sampleCount >= MAX_SAMPLES) {
        isRunning = false;
      }
    }
  }
}