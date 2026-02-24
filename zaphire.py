import serial
import time
import numpy as np
from scipy import signal

# --- CONSTANTES PHYSIQUES ET MATÉRIELLES ---
RHO = 1.204  # Densité de l'air en kg/m3 aux conditions standards [cite: 86]
R1 = 0.0139  # Rayon d'entrée ajusté (13.9 mm) suite aux tolérances d'impression 3D [cite: 203]
R2 = 0.0049  # Rayon d'étranglement ajusté (4.9 mm) [cite: 203]
A1 = np.pi * (R1 ** 2)  # Aire de la section 1 [cite: 86]
A2 = np.pi * (R2 ** 2)  # Aire de la section 2 [cite: 86]
MU = 1.81e-5  # Viscosité dynamique de l'air en kg/(m*s)
RUGOSITE = 10e-6  # Rugosité estimée de la paroi imprimée (10 µm)
D_ENTREE = R1 * 2  # Diamètre d'entrée
FACTEUR_ECHELLE = 1.47  # Facteur pour compenser le pont diviseur de tension [cite: 127, 131]
N_SAMPLES = 10000  # 10 secondes d'enregistrement à 1000 Hz [cite: 116, 125]
FE = 1000.0  # Fréquence d'échantillonnage en Hz [cite: 125]


def acquerir_donnees_serie(port='/dev/ttyUSB0', baudrate=115200):
    """
    Établit la connexion avec l'ESP32, envoie un signal de départ,
    et collecte exactement N_SAMPLES données brutes.
    Gère les erreurs de déconnexion et de port matériel.
    """
    print(f"\nTentative de connexion au capteur sur le port {port} à {baudrate} bauds...")
    donnees_brutes = []

    try:
        with serial.Serial(port, baudrate, timeout=2) as ser:
            time.sleep(2)
            print("✅ Capteur détecté ! Envoi du signal d'initialisation...")
            ser.write(b'START\n')
            print(f"Acquisition de {N_SAMPLES} échantillons en cours (~10 secondes)...")

            while len(donnees_brutes) < N_SAMPLES:
                if ser.in_waiting > 0:
                    ligne = ser.readline().decode('utf-8').strip()
                    if ligne:
                        try:
                            valeur = float(ligne)
                            donnees_brutes.append(valeur)
                        except ValueError:
                            pass # On ignore silencieusement les artefacts

            ser.write(b'STOP\n')
            print("Acquisition terminée avec succès.")

    except serial.SerialException as e:
        print("\n❌ [ERREUR CRITIQUE] Impossible de communiquer avec le capteur !")
        print(f"Détails techniques : {e}")
        print("\n--- Pistes de résolution ---")
        print("1. Vérifiez que la carte ESP32 est bien branchée en USB.")
        print("2. Vérifiez le nom du port. Sous Linux, essayez port='/dev/ttyACM0' ou port='/dev/ttyUSB1'.")
        print("3. Assurez-vous d'avoir redémarré votre PC après avoir configuré le groupe 'dialout'.")
        print("----------------------------\n")
        return None  # On retourne None pour signaler l'échec

    return np.array(donnees_brutes)


def convertir_en_debit(donnees_adc):
    """
    Calibre autour de zéro, convertit la tension en pression (Pa),
    puis en débit volumétrique (m^3/s) via Bernoulli.
    """
    # 1. Conversion ADC (0-4095) vers Tension avec mise à l'échelle
    tension = (donnees_adc / 4095.0) * 3.3 * FACTEUR_ECHELLE
    
    # 2. Calibrage (soustraction de la moyenne des 100 premiers points de repos)
    tension_calibree = tension - np.mean(tension[:100])
    tension_calibree[tension_calibree < 0] = 0 
    
    # 3. Conversion Tension -> Pression Différentielle (Pa)
    # Le MXP5010DP mesure de 0 à 10 kPa (10 000 Pa).
    # Sensibilité type : 45 mV / kPa (avec une alimentation 5V)
    # P (kPa) = (Vout / Vs - 0.04) / 0.09
    pression_kpa = (tension_calibree / 5.0 - 0.04) / 0.09
    pression_pa = pression_kpa * 1000.0
    pression_pa[pression_pa < 0] = 0 # Sécurité
    
    # 4. Équation de Bernoulli pour le débit (Q)
    # Q = A1 * sqrt( (2 * dp) / (RHO * (1 - (A1/A2)**2)) )
    denominateur = RHO * (1 - (A1 / A2)**2)
    # On utilise np.abs pour éviter les erreurs de racine carrée sur des micro-fluctuations négatives
    debit = A1 * np.sqrt((2 * np.abs(pression_pa)) / np.abs(denominateur))
    
    return debit


def compenser_perte_charge(debit):
    """
    Calcule le nombre de Reynolds et ajuste le débit pour compenser
    la perte de pression due à la friction.
    """
    debit_compense = np.zeros_like(debit)
    rugosite_relative = RUGOSITE / D_ENTREE
    
    for i, q in enumerate(debit):
        if q <= 0:
            continue
            
        # Vitesse de l'air (m/s)
        v = q / A1
        
        # Nombre de Reynolds
        Re = (RHO * v * D_ENTREE) / MU
        
        # Facteur de friction via Colebrook
        f = resoudre_colebrook(Re, rugosite_relative)
        
        # Calcul de la chute de pression supplémentaire (Darcy-Weisbach)
        # Supposons une longueur effective L de la section convergente (ex: 0.045 m)
        L = 0.045 
        dp_friction = f * (L / D_ENTREE) * (RHO * v**2) / 2.0
        
        # Le vrai débit est légèrement supérieur car la pression mesurée 
        # a subi cette perte. On rajoute cette pression et on recalcule.
        # (Version simplifiée de l'intégration pour maintenir de bonnes performances)
        pression_initiale_estimee = (q / A1)**2 * (RHO * np.abs(1 - (A1 / A2)**2)) / 2.0
        pression_reelle = pression_initiale_estimee + dp_friction
        
        debit_compense[i] = A1 * np.sqrt((2 * pression_reelle) / np.abs(RHO * (1 - (A1 / A2)**2)))
        
    return debit_compense


def calculer_frequence_coupure_fft(signal_debit, fs):
    """
    Calcule la fréquence de coupure dynamique en utilisant la FFT.
    Trouve la fréquence qui contient 90% de la puissance du signal.
    """
    N = len(signal_debit)
    
    # 1. Centrer le signal pour éviter un pic continu (DC) énorme à 0 Hz
    signal_centre = signal_debit - np.mean(signal_debit)
    
    # 2. Calcul de la FFT
    fft_valeurs = np.fft.fft(signal_centre)
    fft_freqs = np.fft.fftfreq(N, 1/fs)
    
    # 3. Prendre seulement le spectre unilatéral (fréquences positives)
    demi_N = N // 2
    fft_valeurs = fft_valeurs[:demi_N]
    fft_freqs = fft_freqs[:demi_N]
    
    # 4. Calcul de la puissance (amplitude au carré)
    puissance = np.abs(fft_valeurs) ** 2
    
    # 5. Puissance cumulative
    puissance_cumulative = np.cumsum(puissance)
    puissance_totale = puissance_cumulative[-1]
    
    # 6. Trouver l'indice où la puissance cumulative atteint 90%
    seuil_90 = 0.90 * puissance_totale
    indice_90 = np.argmax(puissance_cumulative >= seuil_90)
    
    fc_dynamique = fft_freqs[indice_90]
    
    # Sécurité : limiter la fréquence de coupure à des valeurs physiologiques
    if fc_dynamique < 1.0:
        fc_dynamique = 1.0
    elif fc_dynamique > 50.0:
        fc_dynamique = 50.0
        
    return fc_dynamique


def filtrer_signal(debit):
    """
    Identifie dynamiquement la fréquence de coupure via FFT,
    puis applique un filtre passe-bas de Butterworth du 4ème ordre.
    """
    # 1. Identifier la fréquence de coupure dynamique (90% de la puissance via FFT)
    fc = calculer_frequence_coupure_fft(debit, FE)
    print(f"Fréquence de coupure calculée via FFT : {fc:.2f} Hz")
    
    # 2. Appliquer le filtre de Butterworth avec filtrage à phase nulle
    sos = signal.butter(4, fc, 'low', fs=FE, output='sos')
    debit_filtre = signal.sosfiltfilt(sos, debit) 
    
    return debit_filtre


def extraire_parametres_spirometrie(debit_filtre):
    """
    Identifie le PEF et intègre le signal pour trouver FVC et FEV1.
    Convertit les résultats finaux en Litres (L) et Litres par minute (L/min).
    """
    # 1. Identifier le Débit Expiratoire de Pointe (PEF) en m^3/s
    pef = np.max(debit_filtre)
    index_pef = np.argmax(debit_filtre)

    # 2. Définir la durée de l'exhalation (seuil de 4% du maximum)
    seuil = 0.04 * pef

    # Recherche de l'index de début en itérant vers l'arrière
    index_debut = index_pef
    while index_debut > 0 and debit_filtre[index_debut] > seuil:
        index_debut -= 1

    # Recherche de l'index de fin en itérant vers l'avant
    index_fin = index_pef
    while index_fin < len(debit_filtre) - 1 and debit_filtre[index_fin] > seuil:
        index_fin += 1

    # 3. Calculer la Capacité Vitale Forcée (FVC) via la méthode des trapèzes
    # CORRECTION ICI : np.trapz devient np.trapezoid
    fvc_m3 = np.trapezoid(debit_filtre[index_debut:index_fin], dx=1 / FE)

    # 4. Calculer le Volume Expiratoire Maximal Seconde (FEV1)
    index_1s = min(index_debut + int(1.0 * FE), index_fin)
    # CORRECTION ICI : np.trapz devient np.trapezoid
    fev1_m3 = np.trapezoid(debit_filtre[index_debut:index_1s], dx=1 / FE)

    # 5. Conversions d'unités (m^3 vers Litres, et m^3/s vers L/min)
    fvc_litres = fvc_m3 * 1000.0
    fev1_litres = fev1_m3 * 1000.0
    pef_lmin = pef * 60000.0

    return fvc_litres, fev1_litres, pef_lmin

def resoudre_colebrook(Re, rugosite_relative):
    """
    Résout l'équation de Colebrook-White pour trouver le facteur de friction (f)
    en utilisant la méthode de dichotomie (bissection).
    """
    if Re < 2000:
        return 64.0 / Re  # Flux laminaire
        
    # Limites pour la bissection
    f_min = 0.008
    f_max = 0.1
    tolerance = 1e-6
    
    for _ in range(100):
        f_mid = (f_min + f_max) / 2.0
        # Équation : 1/sqrt(f) + 2*log10( (e/D)/3.7 + 2.51/(Re*sqrt(f)) ) = 0
        terme_gauche = 1.0 / np.sqrt(f_mid)
        terme_droite = -2.0 * np.log10((rugosite_relative / 3.7) + (2.51 / (Re * np.sqrt(f_mid))))
        
        erreur = terme_gauche - terme_droite
        
        if abs(erreur) < tolerance:
            return f_mid
            
        if erreur > 0:
            f_min = f_mid
        else:
            f_max = f_mid
            
    return (f_min + f_max) / 2.0



def generer_donnees_simulees():
    """
    Génère 10 000 points de données simulant une expiration physiologique
    pour tester l'algorithme sans le capteur physique MXP5010DP.
    """
    print("Génération d'un signal physiologique simulé...")
    donnees = np.ones(N_SAMPLES) * 250  # Ligne de base (bruit de fond de l'ADC)
    
    # Ajout d'un bruit blanc aléatoire pour simuler les imperfections du capteur
    bruit = np.random.normal(0, 5, N_SAMPLES)
    donnees += bruit
    
    # Simulation d'un souffle qui commence à la 2ème seconde (index 2000)
    debut_souffle = 2000
    duree_souffle = 5000  # Le souffle dure 5 secondes
    
    for i in range(duree_souffle):
        t = i / FE  # Temps écoulé en secondes depuis le début du souffle
        # Modèle mathématique d'un souffle : P(t) = A * t * exp(-k*t)
        # Ces constantes génèrent une courbe qui ressemble à un vrai signal de spirométrie
        valeur_souffle = 12000 * t * np.exp(-4 * t) 
        
        # On ajoute le souffle à la ligne de base
        if debut_souffle + i < N_SAMPLES:
            donnees[debut_souffle + i] += valeur_souffle
            
    # L'ADC de l'ESP32 est limité entre 0 et 4095
    donnees = np.clip(donnees, 0, 4095)
    
    return donnees


def calculer_fvc_theorique(age, taille_cm, sexe):
    """
    Calcule une estimation simplifiée du FVC théorique attendu (en Litres).
    (Basé sur des équations spirométriques standards simplifiées de l'ERS/ATS)
    """
    taille_m = taille_cm / 100.0
    if sexe.lower() == 'homme':
        # Formule simplifiée pour homme
        fvc_pred = (5.76 * taille_m) - (0.026 * age) - 4.34
    elif sexe.lower() == 'femme':
        # Formule simplifiée pour femme
        fvc_pred = (4.43 * taille_m) - (0.026 * age) - 2.89
    else:
        fvc_pred = None
    
    return fvc_pred



def diagnostiquer_patient(fvc, fev1, fvc_pred=None):
    """
    Génère un diagnostic clinique préliminaire basé sur l'indice de Tiffeneau
    et la comparaison avec le FVC théorique.
    """
    if fvc <= 0:
        return "Erreur : Capacité Vitale (FVC) invalide."

    ratio = (fev1 / fvc) * 100

    print("\n--- ANALYSE CLINIQUE ---")
    print(f"Indice de Tiffeneau (FEV1/FVC) : {ratio:.1f} %")
    
    # 1. Test de l'obstruction
    if ratio < 70.0:
        diagnostic = "⚠️ PROFIL OBSTRUCTIF SUSPECTÉ (ex: Asthme, BPCO)."
        # On peut affiner la sévérité avec le FEV1 par rapport à un FEV1 prédit
    else:
        # 2. Test de la restriction (nécessite la valeur théorique)
        if fvc_pred:
            pourcentage_fvc = (fvc / fvc_pred) * 100
            print(f"FVC mesuré vs théorique : {pourcentage_fvc:.1f} %")
            
            if pourcentage_fvc < 80.0:
                diagnostic = "⚠️ PROFIL RESTRICTIF SUSPECTÉ (Capacité pulmonaire réduite)."
            else:
                diagnostic = "✅ PROFIL NORMAL (Aucune anomalie respiratoire majeure détectée)."
        else:
            diagnostic = "❕ PROFIL NORMAL/INDÉTERMINÉ (Ratio normal, mais FVC théorique manquant pour écarter une restriction)."

    return diagnostic


def calculer_fvc_theorique(age, taille_cm, sexe):
    """
    Calcule une estimation du FVC théorique attendu (en Litres).
    Basé sur des équations spirométriques de référence.
    """
    taille_m = taille_cm / 100.0
    if sexe.lower() == 'homme':
        fvc_pred = (5.76 * taille_m) - (0.026 * age) - 4.34
    elif sexe.lower() == 'femme':
        fvc_pred = (4.43 * taille_m) - (0.026 * age) - 2.89
    else:
        fvc_pred = None
    
    return fvc_pred

def diagnostiquer_patient(fvc, fev1, fvc_pred=None):
    """
    Génère un diagnostic clinique préliminaire basé sur l'indice de Tiffeneau
    et la comparaison avec le FVC théorique.
    """
    if fvc <= 0:
        return "Erreur : Capacité Vitale (FVC) invalide."

    ratio = (fev1 / fvc) * 100

    print("\n--- ANALYSE CLINIQUE PRÉLIMINAIRE ---")
    
    # 1. Test du syndrome obstructif
    if ratio < 70.0:
        diagnostic = "⚠️ PROFIL OBSTRUCTIF SUSPECTÉ (ex: Asthme, BPCO). Le ratio FEV1/FVC est inférieur à 70%."
    else:
        # 2. Test du syndrome restrictif (nécessite la valeur théorique)
        if fvc_pred:
            pourcentage_fvc = (fvc / fvc_pred) * 100
            print(f"FVC mesuré vs théorique : {pourcentage_fvc:.1f} %")
            
            if pourcentage_fvc < 80.0:
                diagnostic = "⚠️ PROFIL RESTRICTIF SUSPECTÉ. Ratio normal, mais Capacité Pulmonaire (FVC) réduite."
            else:
                diagnostic = "✅ PROFIL NORMAL. Aucune anomalie respiratoire majeure détectée."
        else:
            diagnostic = "❕ PROFIL INDÉTERMINÉ. Ratio normal, mais manque le FVC théorique pour écarter une restriction."

    return diagnostic


def demander_informations_patient():
    """
    Demande à l'utilisateur de saisir les paramètres physiologiques du patient
    avec des vérifications de sécurité pour éviter les erreurs de frappe.
    """
    print("\n--- INFORMATIONS DU PATIENT ---")
    
    while True:
        try:
            age = int(input("Entrez l'âge du patient (en années) : "))
            if age > 0:
                break
            else:
                print("L'âge doit être supérieur à 0.")
        except ValueError:
            print("Erreur : Veuillez entrer un nombre entier valide.")
            
    while True:
        try:
            taille_cm = float(input("Entrez la taille du patient (en cm) : "))
            if taille_cm > 50:
                break
            else:
                print("La taille semble invalide (doit être > 50 cm).")
        except ValueError:
            print("Erreur : Veuillez entrer un nombre valide.")
            
    while True:
        sexe = input("Entrez le sexe biologique du patient (homme/femme) : ").strip().lower()
        if sexe in ['homme', 'femme']:
            break
        print("Erreur : Veuillez taper 'homme' ou 'femme'.")
        
    return age, taille_cm, sexe

# --- FLUX D'EXÉCUTION PRINCIPAL ---
if __name__ == '__main__':
    print("Démarrage du système Spiromètre...")
    
    # 1. Saisie des données du patient avant de lancer l'acquisition
    age_patient, taille_patient, sexe_patient = demander_informations_patient()
    
    # 2. Acquisition (simulée ou matérielle)
    donnees = acquerir_donnees_serie(port='/dev/ttyUSB0')
    # donnees = generer_donnees_simulees() # pour les tests sans matériel physique
    
    if np.any(donnees):
        print("\nTraitement des données en cours...")
        debit = convertir_en_debit(donnees)
        debit_compense = compenser_perte_charge(debit)
        debit_filtre = filtrer_signal(debit_compense)
        fvc, fev1, pef = extraire_parametres_spirometrie(debit_filtre)
        
        print("\n--- RÉSULTATS DE LA SPIROMÉTRIE ---")
        print(f"FVC (Capacité Vitale Forcée) : {fvc:.2f} L")
        print(f"FEV1 (Volume Expiré en 1s)   : {fev1:.2f} L")
        print(f"PEF (Débit de Pointe)        : {pef:.2f} L/min")
        print(f"Ratio FEV1/FVC               : {(fev1/fvc)*100:.1f} %" if fvc > 0 else "Ratio FEV1/FVC: N/A")
        
        # 3. Diagnostic clinique avec les variables saisies
        fvc_theorique = calculer_fvc_theorique(age=age_patient, taille_cm=taille_patient, sexe=sexe_patient)
        resultat = diagnostiquer_patient(fvc, fev1, fvc_pred=fvc_theorique)
        print(resultat)