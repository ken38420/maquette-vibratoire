// Maquette vibratoire en C++ : même modèle et mêmes données que maquette.py.
// Caisson (m1) sur silentblocs + carte (m2), base secouée par une DSP aléatoire.
//
// Compiler : g++ -O2 -o maquette.exe maquette.cpp
// Lancer   : maquette.exe

#include <cmath>
#include <complex>
#include <cstdio>

using cplx = std::complex<double>;
const double PI = 3.14159265358979;
const double g = 9.81;

// ---------- Données d'entrée (valeurs illustratives) ----------
const double m1 = 25.0, f1 = 25.0, z1 = 0.10;   // caisson : masse [kg], fréquence [Hz], amortissement
const double m2 = 0.25, f2 = 200.0, z2 = 0.03;  // carte
const double profil_f[4] = {20, 80, 350, 2000}; // profil d'entrée : fréquences [Hz]
const double profil_S[4] = {0.01, 0.04, 0.04, 0.007};  // DSP [g²/Hz]
const double jeu_mm = 5.0;                      // débattement disponible [mm]
const double B = 160.0, h = 2.0, L = 25.0;      // carte pour Steinberg [mm]
const double FACTEUR_CENTRE = 16 / (PI * PI);   // flexion au centre de la plaque / flexion de la masse

// DSP d'entrée : droite entre deux points du profil sur un graphe log-log
double dsp(double f) {
    for (int i = 0; i < 3; i++)
        if (f >= profil_f[i] && f <= profil_f[i + 1]) {
            double pente = std::log(profil_S[i + 1] / profil_S[i]) / std::log(profil_f[i + 1] / profil_f[i]);
            return profil_S[i] * std::pow(f / profil_f[i], pente);
        }
    return 0.0;
}

int main() {
    // Raideurs et amortisseurs
    double k1 = (m1 + m2) * std::pow(2 * PI * f1, 2), k2 = m2 * std::pow(2 * PI * f2, 2);
    double c1 = 2 * z1 * std::sqrt(k1 * (m1 + m2)), c2 = 2 * z2 * std::sqrt(k2 * m2);

    // Modes propres : det(K - λM) = 0 donne  a λ² + b λ + c = 0
    double a = m1 * m2, b = -(m1 * k2 + m2 * (k1 + k2)), c = (k1 + k2) * k2 - k2 * k2;
    double delta = std::sqrt(b * b - 4 * a * c);
    double fm1 = std::sqrt((-b - delta) / (2 * a)) / (2 * PI);
    double fm2 = std::sqrt((-b + delta) / (2 * a)) / (2 * PI);

    // Balayage en fréquence (grille logarithmique) et intégration des DSP par trapèzes
    const int n = 20000;
    double aire_in = 0, aire_caisson = 0, aire_carte = 0, aire_debat = 0, aire_drel = 0;
    double f_prec = 0, s_prec[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < n; i++) {
        double f = std::exp(std::log(1.0) + (std::log(3000.0) - std::log(1.0)) * i / (n - 1));
        double w = 2 * PI * f, S = dsp(f);

        // Système 2x2 : (K - ω²M + iωC) z = -M·1, résolu par la règle de Cramer
        cplx A11(k1 + k2 - w * w * m1, w * (c1 + c2)), A12(-k2, -w * c2);
        cplx A22(k2 - w * w * m2, w * c2), det = A11 * A22 - A12 * A12;
        cplx za = (-m1 * A22 + A12 * m2) / det;     // caisson / base
        cplx zb = (-A11 * m2 + m1 * A12) / det;     // carte / base
        cplx T1 = 1.0 - w * w * za, T2 = 1.0 - w * w * zb;

        double s[5] = {S, std::norm(T1) * S, std::norm(T2) * S,
                       std::norm(za) * S * g * g, std::norm(zb - za) * S * g * g};
        if (i > 0) {
            double df = f - f_prec;
            aire_in += 0.5 * (s[0] + s_prec[0]) * df;
            aire_caisson += 0.5 * (s[1] + s_prec[1]) * df;
            aire_carte += 0.5 * (s[2] + s_prec[2]) * df;
            aire_debat += 0.5 * (s[3] + s_prec[3]) * df;
            aire_drel += 0.5 * (s[4] + s_prec[4]) * df;
        }
        f_prec = f;
        for (int j = 0; j < 5; j++) s_prec[j] = s[j];
    }
    double debat = 4 * std::sqrt(aire_debat) * 1000;                       // débattement 4σ [mm]
    double drel = FACTEUR_CENTRE * 3 * std::sqrt(aire_drel) * 1000;        // flexion au centre 3σ [mm]

    // Critère de Steinberg (formule en pouces)
    double po = 25.4, z_adm = 0.00022 * (B / po) / ((h / po) * std::sqrt(L / po)) * po;

    std::printf("Modes propres   : %.1f Hz et %.1f Hz\n", fm1, fm2);
    std::printf("Entree          : %.2f Grms\n", std::sqrt(aire_in));
    std::printf("Caisson         : %.2f Grms\n", std::sqrt(aire_caisson));
    std::printf("Carte           : %.2f Grms\n", std::sqrt(aire_carte));
    std::printf("Marges de securite (MS >= 0 : conforme)\n");
    std::printf("  Separation des modes : %.2f\n", (fm2 / fm1) / 2.0 - 1);
    std::printf("  Debattement 4s       : %.3f mm / %.1f mm -> MS = %.2f\n", debat, jeu_mm, jeu_mm / debat - 1);
    std::printf("  Carte 3s (Steinberg) : %.3f mm / %.3f mm -> MS = %.2f\n", drel, z_adm, z_adm / drel - 1);
    return 0;
}
