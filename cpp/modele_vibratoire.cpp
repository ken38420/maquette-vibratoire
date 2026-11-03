// Maquette 2 ddl : caisson sur silentblocs + carte électronique — version C++.
//
// Même modèle et mêmes équations que src/maquette_vibratoire/modele.py :
//   M z'' + C z' + K z = -M·1·ÿ           (z = x - y, déplacement relatif à la base)
//   Z(ω) = (K - ω²M + iωC)⁻¹ (-M·1)       → H_rel [m / (m/s²)]
//   T(ω) = 1 - ω² H_rel                   (transmissibilité en accélération absolue)
//   Grms = sqrt(∫ |T|² S_in df)
//
// Paramètres : mêmes valeurs ILLUSTRATIVES que config.json (lues en ligne de commande
// pour rester sans dépendance : pas de bibliothèque JSON ni d'algèbre linéaire).
// Unités : SI en interne, DSP en g²/Hz.
//
// Compilation : g++ -std=c++17 -O2 -Wall -Wextra -o modele_vibratoire modele_vibratoire.cpp

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double G_SI = 9.80665;  // m/s², valeur normalisée (identique à scipy.constants.g)

using Mat2 = std::array<std::array<double, 2>, 2>;
using cplx = std::complex<double>;

struct Systeme {
    Mat2 M, K, C;
};

struct Modes {
    std::array<double, 2> f_Hz;
    Mat2 phi;  // phi[ddl][mode], normalisé max |φ| = 1, composante caisson > 0
};

struct DspTrapeze {
    double f_min, f_a, f_b, f_max, niveau_g2_Hz, pente_montee_dB_oct, pente_descente_dB_oct;
};

// Même convention que matrices_2ddl() en Python : k1 fixé par (m1+m2) à f_susp,
// k2 par m2 à f_carte, c_i = 2 ζ_i sqrt(k_i m_i).
Systeme matrices_2ddl(double m1_kg, double m2_kg, double f_susp_Hz, double zeta_susp,
                      double f_carte_Hz, double zeta_carte) {
    if (m1_kg <= 0 || m2_kg <= 0 || f_susp_Hz <= 0 || f_carte_Hz <= 0)
        throw std::invalid_argument("masses et fréquences doivent être > 0");
    if (!(zeta_susp > 0 && zeta_susp < 1) || !(zeta_carte > 0 && zeta_carte < 1))
        throw std::invalid_argument("amortissements doivent être dans ]0, 1[");

    const double m_susp = m1_kg + m2_kg;
    const double k1 = m_susp * std::pow(2 * PI * f_susp_Hz, 2);
    const double k2 = m2_kg * std::pow(2 * PI * f_carte_Hz, 2);
    const double c1 = 2 * zeta_susp * std::sqrt(k1 * m_susp);
    const double c2 = 2 * zeta_carte * std::sqrt(k2 * m2_kg);

    Systeme s;
    s.M = {{{m1_kg, 0.0}, {0.0, m2_kg}}};
    s.K = {{{k1 + k2, -k2}, {-k2, k2}}};
    s.C = {{{c1 + c2, -c2}, {-c2, c2}}};
    return s;
}

// Problème aux valeurs propres généralisé 2x2 (M diagonale) : det(K - λM) = 0
//   m1 m2 λ² - (m1 k22 + m2 k11) λ + (k11 k22 - k12²) = 0
// Racines par la forme stable (évite l'annulation catastrophique de -b ± sqrt(Δ)).
Modes modes_propres(const Systeme& s) {
    const double m1 = s.M[0][0], m2 = s.M[1][1];
    const double k11 = s.K[0][0], k12 = s.K[0][1], k22 = s.K[1][1];
    const double a = m1 * m2;
    const double b = -(m1 * k22 + m2 * k11);
    const double c = k11 * k22 - k12 * k12;
    const double delta = b * b - 4 * a * c;
    if (delta < 0) throw std::runtime_error("discriminant négatif : matrices non symétriques définies ?");
    const double q = -0.5 * (b + std::copysign(std::sqrt(delta), b));
    double lam1 = c / q, lam2 = q / a;
    if (lam1 > lam2) std::swap(lam1, lam2);
    if (lam1 <= 0) throw std::runtime_error("valeur propre non positive");

    Modes m;
    const std::array<double, 2> lam = {lam1, lam2};
    for (int j = 0; j < 2; ++j) {
        m.f_Hz[j] = std::sqrt(lam[j]) / (2 * PI);
        // (K - λM) φ = 0, 1re ligne : (k11 - λ m1) φ1 + k12 φ2 = 0 → φ = [-k12, k11 - λ m1]
        double p1 = -k12, p2 = k11 - lam[j] * m1;
        const double norme = std::max(std::fabs(p1), std::fabs(p2));
        p1 /= norme;
        p2 /= norme;
        if (p1 < 0) { p1 = -p1; p2 = -p2; }
        m.phi[0][j] = p1;
        m.phi[1][j] = p2;
    }
    return m;
}

// H_rel(ω) : résolution du système 2x2 complexe par la règle de Cramer.
std::array<cplx, 2> frf_excitation_base(const Systeme& s, double f_Hz) {
    const double w = 2 * PI * f_Hz;
    cplx A[2][2];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            A[i][j] = cplx(s.K[i][j] - w * w * s.M[i][j], w * s.C[i][j]);
    const cplx r0 = -(s.M[0][0] + s.M[0][1]);  // -M·1
    const cplx r1 = -(s.M[1][0] + s.M[1][1]);
    const cplx det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (std::abs(det) == 0.0) throw std::runtime_error("matrice dynamique singulière");
    return {(r0 * A[1][1] - A[0][1] * r1) / det, (A[0][0] * r1 - r0 * A[1][0]) / det};
}

double dsp_trapeze(double f, const DspTrapeze& d) {
    const double n_montee = d.pente_montee_dB_oct / (10 * std::log10(2.0));
    const double n_descente = d.pente_descente_dB_oct / (10 * std::log10(2.0));
    if (f < d.f_min || f > d.f_max) return 0.0;
    if (f < d.f_a) return d.niveau_g2_Hz * std::pow(f / d.f_a, n_montee);
    if (f <= d.f_b) return d.niveau_g2_Hz;
    return d.niveau_g2_Hz * std::pow(f / d.f_b, n_descente);
}

std::vector<double> grille_log(double f_min, double f_max, int n) {
    std::vector<double> f(n);
    const double lmin = std::log(f_min), lmax = std::log(f_max);
    for (int i = 0; i < n; ++i) f[i] = std::exp(lmin + (lmax - lmin) * i / (n - 1));
    return f;
}

double trapeze(const std::vector<double>& y, const std::vector<double>& x) {
    if (x.size() != y.size()) throw std::logic_error("tailles différentes");
    double s = 0.0;
    for (size_t i = 1; i < x.size(); ++i) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    return s;
}

struct Reponse {
    double grms_entree, grms_caisson, grms_carte, drel_rms_mm, T_max_caisson, T_max_carte;
};

Reponse analyse_aleatoire(const Systeme& s, const std::vector<double>& f, const DspTrapeze& d) {
    std::vector<double> S_in(f.size()), S1(f.size()), S2(f.size()), Sd(f.size());
    Reponse r{};
    for (size_t i = 0; i < f.size(); ++i) {
        const auto H = frf_excitation_base(s, f[i]);
        const double w2 = std::pow(2 * PI * f[i], 2);
        const cplx T1 = 1.0 - w2 * H[0];
        const cplx T2 = 1.0 - w2 * H[1];
        S_in[i] = dsp_trapeze(f[i], d);
        S1[i] = std::norm(T1) * S_in[i];
        S2[i] = std::norm(T2) * S_in[i];
        Sd[i] = std::norm(H[1] - H[0]) * S_in[i] * G_SI * G_SI;  // m²/Hz
        r.T_max_caisson = std::max(r.T_max_caisson, std::abs(T1));
        r.T_max_carte = std::max(r.T_max_carte, std::abs(T2));
    }
    r.grms_entree = std::sqrt(trapeze(S_in, f));
    r.grms_caisson = std::sqrt(trapeze(S1, f));
    r.grms_carte = std::sqrt(trapeze(S2, f));
    r.drel_rms_mm = std::sqrt(trapeze(Sd, f)) * 1e3;
    return r;
}

// Validation : 1 ddl sous bruit blanc, Grms² = S0 π fn (1+4ζ²)/(4ζ) (exact).
// Transmissibilité 1 ddl en forme adimensionnelle r = f/fn : T = (1 + 2iζr) / (1 - r² + 2iζr).
// Teste la grille logarithmique et l'intégration trapèze utilisées par le modèle 2 ddl.
void valider_1ddl_bruit_blanc() {
    const double fn = 100.0, zeta = 0.05, S0 = 0.01;
    const auto f = grille_log(0.01, 1e5, 200000);
    std::vector<double> S(f.size());
    for (size_t i = 0; i < f.size(); ++i) {
        const double r = f[i] / fn;
        const cplx T = cplx(1.0, 2 * zeta * r) / cplx(1.0 - r * r, 2 * zeta * r);
        S[i] = std::norm(T) * S0;
    }
    const double num = std::sqrt(trapeze(S, f));
    const double ana = std::sqrt(S0 * PI * fn * (1 + 4 * zeta * zeta) / (4 * zeta));
    const double ecart = std::fabs(num - ana) / ana;
    if (ecart > 5e-3) throw std::runtime_error("validation 1 ddl échouée");
    std::printf("  1 ddl bruit blanc : Grms num %.4f / analytique %.4f (écart %.1e) : OK\n", num, ana, ecart);
}

double arg(char** argv, int i) {
    char* fin = nullptr;
    const double v = std::strtod(argv[i], &fin);
    if (fin == argv[i] || *fin != '\0') throw std::invalid_argument(std::string("nombre invalide : ") + argv[i]);
    return v;
}

}  // namespace

int calcul(int argc, char** argv);

int main(int argc, char** argv) {
    // Toute erreur arrête le programme avec un message explicite et un code de sortie non nul.
    try {
        return calcul(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ERREUR : %s\n", e.what());
        return 1;
    }
}

int calcul(int argc, char** argv) {
    // Ordre : m1 m2 f_susp zeta_susp f_carte zeta_carte
    if (argc != 7) {
        std::fprintf(stderr,
                     "Usage : %s m1_kg m2_kg f_susp_Hz zeta_susp f_carte_Hz zeta_carte\n"
                     "Exemple (config.json) : %s 25 0.3 25 0.10 200 0.03\n",
                     argv[0], argv[0]);
        return 1;
    }
    const Systeme s = matrices_2ddl(arg(argv, 1), arg(argv, 2), arg(argv, 3), arg(argv, 4), arg(argv, 5), arg(argv, 6));

    // Profil de DSP et grille : mêmes valeurs que config.json (dsp_entree, numerique).
    const DspTrapeze dsp{20.0, 80.0, 350.0, 2000.0, 0.04, 3.0, -3.0};
    const auto f = grille_log(1.0, 3000.0, 40000);

    const Modes m = modes_propres(s);
    const Reponse r = analyse_aleatoire(s, f, dsp);

    std::printf("=== Modes propres (système couplé, non amorti)\n");
    for (int j = 0; j < 2; ++j)
        std::printf("  Mode %d : %7.2f Hz   déformée caisson/carte = %+.3f / %+.3f\n", j + 1, m.f_Hz[j],
                    m.phi[0][j], m.phi[1][j]);
    std::printf("=== Transmissibilité maximale\n  Caisson : %.1f   Carte : %.1f\n", r.T_max_caisson, r.T_max_carte);
    std::printf("=== Réponse aléatoire\n");
    std::printf("  Entrée  : %.2f Grms\n  Caisson : %.2f Grms\n  Carte   : %.2f Grms\n", r.grms_entree,
                r.grms_caisson, r.grms_carte);
    std::printf("  Déplacement relatif carte 3σ : %.3f mm\n", 3 * r.drel_rms_mm);
    std::printf("=== Validations\n");
    valider_1ddl_bruit_blanc();
    return 0;
}
