#include "UAR.h"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Model ARX

ModelARX::ModelARX() : gen(std::random_device{}()) {
    m_A = {-0.5, 0.0, 0.0};
    m_B = {0.5, 0.0, 0.0};
    m_k = 1;
    // Inicjalizacja buforów zerami
    m_historia_u.resize(20, 0.0);
    m_historia_y.resize(20, 0.0);
}

void ModelARX::setParams(const std::vector<double>& wsp_a, const std::vector<double>& wsp_b, unsigned int opoznienie_k) {
    m_A = wsp_a;
    m_B = wsp_b;
    m_k = (opoznienie_k < 1) ? 1 : opoznienie_k;

    // Rozszerz bufory jeśli potrzeba (z zapasem)
    size_t required = m_k + m_B.size() + m_A.size() + 5;
    if (m_historia_u.size() < required) {
        m_historia_u.resize(required, 0.0);
        m_historia_y.resize(required, 0.0);
    }
}

void ModelARX::setLimity(double minU, double maxU, double minY, double maxY, bool wlaczone) {
    m_minU = minU; m_maxU = maxU;
    m_minY = minY; m_maxY = maxY;

    m_ogranicz_u = wlaczone;
    m_ogranicz_y = wlaczone;
}

void ModelARX::setSzum(double std_dev) {
    m_szum = std_dev;
    if (m_szum > 0.0) dist = std::normal_distribution<double>(0.0, m_szum);
}

double ModelARX::symuluj(double u_raw) {
    // Nasycenie wejścia przed obliczeniami
    double u = u_raw;
    if (m_ogranicz_u) {
        if (u > m_maxU) u = m_maxU;
        if (u < m_minU) u = m_minU;
    }

    // Aktualizacja historii sterowania
    m_historia_u.push_front(u);
    if (m_historia_u.size() > 100) m_historia_u.pop_back();

    // Obliczenia ARX
    double y_temp = 0.0;

    // Część od sterowania (B)
    for (size_t i = 0; i < m_B.size(); ++i) {
        size_t idx = m_k + i;
        if (idx < m_historia_u.size()) {
            y_temp += m_B[i] * m_historia_u[idx];
        }
    }

    // Część od wyjścia (A) - odejmujemy
    for (size_t i = 0; i < m_A.size(); ++i) {
        size_t idx = i;
        if (idx < m_historia_y.size()) {
            y_temp -= m_A[i] * m_historia_y[idx];
        }
    }

    // Dodanie szumu
    if (m_szum > 0.0001) {
        y_temp += dist(gen);
    }

    // Nasycenie wyjścia
    if (m_ogranicz_y) {
        if (y_temp > m_maxY) y_temp = m_maxY;
        if (y_temp < m_minY) y_temp = m_minY;
    }

    // Aktualizacja historii wyjścia
    m_historia_y.push_front(y_temp);
    if (m_historia_y.size() > 100) m_historia_y.pop_back();

    return y_temp;
}

void ModelARX::reset() {
    std::fill(m_historia_u.begin(), m_historia_u.end(), 0.0);
    std::fill(m_historia_y.begin(), m_historia_y.end(), 0.0);
}

// Regulator PID

RegulatorPID::RegulatorPID() {}

void RegulatorPID::setNastawy(double k, double Ti, double Td, LiczCalk tryb) {
    m_k = k;
    m_Ti = Ti;
    m_Td = Td;

    if (tryb != m_liczCalk)
    {
        if (tryb == LiczCalk::Wew)
            m_suma_e = m_suma_e / m_Ti * 1.0;
        else
            m_suma_e = m_suma_e * m_Ti * 1.0;
    }


    m_liczCalk = tryb;
}


double RegulatorPID::symuluj(double e) {
    // Proporcjonalna
    m_u_P = m_k * e * 1.0;

    double I_temp = 0.0;

    // Całkująca
    if (m_Ti == 0.0) {
        m_u_I = 0.0; // Człon wyłączony
    } else {
        if (m_liczCalk == LiczCalk::Wew) { // Stała PRZED sumą
            m_suma_e += e / m_Ti;
            I_temp = m_suma_e * 1.0;
        } else { // Stała W sumie (Pod całką)
            m_suma_e += e;  // / m_Ti;
            I_temp = m_suma_e / m_Ti * 1.0;  // * m_k;
        }
    }

    m_u_I = m_k * I_temp;

    // Różniczkująca
    m_u_D = m_Td * (e - m_prev_e) * m_k;
    m_prev_e = e;

    return m_u_P + m_u_I + m_u_D;
}

void RegulatorPID::reset() {
    m_u_P = m_u_I = m_u_D = 0.0;
    m_suma_e = 0.0;
    m_prev_e = 0.0;
}

void RegulatorPID::resetMemory() {
    m_suma_e = 0.0;
    m_prev_e = 0.0;

    m_u_I = 0.0;
}

// Generator

GeneratorWartosci::GeneratorWartosci() {}

void GeneratorWartosci::setParams(TrybGen tryb, double okres_rzecz, double ampl, double off, double wyp, int interwal_ms) {
    m_tryb = tryb;
    m_T_RZ = okres_rzecz;
    m_A = ampl;
    m_S = off;
    m_p = wyp;
    m_T_T = interwal_ms;
    aktualizujT();
}

void GeneratorWartosci::aktualizujT() {
    // Przeliczenie okresu w sekundach na okres w próbkach
    if (m_T_T <= 0) m_T_T = 200;
    double probek_na_sekunde = m_T_RZ * 1000.0 / m_T_T;
    m_T_probki = static_cast<int>(probek_na_sekunde);
    if (m_T_probki < 1) m_T_probki = 1;
}

double GeneratorWartosci::generuj() {
    // Obliczenie fazy sygnału
    int faza = m_i % m_T_probki;
    double val = 0.0;

    if (m_tryb == TrybGen::Sin) {
        double ratio = (double)faza / (double)m_T_probki;
        val = m_A * std::sin(ratio * 2.0 * M_PI) + m_S;
    } else {
        // Prostokąt
        double limit = m_p * m_T_probki;
        if (faza < limit) val = m_A + m_S;
        else val = -m_A + m_S;
    }

    m_w_i = val;
    m_i++;
    return val;
}

void GeneratorWartosci::reset() {
    m_i = 0;
    m_w_i = 0.0;
}

// UAR

ProstyUAR::ProstyUAR() {}

double ProstyUAR::symuluj() {
    // Wyznacz wartość zadaną
    double w = m_genWart.generuj();

    // Oblicz uchyb (wartość zadana - poprzednie wyjście obiektu)
    // Sprzężenie zwrotne bierze y z poprzedniego kroku
    m_e_i = w - m_y_i;

    // Regulator PID
    m_u_i = m_PID.symuluj(m_e_i);

    // ARX
    m_y_i = m_ARX.symuluj(m_u_i);

    return m_y_i;
}

void ProstyUAR::reset() {
    m_ARX.reset();
    m_PID.reset();
    m_genWart.reset();
    m_y_i = 0.0;
    m_e_i = 0.0;
    m_u_i = 0.0;
}

void ProstyUAR::resetPID()  {
    m_PID.resetMemory();
}
