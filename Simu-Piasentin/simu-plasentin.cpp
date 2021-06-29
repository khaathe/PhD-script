// Author : Nicola Piasentin
// Master Thesis Project
// The control of acidity intumour cells : a biophysical model
// GSL libraries needed

/*
    Compile command
    The program need to link the gsl library during compilation
    g++ simu-plasentin.cpp -o simu.exe `gsl-config  --cflags --libs`
    
    Or add those flag :
    -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

ofstream outData, outDatapH;

double m_CO2_C_old, m_H_C_old, m_HCO3_C_old, m_CO2_c_old, m_H_c_old, m_HCO3_c_old;
const double Pi = M_PI;

//sensors
const double SensO2 = 1.0;
const double SensATP = 1.0;

struct cell_params
{
    double MW_H;
    double MW_CO2;
    double MW_O2;
    double MW_HCO3;
    double MW_AcL;
    double r_C;
    double PM_CO2;
    double gAcL;
    double q_O2;
    double k1;
    double k2;
    double VMAXAcL;
    double K_mAcL;
    double a2cH_slope;
    double a2cH_thr;
    double c2aH_slope;
    double c2aH_thr;
    double VMAXNHE;

    double K_mNHE;
    double a;
    double l_NHE;
    double pH0_NHE;
    double VMAXTHCO3;
    double K_mTHCO3;
    double l_THCO3;
    double pHe0_THCO3;
    double g_THCO3;
    double pHi0_THCO3;
    double VMAXCA9;
    double K_mCA9;
    double d_CA9;
    double V_c;
    double dt;
};

int cell(const gsl_vector *x, void *params, gsl_vector *f)
{
    double MW_H = ((struct cell_params *)params)->MW_H;
    double MW_CO2 = ((struct cell_params *)params)->MW_CO2;
    double MW_O2 = ((struct cell_params *)params)->MW_O2;
    double MW_HCO3 = ((struct cell_params *)params)->MW_HCO3;
    double MW_AcL = ((struct cell_params *)params)->MW_AcL;
    double r_C = ((struct cell_params *)params)->r_C;
    double PM_CO2 = ((struct cell_params *)params)->PM_CO2;
    double gAcL = ((struct cell_params *)params)->gAcL;
    double q_O2 = ((struct cell_params *)params)->q_O2;
    double k1 = ((struct cell_params *)params)->k1;
    double k2 = ((struct cell_params *)params)->k2;
    double VMAXAcL = ((struct cell_params *)params)->VMAXAcL;
    double K_mAcL = ((struct cell_params *)params)->K_mAcL;
    double a2cH_slope = ((struct cell_params *)params)->a2cH_slope;
    double a2cH_thr = ((struct cell_params *)params)->a2cH_thr;
    double c2aH_slope = ((struct cell_params *)params)->c2aH_slope;
    double c2aH_thr = ((struct cell_params *)params)->c2aH_thr;
    double VMAXNHE = ((struct cell_params *)params)->VMAXNHE;
    double K_mNHE = ((struct cell_params *)params)->K_mNHE;
    double a = ((struct cell_params *)params)->a;
    double l_NHE = ((struct cell_params *)params)->l_NHE;
    double pH0_NHE = ((struct cell_params *)params)->pH0_NHE;
    double VMAXTHCO3 = ((struct cell_params *)params)->VMAXTHCO3;
    double K_mTHCO3 = ((struct cell_params *)params)->K_mTHCO3;
    double l_THCO3 = ((struct cell_params *)params)->l_THCO3;
    double pHe0_THCO3 = ((struct cell_params *)params)->pHe0_THCO3;
    double g_THCO3 = ((struct cell_params *)params)->g_THCO3;
    double pHi0_THCO3 = ((struct cell_params *)params)->pHi0_THCO3;
    double VMAXCA9 = ((struct cell_params *)params)->VMAXCA9;
    double K_mCA9 = ((struct cell_params *)params)->K_mCA9;
    double d_CA9 = ((struct cell_params *)params)->d_CA9;
    double V_c = ((struct cell_params *)params)->V_c;
    double dt = ((struct cell_params *)params)->dt;

    const double m_CO2_C = gsl_vector_get(x, 0);
    const double m_H_C = gsl_vector_get(x, 1);
    const double m_HCO3_C = gsl_vector_get(x, 2);
    const double m_H_c = gsl_vector_get(x, 3);
    const double m_HCO3_c = gsl_vector_get(x, 4);

    //intracellular carbondioxide dynamics

    const double y0 = m_CO2_C-m_CO2_C_old-dt * (
        //internalrate
        SensO2 * q_O2 * MW_CO2 / MW_O2
        //chemicalequilibrium
        -k1 * m_CO2_C + k2 * m_H_C * m_HCO3_C * 1000 * MW_CO2 / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H * MW_HCO3)
        //diffusion
        + PM_CO2 * (m_CO2_c_old / V_c-m_CO2_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0))) * (4.0 * Pi * pow(r_C, 2.0))
    );

    //intracellular hydrogen dynamics
    const double y1 = m_H_C - m_H_C_old - dt * (
        //internal rate
        SensATP * gAcL * MW_H / MW_AcL
        //chemical equilibrium
        + k1 * m_CO2_C * MW_H / MW_CO2 - k2 * m_H_C * m_HCO3_C * 1000 / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_HCO3)
        //nu_MCT_in->out
        - (2.0 - tanh(c2aH_slope * (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * Pi * pow(r_C, 2.0)) * m_H_C / ((4.0 / 3.0 * Pi * pow(r_C, 3.0)) * K_mAcL * MW_H / MW_AcL + m_H_C)
        //nu_MCT_out->in
        + (2.0 - tanh(a2cH_slope * (-log10(1000 * m_H_c / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * Pi * pow(r_C, 2.0)) * m_H_c / (V_c * K_mAcL * MW_H / MW_AcL + m_H_c)
        //nu_NHE_in->out
        - SensO2 * SensATP * 0.5 * (1.0 + ((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (4.0 * Pi * pow(r_C, 2.0)) * pow(m_H_C, a) / (pow(4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H * K_mNHE / 1000, a) + pow(m_H_C, a))
    );

    //intracellular bicarbonate ions dynamics
    const double y2 = m_HCO3_C - m_HCO3_C_old - dt * (
        //chemical equilibrium
        k1 * MW_HCO3 / MW_CO2 * m_CO2_C - k2 * m_H_C * m_HCO3_C * 1000 / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)
        //nu_THCO3_out->in
        + SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * Pi * pow(r_C, 2.0)) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
    );

    //extracellular hydrogen dynamics
    const double y3 = m_H_c - m_H_c_old - dt * (
        //chemical equilibrium
        k1 * m_CO2_c_old * MW_H / MW_CO2 - k2 * m_H_c * m_HCO3_c * 1000 / (V_c * MW_HCO3)
        //nu_MCT_in->out
        + (2.0 - tanh(c2aH_slope * (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H))) - c2aH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * Pi * pow(r_C, 2.0)) * m_H_C / ((4.0 / 3.0 * Pi * pow(r_C, 3.0)) * K_mAcL * MW_H / MW_AcL + m_H_C)
        //nu_MCT_out->in
        - (2.0 - tanh(a2cH_slope * (-log10(1000 * m_H_c / (V_c * MW_H))) - a2cH_thr)) * VMAXAcL * MW_H / MW_AcL 
        * (4.0 * Pi * pow(r_C, 2.0)) * m_H_c / (V_c * K_mAcL * MW_H / MW_AcL + m_H_c) 
        //nu_NHE_in->out 
        + SensATP * SensO2 * 0.5 * (1.0 + ((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE) / (l_NHE + abs((-log10(1000 * m_H_c / (V_c * MW_H))) - pH0_NHE))) 
        * VMAXNHE * (4.0 * Pi * pow(r_C, 2.0)) * pow(m_H_C, a) / (pow(4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H * K_mNHE / 1000, a) + pow(m_H_C, a)) 
        //nu_CA9 
        + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * 4.0 * Pi * pow(r_C, 2.0) * m_CO2_c_old / (V_c * K_mCA9 * MW_CO2 / 1000 + m_CO2_c_old) 
        * MW_H / MW_CO2
    );

    //extracellular bicarbonate ions dynamics
    const double y4 = m_HCO3_c - m_HCO3_c_old - dt * (
        //chemical equilibrium
        k1 * m_CO2_c_old * MW_HCO3 / MW_CO2 - k2 * m_H_c * m_HCO3_c * 1000 / (V_c * MW_H)
        //nu_THCO3_out->in
        - SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) 
        * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * Pi * pow(r_C, 2.0)) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
        //nu_CA9
        + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * 4.0 * Pi * pow(r_C, 2.0) * m_CO2_c_old / (V_c * K_mCA9 * MW_CO2 / 1000 + m_CO2_c_old) * MW_HCO3 / MW_CO2
    );

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);
    gsl_vector_set(f, 2, y2);
    gsl_vector_set(f, 3, y3);
    gsl_vector_set(f, 4, y4);

    return GSL_SUCCESS;
}

int main(void)
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status, time, j, k, perc;
    double dt, pH; //,pH_temp;
    size_t iter = 0;

    // respect units!
    const double MW_H = 1.0;                    // g/mol
    const double MW_CO2 = 44.0;                 // g/mol
    const double MW_O2 = 32.0;                  // g/mol
    const double MW_HCO3 = 61.0;                // g/mol
    const double MW_AcL = 90.1;                 // g/mol
    const double r_C = 6.55;                    // mim
    const double PM_CO2 = 3.2 * pow(10, 4);     // mim/s
    const double gAcL = 3.8 * pow(10,-4);       // pg/s
    const double q_O2 = 3.5 * pow(10,-5);       // pg/s
    const double k1 = 0.144;                    // 1/s
    const double k2 = 1.9 * pow(10, 5);         // 1/(M*s)
    const double VMAXAcL = 9.58 * pow(10,-5);   // pg/(s*mim^2)
    const double K_mAcL = 0.405 * pow(10,-3);   // pg/mim^3
    const double a2cH_slope = 1.5;             // adim
    const double a2cH_thr = 7.0;               // adim
    const double c2aH_slope = 1.5;              // adim
    const double c2aH_thr = 7.0;                // adim
    const double VMAXNHE = 5.15 * pow(10,-7);   // pg/(s*mim^2)
    const double K_mNHE = 0.196 * pow(10,-6);   // pg/mim^3
    const double a = 2.67;                      // adim
    const double l_NHE = 0.076;                 // adim
    const double pH0_NHE = 7.1;                 // adim
    const double VMAXTHCO3 = 2.02 * pow(10,-5); // pg/(s*mim^2)
    const double K_mTHCO3 = 7.38 * pow(10,-3);  // pg/mim^3
    const double l_THCO3 = 1.63;                // adim
    const double pHe0_THCO3 = 6.85;             // adim
    const double g_THCO3 = 4.2;                 // adim
    const double pHi0_THCO3 = 6.90;             // adim
    const double VMAXCA9 = 9.47 * pow(10,-2);   // pg/(s*mim^2)
    const double K_mCA9 = 7.2 * pow(10,-3);     // pg/mim^3
    const double d_CA9 = 7.3;                   // adim
    const double V_c = 1.0 * pow(10, 12);       // mim^3
    const double pKa = -log10(k1 / k2);          // adim
    const double pH_cell = 7.40;                // adim

    const int max_iter = 1000; //maxnumberofiterationsforNewton-Raphson

    // input from the user
    cout << "Time interval:" << endl;
    cin >> dt;
    cout << "Total integration time:" << endl;
    cin >> time;
    cout << "Starting pH : " << endl;
    cin >> pH;

    // starting conditions
    m_CO2_C_old = 5.39 * pow(10.0, -5) * (4.0 / 3.0 * Pi * pow(r_C, 3.0));
    m_H_C_old = pow(10.0, -pH_cell - 3.0) * (4.0 / 3.0 * Pi * pow(r_C, 3.0));
    m_HCO3_C_old = MW_HCO3 / MW_CO2 * m_CO2_C_old * pow(10.0, pH_cell - pKa);

    m_CO2_c_old = 5.39 * pow(10.0,-5) * V_c;
    m_H_c_old = pow(10.0, -pH - 3.0) * V_c;
    m_HCO3_c_old = MW_HCO3 / MW_CO2 * m_CO2_c_old * pow(10.0, pH - pKa);

    // friendly reminder
    cout << endl;
    cout << "***********************************************" << endl;
    cout << endl;
    cout << "Starting parameters" << endl;
    cout << endl;
    cout << "r_C\t(mim): " << r_C << endl;
    cout << "m_CO2_C\t(pg): " << m_CO2_C_old << endl;
    cout << "m_H_C\t(pg): " << m_H_C_old << endl;
    cout << "m_HCO3_C\t(pg): " << m_HCO3_C_old << endl;
    cout << "m_CO2_c\t(pg): " << m_CO2_c_old << endl;
    cout << "m_H_c\t(pg): " << m_H_c_old << endl;
    cout << "m_HCO3_c\t(pg): " << m_HCO3_c_old << endl;
    cout << "V_c\t(mim^3): " << V_c << endl;
    cout << "dt: " << dt << endl;
    cout << "steps: " << time << endl;
    cout << "time: " << time * dt << endl;
    cout << "startingpH: " << pH << endl;
    cout << "sensO2: " << SensO2 << endl;
    cout << "sensATP: " << SensATP << endl;
    cout << endl;
    cout << "***********************************************" << endl;
    cout << endl;
    cout << "Running..." << endl;
    cout << endl;

    //FYI
    perc = 10;
    cout << "-- completed at : 0%" << endl;

    const size_t n = 5;
    struct cell_params cell_p = {
        MW_H, MW_CO2, MW_O2, MW_HCO3, MW_AcL, r_C, PM_CO2, gAcL, q_O2, k1, k2, VMAXAcL, K_mAcL, a2cH_slope,
        a2cH_thr, c2aH_slope, c2aH_thr, VMAXNHE, K_mNHE, a, l_NHE, pH0_NHE, VMAXTHCO3, K_mTHCO3, l_THCO3,
        pHe0_THCO3, g_THCO3, pHi0_THCO3, VMAXCA9, K_mCA9, d_CA9, V_c, dt
    };

    gsl_multiroot_function cell_f = {&cell, n, &cell_p};

    double x_init[n] = { m_CO2_C_old, m_H_C_old, m_HCO3_C_old, m_H_c_old, m_HCO3_c_old }; // starting point
    gsl_vector *x = gsl_vector_alloc(n);

    for (k = 0; k < n; k++)
    {
        gsl_vector_set(x, k, x_init[k]);
    }

    T = gsl_multiroot_fsolver_dnewton; // discrete Newton (discrete Jacobian)
    s = gsl_multiroot_fsolver_alloc(T, n);

    gsl_multiroot_function f = cell_f;

    outData.open("cell_output.txt");      // output masses
    outDatapH.open("cell_output_pH.txt"); // output pH

    gsl_multiroot_fsolver_set(s, &f, x);

    // output on masses file
    outData << "# Starting parameters" << endl;

    outData << endl;
    outData << "r_C\t(mim): " << r_C << endl;
    outData << "# m_CO2_C\t(pg):" << m_CO2_C_old << endl;
    outData << "# m_H_C\t(pg):" << m_H_C_old << endl;
    outData << "# m_HCO3_C\t(pg):" << m_HCO3_C_old << endl;
    outData << "# m_CO2_c\t(pg):" << m_CO2_c_old << endl;
    outData << "# m_H_c\t(pg):" << m_H_c_old << endl;
    outData << "# m_HCO3_c\t(pg):" << m_HCO3_c_old << endl;
    outData << "# V_c\t(mim^3):" << V_c << endl;
    outData << "# dt:" << dt << endl;
    outData << "# steps:" << time << endl;
    outData << "# time:" << time * dt << endl;
    outData << "# startingpH:" << pH << endl;
    outData << "# sensO2:" << SensO2 << endl;
    outData << "# sensATP:" << SensATP << endl;
    outData << endl;
    outData << "# step\tm_CO2_C\tm_H_C\tm_HCO3_C\tm_H_c\tm_HCO3_c" << endl;
    outData << endl;

    // output on pH file
    outDatapH << "# Starting parameters" << endl;
    outDatapH << endl;
    outDatapH << "r_C\t(mim): " << r_C << endl;
    outDatapH << "# m_CO2_C\t(pg):" << m_CO2_C_old << endl;
    outDatapH << "# m_H_C\t(pg):" << m_H_C_old << endl;
    outDatapH << "# m_HCO3_C\t(pg):" << m_HCO3_C_old << endl;
    outDatapH << "# m_CO2_c\t(pg):" << m_CO2_c_old << endl;
    outDatapH << "# m_H_c\t(pg):" << m_H_c_old << endl;
    outDatapH << "# m_HCO3_c\t(pg):" << m_HCO3_c_old << endl;
    outDatapH << "# V_c\t(mim^3):" << V_c << endl;
    outDatapH << "# dt:" << dt << endl;
    outDatapH << "# steps:" << time << endl;
    outDatapH << "# time:" << time * dt << endl;
    outDatapH << "# startingpH:" << pH << endl;
    outDatapH << "# sensO2:" << SensO2 << endl;
    outDatapH << "# sensATP:" << SensATP << endl;
    outDatapH << endl;
    outDatapH << "# step\tpH_C\tpH_c\tpH_CHH\tpH_cHH" << endl;
    outDatapH << endl;

    // gonna need them
    outData << fixed;
    outData << setprecision(15);

    outDatapH << fixed;
    outDatapH << setprecision(15);

    // first term
    outData << "0" << "\t" << m_CO2_C_old << "\t" << m_H_C_old << "\t" << m_HCO3_C_old << "\t" << m_H_c_old << "\t" << m_HCO3_c_old << endl;
    outDatapH << "0" << "\t"
        << -log10(1000 * m_H_C_old / (4.0 / 3.0 * Pi * pow(r_C, 3.0))) << "\t"
        << -log10(1000 * m_H_c_old / V_c) << "\t" << log10((m_HCO3_C_old * MW_CO2 * k2) / (m_CO2_C_old * MW_HCO3 * k1))
        << "\t" << log10((m_HCO3_c_old * MW_CO2 * k2) / (m_CO2_c_old * MW_HCO3 * k1)) << endl;

    // starts the time
    for (j = 1; j < time; j++)
    {
        iter = 0;

        do
        {
            iter++;

            status = gsl_multiroot_fsolver_iterate(s);

            if (status) // check if solver is stuck
                break;

            status = gsl_multiroot_test_residual(s->f, 1e-6);

        } while (status == GSL_CONTINUE && iter < max_iter);

        // lazy...
        if (j * 100 / time == perc)
        {
            cout << "-- completed at:" << perc << "%" << endl;
            perc = perc + 10;
        }

        // output control
        if (j % 100 == 0)
        {

            outData << j * dt << "\t" << gsl_vector_get(s->x, 0) << "\t" << gsl_vector_get(s->x, 1) << "\t" << gsl_vector_get(s->x, 2)
                << "\t" << gsl_vector_get(s->x, 3) << "\t" << gsl_vector_get(s->x, 4) << endl;

            outDatapH << j * dt << "\t"
                << -log10(1000 * gsl_vector_get(s->x, 1) / (4.0 / 3.0 * Pi * pow(r_C, 3.0))) << "\t"
                << -log10(1000 * gsl_vector_get(s->x, 3) / V_c)
                << "\t" << log10((gsl_vector_get(s->x, 2) * MW_CO2 * k2) / (gsl_vector_get(s->x, 0) * MW_HCO3 * k1))
                << "\t" << log10((gsl_vector_get(s->x, 4) * MW_CO2 * k2) / (m_CO2_c_old * MW_HCO3 * k1)) << endl;
        }

        x_init[0] = gsl_vector_get(s->x, 0);
        x_init[1] = gsl_vector_get(s->x, 1);
        x_init[2] = gsl_vector_get(s->x, 2);
        x_init[3] = gsl_vector_get(s->x, 3);
        x_init[4] = gsl_vector_get(s->x, 4);

        m_CO2_C_old = x_init[0];
        m_H_C_old = x_init[1];
        m_HCO3_C_old = x_init[2];
        m_H_c_old = x_init[3];
        m_HCO3_c_old = x_init[4];

        // pH drop at time j*dt
        /*
            if(j == 200000)
            {
                m_H_C_old = m_H_C_old * pow(10.0, 1.0);
                pH_temp = -log10(m_H_C_old * 1000.0/ (4.0 / 3.0 * Pi * pow(r_C,3.0) * MW_H));
                m_HCO3_C_old = MW_HCO3 / MW_CO2 * m_CO2_C_old * pow(10.0, pH_temp - pKa);
                x_init[1] = m_H_C_old;
                x_init[2] = m_HCO3_C_old;
            }
        */

        for (k = 0; k < n; k++)
        {
            gsl_vector_set(x, k, x_init[k]);
        }

        gsl_multiroot_fsolver_set(s, &f, x);
    }

    // closing
    cout << "-- completed at:" << perc << "%" << endl;
    cout << "-- done !" << endl;
    cout << endl;

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    outData.close();
    outDatapH.close();

    //because windows
    //system("pause");

    return 0;
}