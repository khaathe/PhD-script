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
#include<vector>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
using namespace std::placeholders;
typedef std::vector<double> state_type;

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

void cell (state_type const &x, state_type &dxdt, double t, void *params)
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

    const double m_CO2_C = x[0];
    const double m_H_C = x[1];
    const double m_HCO3_C = x[2];
    const double m_H_c = x[3];
    const double m_HCO3_c = x[4];

    //intracellular carbondioxide dynamics

    dxdt[0] = 
        //internalrate
        SensO2 * q_O2 * MW_CO2 / MW_O2
        //chemicalequilibrium
        -k1 * m_CO2_C + k2 * m_H_C * m_HCO3_C * 1000 * MW_CO2 / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H * MW_HCO3)
        //diffusion
        + PM_CO2 * (m_CO2_c_old / V_c-m_CO2_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0))) * (4.0 * Pi * pow(r_C, 2.0))
    ;

    //intracellular hydrogen dynamics
    dxdt[1] = 
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
    ;

    //intracellular bicarbonate ions dynamics
    dxdt[2] = 
        //chemical equilibrium
        k1 * MW_HCO3 / MW_CO2 * m_CO2_C - k2 * m_H_C * m_HCO3_C * 1000 / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)
        //nu_THCO3_out->in
        + SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * Pi * pow(r_C, 2.0)) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
    ;

    //extracellular hydrogen dynamics
    dxdt[3] = 
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
    ;

    //extracellular bicarbonate ions dynamics
    dxdt[4] = 
        //chemical equilibrium
        k1 * m_CO2_c_old * MW_HCO3 / MW_CO2 - k2 * m_H_c * m_HCO3_c * 1000 / (V_c * MW_H)
        //nu_THCO3_out->in
        - SensATP * (0.5) * (1.0 + tanh(l_THCO3 * (-log10(1000 * m_H_c / (V_c * MW_H)) - pHe0_THCO3))) 
        * (0.5) * (1.0 + tanh(g_THCO3 * (pHi0_THCO3 - (-log10(1000 * m_H_C / (4.0 / 3.0 * Pi * pow(r_C, 3.0) * MW_H)))))) 
        * VMAXTHCO3 * (4.0 * Pi * pow(r_C, 2.0)) * m_HCO3_c / (V_c * K_mTHCO3 * MW_HCO3 / 1000 + m_HCO3_c)
        //nu_CA9
        + (3.0 + 2.0 * tanh(-d_CA9 * SensO2)) * VMAXCA9 * 4.0 * Pi * pow(r_C, 2.0) * m_CO2_c_old / (V_c * K_mCA9 * MW_CO2 / 1000 + m_CO2_c_old) * MW_HCO3 / MW_CO2
    ;
}

int main(void)
{
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

    outData.open("cell_output.txt");      // output masses
    outDatapH.open("cell_output_pH.txt"); // output pH

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
    
    double x_init[n] = { m_CO2_C_old, m_H_C_old, m_HCO3_C_old, m_H_c_old, m_HCO3_c_old }; // starting point
    // starts the time
    for (j = 1; j < time; j++)
    {
        state_type x;
        x.resize(n);
        copy(x_init, x);
        double t0 = 0.0;
        double end_time = dt;
        double step_size = dt/10;
        auto ode = std::bind(cell, _1, _2, _3, &cell_p);

        // Solve ODEs
        typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
        typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
        double abs_err = 1.0e-10, rel_err = 1.0e-6, a_x = 1.0, a_dxdt = 1.0;
        controlled_stepper_type controlled_stepper( default_error_checker<double, range_algebra, default_operations>(abs_err, rel_err, a_x, a_dxdt) );
        integrate_adaptive( controlled_stepper, ode, x, t0, end_time, step_size);

        // lazy...
        if (j * 100 / time == perc)
        {
            cout << "-- completed at:" << perc << "%" << endl;
            perc = perc + 10;
        }

        // output control
        if (j % 100 == 0)
        {

            outData << j * dt << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << x[4] << endl;

            outDatapH << j * dt << "\t"
                << -log10(1000 * x[1] / (4.0 / 3.0 * Pi * pow(r_C, 3.0))) << "\t"
                << -log10(1000 * x[3] / V_c)
                << "\t" << log10( ( x[2] * MW_CO2 * k2) / ( x[0] * MW_HCO3 * k1) )
                << "\t" << log10((x[4] * MW_CO2 * k2) / (m_CO2_c_old * MW_HCO3 * k1)) << endl;
        }

        copy(x, x_init);
        
        m_CO2_C_old = x_init[0];
        m_H_C_old = x_init[1];
        m_HCO3_C_old = x_init[2];
        m_H_c_old = x_init[3];
        m_HCO3_c_old = x_init[4];
    }

    // closing
    cout << "-- completed at:" << perc << "%" << endl;
    cout << "-- done !" << endl;
    cout << endl;

    outData.close();
    outDatapH.close();

    //because windows
    //system("pause");

    return 0;
}