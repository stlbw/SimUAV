#include <cmath>
/**
 * Computes the air density at a certain altitude using a simplified mathematical model
 * @param h
 * @return
 */

struct aero_condition{
    double rho_SL = 1.225;
    double T_0 = 288.15;
    double p_sl = 101325;
    double rho_h = 0;
    double T_h = 0;
    double p_h = 0;
    double altitude = 0;
    int input;
};




aero_condition AtmosphereChoice (aero_condition hb) {

    cout << "\n\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ISA ATMOSPHERE\\\\\\\\\\\\\\\\\\\\\\\\";
    cout << "\nThe simulation refers to the ISA atmospheric model with the following values for altitude h = 0 m (Sea Level):\n\n";
    cout << "\tPressure: \t\tP = " << hb.p_sl << " Pa" << endl;
    cout << "\tTemperature \t\tT = " << hb.T_0 <<  " K" << endl;
    cout << "\tDensity: \t\trho = " << hb.rho_SL <<  " kg/m^3" << endl;
    //cout << "\tSound Speed: \ta = %g m/s\n" << db.speed_sl;

    int input;
    double pressure, temperature, density;

    cout << ("\nIf you do not wish to proceed with the above parameters, you can change them by re-entering them manually or choosing a different altitude.\n");
    do {
        cout << "\nPress:\n";
        cout << "\t1 - if you wish to proceed with the above values\n";
        cout << "\t2 - if you wish to change the initial values to h = 0\n";
        cout << "\t3 - if you wish to enter values manually at a specific altitude (altitude is computed automatically from these values)\n";

        cin >> input;

        switch(input)
        {
            case 1:
                break;
            case 2:
                cout << "Enter a pressure value (h = 0) in Pa: ";
                cin >> pressure;
                if (pressure < 0) {
                    do {
                        cout << "[!]WARNING enter a positive number\n";
                        cout << "Enter a pressure value (h = 0) in Pa: ";
                        cin >> pressure;
                    } while (pressure < 0);
                }
                hb.p_sl = pressure;

                cout << "Enter a temperature value (h=0) in K: ";
                cin >> temperature;
                hb.T_0 = temperature;

                cout << "Enter a density value (h=0) in kg/m^3: ";
                cin >> density;
                if (density < 0) {
                    do {
                        cout << "[!]WARNING enter a positive number\n";
                        cout << "Enter a density value (h=0) in kg/m^3: ";
                        cin >> density;
                    } while (density < 0);
                }
                hb.rho_SL = density;
                break;
            case 3:
                cout <<"Enter a pressure value in Pa: ";
                cin >> pressure;
                if (pressure < 0) {
                    do {
                        cout << "[!]WARNING enter a positive number\n";
                        cout << "Enter a pressure value in Pa: ";
                        cin >> pressure;
                    } while (pressure < 0);
                }
                hb.p_h = pressure;

                cout <<"Enter a temperature value in K: ";
                cin >> temperature;
                hb.T_h = temperature;

                cout << "Enter a density value in kg/m^3: ";
                cin >> density;
                if (density < 0) {
                    do {
                        cout << "[!]WARNING enter a positive number\n";
                        cout << "Enter a density value in kg/m^3: ";
                        cin >> density;
                    } while (density < 0);
                }
                hb.rho_h = density;

                double grad = 0.0065;
                hb.altitude = (hb.T_0 - hb.T_h) / grad;
                break;
        }
    }
    while(input!=1 && input!=2 && input!=3);

    hb.input = input;

    return hb;
}



/*double computeDensity(double h) {
    double rho_SL = 1.225;
    double grad = 0.0065;
    double m = 4.2561;
    double T_0 = 288.15;
    double ans = ((T_0 - grad * h) / T_0);
    double rho = rho_SL * pow(ans, m);

    return rho;
}*/

double AtmosphereCalc (int input, aero_condition hb, double h_ref){
    double R = 287.05;
    double grad = 0.0065;
    double m = 4.2561;

    switch (input)
    {
        case 1:
            hb.T_h = hb.T_0 - grad * h_ref;
            hb.p_h = hb.p_sl * pow((hb.T_h / hb.T_0), m + 1);
            hb.rho_h = hb.p_h / (R * hb.T_h);
            break;
        case 2:
            hb.T_h = hb.T_0 - grad * h_ref;
            hb.p_h = hb.p_sl * pow((hb.T_h / hb.T_0), m + 1);
            hb.rho_h = hb.p_h / (R * hb.T_h);
            break;
        case 3:
            break;
    }

    double rho = hb.rho_h;

    return rho;
}