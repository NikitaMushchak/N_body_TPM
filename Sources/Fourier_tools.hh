#pragma once
#include "ai.hh"




// Cooley-Tukey FFT

/**/

/**/;
 inline void conjugate2(std::vector< std::vector<double> > &vector){
    std::size_t length = vector.size();

    for(std::size_t i = 0; i < length; ++i){
        vector[i][1] = -vector[i][1];
    }
}

 inline void fft2(std::vector< std::vector<double> > &vector){
    std::size_t k = vector.size();
    std::size_t j = 0;
    std::size_t n = 0;

    const double length = (double) k;

    double thetaT = 3.14159265358979 / length;
    double swap0 = 0.;
    double swap1 = 0.;
    double T0 = 1.;
    double T1 = 0.;
    double phiT0 = cos(thetaT);
    double phiT1 = -sin(thetaT);

    while(k > 1){
        n = k;

        k >>= 1;

        swap0 = phiT0;
        swap1 = phiT1;

        phiT0 = swap0 * swap0 - swap1 * swap1;
        phiT1 = 2. * swap0 * swap1;

        T0 = 1.;
        T1 = 0.;

        for(std::size_t l = 0; l < k; ++l)
        {
            for(std::size_t i = l; i < length; i += n)
            {
                j = i + k;

                swap0 = vector[i][0] - vector[j][0];
                swap1 = vector[i][1] - vector[j][1];

                vector[i][0] += vector[j][0];
                vector[i][1] += vector[j][1];

                vector[j][0] = swap0 * T0 - swap1 * T1;
                vector[j][1] = swap0 * T1 + swap1 * T0;
            }

            swap0 = T0;

            T0 = swap0 * phiT0 - T1 * phiT1;
            T1 = swap0 * phiT1 + T1 * phiT0;
        }
    }

    std::size_t m = (std::size_t) log2(length);

    for(std::size_t i = 0; i < length; ++i){
        j = i;

        j = (((j & 0xaaaaaaaa) >> 1) | ((j & 0x55555555) << 1));
        j = (((j & 0xcccccccc) >> 2) | ((j & 0x33333333) << 2));
        j = (((j & 0xf0f0f0f0) >> 4) | ((j & 0x0f0f0f0f) << 4));
        j = (((j & 0xff00ff00) >> 8) | ((j & 0x00ff00ff) << 8));
        j = ((j >> 16) | (j << 16)) >> (32 - m);

        if(j > i){
            swap0 = vector[i][0];
            swap1 = vector[i][1];

            vector[i][0] = vector[j][0];
            vector[i][1] = vector[j][1];

            vector[j][0] = swap0;
            vector[j][1] = swap1;
        }
    }
}

 inline void ifft2(std::vector< std::vector<double> > &vector){
    const double length = (double) vector.size();

    conjugate2(vector);

    fft2(vector);
		conjugate2(vector);

    for(std::size_t i = 0; i < length; ++i){
        vector[i][0] /= length;
				vector[i][1] /= length;
    }
}
inline  void fft2_Real(std::vector< std::vector<double> > &vector) {

    std::size_t N = vector.size();
    std::size_t N_2 = N / 2;

    double thetaT2 = 3.14159265358979 / N_2;
    double swap0 = 0.;
    double U0 = 1.;
    double U1 = 0.;
    double phiT0_2 = cos(thetaT2);
    double phiT1_2 = sin(thetaT2);

    //std::vector<double> vector0, vector1;
    std::vector< std::vector<double> > vector0;
    vector0.resize(N_2);

    std::size_t i;
    for (i = 0; i < N_2; ++i)
    {
        vector0[i].resize(2);
        vector0[i][0] = vector[2 * i][0] / 2.0;
        vector0[i][1] = vector[2 * i + 1][0] / 2.0;
    }

    fft2(vector0);

    for (i = 0; i < N_2; ++i) {

        std::size_t k = (N_2 - i) % N_2;

        double Ax = vector0[i][0];
        double Ay = vector0[i][1];
        double Bx = vector0[k][0];
        double By = vector0[k][1];

        double Cx = Ax + Bx;
        double Cy = Ay - By;
        double Dx = Ay + By;
        double Dy = Ax - Bx;

        Ax = U0 * Dx + U1 * Dy;
        Ay = U1 * Dx - U0 * Dy;

        vector[i][0] = Ax + Cx;
        vector[i][1] = Ay + Cy;
        vector[i + N_2][0] = Cx - Ax;
        vector[i + N_2][1] = Cy - Ay;

        swap0 = U0;
        U0 = swap0 * phiT0_2 - U1 * phiT1_2;
        U1 = swap0 * phiT1_2 + U1 * phiT0_2;
    }
}

// Îáðàòíîå FFT ñ âåùåñòâåííûì âûõîäîì

inline void ifft2_Real(std::vector< std::vector<double> > &vector) {

    std::size_t N = vector.size();
    std::size_t N_2 = N / 2;

    double thetaT2 = 3.14159265358979 / N_2;
    double swap0 = 0.;
    double U0 = 1.;
    double U1 = 0.;
    double phiT0_2 = cos(thetaT2);
    double phiT1_2 = sin(thetaT2);

    //std::vector<double> vector0, vector1;
    std::vector< std::vector<double> > vector0;
    vector0.resize(N_2);

    std::size_t i;
    for (i = 0; i < N_2; ++i) {

        double X1_Re = vector[i][0];
        double X1_Im = vector[i][1];
        double X2_Re = vector[i + N_2][0];
        double X2_Im = vector[i + N_2][1];

        double Y0_Re = X1_Re + X2_Re;
        double Y0_Im = X1_Im + X2_Im;
        double Y1_Re = U0 * (X1_Re - X2_Re) - U1 * (X1_Im - X2_Im);
        double Y1_Im = U1 * (X1_Re - X2_Re) + U0 * (X1_Im - X2_Im);

        vector0[i].resize(2);
        vector0[i][0] = (Y0_Re - Y1_Im) / 2.0;
        vector0[i][1] = (Y0_Im + Y1_Re) / 2.0;

        swap0 = U0;
        U0 = swap0 * phiT0_2 - U1 * phiT1_2;
        U1 = swap0 * phiT1_2 + U1 * phiT0_2;
    }

    ifft2(vector0);

    for (i = 0; i < N_2; ++i) {

        vector[2 * i][0] = vector0[i][0];
        vector[2 * i + 1][0] = vector0[i][1];
        vector[i][1] = vector[i + N_2][1] = 0;
    }
}
