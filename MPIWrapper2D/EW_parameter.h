#ifndef EW_PARAMETER_H
#define EW_PARAMETER_H

#include "eigen/Eigen/Eigen"
#include "./ParaSite/ParaSite.h"
#include "EW_helper.h"

namespace Electroweak{
    typedef double Real;
    typedef float RealSave;
    typedef std::complex<Real> Cmplx;

    typedef Eigen::Matrix<Cmplx, 2, 1> SU2vector;
    typedef Eigen::Matrix<Cmplx, 2, 2> SU2matrix;
    typedef Cmplx U1matrix;

    //---math constants---
    const Real ERROR = 1e-10;
    const Real EVO_ERROR = 1e-10;
	constexpr Cmplx I(0.0, 1.0);
	constexpr Real PI = 3.14159265358979;

    //---constants for electroweak theory
    constexpr Real g = 0.65;
    constexpr Real sin2w = 0.22;
    const Real sinw = sqrt(sin2w);
    constexpr Real cos2w = 1.0 - sin2w;
    const Real cosw = sqrt(cos2w);

    const Real gp = g*sinw / cosw; //g'=g*tanw
    const Real gz = sqrt(g*g + gp*gp);
    const Real ge = g*gp / gz;
    const Real damping = std::stod(ReadConfigIni("BUBBLE_DAMPING")); //Phi field damping
    const Real DAMPING_HIGGS = damping;
    const Real DAMPING_GAUGE = 0; //such assignment does not spoil gauge invariance.
    constexpr int FermiFamily = 1;

    constexpr int DIM = 3;
    constexpr int SU2DIM = 3;
    constexpr int CYCLE = 2;
    constexpr int halo = 0;

    const ParaSite::IndexType N = std::stoi(ReadConfigIni("EWLAT_CUBE_NSIZE"));
    const ParaSite::IndexType nSize[DIM] = { N,N,N }; //evolution only works for Nx=Ny=Nz yet.
    const int startStep = std::stoi(ReadConfigIni("EWLAT_START_STEP")); //isContinue: if Continue, startStep = Last time Ntimesteps.
    const int Ntimesteps = std::stoi(ReadConfigIni("EWLAT_END_STEP")); //### running steps ###
    const int isCalc = std::stoi(ReadConfigIni("EWLAT_IS_CALC")); //calculate related quantities every *isCalc* steps
    const int isSave = std::stoi(ReadConfigIni("EWLAT_IS_SAVE"));//save every *isSave* steps, must be an integer multiple of *isCalc*

    const int SAMPLE_NUMBER = std::stoi(ReadConfigIni("SAMPLE_NUMBER"));
    const Real DX = std::stod(ReadConfigIni("EWLAT_DX"));
    const Real DX2 = DX*DX;
    const Real DX3 = DX*DX*DX;
	const Real DX4 = DX * DX*DX*DX;
    const Real DT = DX / 4;

    const Real CENTER_POS[DIM] = { nSize[0] / 2 - 0.5 , nSize[1] / 2 - 0.5 , nSize[2] / 2 - 0.5 };

    const Real lambda = std::stod(ReadConfigIni("lambda"));//lambda=0.129
    const Real v = 6.0;
    const Real v2 = v*v;

    //---auxilliary constants---
    const Real UNITY_E = 2.0 / (gp*DX*DT);
    const Real UNITY_F = 1.0 / (g*DX*DT);
    const Real UNITY_EB = UNITY_E*DT / DX; //==2.0/(gp*DX*DX)
    const Real UNITY_FB = UNITY_F*DT / DX; //==1.0/(g*DX*DX)

    //--boundary constants--
    enum LATTICE_BOUNDARY_TYPE {
        INSIDE,
        FACE,
        EDGE,
        CORNER,
        OUTSIDE,
        UNDETERMINED
    };

    //---Pauli matrices---
    const SU2matrix Ident = SU2matrix::Identity();
    const SU2matrix Pauli[SU2DIM] = {
        (SU2matrix() << 0., 1., 1., 0.).finished() ,
        (SU2matrix() << 0., -I, I, 0.).finished() ,
        (SU2matrix() << 1., 0., 0., -1.).finished()
    };
    const SU2matrix iPauli[SU2DIM] = {
        (SU2matrix() << 0., I, I, 0.).finished() ,
        (SU2matrix() << 0., 1., -1., 0.).finished() ,
        (SU2matrix() << I, 0., 0., -I).finished()
    };

    constexpr int X_AXIS = 0;
    constexpr int Y_AXIS = 1;
    constexpr int Z_AXIS = 2;

}
#endif