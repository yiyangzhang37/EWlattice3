#include "EW_parameter.h"
#include "EW_functions.h"
using namespace Eigen;

//---Functions for Calculation of Physical Quantities---
//01/23/2017 The updatedTimeStep and nowTimeStep are modified, but these functions are not modified yet.
void CalcEnergyDensity_KS(Site x,
	Field<RealSave>& saveData, int saveCol,
	Field<SU2vector>& phi, Field<SU2vector>& pi,
	Field<SU2matrix>& U, Field<SU2matrix>& F,
	Field<U1matrix>& V, Field<U1matrix>& E,
	const int calcTimeStep,
	Real& totalEnergy,
	Real* arrVar) {
	//in this function, the energy is calculated as the current step, not the updated step.
	//The problem of this function is, when initial condition is not accurate, the calculated energy density may be not positive-definite.
	int nowTime = calcTimeStep % CYCLE;
	int preTime;
	if (CYCLE == 2) {
		preTime = (calcTimeStep + 1) % CYCLE;//ONLY true for CYCLE=2
	}
	else {
		preTime = nowTime; //true for CYCLE=1
	}
	totalEnergy = 0;
	for (x.first(); x.test(); x.next()) {
		// ###NOT including sites on boundaries###
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			Real eKinetic = 0, eGrad = 0, ePotential = 0, eU1 = 0;
			Cmplx eSU2(0, 0);
			eKinetic = (0.5*(pi(x, nowTime) + pi(x, preTime))).squaredNorm();
			ePotential = Vphi(x, phi, nowTime);
			for (int i = 0; i<DIM; ++i) {
				eGrad += 1.0 / (DX*DX)*(U(x, i, nowTime)*V(x, i, nowTime)*phi(x + i, nowTime) - phi(x, nowTime)).squaredNorm();
				eU1 += 4.0 / (pow(DT*DX*gp, 2)) * (1.0 - 1.0 / UNITY_E * 0.5*(E(x, i, nowTime) + E(x, i, preTime)).real());
				eSU2 += 2.0 / (pow(DT*DX*g, 2)) * (2.0 - 1.0 / UNITY_F * 0.5*(F(x, i, nowTime) + F(x, i, preTime)).trace());
				for (int j = 0; j < DIM; ++j) {
					eU1 += 2.0 / (pow(DX, 4)*pow(gp, 2)) * (1.0 - Plaquette(x, V, i, j, nowTime).real());
					eSU2 += 1.0 / (pow(DX, 4)*pow(g, 2)) * (2.0 - Plaquette(x, U, i, j, nowTime).trace());
				}
			}
#ifdef TEST_MODE
			if (eSU2.imag()>ERROR) cout << "SU2 FIELD ERROR:ENERGY." << endl;
#endif
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = static_cast<RealSave> (DX*DX*DX * (eKinetic + eGrad + ePotential + eU1 + eSU2.real()));
			totalEnergy += eKinetic + ePotential + eGrad + eU1 + eSU2.real();
			if (arrVar != NULL) {
				arrVar[0] += eKinetic;
				arrVar[1] += eGrad;
				arrVar[2] += ePotential;
				arrVar[3] += eU1;
				arrVar[4] += eSU2.real();
			}
		}
	}
	totalEnergy *= DX3;
	if (arrVar != NULL) {
		arrVar[0] *= DX3;
		arrVar[1] *= DX3;
		arrVar[2] *= DX3;
		arrVar[3] *= DX3;
		arrVar[4] *= DX3;
	}
	return;
}

[[deprecated]] void CalcEnergyDensityComponents_KS(Site x,
	Field<RealSave>& saveData, int saveCol,
	Field<SU2vector>& phi, Field<SU2vector>& pi,
	Field<SU2matrix>& U, Field<SU2matrix>& F,
	Field<U1matrix>& V, Field<U1matrix>& E,
	const int component,
	const int calcTimeStep,
	Real& totalEnergyComp,
	Real* arrVar) {
	//in this function, the energy is calculated as the current step, not the updated step.
	//The problem of this function is, when initial condition is not accurate, the calculated energy density may be not positive-definite.
	//component: 0==Higgs Kinetic; 1==Higgs Gradient; 2==Higgs potential; 3==U1; 4==SU2.
	int updateTime = calcTimeStep % CYCLE;
	int nowTime = (calcTimeStep + 1) % CYCLE;//ONLY true for CYCLE=2
	totalEnergyComp = 0;
	for (x.first(); x.test(); x.next()) {
		int isBoundary = 0;
		for (int j = 0; j<DIM; ++j) {
			if (x.coord(j) == 0 || x.coord(j) == nSize[j] - 1) { isBoundary = 1; break; }
		}
		// ###NOT including sites on boundaries###
		if (isBoundary > 0) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			Real eKinetic = 0, eGrad = 0, ePotential = 0, eU1 = 0;
			Cmplx eSU2(0, 0);
			if (component == 0) eKinetic = (0.5*(pi(x, updateTime) + pi(x, nowTime))).squaredNorm();
			if (component == 2) ePotential = Vphi(x, phi, updateTime);
			for (int i = 0; i<DIM; ++i) {
				if (component == 1) eGrad += 1.0 / (DX*DX)*(U(x, i, updateTime)*V(x, i, updateTime)*phi(x + i, updateTime) - phi(x, updateTime)).squaredNorm();
				if (component == 3) eU1 += 4.0 / (pow(DT*DX*gp, 2)) * (1.0 - 1.0 / UNITY_E * 0.5*(E(x, i, updateTime) + E(x, i, nowTime)).real());
				if (component == 4) eSU2 += 2.0 / (pow(DT*DX*g, 2)) * (2.0 - 1.0 / UNITY_F * 0.5*(F(x, i, updateTime) + F(x, i, nowTime)).trace());
				for (int j = 0; j < DIM; ++j) {
					if (component == 3) eU1 += 2.0 / (pow(DX, 4)*pow(gp, 2)) * (1.0 - Plaquette(x, V, i, j, updateTime).real());
					if (component == 4)eSU2 += 1.0 / (pow(DX, 4)*pow(g, 2)) * (2.0 - Plaquette(x, U, i, j, updateTime).trace());
				}
			}
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = static_cast<RealSave> (DX*DX*DX * (eKinetic + eGrad + ePotential + eU1 + eSU2.real()));
			totalEnergyComp += eKinetic + ePotential + eGrad + eU1 + eSU2.real();
		}
	}
	totalEnergyComp *= DX * DX*DX;
	return;
}
#ifdef IMPROVED_EVOLUTION
void CalcEnergyDensity_IMPROVED(Site x,
	Field<RealSave>& saveData, const int saveCol,
	Field<SU2vector>& phi, Field<SU2vector>& pi,
	Field<SU2matrix>& U, Field<SU2matrix>& F,
	Field<U1matrix>& V, Field<U1matrix>& E,
	const int updatedTimeStep,
	Real& totalEnergy,
	Real* arrVar) {
	int updateTime = updatedTimeStep % CYCLE;
	int pastTime = (updatedTimeStep + 1) % CYCLE;
	totalEnergy = 0;
	arrVar[0] = 0;
	arrVar[1] = 0;
	arrVar[2] = 0;
	for (x.first(); x.test(); x.next()) {
		int isBoundary = 0;
		for (int j = 0; j<DIM; ++j) {
			if (x.coord(j) == 0 || x.coord(j) == nSize[j] - 1) { isBoundary = 1; break; }
		}
		// ###NOT including sites on boundaries###
		if (isBoundary > 0) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			Cmplx eScalar(0, 0);
			Real eU1 = 0;
			Cmplx eSU2(0, 0);
			eScalar = 0.5*(pi(x, updateTime).squaredNorm() + pi(x, pastTime).squaredNorm()) + (lambda*pow((v2 - phi(x, updateTime).squaredNorm()), 2));
			for (int i = 0; i<DIM; ++i) {
				eScalar += 1.0 / (DX*DX)*(phi(x, updateTime) - U(x - i, i, updateTime).adjoint()*conj(V(x - i, i, updateTime))*phi(x - i, updateTime)).dot(
					-1. / 12.*U(x, i, updateTime)*V(x, i, updateTime)*phi(x + i, updateTime) + 1.25*phi(x, updateTime)
					- 1.25*U(x - i, i, updateTime).adjoint()*conj(V(x - i, i, updateTime))*phi(x - i, updateTime)
					+ 1. / 12.*U(x - i, i, updateTime).adjoint()*conj(V(x - i, i, updateTime))*U(x - i - i, i, updateTime)*conj(V(x - i - i, i, updateTime))*phi(x - i - i, updateTime));
				//TODO: electric energy
				eU1 += 4.0 / (pow(DT*DX*gp, 2)) * (1.0 - 1.0 / UNITY_E * E(x, i, updateTime).real());
				eSU2 += 2.0 / (pow(DT*DX*g, 2)) * (2.0 - 1.0 / UNITY_F * F(x, i, updateTime).trace());

				for (int j = 0; j<i; ++j) {
					eU1 += 4.0 / (pow(DX, 4)*pow(gp, 2)) * (5. / 3.*(1. - Plaquette(x, V, i, j, updateTime).real()) - 1. / 12.*(1 - PlaquetteRect12(x, V, i, j, updateTime).real()) - 1. / 12.*(1 - PlaquetteRect12(x, V, j, i, updateTime).real()));
					eSU2 += 4.0 / (pow(DX, 4)*pow(g, 2)) * (5. / 3.*(1. - 0.5*Plaquette(x, U, i, j, updateTime).trace()) - 1. / 12.*(1. - 0.5*PlaquetteRect12(x, U, i, j, updateTime).trace()) - 1. / 12.*(1. - 0.5*PlaquetteRect12(x, U, j, i, updateTime).trace()));
				}
			}
#ifdef TEST_MODE
			if (eScalar.imag()>ERROR) cout << "SCALAR FIELD ERROR: IMAGINARY ENERGY." << eScalar.imag() << endl;
			if (eSU2.imag()>ERROR) cout << "SU2 FIELD ERROR:ENERGY." << endl;
#endif
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = static_cast<RealSave> (DX*DX*DX * (eScalar.real() + eU1 + eSU2.real()));
			totalEnergy += DX * DX*DX * (eScalar.real() + eU1 + eSU2.real());
			arrVar[0] += eScalar.real()*DX*DX*DX;
			arrVar[1] += eU1 * DX*DX*DX;
			arrVar[2] += eSU2.real()*DX*DX*DX;
		}
	}
	return;
}
#endif
void CalcChernSimonsNumber(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2matrix>& U, Field<U1matrix>& V,
	const int updatedTimeStep,
	Real& totalCSnumber,
	Real* arrVar) {
	int updateTime = updatedTimeStep % CYCLE;
	totalCSnumber = 0;
	if (arrVar != NULL) arrVar[0] = 0;
	const Real pre_factor = FermiFamily / (32.0*PI*PI) * DX3;
	const Real g3 = 2.0*g*g*g / 3.0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			const Real SU2W[DIM][SU2DIM] = {
				{ SU2_W_SITE(x, U, X_AXIS, 0, updateTime), SU2_W_SITE(x, U, X_AXIS, 1, updateTime), SU2_W_SITE(x, U, X_AXIS, 2, updateTime) },
			{ SU2_W_SITE(x, U, Y_AXIS, 0, updateTime), SU2_W_SITE(x, U, Y_AXIS, 1, updateTime), SU2_W_SITE(x, U, Y_AXIS, 2, updateTime) },
			{ SU2_W_SITE(x, U, Z_AXIS, 0, updateTime), SU2_W_SITE(x, U, Z_AXIS, 1, updateTime), SU2_W_SITE(x, U, Z_AXIS, 2, updateTime) } };
			const Real U1Y[DIM] = {
				U1_Y_SITE(x, V, X_AXIS, updateTime),
				U1_Y_SITE(x, V, Y_AXIS, updateTime),
				U1_Y_SITE(x, V, Z_AXIS, updateTime)
			};
			Real U1B[DIM][DIM];
			Real SU2B[DIM][DIM][SU2DIM];
			for (int i = 0; i < DIM; ++i) {
				for (int j = i; j < DIM; ++j) {
					U1B[i][j] = U1_B_SITE(x, V, i, j, updateTime);
					U1B[j][i] = -U1B[i][j];
					for (int a = 0; a < SU2DIM; ++a) {
						SU2B[i][j][a] = SU2_B_SITE(x, U, i, j, a, updateTime);
						SU2B[j][i][a] = -SU2B[i][j][a];
					}
				}
			}
			for (int i = 0; i<DIM; ++i) {
				int j = (i + 1) % DIM;
				int k = (i + 2) % DIM;
				tempSave -= gp * gp*(U1B[i][j] * U1Y[k] - U1B[i][k] * U1Y[j]);
				for (int a = 0; a<SU2DIM; ++a) {
					int b = (a + 1) % DIM;
					int c = (a + 2) % DIM;
					tempSave += g * g*(SU2B[i][j][a] * SU2W[k][a] - SU2B[i][k][a] * SU2W[j][a]);
					if (arrVar != NULL) {
						Real localVar = -g3 * SU2W[i][a] * (SU2W[j][b] * SU2W[k][c] - SU2W[j][c] * SU2W[k][b]);
						arrVar[0] += localVar;
						tempSave += localVar;
					}
					else {
						tempSave -= g3 * SU2W[i][a] * (SU2W[j][b] * SU2W[k][c] - SU2W[j][c] * SU2W[k][b]);
					}
					//A modification here 03/03/16. should not change the result
					//tempSave += g*g * ((SU2_B_SITE(x, U, i, j, a, updateTime)*SU2_W_SITE(x, U, k, a, updateTime) + SU2_B_SITE(x, U, k, i, a, updateTime)*SU2_W_SITE(x, U, j, a, updateTime)) - 2.0*g / 3.0*(SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, b, updateTime)*SU2_W_SITE(x, U, k, c, updateTime) - SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, c, updateTime)*SU2_W_SITE(x, U, k, b, updateTime)));
				}
			}
			tempSave *= pre_factor;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
			totalCSnumber += tempSave;
		}
	}
	if (arrVar != NULL) arrVar[0] *= pre_factor;
	return;
}

//backup - IMPROVEMENT-07032017
/*
void CalcChernSimonsNumber(Site x,
Field<RealSave>& saveinfo, const int saveCol,
Field<SU2matrix>& U, Field<U1matrix>& V,
const int updatedTimeStep,
Real& totalCSnumber,
Real* arrVar) {
int updateTime = updatedTimeStep % CYCLE;
totalCSnumber = 0;
if (arrVar != NULL) arrVar[0] = 0;
const Real pre_factor = FermiFamily / (32.0*PI*PI) * DX3;
const Real g3 = 2.0*g*g*g / 3.0;
for (x.first(); x.test(); x.next()) {
Real tempSave = 0;
if (isBoundary(x)) {
if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
continue;
}
else {
for (int i = 0; i<DIM; ++i) {
int j = (i + 1) % DIM;
int k = (i + 2) % DIM;
tempSave -= gp*gp *(U1_B_SITE(x, V, i, j, updateTime)*U1_Y_SITE(x, V, k, updateTime) - U1_B_SITE(x, V, i, k, updateTime)*U1_Y_SITE(x, V, j, updateTime));
for (int a = 0; a<SU2DIM; ++a) {
int b = (a + 1) % DIM;
int c = (a + 2) % DIM;
tempSave += g*g * (SU2_B_SITE(x, U, i, j, a, updateTime)*SU2_W_SITE(x, U, k, a, updateTime) - SU2_B_SITE(x, U, i, k, a, updateTime)*SU2_W_SITE(x, U, j, a, updateTime));
if (arrVar != NULL) {
Real localVar = -g3*(SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, b, updateTime)*SU2_W_SITE(x, U, k, c, updateTime) - SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, c, updateTime)*SU2_W_SITE(x, U, k, b, updateTime));
arrVar[0] += localVar;
tempSave += localVar;
}
else {
tempSave -= g3*SU2_W_SITE(x, U, i, a, updateTime)*(SU2_W_SITE(x, U, j, b, updateTime)*SU2_W_SITE(x, U, k, c, updateTime) - SU2_W_SITE(x, U, j, c, updateTime)*SU2_W_SITE(x, U, k, b, updateTime));
}
//A modification here 03/03/16. should not change the result
//tempSave += g*g * ((SU2_B_SITE(x, U, i, j, a, updateTime)*SU2_W_SITE(x, U, k, a, updateTime) + SU2_B_SITE(x, U, k, i, a, updateTime)*SU2_W_SITE(x, U, j, a, updateTime)) - 2.0*g / 3.0*(SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, b, updateTime)*SU2_W_SITE(x, U, k, c, updateTime) - SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, c, updateTime)*SU2_W_SITE(x, U, k, b, updateTime)));
}
}
tempSave *= pre_factor;
if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
totalCSnumber += tempSave;
}
}
if (arrVar != NULL) arrVar[0] *= pre_factor;
return;
}
*/

void CalcTopoDegree_Gauge(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2matrix>& U,
	const int updatedTimeStep,
	Real& totalTopoDegree,
	Real* arrVar) {
	int updateTime = updatedTimeStep % CYCLE;
	totalTopoDegree = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			for (int i = 0; i<DIM; ++i) {
				int j = (i + 1) % DIM;
				int k = (i + 2) % DIM;
				for (int a = 0; a<SU2DIM; ++a) {
					int b = (a + 1) % DIM;
					int c = (a + 2) % DIM;
					tempSave -= 2.0*g*g*g / 3.0*(SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, b, updateTime)*SU2_W_SITE(x, U, k, c, updateTime) - SU2_W_SITE(x, U, i, a, updateTime)*SU2_W_SITE(x, U, j, c, updateTime)*SU2_W_SITE(x, U, k, b, updateTime));
				}
			}
			tempSave *= FermiFamily / (32.0*PI*PI) * pow(DX, 3);
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
			totalTopoDegree += tempSave;
		}
	}
	return;
}

void CalcTopoDegree_Higgs_2(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi,
	const int updatedTimeStep,
	Real& totalTopoDegree,
	Real* arrVar) {
	int updateTime = updatedTimeStep % CYCLE;
	totalTopoDegree = 0;
	for (x.first(); x.test(); x.next()) {
		Cmplx tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			for (int i = 0; i < DIM; ++i) {
				int j = (i + 1) % DIM;
				int k = (i + 2) % DIM;
				tempSave -= 1.0 / (24.0 * DX3) * (SU2invphi(x, updateTime, phi)*(SU2phi(x + i, updateTime, phi) - SU2phi(x - i, updateTime, phi)) *
					(SU2invphi(x, updateTime, phi)*(SU2phi(x + j, updateTime, phi) - SU2phi(x - j, updateTime, phi)) * SU2invphi(x, updateTime, phi)*(SU2phi(x + k, updateTime, phi) - SU2phi(x - k, updateTime, phi))
						- SU2invphi(x, updateTime, phi)*(SU2phi(x + k, updateTime, phi) - SU2phi(x - k, updateTime, phi)) * SU2invphi(x, updateTime, phi)*(SU2phi(x + j, updateTime, phi) - SU2phi(x - j, updateTime, phi)))).trace();
			}
			if (tempSave.imag() > ERROR) {
				COUT << "Higgs Topological Degree error: complex value." << endl;
				exit(EXIT_FAILURE);
			}
			tempSave *= FermiFamily / (8.0*PI*PI) * DX3;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave.real());
			if (!std::isnan(tempSave.real())) totalTopoDegree += tempSave.real();
		}
	}
	return;
}

void CalcHelicity(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi,
	Field<SU2matrix>& U,
	Field<U1matrix>& V,
	const int updatedTimeStep,
	Real& totalHelicity,
	const int isScalarTermIncluded,
	Real* arrVar) {
	//In this function, ALL boundary sites are neglected (therefore set to 0)
	int updateTime = updatedTimeStep % CYCLE;
	totalHelicity = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			for (int i = 0; i < DIM; ++i) {
				int j = (i + 1) % DIM;
				int k = (i + 2) % DIM;
				tempSave += EM_A_SITE(x, updateTime, i, phi, U, V) *EM_F_SITE(x, updateTime, j, k, phi, U, V);
			}
			tempSave *= DX3;
			totalHelicity += tempSave;
		}
		if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
	}
	return;
}

[[deprecated]] void CalcHelicityCurlA(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi,
	Field<SU2matrix>& U,
	Field<U1matrix>& V,
	const int updatedTimeStep,
	Real& totalHelicity,
	Real* arrVar) {
	//In this function, ALL boundary sites are neglected (therefore set to 0)
	int updateTime = updatedTimeStep % CYCLE;
	totalHelicity = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		int isBoundary = 0;
		for (int i = 0; i<DIM; ++i) {
			if (x.coord(i) == 0 || x.coord(i) == nSize[i] - 1) { ++isBoundary; break; }
		}
		if (isBoundary) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			for (int i = 0; i<DIM; ++i) {
				int j = (i + 1) % DIM;
				int k = (i + 2) % DIM;
				tempSave += EM_A_SITE(x, updateTime, i, phi, U, V)*
					((EM_A_SITE(x + j, updateTime, k, phi, U, V) - EM_A_SITE(x - j, updateTime, k, phi, U, V))
						- (EM_A_SITE(x + k, updateTime, j, phi, U, V) - EM_A_SITE(x - k, updateTime, j, phi, U, V)));
			}
			tempSave *= 0.5 / DX * pow(DX, 3);
			totalHelicity += tempSave;
		}
		if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
	}
	return;
}

void CalcEMmagneticEnergy(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi,
	Field<SU2matrix>& U,
	Field<U1matrix>& V,
	const int updatedTimeStep,
	Real& totalEMenergy,
	const int isScalarTermIncluded,
	Real* arrVar) {
	//In this function, ALL boundary sites are neglected (therefore set to 0)
	int updateTime = updatedTimeStep % CYCLE;
	totalEMenergy = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			tempSave += 0.5*(pow(EM_F_SITE(x, updateTime, 0, 1, phi, U, V), 2)
				+ pow(EM_F_SITE(x, updateTime, 1, 2, phi, U, V), 2)
				+ pow(EM_F_SITE(x, updateTime, 0, 2, phi, U, V), 2));
			tempSave *= DX3;
			totalEMenergy += tempSave;
		}
		if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
	}
	return;
}

void CalcScalarToDirectedVEV(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi, const int updatedTimeStep, Real& avgDev) {
	int updateTime = updatedTimeStep % CYCLE;
	avgDev = 0;
	SU2vector vev;
	vev << 0.0, v;
	Real siteDev = 0;
	for (x.first(); x.test(); x.next()) {
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			siteDev = (phi(x, updateTime) - vev).squaredNorm() / vev.squaredNorm();
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(siteDev);
			avgDev += siteDev;
		}
	}
	avgDev /= (nSize[0] - 2)*(nSize[1] - 2)*(nSize[2] - 2); //This value has no physical meaning, just a man-made measure
	return;
}

void CalcHiggsMag2(LATfield2::Site x,
	LATfield2::Field<RealSave>& saveinfo, const int saveCol,
	LATfield2::Field<SU2vector>& phi,
	const int calcTimeStep,
	Real& totalOutput) {
	const int T = calcTimeStep % CYCLE;
	totalOutput = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			tempSave = phi(x, T).squaredNorm();
			//if (tempSave < VECN_THRESHOLD) totalOutput += 1.0;
			totalOutput += tempSave;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
		}
	}
	return;
}

void CalcScalarDirectionComponent(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi, const int direction, const int calcTimeStep, Real& avgComp) {
	int calcTime = calcTimeStep % CYCLE;
	avgComp = 0;
	for (x.first(); x.test(); x.next()) {
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(vecN(phi(x, calcTime), direction));
		}
	}
	return;
}

void CalcEnergyDensity_Static(LATfield2::Site x,
	Field<RealSave>& saveData, const int saveCol,
	Field<SU2vector>& phi,
	Field<SU2matrix>& U,
	Field<U1matrix>& V,
	const int calcTimeStep,
	Real& totalEnergy,
	Real* arrVar) {
	//in this function, the energy is calculated as the current step, not the updated step.
	//The problem of this function is, when initial condition is not accurate, the calculated energy density may be not positive-definite.
	int nowTime = calcTimeStep % CYCLE;
	totalEnergy = 0;
	for (x.first(); x.test(); x.next()) {
		int isBoundary = 0;
		for (int j = 0; j<DIM; ++j) {
			if (x.coord(j) == 0 || x.coord(j) == nSize[j] - 1) { isBoundary = 1; break; }
		}
		// ###NOT including sites on boundaries###
		if (isBoundary > 0) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			Real eKinetic = 0, eGrad = 0, ePotential = 0, eU1 = 0;
			Cmplx eSU2(0, 0);
			ePotential = Vphi(x, phi, nowTime);
			for (int i = 0; i<DIM; ++i) {
				eGrad += 1.0 / (DX*DX)*(U(x, i, nowTime)*V(x, i, nowTime)*phi(x + i, nowTime) - phi(x, nowTime)).squaredNorm();
				for (int j = 0; j < DIM; ++j) {
					eU1 += 2.0 / (pow(DX, 4)*pow(gp, 2)) * (1.0 - Plaquette(x, V, i, j, nowTime).real());
					eSU2 += 1.0 / (pow(DX, 4)*pow(g, 2)) * (2.0 - Plaquette(x, U, i, j, nowTime).trace());
				}
			}
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = static_cast<RealSave> (DX*DX*DX * (eKinetic + eGrad + ePotential + eU1 + eSU2.real()));
			//saveData(x, saveCol) = CovDphi(x, nowTime, 2, phi, U, V).norm() < CONSTRAINTS_THRESHOLD_D ? 1 : 0;
			//saveData(x, saveCol) = phi(x,nowTime).norm() >= CONSTRAINTS_THRESHOLD ? 1 : 0;
			totalEnergy += eKinetic + ePotential + eGrad + eU1 + eSU2.real();
			if (arrVar != NULL) {
				arrVar[0] += eKinetic;
				arrVar[1] += eGrad;
				arrVar[2] += ePotential;
				arrVar[3] += eU1;
				arrVar[4] += eSU2.real();
			}
		}
	}
	totalEnergy *= DX * DX*DX;
	if (arrVar != NULL) {
		arrVar[0] *= DX * DX*DX;
		arrVar[1] *= DX * DX*DX;
		arrVar[2] *= DX * DX*DX;
		arrVar[3] *= DX * DX*DX;
		arrVar[4] *= DX * DX*DX;
	}
	return;
}

void CalcEMmagneticField(LATfield2::Site x,
	LATfield2::Field<RealSave>& saveData, const int saveCol,
	const int dir,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const int calcTimeStep,
	Real& totalMagneticEnergy) {
	int nowTime = calcTimeStep % CYCLE;
	totalMagneticEnergy = 0;
	for (x.first(); x.test(); x.next()) {
		// ###NOT including sites on boundaries###
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			const int i = (dir + 1) % DIM;
			const int j = (dir + 2) % DIM;
			saveData(x, saveCol) = -EM_F_SITE(x, nowTime, i, j, phi, U, V);
		}
	}
	return;
}

void CalcEMelectricField(LATfield2::Site x,
	LATfield2::Field<RealSave>& saveData, const int saveCol,
	const int dir,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const int calcTimeStep,
	Real& total) {
	int nowTime = calcTimeStep % CYCLE;
	total = 0;
	for (x.first(); x.test(); x.next()) {
		// ###NOT including sites on boundaries###
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			const int i = (dir + 1) % DIM;
			const int j = (dir + 2) % DIM;
			saveData(x, saveCol) = EM_E_SITE(x, nowTime, dir, phi, U, V, pi, F, E);
		}
	}
	return;
}

void CalcZmagneticField(LATfield2::Site x,
	LATfield2::Field<RealSave>& saveData, const int saveCol,
	const int dir,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const int calcTimeStep,
	Real& totalMagneticEnergy) {
	int nowTime = calcTimeStep % CYCLE;
	totalMagneticEnergy = 0;
	for (x.first(); x.test(); x.next()) {
		// ###NOT including sites on boundaries###
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveData(x, saveCol) = 0;
		}
		else {
			const int i = (dir + 1) % DIM;
			const int j = (dir + 2) % DIM;
			saveData(x, saveCol) = -Z_F_SITE(x, nowTime, i, j, phi, U, V);
		}
	}
	return;
}

void CalcMinHiggs2(LATfield2::Site x,
	LATfield2::Field<RealSave>& saveData, const int saveCol,
	LATfield2::Field<SU2vector>& phi,
	const int calcTimeStep,
	Real& min) {
	int nowTime = calcTimeStep % CYCLE;
	min = v2;
	for (x.first(); x.test(); x.next()) {
		if (isBoundary(x)) {
			;
		}
		else {
			const auto phi2 = phi(x, nowTime).squaredNorm();
			min = phi2 < min ? phi2 : min;
		}
	}
	return;
}

void TestFields(Site x,
	Field<SU2vector>& phi, Field<SU2vector>& pi,
	Field<SU2matrix>& U, Field<SU2matrix>& F,
	Field<U1matrix>& V, Field<U1matrix>& E,
	const int updatedTimeStep) {
	//in this function, the energy is calculated as the current step, not the updated step.
	//The problem of this function is, when initial condition is not accurate, the calculated energy density may be not positive-definite.
	int nowTime = updatedTimeStep % CYCLE;
	std::cout << "Test Fields." << std::endl;
	for (x.first(); x.test(); x.next()) {
		for (int i = 0; i < DIM; ++i) {
			if (abs(E(x, i, nowTime).imag())>0) {
				COUT << "( " << x.coord(0) << " , " << x.coord(1) << " , " << x.coord(2) << " )" << endl;
			}
		}
	}
	return;
}

void ErrorInfo(LATfield2::Site x,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const int nowTime,
	const int direction) {
	cerr << "Lattice Coordinate: (" << x.coord(0) << "," << x.coord(1) << "," << x.coord(2) << ")" << endl;
	cerr << "Physical Coordinate: (" << rx(x, 0) << "," << rx(x, 1) << "," << rx(x, 2) << ")" << endl;
	cerr << "Higgs: " << endl;
	cerr << phi(x, nowTime) << endl;
	if (direction != -1) {
		Site y = x + direction;
		cerr << "Connecting Lattice Coordinate: (" << y.coord(0) << "," << y.coord(1) << "," << y.coord(2) << ")" << endl;
		cerr << "Connecting Physical Coordinate: (" << rx(y, 0) << "," << rx(y, 1) << "," << rx(y, 2) << ")" << endl;
	}
	cerr << "SU2 Field: " << endl;
	cerr << "Direction: " << direction << endl;
	cerr << U(x, direction, nowTime) << endl;
	cerr << "W_i^1: " << SU2_W(x, U, direction, 0, nowTime) << endl;
	cerr << "W_i^2: " << SU2_W(x, U, direction, 1, nowTime) << endl;
	cerr << "W_i^3: " << SU2_W(x, U, direction, 2, nowTime) << endl;
	return;
}

void CalcHelicityFlux(Site x,
	Field<RealSave>& saveinfo, const int saveCol,
	Field<SU2vector>& phi, Field<SU2vector>& pi,
	Field<SU2matrix>& U, Field<SU2matrix>& F,
	Field<U1matrix>& V, Field<U1matrix>& E,
	const int calcTimeStep,
	Real& totalHelicityFlux,
	const int isScalarTermIncluded,
	Real* arrVar) {
	//The flux is calculated on the layer next to the boundary layer.
	int calcTime = calcTimeStep % CYCLE;
	totalHelicityFlux = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		int isCalcLayer = 0;
		int surfDir[DIM] = { 0,0,0 };
		for (int i = 0; i<DIM; ++i) {
			if (x.coord(i) == 1) { ++isCalcLayer; surfDir[i] = -1; }
			if (x.coord(i) == nSize[i] - 2) { ++isCalcLayer; surfDir[i] = 1; }
		}
		if (isCalcLayer == 1) {
			for (int i = 0; i<DIM; ++i) {
				if (surfDir[i] == 0) continue;
				else {
					int j = (i + 1) % DIM;
					int k = (i + 2) % DIM;
					tempSave += surfDir[i] * (EM_A_SITE_EW(x, calcTime, j, phi, U, V)*EM_E_SITE_EW(x, calcTime, k, phi, U, V, pi, F, E, isScalarTermIncluded) - EM_A_SITE_EW(x, calcTime, k, phi, U, V)*EM_E_SITE_EW(x, calcTime, j, phi, U, V, pi, F, E, isScalarTermIncluded));
				}
			}
			tempSave *= DX * DX * DT;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
			totalHelicityFlux += tempSave;
		}
	}
	return;
}

void CalcScalarQuantity(
	const Real(*scaSingleCell)(const LATfield2::Site& x, const int T,
		LATfield2::Field<SU2vector>& phi,
		LATfield2::Field<SU2matrix>& U,
		LATfield2::Field<U1matrix>& V),
	LATfield2::Site x,
	LATfield2::Field<RealSave>& saveinfo, const int saveCol,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const int calcTimeStep,
	Real& totalOutput) {
	const int T = calcTimeStep % CYCLE;
	totalOutput = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			tempSave = scaSingleCell(x, T, phi, U, V);
			totalOutput += tempSave;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
		}
	}
	totalOutput *= DX3;
	return;
}

void CalcScalarQuantity(
	const Real(*scaSingleCell)(const LATfield2::Site& x, const int T,
		LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
		LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
		LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E),
	LATfield2::Site x,
	LATfield2::Field<RealSave>& saveinfo, const int saveCol,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const int calcTimeStep,
	Real& totalOutput) {
	const int T = calcTimeStep % CYCLE;
	totalOutput = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			tempSave = scaSingleCell(x, T, phi, pi, U, F, V, E);
			totalOutput += tempSave;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
		}
	}
	totalOutput *= DX3;
	return;
}

void CalcVectorQuantity(
	const Real(*vecSingleCell)(const LATfield2::Site& x, const int T, const int dir, LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
		LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
		LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E),
	LATfield2::Site x,
	LATfield2::Field<RealSave>& saveinfo, const int saveCol,
	const int dir,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const int calcTimeStep,
	Real& totalOutput) {
	const int T = calcTimeStep % CYCLE;
	totalOutput = 0;
	for (x.first(); x.test(); x.next()) {
		Real tempSave = 0;
		if (isBoundary(x)) {
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = 0.;
			continue;
		}
		else {
			tempSave = vecSingleCell(x, T, dir, phi, pi, U, F, V, E);
			totalOutput += tempSave;
			if (saveCol != SaveInfo::UNSAVE) saveinfo(x, saveCol) = static_cast<RealSave>(tempSave);
		}
	}
	totalOutput *= DX3;
	return;
}

//---Single-cell quantity---
//Gauss Constraints
const Real GaussU1(const LATfield2::Site& x, const int cycleTime,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<U1matrix>& E) {
	if (isBoundary(x)) { return 0.0; }
	else {
		Real u1term = 0;
		for (int i = 0; i < DIM; ++i) {
			u1term += E(x, i, cycleTime).imag() - E(x - i, i, cycleTime).imag();
		}
		return u1term / DX - gp * (pi(x, cycleTime).dot(phi(x, cycleTime))).imag();
	}
}
const Real GaussSU2(const LATfield2::Site& x, const int cycleTime,
	const int su2dir,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F) {
	if (isBoundary(x)) { return 0.0; }
	else {
		Real su2term = 0;
		for (int i = 0; i < DIM; ++i) {
			su2term += (I*Pauli[su2dir] * (F(x, i, cycleTime) - U(x - i, i, cycleTime).adjoint()*F(x - i, i, cycleTime)*U(x - i, i, cycleTime))).trace().real();
		}
		return su2term / DX - g * (pi(x, cycleTime).dot(I*Pauli[su2dir] * phi(x, cycleTime))).real();
	}
}
const Real HamiltonianGC(const LATfield2::Site& x, const int cycleTime,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return 0.5*(pow(GaussU1(x, cycleTime, phi, pi, E), 2) + pow(GaussSU2(x, cycleTime, 0, phi, pi, U, F), 2) + pow(GaussSU2(x, cycleTime, 1, phi, pi, U, F), 2) + pow(GaussSU2(x, cycleTime, 2, phi, pi, U, F), 2));
}
const Real EMflux_SingleCell(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V) {
	const int&& dirX = 0;
	const int&& dirY = 1;
	const int&& dirZ = 2;
	const Real&& fluxZ = EM_F_PLAQ(x + dirZ, T, dirX, dirY, phi, U, V) - EM_F_PLAQ(x, T, dirX, dirY, phi, U, V);
	const Real&& fluxX = EM_F_PLAQ(x + dirX, T, dirY, dirZ, phi, U, V) - EM_F_PLAQ(x, T, dirY, dirZ, phi, U, V);
	const Real&& fluxY = EM_F_PLAQ(x + dirY, T, dirZ, dirX, phi, U, V) - EM_F_PLAQ(x, T, dirZ, dirX, phi, U, V);
	return (fluxX + fluxY + fluxZ)*DX2;
}
const Real EM_electricEnergy(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return 0.5*(
		pow(EM_E_SITE(x, T, 0, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE(x, T, 1, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE(x, T, 2, phi, U, V, pi, F, E), 2));
}
const Real EM_electricEnergy_noScalar(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return 0.5*(
		pow(EM_E_SITE_noScalar(x, T, 0, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE_noScalar(x, T, 1, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE_noScalar(x, T, 2, phi, U, V, pi, F, E), 2));
}
const Real EM_electricEnergy_pureScalar(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return 0.5*(
		pow(EM_E_SITE_pureScalar(x, T, 0, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE_pureScalar(x, T, 1, phi, U, V, pi, F, E), 2) +
		pow(EM_E_SITE_pureScalar(x, T, 2, phi, U, V, pi, F, E), 2));
}
const Real EM_magneticEnergy(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return  0.5*(pow(EM_F_SITE(x, T, 0, 1, phi, U, V), 2)
		+ pow(EM_F_SITE(x, T, 1, 2, phi, U, V), 2)
		+ pow(EM_F_SITE(x, T, 0, 2, phi, U, V), 2));
}
const Real EM_magneticEnergy_noScalar(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return  0.5*(pow(EM_F_SITE_VAC(x, T, 0, 1, phi, U, V), 2)
		+ pow(EM_F_SITE_VAC(x, T, 1, 2, phi, U, V), 2)
		+ pow(EM_F_SITE_VAC(x, T, 0, 2, phi, U, V), 2));
}
const Real EM_magneticEnergy_pureScalar(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	return  0.5*(pow(EM_F_SITE_pureScalar(x, T, 0, 1, phi, U, V), 2)
		+ pow(EM_F_SITE_pureScalar(x, T, 1, 2, phi, U, V), 2)
		+ pow(EM_F_SITE_pureScalar(x, T, 0, 2, phi, U, V), 2));
}

const Real MomentumDensity(const LATfield2::Site& x, const int T, const int dir,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	const int Tf = (T + 1) % CYCLE;
	Real scalarPart = (pi(x, T) + pi(x, Tf)).dot(CovDphi_Sym(x, T, dir, phi, U, V)).real();
	Real gaugePart = 0;
	for (int k = 0; k < DIM; ++k) {
		Real su2GaugePart = 0;
		for (int a = 0; a < SU2DIM; ++a) {
			su2GaugePart += 0.5*(SU2_E_SITE(x, F, k, a, T) + SU2_E_SITE(x, F, k, a, Tf))*SU2_B_SITE(x, U, k, dir, a, T);
		}
		Real u1GaugePart = 0.5*(U1_E_SITE(x, E, k, T) + U1_E_SITE(x, E, k, Tf))*U1_B_SITE(x, V, k, dir, T);
		gaugePart += u1GaugePart + su2GaugePart;
	}
	return scalarPart + gaugePart;
}
const Real AngularMomentumDensity(const LATfield2::Site& x, const int T, const int dir,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	const int i = (dir + 1) % DIM;
	const int j = (dir + 2) % DIM;
	return rx(x, i)*MomentumDensity(x, T, j, phi, pi, U, F, V, E) - rx(x, j)*MomentumDensity(x, T, i, phi, pi, U, F, V, E);
}

const Real HiggsWinding(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V) {
	Cmplx hw = 0;
	const SU2matrix Uphi_p[DIM] = {
		SU2invphi(x + X_AXIS,T,phi),
		SU2invphi(x + Y_AXIS,T,phi),
		SU2invphi(x + Z_AXIS,T,phi) };
	const SU2matrix Uphi_m[DIM] = {
		SU2invphi(x - X_AXIS,T,phi),
		SU2invphi(x - Y_AXIS,T,phi),
		SU2invphi(x - Z_AXIS,T,phi) };
	const SU2matrix&& Uphi = SU2invphi(x, T, phi);
	for (int i = 0; i < DIM; ++i) {
		const int j = (i + 1) % DIM;
		const int k = (i + 2) % DIM;
		hw += ((Uphi_p[i] - Uphi_m[i]).adjoint()*(
			(Uphi_p[j] - Uphi_m[j])*(Uphi_p[k] - Uphi_m[k]).adjoint() -
			(Uphi_p[k] - Uphi_m[k])*(Uphi_p[j] - Uphi_m[j]).adjoint()
			)*Uphi).trace();
	}
	hw *= 1.0 / (24.0*8.0*PI*PI*DX3);
	if (std::isnan(hw.real())) {
		return 0.0;
	}
	else {
		//somehow there needs to be a minus sign.
		return hw.real();
	}
}

const Real Higgs2VEV2(const LATfield2::Site& x, const int T,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V) {
	return (phi(x, T).squaredNorm() - VEV2()) / VEV2();
}

//---Miscellaneous---
Real radius(const LATfield2::Site& x, const Real dx, const Real* center) {
	const Real&& rad = pow((x.coord(0) - center[0]), 2) + pow((x.coord(1) - center[1]), 2) + pow((x.coord(2) - center[2]), 2);
	return sqrt(rad)*dx;
}

Real EM_A_SITE_EW(const LATfield2::Site& x, const int cycleTime, const int dir_i,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	Cmplx Apotential(0, 0);
	for (int a = 0; a < SU2DIM; ++a) {
		Apotential += vecN(phi(x, cycleTime), a)*SU2_W_SITE(x, U, dir_i, a, cycleTime);
	}
	Apotential *= sinw;
	Apotential += cosw * U1_Y_SITE(x, V, dir_i, cycleTime);
#ifdef TEST_MODE
	if (abs(Apotential.imag()) > ERROR) cerr << "EM POTENTIAL ERROR: NON-REAL VALUE." << endl;
#endif
	return Apotential.real();
}

//Real EM_A_SITE_EW_V2(const LATfield2::Site& x, const int cycleTime, const int dir_i,
//	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
//	Cmplx Apotential(0, 0);
//	for (int a = 0; a < SU2DIM; ++a) {
//		Apotential += vecN_VAC(phi(x, cycleTime), a)*SU2_W_SITE(x, U, dir_i, a, cycleTime);
//	}
//	Apotential *= sinw;
//	Apotential += cosw*U1_Y_SITE(x, V, dir_i, cycleTime);
//#ifdef TEST_MODE
//	if (abs(Apotential.imag()) > ERROR) cerr << "EM POTENTIAL ERROR: NON-REAL VALUE." << endl;
//#endif
//	return Apotential.real();
//}

Real EM_A_SITE_HOPF(const LATfield2::Site& x, const int cycleTime, const int dir_i,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	Cmplx Apotential(0, 0);
	for (int a = 0; a < SU2DIM; ++a) {
		Apotential += -phi(x, cycleTime).dot(Pauli[a] * phi(x, cycleTime)) / phi(x, cycleTime).squaredNorm()*0.5*(SU2_W(x, U, dir_i, a, cycleTime) + SU2_W(x - dir_i, U, dir_i, a, cycleTime));
	}
	Apotential *= sinw;
	Apotential += cosw * 0.5*(U1_Y(x, V, dir_i, cycleTime) + U1_Y(x - dir_i, V, dir_i, cycleTime));
	Apotential += I / DX * sinw / g * (phi(x + dir_i, cycleTime) - phi(x - dir_i, cycleTime)).dot(phi(x, cycleTime)); //Hopf curvature induced potential.
#ifdef TEST_MODE
	if (abs(Apotential.imag()) > ERROR) cerr << "EM POTENTIAL ERROR: NON-REAL VALUE." << Apotential.imag() << endl;
#endif
	return Apotential.real();
}

Real EM_A_LINK_SU2INV(const LATfield2::Site& x, const int T, const int dir_i,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	Cmplx Apotential(0, 0);
	Apotential = I * sinw / g * (CovDphi(x, T, dir_i, phi, U, V).dot(phi(x, T)) - phi(x, T).dot(CovDphi(x, T, dir_i, phi, U, V)));
	Apotential += U1_Y(x, V, dir_i, T) / cosw;
	return Apotential.real();
}

Real EM_A_SITE_SU2INV(const LATfield2::Site& x, const int T, const int dir_i,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	return 0.5*(EM_A_LINK_SU2INV(x, T, dir_i, phi, U, V) + EM_A_LINK_SU2INV(x - dir_i, T, dir_i, phi, U, V));
}

Real EM_F_SITE_TanmayEW(const LATfield2::Site& x, const int nowTime,
	const int dir_i, const int dir_j,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V,
	const int isScalarTermIncluded) {
	Cmplx Fstrength(0, 0);
	for (int a = 0; a<SU2DIM; ++a) {
		Fstrength += vecN(phi(x, nowTime), a) * SU2_B_SITE(x, U, dir_i, dir_j, a, nowTime);
	}
	Fstrength *= sinw;
	Fstrength += cosw * U1_B_SITE(x, V, dir_i, dir_j, nowTime);
	if (isScalarTermIncluded) {
		Fstrength += (-2.0) * I / g * sinw * (CovDphi_Sym(x, nowTime, dir_i, phi, U, V).dot(CovDphi_Sym(x, nowTime, dir_j, phi, U, V)) - CovDphi_Sym(x, nowTime, dir_j, phi, U, V).dot(CovDphi_Sym(x, nowTime, dir_i, phi, U, V))) / phiSqNorm(phi(x, nowTime));
	}
#ifdef TEST_MODE
	if (abs(Fstrength.imag()) > ERROR) cerr << "EM FIELD STRENGTH ERROR: NON-REAL VALUE." << endl;
#endif
	return Fstrength.real();
}

Real EM_E_SITE_EW(const LATfield2::Site& x, const int cycleTime, const int dir_i,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V,
	Field<SU2vector>& pi, Field<SU2matrix>& F, Field<U1matrix>& E,
	const int isScalarCurrentIncluded) {
	//This is a little inconsistent, bub EM electric field is defined at t moments, for convenience.
	int nowTime = cycleTime;
	int preTime = (cycleTime + 1) % CYCLE;
	Cmplx efield = 0.;
	for (int a = 0; a<SU2DIM; ++a) {
		efield += 0.25* vecN(phi(x, nowTime), a) * (-iPauli[a] * (F(x, dir_i, nowTime) + F(x, dir_i, preTime) + F(x - dir_i, dir_i, nowTime) + F(x - dir_i, dir_i, preTime))).trace();
	}
	efield *= sinw;
	efield += 0.25*cosw * (E(x, dir_i, nowTime) + E(x - dir_i, dir_i, nowTime) + E(x, dir_i, preTime) + E(x - dir_i, dir_i, preTime)).imag();
	if (isScalarCurrentIncluded) {
		efield += (-I)*2.0 * sinw / g * (pi(x, nowTime).dot(CovDphi_Sym(x, nowTime, dir_i, phi, U, V)) - CovDphi_Sym(x, nowTime, dir_i, phi, U, V).dot(pi(x, nowTime))) / phiSqNorm(phi(x, nowTime));
	}
	return efield.real();
}

Real EM_F_PLAQ_TanmayEW(const LATfield2::Site& x,
	const int T, const int dir_i, const int dir_j,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2matrix>& U, LATfield2::Field<U1matrix>& V,
	const int isScalarCurrentIncluded) {
	const SU2vector cPhi = plaqPhi(x, dir_i, dir_j, phi, T);
	Cmplx Fstrength(0, 0);
	for (int a = 0; a<SU2DIM; ++a) {
		Fstrength += vecN(cPhi, a) * SU2_B(x, U, dir_i, dir_j, a, T);
	}
	Fstrength *= sinw;
	Fstrength += cosw * U1_B(x, V, dir_i, dir_j, T);
	if (isScalarCurrentIncluded) {
		Fstrength -= 2.0 * I / g * sinw * (CovDphi_Plaq(x, T, dir_i, dir_j, phi, U, V).dot(CovDphi_Plaq(x, T, dir_j, dir_i, phi, U, V)) - CovDphi_Plaq(x, T, dir_j, dir_i, phi, U, V).dot(CovDphi_Plaq(x, T, dir_i, dir_j, phi, U, V))) / phiSqNorm(cPhi);
	}
	return Fstrength.real();
}

//Real EM_F_SITE_TanmayEW_V2(const LATfield2::Site& x, const int T, const int dir_i, const int dir_j,
//	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2matrix>& U, LATfield2::Field<U1matrix>& V,
//	const int isScalarCurrentIncluded) {
//	Cmplx Fstrength(0, 0);
//	for (int a = 0; a<SU2DIM; ++a) {
//		Fstrength += vecN_VAC(phi(x, T), a) * SU2_B_SITE(x, U, dir_i, dir_j, a, T);
//	}
//	Fstrength *= sinw;
//	Fstrength += cosw*U1_B_SITE(x, V, dir_i, dir_j, T);
//	if (isScalarCurrentIncluded) {
//		Fstrength += (-2.0) * I / g * sinw * (CovDphi_Sym(x, T, dir_i, phi, U, V).dot(CovDphi_Sym(x, T, dir_j, phi, U, V)) - CovDphi_Sym(x, T, dir_j, phi, U, V).dot(CovDphi_Sym(x, T, dir_i, phi, U, V))) / v2;
//	}
//#ifdef TEST_MODE
//	if (abs(Fstrength.imag()) > ERROR) cerr << "EM FIELD STRENGTH ERROR: NON-REAL VALUE." << endl;
//#endif
//	return Fstrength.real();
//}

Real EM_F_SITE_SU2INV(const LATfield2::Site& x, const int T,
	const int dir_i, const int dir_j,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	return ((EM_A_LINK_SU2INV(x, T, dir_j, phi, U, V) - EM_A_LINK_SU2INV(x - dir_i, T, dir_j, phi, U, V)) - (EM_A_LINK_SU2INV(x, T, dir_i, phi, U, V) - EM_A_LINK_SU2INV(x - dir_j, T, dir_i, phi, U, V))) / DX;
}

Real Z_F_SITE_VAC(const LATfield2::Site& x, const int nowTime,
	const int dir_i, const int dir_j,
	Field<SU2vector>& phi, Field<SU2matrix>& U, Field<U1matrix>& V) {
	Cmplx Fstrength(0, 0);
	for (int a = 0; a<SU2DIM; ++a) {
		Fstrength += vecN(phi(x, nowTime), a) * SU2_B_SITE(x, U, dir_i, dir_j, a, nowTime);
	}
	Fstrength *= cosw;
	Fstrength -= sinw * U1_B_SITE(x, V, dir_i, dir_j, nowTime);
	return Fstrength.real();
}


//derivatives

//---Pauli Matrices---
SU2matrix PauliSph(LATfield2::Site x, const int su2dir, const int linkDir) {
	Real cart[DIM];
	Real sph[DIM];
	for (int i = 0; i < DIM; ++i) {
		cart[i] = rx(x, i);
		if (i == linkDir) cart[i] += 0.5*DX;
	}
	Cart2Sph(cart, sph);
	Real& theta = sph[1];
	Real& psi = sph[2];
	SU2matrix pauli;
	if (su2dir == 0) {
		pauli(0, 0) = cos(theta);
		pauli(0, 1) = sin(theta)*(cos(psi) - I * sin(psi));
		pauli(1, 0) = conj(pauli(0, 1));
		pauli(1, 1) = -pauli(0, 0);
	}
	if (su2dir == 1) {
		pauli(0, 0) = sin(theta);
		pauli(0, 1) = cos(theta)*(cos(psi) - I * sin(psi));
		pauli(1, 0) = conj(pauli(0, 1));
		pauli(1, 1) = -pauli(0, 0);
	}
	if (su2dir == 2) {
		pauli(0, 0) = 0.;
		pauli(0, 1) = -I * (cos(psi) - I * sin(psi));
		pauli(1, 0) = conj(pauli(0, 1));
		pauli(1, 1) = -pauli(0, 0);
	}
	return pauli;
}

//---File functions---
std::string make_datafile_name(std::string name, int num)
{
	std::ostringstream file;
	if (num == -1) {
		file << name << ".dat";
	}
	else {
		file << name << "_" << std::setw(3) << std::setfill('0') << num << ".dat";
	}
	return file.str();
}
std::string make_datfast_name(std::string name, int num)
{
	std::ostringstream file;
	if (num == -1) {
		file << name << ".datfast";
	}
	else {
		file << name << "_" << std::setw(3) << std::setfill('0') << num << ".datfast";
	}
	return file.str();
}

void outputFields(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const int cycleTime, const std::string id, int isSaveDiskSpace) {
	if (isSaveDiskSpace) {
		Site x(lat);
		Field<Real> opU(lat, DIM*SU2DIM, CYCLE); //0-2: x direction; 3-5: y direction; 6-8: z direction.
		Field<Real> opF(lat, DIM*SU2DIM, CYCLE);
		Field<Real> opV(lat, DIM, CYCLE);
		Field<Real> opE(lat, DIM, CYCLE);
		for (x.first(); x.test(); x.next()) {
			for (int cyc = 0; cyc < CYCLE; ++cyc) {
				for (int i = 0; i < DIM; ++i) {
					opV(x, i, cyc) = V(x, i, cyc).imag();
					opE(x, i, cyc) = E(x, i, cyc).imag();
					opU(x, SU2DIM*i + 0, cyc) = (I*Pauli[0] * U(x, i, cyc)).trace().real();
					opU(x, SU2DIM*i + 1, cyc) = (I*Pauli[1] * U(x, i, cyc)).trace().real();
					opU(x, SU2DIM*i + 2, cyc) = (I*Pauli[2] * U(x, i, cyc)).trace().real();
					opF(x, SU2DIM*i + 0, cyc) = (I*Pauli[0] * F(x, i, cyc)).trace().real();
					opF(x, SU2DIM*i + 1, cyc) = (I*Pauli[1] * F(x, i, cyc)).trace().real();
					opF(x, SU2DIM*i + 2, cyc) = (I*Pauli[2] * F(x, i, cyc)).trace().real();
					//TODO
				}
			}
		}
		phi.fastwrite(make_datafile_name("sphi_", cycleTime)); //s for small file size.
		pi.fastwrite(make_datafile_name("spi_", cycleTime));
		opU.fastwrite(make_datafile_name("sU_", cycleTime));
		opV.fastwrite(make_datafile_name("sV_", cycleTime));
		opF.fastwrite(make_datafile_name("sF_", cycleTime));
		opE.fastwrite(make_datafile_name("sE_", cycleTime));
	}
	else {
		phi.fastwrite(make_datafile_name(id + "phi_", cycleTime));
		pi.fastwrite(make_datafile_name(id + "pi_", cycleTime));
		U.fastwrite(make_datafile_name(id + "U_", cycleTime));
		V.fastwrite(make_datafile_name(id + "V_", cycleTime));
		F.fastwrite(make_datafile_name(id + "F_", cycleTime));
		E.fastwrite(make_datafile_name(id + "E_", cycleTime));
	}
	return;
}

void inputFields(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const int cycleTime, const std::string& prefix, const int isSaveDiskSpace) {
	if (isSaveDiskSpace) {
		Site x(lat);
		Field<Real> inU(lat, DIM*SU2DIM, CYCLE); //0-2: x direction; 3-5: y direction; 6-8: z direction.
		Field<Real> inF(lat, DIM*SU2DIM, CYCLE);
		Field<Real> inV(lat, DIM, CYCLE);
		Field<Real> inE(lat, DIM, CYCLE);
		phi.read(make_datfast_name("sort_sphi_", cycleTime));
		pi.read(make_datfast_name("sort_spi_", cycleTime));
		inU.read(make_datfast_name("sort_sU_", cycleTime));
		inF.read(make_datfast_name("sort_sF_", cycleTime));
		inV.read(make_datfast_name("sort_sV_", cycleTime));
		inE.read(make_datfast_name("sort_sE_", cycleTime));
		//TODO
	}
	else {
		phi.read(make_datfast_name(prefix + "phi_", cycleTime));
		pi.read(make_datfast_name(prefix + "pi_", cycleTime));
		U.read(make_datfast_name(prefix + "U_", cycleTime));
		V.read(make_datfast_name(prefix + "V_", cycleTime));
		F.read(make_datfast_name(prefix + "F_", cycleTime));
		E.read(make_datfast_name(prefix + "E_", cycleTime));
	}
	phi.updateHalo();
	pi.updateHalo();
	U.updateHalo();
	V.updateHalo();
	F.updateHalo();
	E.updateHalo();
	return;
}

void SimpleInputFields(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E) {
	phi.read(ReadConfigIni("FILE_PATH_PHI"));
	pi.read(ReadConfigIni("FILE_PATH_PI"));
	U.read(ReadConfigIni("FILE_PATH_U"));
	V.read(ReadConfigIni("FILE_PATH_V"));
	F.read(ReadConfigIni("FILE_PATH_F"));
	E.read(ReadConfigIni("FILE_PATH_E"));
	phi.updateHalo();
	pi.updateHalo();
	U.updateHalo();
	V.updateHalo();
	F.updateHalo();
	E.updateHalo();
	return;
}
void SimpleInputFields_phiUV(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V) {
	phi.read(ReadConfigIni("FILE_PATH_PHI"));
	U.read(ReadConfigIni("FILE_PATH_U"));
	V.read(ReadConfigIni("FILE_PATH_V"));
	phi.updateHalo();
	U.updateHalo();
	V.updateHalo();
	return;
}
void SimpleInputFields_piEF(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& E) {
	pi.read(ReadConfigIni("FILE_PATH_PI"));
	F.read(ReadConfigIni("FILE_PATH_F"));
	E.read(ReadConfigIni("FILE_PATH_E"));
	pi.updateHalo();
	F.updateHalo();
	E.updateHalo();
	return;
}
void InputFieldsToMultiCycle(LATfield2::Lattice& lat,
	LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const std::string&& inikey_phi,
	const std::string&& inikey_U,
	const std::string&& inikey_V) {
	LATfield2::Site x(lat);
	LATfield2::Field<SU2vector> read_phi(lat, 1);
	LATfield2::Field<SU2matrix> read_U(lat, DIM, 1);
	LATfield2::Field<U1matrix> read_V(lat, DIM, 1);
	read_phi.read(ReadConfigIni(inikey_phi));
	read_U.read(ReadConfigIni(inikey_U));
	read_V.read(ReadConfigIni(inikey_V));
	for (x.first(); x.test(); x.next()) {
		for (int c = 0; c < CYCLE; ++c) {
			phi(x, c) = read_phi(x, 0);
			for (int i = 0; i < DIM; ++i) {
				U(x, i, c) = read_U(x, i, 0);
				V(x, i, c) = read_V(x, i, 0);
			}
		}
	}
	read_phi.~Field<SU2vector>();
	read_U.~Field<SU2matrix>();
	read_V.~Field<U1matrix>();
	phi.updateHalo();
	U.updateHalo();
	V.updateHalo();
	return;
}
void ParallelOutputFields(LATfield2::Field<SU2vector>& phi, LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& U, LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& V, LATfield2::Field<U1matrix>& E,
	const std::string& prefix, const int time) {
	phi.fastwrite(make_datafile_name(prefix + "phi", time));
	U.fastwrite(make_datafile_name(prefix + "U", time));
	V.fastwrite(make_datafile_name(prefix + "V", time));
	pi.fastwrite(make_datafile_name(prefix + "pi", time));
	F.fastwrite(make_datafile_name(prefix + "F", time));
	E.fastwrite(make_datafile_name(prefix + "E", time));
	return;
}
void ParallelOutputFields_phiUV(LATfield2::Field<SU2vector>& phi,
	LATfield2::Field<SU2matrix>& U,
	LATfield2::Field<U1matrix>& V,
	const std::string& prefix, const int time) {
	phi.fastwrite(make_datafile_name(prefix + "phi", time));
	U.fastwrite(make_datafile_name(prefix + "U", time));
	V.fastwrite(make_datafile_name(prefix + "V", time));
	return;
}
void ParallelOutputFields_piEF(LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<SU2matrix>& F,
	LATfield2::Field<U1matrix>& E,
	const std::string& prefix, const int time) {
	pi.fastwrite(make_datafile_name(prefix + "pi", time));
	F.fastwrite(make_datafile_name(prefix + "F", time));
	E.fastwrite(make_datafile_name(prefix + "E", time));
	return;
}

void InitializeMomentumField(LATfield2::Site x,
	LATfield2::Field<SU2vector>& pi,
	LATfield2::Field<U1matrix>& E,
	LATfield2::Field<SU2matrix>& F,
	const int initTimeStep) {
	const int nowTime = initTimeStep % CYCLE;
	for (x.first(); x.test(); x.next()) {
		pi(x, nowTime) << 0, 0;
		E(x, 0, nowTime) = UNITY_E * U1matrix(1, 0);
		E(x, 1, nowTime) = UNITY_E * U1matrix(1, 0);
		E(x, 2, nowTime) = UNITY_E * U1matrix(1, 0);
		F(x, 0, nowTime) = UNITY_F * Ident;
		F(x, 1, nowTime) = UNITY_F * Ident;
		F(x, 2, nowTime) = UNITY_F * Ident;
	}
	pi.updateHalo();
	E.updateHalo();
	F.updateHalo();
	return;
}
