#ifndef EW_FUNCTIONS_H
#define EW_FUNCTIONS_H

#include <vector>
#include "EW_parameter.h"
#include <fstream>
#include <string>
#include <system_error>
#include <iomanip>

namespace Electroweak {
	using namespace ParaSite;

	//#define EM_F_SITE EM_F_SITE_TanmayEW//choose a way to calculate EM field. //Tanmay's
#define EM_F_SITE EM_F_SITE_TanmayEWFull
#define EM_A_SITE EM_A_SITE_EW
#define EM_E_SITE EM_E_SITE_EWFull
#define Z_F_SITE Z_F_SITE_VAC
#define EM_F_PLAQ EM_F_PLAQ_VAC
#define vecN vecN_VAC
#define phiSqNorm phiSqNorm_VAC

	//---class---

	class BoundarySite {
	public:
		ParaSite::Site<DIM> bs;
		std::vector<int> n;
		//Constructor
		BoundarySite();

		//functions
		int bdPos();
		int bdGaugeFieldPos(int direction);
		std::vector<Real> bdGaugeFieldNormalVector(int direcion); //Normalized
	};

	template<int D>
	const SU2vector DVDphi_ZeroT(const Site<D>& x, const Field<SU2vector, D>& phi, const int cycleTime) {
		return 2.0 * lambda*(phi(x, cycleTime).squaredNorm() - v2)*phi(x, cycleTime);
	}

	template<int D>
	const Real Vphi_ZeroT(const Site<D>& x, const Field<SU2vector, D>& phi, const int cycleTime) {
		return  lambda * pow((phi(x, cycleTime).squaredNorm() - v2), 2);
	}

	inline const Real VEV2_ZeroT() {
		return v2;
	}


#define DVDphi DVDphi_ZeroT<DIM>
#define Vphi Vphi_ZeroT<DIM>
#define VEV2 VEV2_ZeroT<DIM>



	//---auxilliary inline functions---
	template<int D>
	inline Real rx(const Site<D>& x, const int dir_j, const Real dx = DX, const Real* center_position = CENTER_POS) {
		return (x.coord(dir_j) - CENTER_POS[dir_j])*dx;
	}
	template<int D>
	inline Real rgx(const Site<D>& x, const int dir_j, const Real dx = DX, const Real* center_position = CENTER_POS) {
		return (x.coord(dir_j) - CENTER_POS[dir_j] + 0.5)*dx;
	}
	template<int D>
	Real radius(const Site<D>& x, const Real dx = DX, const Real* center = CENTER_POS) {
		auto rad = pow((x.coord(0) - center[0]), 2) + pow((x.coord(1) - center[1]), 2) + pow((x.coord(2) - center[2]), 2);
		return sqrt(rad)*dx;
	}
	template<int D>
	inline Real distance(const Site<D>& x, const Site<D>& y, const Real dx = DX) {
		return sqrt(pow((x.coord(0) - y.coord(0)), 2) + pow((x.coord(1) - y.coord(1)), 2) + pow((x.coord(2) - y.coord(2)), 2)) * dx;
	}
	template<int D>
	inline Real pbc_distance(const Site<D>& x, const Site<D>& y, const Real dx = DX) {
		auto x1_lat = std::abs(x.coord(0) - y.coord(0));
		auto x2_lat = std::abs(x.coord(1) - y.coord(1));
		auto x3_lat = std::abs(x.coord(2) - y.coord(2));
		auto x1_pbc = std::min(x1_lat, nSize[0] - x1_lat);
		auto x2_pbc = std::min(x2_lat, nSize[1] - x2_lat);
		auto x3_pbc = std::min(x3_lat, nSize[2] - x3_lat);
		return sqrt(x1_pbc*x1_pbc + x2_pbc * x2_pbc + x3_pbc * x3_pbc)*dx;
	}

	inline void Cart2Sph(const Real* cart, Real* sph) {
		sph[0] = sqrt(cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2]);
		if (std::abs(sph[0] - std::abs(cart[2])) < ERROR) { //on the z-axis
			sph[2] = 0;
			if (sph[0] < ERROR) { //at origin
				sph[1] = 0;
			}
			else {
				sph[1] = acos(cart[2] / sph[0]);
			}
		}
		else {
			sph[1] = acos(cart[2] / sph[0]);
			sph[2] = atan2(cart[1], cart[0]);
		}
	}
	inline Real arctan(const Real x) {
		Real y = atan(x);
		if (y >= 0) { return y; }
		else { return (y + PI); }
	}

	//Pauli matrices
	template<int D>
	SU2matrix PauliSph(const Site<D>& x, const int su2dir, const int linkDir = -1) {
		Real cart[D];
		Real sph[D];
		for (int i = 0; i < D; ++i) {
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

	//Plaquette
	template<int D>
	inline U1matrix Plaquette(const Site<D>& x, const Field<U1matrix, D>& v, const int i, const int j, const int cycleTime) {
		//return v(x, j, cycleTime) * v(x + j, i, cycleTime) * conj(v(x + i, j, cycleTime)) * conj(v(x, i, cycleTime));
		return v(x, j, cycleTime) * v(x + j, i, cycleTime) * conj(v(x, i, cycleTime)*v(x + i, j, cycleTime));
	}
	template<int D>
	inline SU2matrix Plaquette(const Site<D>& x, const Field<SU2matrix, D>& u, const int i, const int j, const int cycleTime) {
		//return u(x, j, cycleTime) * u(x + j, i, cycleTime) * u(x + i, j, cycleTime).adjoint() * u(x, i, cycleTime).adjoint();
		return u(x, j, cycleTime) * u(x + j, i, cycleTime) * (u(x, i, cycleTime)*u(x + i, j, cycleTime)).adjoint();
	}
	template<int D>
	inline const U1matrix Plaquette_SITE(const Site<D>& x, const Field<U1matrix, D>& v, const int i, const int j, const int T) {
		return v(x - i - j, j, T)*v(x - i, j, T)*v(x - i + j, i, T)*v(x + j, i, T)*conj(v(x - i - j, i, T)*v(x - j, i, T)*v(x + i - j, j, T)*v(x + i, j, T));
	}
	template<int D>
	inline const SU2matrix Plaquette_SITE(const Site<D>& x, const Field<SU2matrix, D>& v, const int i, const int j, const int T) {
		return v(x - i - j, j, T)*v(x - i, j, T)*v(x - i + j, i, T)*v(x + j, i, T)*(v(x - i - j, i, T)*v(x - j, i, T)*v(x + i - j, j, T)*v(x + i, j, T)).adjoint();
	}
	template<int D>
	inline U1matrix PlaquetteRect12(const Site<D>& x, const Field<U1matrix, D>& v, const int i, const int j, const int cycleTime) {
		return v(x, j, cycleTime)*v(x + j, j, cycleTime)*v(x + j + j, i, cycleTime)*conj(v(x + i + j, j, cycleTime))*conj(v(x + i, j, cycleTime))*conj(v(x, i, cycleTime));
	}
	template<int D>
	inline SU2matrix PlaquetteRect12(const Site<D>& x, const Field<SU2matrix, D>& u, const int i, const int j, const int cycleTime) {
		return u(x, j, cycleTime)*u(x + j, j, cycleTime)*u(x + j + j, i, cycleTime)*u(x + i + j, j, cycleTime).adjoint()*u(x + i, j, cycleTime).adjoint()*u(x, i, cycleTime).adjoint();
	}
	template<int D>
	inline U1matrix PlaquetteSum4(const Site<D>& x, const Field<U1matrix, D>& v, const int i, const int cycleTime) {
		//sum of 4 square plaquettes starting at (x,t) and orienting in U(x,i)
		U1matrix sum(0, 0);
		for (int j = 0; j < D; ++j) {
			if (j != i) {
				//sum += v(x, i, cycleTime)*v(x + i, j, cycleTime)*conj(v(x + j, i, cycleTime))*conj(v(x, j, cycleTime))
				//	+ v(x, i, cycleTime)*conj(v(x + i - j, j, cycleTime))*conj(v(x - j, i, cycleTime))*v(x - j, j, cycleTime);
				sum += v(x, i, cycleTime)*v(x + i, j, cycleTime)*conj(v(x, j, cycleTime)*v(x + j, i, cycleTime))
					+ v(x, i, cycleTime)*conj(v(x - j, i, cycleTime)*v(x + i - j, j, cycleTime))*v(x - j, j, cycleTime);
			}
		}
		return sum;
	}
	template<int D>
	inline SU2matrix PlaquetteSum4(const Site<D>& x, const Field<SU2matrix, D>& v, const int i, const int cycleTime) {
		SU2matrix sum = SU2matrix::Zero();
		for (int j = 0; j < D; ++j) {
			if (j != i) {
				//sum += v(x, i, cycleTime)*v(x + i, j, cycleTime)*v(x + j, i, cycleTime).adjoint()*v(x, j, cycleTime).adjoint()
				//	+ v(x, i, cycleTime)*v(x + i - j, j, cycleTime).adjoint()*v(x - j, i, cycleTime).adjoint()*v(x - j, j, cycleTime);
				sum += v(x, i, cycleTime)*v(x + i, j, cycleTime)*(v(x, j, cycleTime)*v(x + j, i, cycleTime)).adjoint()
					+ v(x, i, cycleTime)*(v(x - j, i, cycleTime)*v(x + i - j, j, cycleTime)).adjoint()*v(x - j, j, cycleTime);
			}
		}
		return sum;
	}
	template<int D>
	inline U1matrix PlaquetteRectSum12(const Site<D>& x, const Field<U1matrix, D>& v, const int i, const int cycleTime) {
		//sum of 12 rectangular(1,2) plaquettes starting at (x,t) and orienting in U(x,i)
		U1matrix sum(0, 0);
		for (int j = 0; j < D; ++j) {
			if (j != i) {
				sum += v(x, i, cycleTime)*v(x + i, i, cycleTime)*v(x + i + i, j, cycleTime)*conj(v(x + j + i, i, cycleTime))*conj(v(x + j, i, cycleTime))*conj(v(x, j, cycleTime))
					+ v(x, i, cycleTime)*v(x + i, j, cycleTime)*conj(v(x + j, i, cycleTime))*conj(v(x - i + j, i, cycleTime))*conj(v(x - i, j, cycleTime))*v(x - i, i, cycleTime)
					+ v(x, i, cycleTime)*v(x + i, i, cycleTime)*conj(v(x + i + i - j, j, cycleTime))*conj(v(x + i - j, i, cycleTime))*conj(v(x - j, i, cycleTime))*v(x - j, j, cycleTime)
					+ v(x, i, cycleTime)*conj(v(x + i - j, j, cycleTime))*conj(v(x - j, i, cycleTime))*conj(v(x - i - j, i, cycleTime))*v(x - i - j, j, cycleTime)*v(x - i, i, cycleTime)
					+ v(x, i, cycleTime)*v(x + i, j, cycleTime)*v(x + i + j, j, cycleTime)*conj(v(x + j + j, i, cycleTime))*conj(v(x + j, j, cycleTime))*conj(v(x, j, cycleTime))
					+ v(x, i, cycleTime)*conj(v(x + i - j, j, cycleTime))*conj(v(x + i - j - j, j, cycleTime))*conj(v(x - j - j, i, cycleTime))*v(x - j - j, j, cycleTime)*v(x - j, j, cycleTime);
				//check on Monday��checked.
			}
		}
		return sum;
	}
	template<int D>
	inline SU2matrix PlaquetteRectSum12(const Site<D>& x, const Field<SU2matrix, D>& v, const int i, const int cycleTime) {
		//sum of 12 rectangular(1,2) plaquettes starting at (x,t) and orienting in U(x,i)
		SU2matrix sum = SU2matrix::Zero();
		for (int j = 0; j < D; ++j) {
			if (j != i) {
				sum += v(x, i, cycleTime)*v(x + i, i, cycleTime)*v(x + i + i, j, cycleTime)*v(x + j + i, i, cycleTime).adjoint()*v(x + j, i, cycleTime).adjoint()*v(x, j, cycleTime).adjoint()
					+ v(x, i, cycleTime)*v(x + i, j, cycleTime)*v(x + j, i, cycleTime).adjoint()*v(x - i + j, i, cycleTime).adjoint()*v(x - i, j, cycleTime).adjoint()*v(x - i, i, cycleTime)
					+ v(x, i, cycleTime)*v(x + i, i, cycleTime)*v(x + i + i - j, j, cycleTime).adjoint()*v(x + i - j, i, cycleTime).adjoint()*v(x - j, i, cycleTime).adjoint()*v(x - j, j, cycleTime)
					+ v(x, i, cycleTime)*v(x + i - j, j, cycleTime).adjoint()*v(x - j, i, cycleTime).adjoint()*v(x - i - j, i, cycleTime).adjoint()*v(x - i - j, j, cycleTime)*v(x - i, i, cycleTime)
					+ v(x, i, cycleTime)*v(x + i, j, cycleTime)*v(x + i + j, j, cycleTime)*v(x + j + j, i, cycleTime).adjoint()*v(x + j, j, cycleTime).adjoint()*v(x, j, cycleTime).adjoint()
					+ v(x, i, cycleTime)*v(x + i - j, j, cycleTime).adjoint()*v(x + i - j - j, j, cycleTime).adjoint()*v(x - j - j, i, cycleTime).adjoint()*v(x - j - j, j, cycleTime)*v(x - j, j, cycleTime);
				//check on Monday��checked.
			}
		}
		return sum;
	}

	//derivatives
	template<int D>
	inline const SU2vector CovDphi_Sym(const Site<D>& x, 
		const int cycleTime,
		const int component,
		const Field<SU2vector, D>& phi, 
		const Field<SU2matrix, D>& U, 
		const Field<U1matrix, D>& V) {
		return 0.5 / DX * (U(x, component, cycleTime)*V(x, component, cycleTime)*phi(x + component, cycleTime) - U(x - component, component, cycleTime).adjoint()*conj(V(x - component, component, cycleTime))*phi(x - component, cycleTime));
	}
	
	template<int D>
	inline const SU2vector CovDphi(const Site<D>& x, 
		const int cycleTime,
		const int component,
		const Field<SU2vector, D>& phi, 
		const Field<SU2matrix, D>& U, 
		const Field<U1matrix, D>& V) {
		return (U(x, component, cycleTime)*V(x, component, cycleTime)*phi(x + component, cycleTime) - phi(x, cycleTime)) / DX;
	}
	
	template<int D>
	inline const SU2vector CovDphi_Plaq(const Site<D>& x, 
		const int T,
		const int dir_dx, const int dir_plaq,
		const Field<SU2vector, D>& phi, 
		const Field<SU2matrix, D>& U, 
		const Field<U1matrix, D>& V) {
		return 0.5*(CovDphi(x, T, dir_dx, phi, U, V) + CovDphi(x + dir_plaq, T, dir_dx, phi, U, V));
	}
	
	template<int D>
	inline const SU2vector dphidx_Site(const Site<D>& x, const int T,
		const int ddir, const Field<SU2vector, D>& phi) {
		return (phi(x + ddir, T) - phi(x - ddir, T)) / (2.0*DX);
	}
	
	template<int D>
	const Real dWdx_Link(const Site<D>& x, const int T,
		const int ddir, const int vdir, const int su2a,
		const Field<SU2matrix, D>& U) {
		if (ddir == vdir) {
			if (x.coord(vdir) == 0) {
				return (SU2_W(x + ddir, U, vdir, su2a, T) - SU2_W(x, U, vdir, su2a, T)) / DX;
			}
			else if (x.coord(vdir) == nSize[vdir] - 2) {
				return (SU2_W(x, U, vdir, su2a, T) - SU2_W(x - ddir, U, vdir, su2a, T)) / DX;
			}
			else {
				return (SU2_W(x + ddir, U, vdir, su2a, T) - SU2_W(x - ddir, U, vdir, su2a, T)) / (2.0*DX);
			}
		}
		else {
			return (SU2_W(x + ddir, U, vdir, su2a, T) - SU2_W(x - ddir, U, vdir, su2a, T)) / (2.0*DX);
		}
	}
	
	template<int D>
	const Real dYdx_Link(const Site<D>& x, const int T,
		const int ddir, const int vdir,
		const Field<U1matrix, D>& V) {
		if (ddir == vdir) {
			if (x.coord(vdir) == 0) {
				return (U1_Y(x + ddir, V, vdir, T) - U1_Y(x, V, vdir, T)) / DX;
			}
			else if (x.coord(vdir) == nSize[vdir] - 2) {
				return (U1_Y(x, V, vdir, T) - U1_Y(x - ddir, V, vdir, T)) / DX;
			}
			else {
				return (U1_Y(x + ddir, V, vdir, T) - U1_Y(x - ddir, V, vdir, T)) / (2.0*DX);
			}
		}
		else {
			return (U1_Y(x + ddir, V, vdir, T) - U1_Y(x - ddir, V, vdir, T)) / (2.0*DX);
		}
	}

	//Higgs fields
	template<int D>
	inline SU2vector unitPhi(const Site<D>& x, const Field<SU2vector, D>& phi, const int cycleTime) {
		if (phi(x, cycleTime).norm() > ERROR) return phi(x, cycleTime) / phi(x, cycleTime).norm();
		else return (SU2vector() << 0, 0).finished();
	}

	inline SU2vector unitPhi(const SU2vector& phi) {
		if (phi.norm() > ERROR) return phi / phi.norm();
		else return (SU2vector() << 0, 0).finished();
	}

	template<int D>
	inline SU2vector linkPhi(const Site<D>& x, const int dir, Field<SU2vector, D>& phi, const int cycleTime) {
		return 0.5*(phi(x, cycleTime) + phi(x + dir, cycleTime));
	}

	template<int D>
	inline SU2vector plaqPhi(const Site<D>& x, const int dir_i, const int dir_j, Field<SU2vector, D>& phi, const int T) {
		return 0.25*(phi(x, T) + phi(x + dir_i, T) + phi(x + dir_j, T) + phi(x + dir_i + dir_j, T));
	}

	inline Real vecN_UNITY(const SU2vector& phi, const int su2dir) {
		return -phi.dot(Pauli[su2dir] * phi).real() / phi.squaredNorm();
	}

	inline Real vecN_VAC(const SU2vector& phi, const int su2dir) {
		return -phi.dot(Pauli[su2dir] * phi).real() / v2;
	}

	inline Real phiSqNorm_UNITY(const SU2vector& phi) {
		return phi.squaredNorm();
	}
	inline Real phiSqNorm_VAC(const SU2vector& phi) {
		return v2;
	}

	template<int D>
	inline SU2matrix SU2phi(const Site<D>& x, 
		const int cycleTime,
		const Field<SU2vector, D>& phi) {
		return 1.0 / phi(x, cycleTime).norm()*(SU2matrix() << conj(phi(x, cycleTime)(1)), phi(x, cycleTime)(0), -conj(phi(x, cycleTime)(0)), phi(x, cycleTime)(1)).finished();
	}

	template<int D>
	inline SU2matrix SU2invphi(const Site<D>& x, 
		const int cycleTime,
		const Field<SU2vector, D>& phi) {
		return 1.0 / phi(x, cycleTime).norm()*(SU2matrix() << phi(x, cycleTime)(1), -phi(x, cycleTime)(0), conj(phi(x, cycleTime)(0)), conj(phi(x, cycleTime)(1))).finished();
	}

	//Link fields to Gauge fields or field strength
	template<int D>
	inline const Real SU2_E(const Site<D>& x, const Field<SU2matrix, D>& F, const int dir_i, const int su2_a, const int cycleTime) {
		return (-iPauli[su2_a] * F(x, dir_i, cycleTime)).trace().real();
	}
	template<int D>
	inline const Real SU2_B(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int dir_j, const int su2_a, const int cycleTime) {
		return UNITY_FB * (-iPauli[su2_a] * Plaquette(x, U, dir_i, dir_j, cycleTime)).trace().real();
	}
	template<int D>
	inline const Real SU2_W(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int su2_a, const int cycleTime) {
		return 1.0 / (g*DX)*(iPauli[su2_a] * U(x, dir_i, cycleTime)).trace().real();
	}
	inline const Real SU2_W(const SU2matrix& U, const int su2a) {
		return 1.0 / (g*DX)*(iPauli[su2a] * U).trace().real();
	}
	template<int D>
	inline Real U1_E(const Site<D>& x, const Field<U1matrix, D>& E, const int dir_i, const int cycleTime) {
		return E(x, dir_i, cycleTime).imag();
	}
	template<int D>
	inline Real U1_B(const Site<D>& x, const Field<U1matrix, D>& V, const int dir_i, const int dir_j, const int cycleTime) {
		return UNITY_EB * Plaquette(x, V, dir_i, dir_j, cycleTime).imag();
	}
	template<int D>
	inline const Real U1_Y(const Site<D>& x, const Field<U1matrix, D>& V, const int dir_i, const int cycleTime) {
		return -2.0 / (gp*DX)*V(x, dir_i, cycleTime).imag();
	}
	inline const Real U1_Y(const U1matrix& V) {
		return -2.0 / (gp*DX)*V.imag();
	}
	template<int D>
	inline Real SU2_Ucomp(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int su2_a, const int cycleTime) {
		return (iPauli[su2_a] * U(x, dir_i, cycleTime)).trace().real();
	}
	template<int D>
	inline Real SU2_E_SITE(const Site<D>& x, const Field<SU2matrix, D>& F, const int dir_i, const int su2_a, const int cycleTime) {
		return 0.5*(-iPauli[su2_a] * (F(x, dir_i, cycleTime) + F(x - dir_i, dir_i, cycleTime))).trace().real();
	}
	//backup - IMPROVEMENT-07032017
	//inline Real SU2_B_SITE(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int dir_j, const int su2_a, const int cycleTime) {
	//	return 0.25*UNITY_FB*(-iPauli[su2_a] * (Plaquette(x, U, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_i, U, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_j, U, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_i-dir_j, U, dir_i, dir_j, cycleTime))).trace().real();
	//}
	template<int D>
	inline Real SU2_B_SITE(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int dir_j, const int su2_a, const int cycleTime) {
		return 0.25*UNITY_FB*(-iPauli[su2_a] * Plaquette_SITE(x, U, dir_i, dir_j, cycleTime)).trace().real();
	}
	template<int D>
	inline Real SU2_W_SITE(const Site<D>& x, const Field<SU2matrix, D>& U, const int dir_i, const int su2_a, const int cycleTime) {
		return 0.5 / (g*DX)*(iPauli[su2_a] * (U(x, dir_i, cycleTime) + U(x - dir_i, dir_i, cycleTime))).trace().real();
	}
	template<int D>
	inline Real U1_E_SITE(const Site<D>& x, const Field<U1matrix, D>& E, const int dir_i, const int cycleTime) {
		return 0.5*(E(x - dir_i, dir_i, cycleTime) + E(x, dir_i, cycleTime)).imag();
	}
	//backup - IMPROVEMENT-07032017
	//inline Real U1_B_SITE(const Site<D>& x, const Field<U1matrix, D>& V, const int dir_i, const int dir_j, const int cycleTime) {
	//	return 0.25*UNITY_EB*(Plaquette(x, V, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_i, V, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_j, V, dir_i, dir_j, cycleTime)+ Plaquette(x-dir_i-dir_j, V, dir_i, dir_j, cycleTime)).imag();
	//}
	template<int D>
	inline Real U1_B_SITE(const Site<D>& x, const Field<U1matrix, D>& V, const int dir_i, const int dir_j, const int cycleTime) {
		return 0.25*UNITY_EB*Plaquette_SITE(x, V, dir_i, dir_j, cycleTime).imag();
	}
	template<int D>
	inline Real U1_Y_SITE(const Site<D>& x, const Field<U1matrix, D>& V, const int dir_i, const int cycleTime) {
		return -1.0 / (gp*DX)*(V(x, dir_i, cycleTime) + V(x - dir_i, dir_i, cycleTime)).imag();
	}

	//EM and Z field
	template<int D>
	Real EM_A_SITE_EW(const Site<D>& x, const int cycleTime, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
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
	//Real EM_A_SITE_EW_V2(const Site<D>& x, const int cycleTime, const int dir_i,
	//	LATfield2::Field<SU2vector>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V);

	template<int D>
	Real EM_A_LINK_SU2INV(const Site<D>& x, const int T, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		Cmplx Apotential(0, 0);
		Apotential = I * sinw / g * (CovDphi(x, T, dir_i, phi, U, V).dot(phi(x, T)) - phi(x, T).dot(CovDphi(x, T, dir_i, phi, U, V)));
		Apotential += U1_Y(x, V, dir_i, T) / cosw;
		return Apotential.real();
	}

	template<int D>
	Real EM_A_SITE_SU2INV(const Site<D>& x, const int T, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		return 0.5*(EM_A_LINK_SU2INV(x, T, dir_i, phi, U, V) + EM_A_LINK_SU2INV(x - dir_i, T, dir_i, phi, U, V));
	}

	template<int D>
	Real EM_F_SITE_TanmayEW(const Site<D>& x,
		const int cycleTime, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const int isScalarCurrentIncluded = 1) {
		Cmplx Fstrength(0, 0);
		for (int a = 0; a<SU2DIM; ++a) {
			Fstrength += vecN(phi(x, cycleTime), a) * SU2_B_SITE(x, U, dir_i, dir_j, a, cycleTime);
		}
		Fstrength *= sinw;
		Fstrength += cosw * U1_B_SITE(x, V, dir_i, dir_j, cycleTime);
		if (isScalarCurrentIncluded) {
			Fstrength += (-2.0) * I / g * sinw * (CovDphi_Sym(x, cycleTime, dir_i, phi, U, V).dot(CovDphi_Sym(x, cycleTime, dir_j, phi, U, V)) - CovDphi_Sym(x, cycleTime, dir_j, phi, U, V).dot(CovDphi_Sym(x, cycleTime, dir_i, phi, U, V))) / phiSqNorm(phi(x, cycleTime));
		}
#ifdef TEST_MODE
		if (abs(Fstrength.imag()) > ERROR) cerr << "EM FIELD STRENGTH ERROR: NON-REAL VALUE." << endl;
#endif
		return Fstrength.real();
	}
	
	template<int D>
	inline Real EM_F_SITE_VAC(const Site<D>& x,
		const int cycleTime, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		//forcefully ignore the Higgs term.
		return EM_F_SITE_TanmayEW(x, cycleTime, dir_i, dir_j, phi, U, V, 0);
	}
	template<int D>
	inline Real EM_F_SITE_TanmayEWFull(const Site<D>& x,
		const int cycleTime, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		//forcefully include the Higgs term.
		return EM_F_SITE_TanmayEW(x, cycleTime, dir_i, dir_j, phi, U, V, 1);
	}
	template<int D>
	inline Real EM_F_SITE_pureScalar(const Site<D>& x,
		const int cycleTime, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		const auto f = (-2.0) * I / g * sinw * (CovDphi_Sym(x, cycleTime, dir_i, phi, U, V).dot(CovDphi_Sym(x, cycleTime, dir_j, phi, U, V)) - CovDphi_Sym(x, cycleTime, dir_j, phi, U, V).dot(CovDphi_Sym(x, cycleTime, dir_i, phi, U, V))) / phiSqNorm(phi(x, cycleTime));
		return f.real();
	}
	template<int D>
	Real EM_F_PLAQ_TanmayEW(const Site<D>& x,
		const int T, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const int isScalarCurrentIncluded = 1) {
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
	template<int D>
	inline Real EM_F_PLAQ_VAC(const Site<D>& x,
		const int T, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		//forcefully ignore the Higgs term.
		return EM_F_PLAQ_TanmayEW(x, T, dir_i, dir_j, phi, U, V, 0);
	}
	template<int D>
	inline Real EM_F_PLAQ_TanmayFull(const Site<D>& x,
		const int T, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		//forcefully ignore the Higgs term.
		return EM_F_PLAQ_TanmayEW(x, T, dir_i, dir_j, phi, U, V, 1);
	}

	template<int D>
	Real EM_F_SITE_SU2INV(const Site<D>& x,
		const int T, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		return ((EM_A_LINK_SU2INV(x, T, dir_j, phi, U, V) - EM_A_LINK_SU2INV(x - dir_i, T, dir_j, phi, U, V)) - (EM_A_LINK_SU2INV(x, T, dir_i, phi, U, V) - EM_A_LINK_SU2INV(x - dir_j, T, dir_i, phi, U, V))) / DX;
	}
	
	template<int D>
	Real Z_F_SITE_VAC(const Site<D>& x, const int cycleTime, const int dir_i, const int dir_j,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V) {
		Cmplx Fstrength(0, 0);
		for (int a = 0; a<SU2DIM; ++a) {
			Fstrength += vecN(phi(x, cycleTime), a) * SU2_B_SITE(x, U, dir_i, dir_j, a, cycleTime);
		}
		Fstrength *= cosw;
		Fstrength -= sinw * U1_B_SITE(x, V, dir_i, dir_j, cycleTime);
		return Fstrength.real();
	}

	template<int D>
	Real EM_E_SITE_EW(const Site<D>& x, const int cycleTime, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const Field<SU2vector, D>& pi, const Field<SU2matrix, D>& F, const Field<U1matrix, D>& E,
		const int isScalarCurrentIncluded = 1) {
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

	template<int D>
	inline Real EM_E_SITE_EWFull(const Site<D>& x, const int cycleTime, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const Field<SU2vector, D>& pi, const Field<SU2matrix, D>& F, const Field<U1matrix, D>& E) {
		return EM_E_SITE_EW(x, cycleTime, dir_i, phi, U, V, pi, F, E, 1);
	}

	template<int D>
	inline Real EM_E_SITE_noScalar(const Site<D>& x, const int cycleTime, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const Field<SU2vector, D>& pi, const Field<SU2matrix, D>& F, const Field<U1matrix, D>& E) {
		return EM_E_SITE_EW(x, cycleTime, dir_i, phi, U, V, pi, F, E, 0);
	}

	template<int D>
	inline Real EM_E_SITE_pureScalar(const Site<D>& x, const int cycleTime, const int dir_i,
		const Field<SU2vector, D>& phi, const Field<SU2matrix, D>& U, const Field<U1matrix, D>& V,
		const Field<SU2vector, D>& pi, const Field<SU2matrix, D>& F, const Field<U1matrix, D>& E) {
		const auto efield = (-I)*2.0 * sinw / g * (pi(x, cycleTime).dot(CovDphi_Sym(x, cycleTime, dir_i, phi, U, V)) - CovDphi_Sym(x, cycleTime, dir_i, phi, U, V).dot(pi(x, cycleTime))) / phiSqNorm(phi(x, cycleTime));
		return efield.real();
	}

	// Gauge fields to links fields
	template<int D>
	int u1field2V(const Site<D>& x, const int dir, const Field<U1matrix, D>& V, const int cycleTime, Real imagV, const int isLargeField = 0) {
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			imagV = sin(imagV);
		}
		Real reFutureV;
		if (pow(imagV, 2) <= 1.0) reFutureV = sqrt(1.0 - pow(imagV, 2));
		else {
			std::cerr << "U1 FIELD ERROR: NON-UNITARY." << std::endl;
			return 0;
		}
		V(x, dir, cycleTime) = U1matrix(reFutureV, imagV);
		return 1;
	}
	
	template<int D>
	int u1field2E(const Site<D>& x, const int dir, const Field<U1matrix, D>& E, const int cycleTime, Real imagE, const int isLargeField = 0) {
		Real imagV = imagE / UNITY_E;
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			imagV = sin(imagV);
		}
		Real reFutureV;
		if (pow(imagV, 2) <= 1.0) reFutureV = sqrt(1.0 - pow(imagV, 2));
		else {
			std::cerr << "U1 FIELD ERROR: NON-UNITARY." << std::endl;
			return 0;
		}
		E(x, dir, cycleTime) = UNITY_E * U1matrix(reFutureV, imagV);
		return 1;
	}

	inline const U1matrix u1field2E(Real imagE, const int isLargeField = 0) {
		Real imagV = imagE / UNITY_E;
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			imagV = sin(imagV);
		}
		Real reFutureV;
		if (pow(imagV, 2) <= 1.0) reFutureV = sqrt(1.0 - pow(imagV, 2));
		else {
			std::cerr << "U1 FIELD ERROR: NON-UNITARY." << std::endl;
			return 0;
		}
		return UNITY_E * U1matrix(reFutureV, imagV);
	}

	template<int D>
	inline int su2field2U(const Site<D>& x, const int dir, const Field<SU2matrix, D>& U, const int cycleTime, Cmplx* su2_Ucomp, const int isLargeField = 0) {
		Real su2_Usum = 0;
		Real su2_U0;
		for (int a = 0; a < SU2DIM; ++a) {
			if (su2_Ucomp[a].imag() > EVO_ERROR) {
				std::cerr << "SU2 FIELD ERROR:COMPLEX VALUE" << std::endl;
				return 0;
			}
			su2_Usum += pow(su2_Ucomp[a].real() / 2.0, 2);
		}
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			if (su2_Usum > ERROR) {
				for (int a = 0; a < SU2DIM; ++a) {
					su2_Ucomp[a] = su2_Ucomp[a] / sqrt(su2_Usum)*std::abs(sin(sqrt(su2_Usum)));
				}
				su2_Usum = pow(sin(sqrt(su2_Usum)), 2);
			}
		}
		if (su2_Usum <= 1.0) su2_U0 = sqrt(1.0 - su2_Usum);
		else {
			std::cerr << "SU2 FIELD ERROR: NON-UNITARY." << std::endl;
			return 0;
		}
		U(x, dir, cycleTime) = su2_U0 * Ident;
		for (int a = 0; a < SU2DIM; ++a) {
			U(x, dir, cycleTime) += -I / 2.0 * Pauli[a] * su2_Ucomp[a].real();
		}
		return 1;
	}

	template<int D>
	inline int su2field2F(const Site<D>& x, const int dir, const Field<SU2matrix, D>& F, const int cycleTime, Cmplx* su2_Fcomp, const int isLargeField = 0) {
		Cmplx su2_Ucomp[SU2DIM] = { su2_Fcomp[0] / UNITY_F, su2_Fcomp[1] / UNITY_F, su2_Fcomp[2] / UNITY_F };
		Real su2_Usum = 0;
		Real su2_U0;
		for (int a = 0; a < SU2DIM; ++a) {
			if (su2_Ucomp[a].imag() > EVO_ERROR) {
				std::cerr << "SU2 FIELD ERROR:COMPLEX VALUE" << std::endl;
				return 0;
			}
			su2_Usum += pow(su2_Ucomp[a].real() / 2.0, 2);
		}
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			if (su2_Usum > ERROR) {
				for (int a = 0; a < SU2DIM; ++a) {
					su2_Ucomp[a] = su2_Ucomp[a] / std::sqrt(su2_Usum)*std::abs(sin(std::sqrt(su2_Usum)));
				}
				su2_Usum = pow(sin(sqrt(su2_Usum)), 2);
			}
		}
		if (su2_Usum <= 1.0) su2_U0 = sqrt(1.0 - su2_Usum);
		else {
			std::cerr << "SU2 FIELD ERROR: NON-UNITARY: su2_Usum = " << su2_Usum << std::endl;
			std::cerr << x.coord(0) << "," << x.coord(1) << "," << x.coord(2) << std::endl;
			std::cerr << su2_Ucomp[0] << "," << su2_Ucomp[1] << "," << su2_Ucomp[2] << std::endl;
			std::cerr << su2_Fcomp[0] << "," << su2_Fcomp[1] << "," << su2_Fcomp[2] << std::endl;
			exit(EXIT_FAILURE);
		}
		F(x, dir, cycleTime) = su2_U0 * Ident;
		for (int a = 0; a < SU2DIM; ++a) {
			F(x, dir, cycleTime) += -I / 2.0 * Pauli[a] * su2_Ucomp[a].real();
		}
		F(x, dir, cycleTime) *= UNITY_F;
		return 1;
	}

	inline const SU2matrix su2field2F(Cmplx* su2_Fcomp, const int isLargeField = 0) {
		Cmplx su2_Ucomp[SU2DIM] = { su2_Fcomp[0] / UNITY_F, su2_Fcomp[1] / UNITY_F, su2_Fcomp[2] / UNITY_F };
		Real su2_Usum = 0;
		Real su2_U0;
		for (int a = 0; a < SU2DIM; ++a) {
			if (su2_Ucomp[a].imag() > EVO_ERROR) {
				std::cerr << "SU2 FIELD ERROR:COMPLEX VALUE" << std::endl;
				exit(EXIT_FAILURE);
			}
			su2_Usum += pow(su2_Ucomp[a].real() / 2.0, 2);
		}
		if (isLargeField) {
			//A dangerous move to remove non-unitary error. (beyond lattice equations)
			if (su2_Usum > ERROR) {
				for (int a = 0; a < SU2DIM; ++a) {
					su2_Ucomp[a] = su2_Ucomp[a] / sqrt(su2_Usum)*std::abs(sin(sqrt(su2_Usum)));
				}
				su2_Usum = pow(sin(sqrt(su2_Usum)), 2);
			}
		}
		if (su2_Usum <= 1.0) su2_U0 = sqrt(1.0 - su2_Usum);
		else {
			std::cerr << "SU2 FIELD ERROR: NON-UNITARY: su2_Usum = " << su2_Usum << std::endl;
			exit(EXIT_FAILURE);
		}
		SU2matrix su2F = su2_U0 * Ident;
		for (int a = 0; a < SU2DIM; ++a) {
			su2F += -I / 2.0 * Pauli[a] * su2_Ucomp[a].real();
		}
		su2F *= UNITY_F;
		return su2F;
	}

	template<int _ = 0>
	SU2matrix su2_W2U(const SU2matrix& gw, const double dx = DX) {
		//This is a new function
		/*gw=sigma^a/(2i)*W^a*/
		assert( std::abs(gw.trace()) < 1e-6 ); // gw should be traceless.
		//compute gW^a*dx (complex)
		Cmplx cwa[SU2DIM] = {(iPauli[0]*gw*dx).trace(), (iPauli[1]*gw*dx).trace(), (iPauli[2]*gw*dx).trace()};
		assert( std::abs(cwa[0].imag()) < 1e-6 );
		assert( std::abs(cwa[1].imag()) < 1e-6 );
		assert( std::abs(cwa[2].imag()) < 1e-6 );
		Real su2_Usum = 0;
		for(auto a = 0; a < SU2DIM; ++a){
			su2_Usum += pow(cwa[a].real(), 2);
		}
		su2_Usum *= 0.25;
		
		assert(su2_Usum <= 1.0);
		Real su2_U0 = sqrt(1.0 - su2_Usum);
		SU2matrix U = su2_U0*Ident;
		for (int a = 0; a < SU2DIM; ++a) {
			U -= 0.5 * iPauli[a] * cwa[a].real();
		}
		return U;
	}

}

#endif
