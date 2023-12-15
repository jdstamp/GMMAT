/*  ttemapit : An R Package for Survival linear Mixed Model Marginal Epistasis test
 *  Copyright (C) 2023  Julian Stamp
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*  Function fitglmm_ai was modified from GMMAT R Package
 *  by Han Chen.
 */

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <bzlib.h>
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
using namespace std;
using namespace arma;
using namespace Rcpp;


typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16;
#endif

// [[Rcpp::export]]
List fitglmm_ai(arma::vec Y, 
                  arma::mat X, 
                  size_t q, 
                  List Phi, 
                  size_t ng, 
                  List group, 
                  arma::vec W, 
                  arma::vec tau, 
                  arma::uvec fixtau)
{
	try {
		size_t n = X.n_rows, p = X.n_cols;
		const size_t q2 = sum(fixtau == 0);
		mat cov(p, p);
		vec alpha(p), eta(n);
		vec diagP = zeros<vec>(n);
		for(size_t i=1; i<=ng; ++i) {
		        uvec group_idx = as<uvec>(group[i-1]) - 1;
			diagP.elem( group_idx ) = tau[i-1] / W.elem( group_idx );
		}
		mat P = diagmat(diagP);
		for(size_t i=1; i<=q; ++i) {
			P = P + tau[i+ng-1] * as<mat>(Phi[i-1]);
		}
		double detSigma = det(P);
		mat Sigma_i = inv_sympd(P);
		mat Sigma_iX = Sigma_i * X;
		cov = inv_sympd(X.t() * Sigma_iX);
		P = Sigma_i - Sigma_iX * cov * Sigma_iX.t();
		vec PY = P * Y;
		alpha = cov * Sigma_iX.t() * Y;
		eta = Y - diagP % (Sigma_i * (Y - X * alpha));
		double logLik = -0.5 * (detSigma - dot(Y, PY));
 		if(q2 > 0) {
      const uvec idxtau = find(fixtau == 0);
      mat AI(q2, q2);
			vec score(q2), PAPY;
			vec APY = PY / W;
			diagP = diagvec(P) / W;
			for(size_t i=0; i<q2; ++i) {
			        if(idxtau[i] < ng) {
				        uvec group_idx = as<uvec>(group[idxtau[i]]) - 1;
					score[i] = dot(APY.elem( group_idx ), PY.elem( group_idx )) - sum(diagP.elem( group_idx ));
					for(size_t j=0; j<=i; ++j) {
					        uvec group_idx2 = as<uvec>(group[idxtau[j]]) - 1;
						AI(i,j) = dot(APY.elem( group_idx ), P.submat( group_idx, group_idx2) * APY.elem( group_idx2 ));
						if(j!=i) {AI(j,i) = AI(i,j);}
					}
				} else {
				        PAPY = P * as<mat>(Phi[idxtau[i]-ng]) * PY;
					score[i] = dot(Y, PAPY) - accu(P % as<mat>(Phi[idxtau[i]-ng]));
					for(size_t j=0; j<=i; ++j) {
					        if(idxtau[j] < ng) {
						        uvec group_idx = as<uvec>(group[idxtau[j]]) - 1;
						        AI(i,j) = dot(APY.elem( group_idx ), PAPY.elem( group_idx ));
							AI(j,i) = AI(i,j);
						} else {
						        AI(i,j) = dot(PY, as<mat>(Phi[idxtau[j]-ng]) * PAPY);
							if(j!=i) {AI(j,i) = AI(i,j);}
						}
					}
				}
			}
			vec Dtau = solve(AI, score);
			return List::create(Named("Dtau") = Dtau, Named("P") = P, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta, Named("logLik") = logLik);
		} else {
			return List::create(Named("Dtau") = R_NilValue, Named("P") = P, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta, Named("logLik") = logLik);
		}
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}


