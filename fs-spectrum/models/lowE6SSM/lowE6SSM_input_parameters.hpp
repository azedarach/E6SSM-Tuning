// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sun 24 Aug 2014 16:10:23

#ifndef lowE6SSM_INPUT_PARAMETERS_H
#define lowE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct lowE6SSM_input_parameters {
   double TanBeta;
   double QQp;
   double QLp;
   double QH1p;
   double QH2p;
   double Qdp;
   double Qup;
   double Qep;
   double QSp;
   double QDxp;
   double QDxbarp;
   double QHpp;
   double QHpbarp;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   double LambdaxInput;
   double MuPrInput;
   double gNInput;
   double vsInput;
   Eigen::Matrix<double,3,3> TYdInput;
   Eigen::Matrix<double,3,3> TYeInput;
   Eigen::Matrix<double,3,3> TKappaInput;
   Eigen::Matrix<double,2,2> TLambda12Input;
   double ALambdaxInput;
   Eigen::Matrix<double,3,3> AYuInput;
   double BMuPrInput;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> me2Input;
   Eigen::Matrix<double,2,2> mH1I2Input;
   Eigen::Matrix<double,2,2> mH2I2Input;
   Eigen::Matrix<double,2,2> msI2Input;
   Eigen::Matrix<double,3,3> mDx2Input;
   Eigen::Matrix<double,3,3> mDxbar2Input;
   double mHp2Input;
   double mHpbar2Input;
   double MassBInput;
   double MassWBInput;
   double MassGInput;
   double MassBpInput;

   lowE6SSM_input_parameters()
      : TanBeta(0), QQp(0), QLp(0), QH1p(0), QH2p(0), Qdp(0), Qup(0), Qep(0), QSp(
   0), QDxp(0), QDxbarp(0), QHpp(0), QHpbarp(0), KappaInput(Eigen::Matrix<
   double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>::Zero()),
   LambdaxInput(0), MuPrInput(0), gNInput(0), vsInput(0), TYdInput(
   Eigen::Matrix<double,3,3>::Zero()), TYeInput(Eigen::Matrix<double,3,3>::Zero
   ()), TKappaInput(Eigen::Matrix<double,3,3>::Zero()), TLambda12Input(
   Eigen::Matrix<double,2,2>::Zero()), TLambdaxInput(0), TYuInput(Eigen::Matrix
   <double,3,3>::Zero()), BMuPrInput(0), mq2Input(Eigen::Matrix<double,3,3>
   ::Zero()), ml2Input(Eigen::Matrix<double,3,3>::Zero()), md2Input(
   Eigen::Matrix<double,3,3>::Zero()), mu2Input(Eigen::Matrix<double,3,3>::Zero
   ()), me2Input(Eigen::Matrix<double,3,3>::Zero()), mH1I2Input(Eigen::Matrix<
   double,2,2>::Zero()), mH2I2Input(Eigen::Matrix<double,2,2>::Zero()),
   msI2Input(Eigen::Matrix<double,2,2>::Zero()), mDx2Input(Eigen::Matrix<double
   ,3,3>::Zero()), mDxbar2Input(Eigen::Matrix<double,3,3>::Zero()), mHp2Input(0
   ), mHpbar2Input(0), MassBInput(0), MassWBInput(0), MassGInput(0),
   MassBpInput(0)

   {}
};

} // namespace flexiblesusy

#endif
