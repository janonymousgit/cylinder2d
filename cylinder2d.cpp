/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2014 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* cylinder2d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR PorousParticleD2Q9Descriptor


// Parameters for the simulation setup
const int N = 1;        // resolution of the model
const int M = 1;        // time discretization refinement
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16.; // max. simulation time in s, SI unit
// const T maxPhysT = 32.; // max. simulation time in s, SI unit
const T L = 0.01/N;     // latticeL
const T mass = 1.;
const T g = 0.;
const T lengthX = 2.2;
const T lengthY = .41+L;
const T centerCylinderX = 0.2;
const T centerCylinderY = 0.2+L/2.;
// const T radiusCylinder = 0.05;
const T radiusCylinder = 0.10;
// T radiusCylinder;
T angle;

/// Stores geometry information in form of material numbers
void prepareGeometry(LBconverter<T> const& converter,
                     SuperGeometry2D<T>& superGeometry)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend(lengthX,lengthY);
  Vector<T,2> center(centerCylinderX,centerCylinderY);
  Vector<T,2> origin;
  IndicatorCircle2D<T> circle(center, radiusCylinder);

  superGeometry.rename(0,2);

  superGeometry.rename(2,1,1,1);

  /// Set material number for inflow
  extend[0] = 2.*L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  superGeometry.rename(2,3,1,inflow);
  /// Set material number for outflow
  origin[0] = lengthX-L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  superGeometry.rename(2,4,1,outflow);
  /// Set material number for cylinder
//  superGeometry.rename(1,5,circle);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice2D<T,DESCRIPTOR>& sLattice,
    LBconverter<T> const& converter,
    Dynamics<T, DESCRIPTOR>& bulkDynamics,
    Dynamics<T, DESCRIPTOR>& designDynamics,
    sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
    sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& offBc,
    SuperGeometry2D<T>& superGeometry)
{

OstreamManager clout(std::cout,"prepareLattice");
clout << "Prepare Lattice ..." << std::endl;

const T omega = converter.getOmega();

/// Material=0 -->do nothing
sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

/// Material=1 -->bulk dynamics
//  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);
sLattice.defineDynamics(superGeometry, 1, &designDynamics);

/// Material=2 -->bounce back
sLattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

/// Material=3 -->bulk dynamics (inflow)
sLattice.defineDynamics(superGeometry, 3, &bulkDynamics);

/// Material=4 -->bulk dynamics (outflow)
sLattice.defineDynamics(superGeometry, 4, &bulkDynamics);

/// Setting of the boundary conditions
sBoundaryCondition.addVelocityBoundary(superGeometry, 3, omega);
sBoundaryCondition.addPressureBoundary(superGeometry, 4, omega);

/// Initial conditions
AnalyticalConst2D<T,T> rhoF(1);
std::vector<T> velocity(2,T(0));
AnalyticalConst2D<T,T> uF(velocity);

// Initialize all values of distribution functions to their local equilibrium
sLattice.defineRhoU(superGeometry, 1, rhoF, uF);
sLattice.iniEquilibrium(superGeometry, 1, rhoF, uF);
sLattice.defineRhoU(superGeometry, 3, rhoF, uF);
sLattice.iniEquilibrium(superGeometry, 3, rhoF, uF);
sLattice.defineRhoU(superGeometry, 4, rhoF, uF);
sLattice.iniEquilibrium(superGeometry, 4, rhoF, uF);

// Make the lattice ready for simulation
sLattice.initialize();

clout << "Prepare Lattice ... OK" << std::endl;
}

class ParticleDynamics {
private:
  SuperLattice2D<T, DESCRIPTOR>& _sLattice;
  LBconverter<T> const& _converter;
  SuperGeometry2D<T>& _superGeometry;

  std::vector<T> _A;

  Vector<T,2> _pos;
  std::vector<T> _vel;
  std::vector<T> _acc;
  T _theta;
  T _omega;
  T _alpha;

public:
  ParticleDynamics(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                   LBconverter<T> const& converter,
                   SuperGeometry2D<T>& superGeometry)
    : _sLattice(sLattice), _converter(converter), _superGeometry(superGeometry),
      _A(3,T()), _pos(2,T()), _vel(2,T()), _acc(2,T()), _theta(0.), _omega(0.), _alpha(0.) {
  };

  void setPos(T xPos, T yPos) {
    _pos[0] = xPos;
    _pos[1] = yPos;
  };

  void setVel(T xVel, T yVel) {
    _vel[0] = xVel;
    _vel[1] = yVel;
  };

  void setTheta(T theta) {
    _theta = theta;
  };
  void setOmega(T omega) {
    _omega = omega;
  };

  Vector<T,2> getPos() {
    return _pos;
  };

  std::vector<T> getVel() {
    return _vel;
  };

  T getTheta() {
    return _theta;
  };

  T getOmega() {
    return _omega;
  };

  T getTorque() {
    return _A[2];
  }

  void computeBoundaryForce(SmoothIndicatorF2D<T,T>& indicator) {
    SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR> f(_sLattice, _superGeometry, indicator, _converter);
    SuperSumIndicator2D<T,DESCRIPTOR> sumF(f, _superGeometry, indicator);

    T F[5]= {0};
    int input[1];
    sumF(F, input);

    /// get particle acceleration through boundary force and gravity (and buoyancy)
    _A[0] = F[0]/mass - g;
    _A[1] = F[1]/mass;

    //circle
    _A[2] = 0.;//F[2] / (0.5*mass*pow(radius,2)); //T=J*a, J=1/2m*r^2 for circle
    //cuboid
//    _A[2] = F[2] / (1./12.*mass*(pow(2*radius,2)+pow(2*radius,2)));  //for cuboid
  };


  void addParticleColl(SmoothIndicatorF2D<T,T>& indicator, const Vector<T,2>& pos2, T delta) {
    T rad = indicator.getRadius();// + delta;//radius+delta;//
    T dist = std::sqrt(std::pow(_pos[0] - pos2[0], 2) + std::pow(_pos[1] - pos2[1], 2));
    T massInv = 1. / mass;
    T e1 = 1e-7;//delta;//
    T e = 1e-7;//delta * delta;//
    if (dist > 2.*rad + delta) {}
    else if (dist > 2.*rad) {
      _A[0] += massInv * (_pos[0] - pos2[0]) * std::pow((2.*rad + delta - dist), 2) /  e;
      _A[1] += massInv * (_pos[1] - pos2[1]) * std::pow((2.*rad + delta - dist), 2) /  e;
    } else {
      _A[0] += massInv * (_pos[0] - pos2[0]) * (2.*rad - dist) /  e1;
      _A[1] += massInv * (_pos[1] - pos2[1]) * (2.*rad - dist) /  e1;
    }
  }

  void addWallColl(SmoothIndicatorF2D<T,T>& indicator, T delta) {
    std::vector<T> dx(2,T());

    T w1 = delta/2.;//1e-7/2.;//
    T w = delta/2.;//1e-7/2.;//

    T rad = indicator.getRadius();// + delta;
    T massInv = 1. / mass;

    dx[0] = lengthX-_converter.getLatticeL() - _pos[0];
    dx[1] = lengthY-_converter.getLatticeL() - _pos[1];

    for (int i=0; i<2; i++) {
      if (dx[i] <= rad) {
        _A[i] += massInv*-dx[i]*(rad - dx[i])/w1;
      }
      if (_pos[i] <= rad) {
        _A[i] += massInv*_pos[i]*(rad - _pos[i])/w1;
      }
      if (dx[i] > rad && dx[i] <= rad+delta) {
        _A[i] += massInv*-dx[i]*std::pow((rad + delta - dx[i]),2)/w;
      }
      if (_pos[i] > rad && _pos[i] <= rad+delta) {
        _A[i] += massInv*_pos[i]*std::pow((rad + delta - _pos[i]),2)/w;
      }
    }
  };

  void addCollisionModel(SmoothIndicatorF2D<T,T>& indicator, const Vector<T,2>& pos2) {
    this->addParticleColl(indicator, pos2, _converter.getLatticeL());
    this->addWallColl(indicator, 1.*_converter.getLatticeL()); //10.*converter
  };

  void eulerIntegration() {
    T time = _converter.physTime();

    for (int i=0; i<2; i++) {
      _vel[i] += _A[i]*time;
      _pos[i] += _vel[i]*time;
    }
    _omega += _A[2]*time;   //angular velocity
    _theta += _omega*time;  //angle
  };

  void verletIntegration() {
    T time = _converter.physTime();

    _pos[0] += _vel[0]*time + (0.5*_acc[0]*time*time);
    _pos[1] += _vel[1]*time + (0.5*_acc[1]*time*time);

    T avgAx = (_acc[0] + _A[0])/2.;
    T avgAy = (_acc[1] + _A[1])/2.;
    _vel[0] += avgAx*time;
    _vel[1] += avgAy*time;

    _acc[0] = _A[0];
    _acc[1] = _A[1];

    _theta += _omega*time + (0.5*_alpha*time*time);
    T avgAlpha = (_alpha + _A[2])/2.;
    _omega += avgAlpha*time;
    _alpha = _A[2];

    return;
  };

  void updateParticleDynamics(std::string name) {
    if (name == "euler") {
      this->eulerIntegration();
    } else if (name == "verlet") {
      this->verletIntegration();
    } else {
      std::cout << "ERROR: no valid integration...use 'euler' or 'verlet'" << std::endl;
    }
  };

  void resetParticleField(SmoothIndicatorF2D<T, T>& indicator) {
    /// Analytical2D functor for particle motion (trans+rot)
    AnalyticalConst2D<T, T> velocity(0., 0.);

    /// defining porosity (additive inverse to one of indicator operator())
    AnalyticalConst2D<T, T> zero(0.);
    AnalyticalConst2D<T, T> one(1.);

    IndicatorCircle2D<T> circleI(_pos, indicator.getRadius() + .2 * _converter.getLatticeL());
    /// and dynamics field (translation and rotation)
    _sLattice.defineExternalField(_superGeometry, circleI,
                                  DESCRIPTOR<T>::ExternalField::velNumerator,
                                  DESCRIPTOR<T>::ExternalField::sizeOfVelNum,
                                  velocity);
    _sLattice.defineExternalField(_superGeometry, circleI,
                                  DESCRIPTOR<T>::ExternalField::velDenominator,
                                  DESCRIPTOR<T>::ExternalField::sizeOfVelDenom,
                                  zero);
    /// defining porous field
    _sLattice.defineExternalField(_superGeometry, circleI,
                                  DESCRIPTOR<T>::ExternalField::porosityIsAt,
                                  1, one);
  }

  void addParticleField(SmoothIndicatorF2D<T, T>& indicator) {
    /// converting in lattice velocities
    std::vector<T> velL(2, T());
    velL[0] = _converter.latticeVelocity(_vel[0]);
    velL[1] = _converter.latticeVelocity(_vel[1]);
    T omegaL = _converter.latticeVelocity(_omega);//_converter.physTime() * _omega; //

    /// Analytical2D functor for particle motion (trans+rot)
    ParticleU2D<T, T> velocity(indicator, velL, omegaL);

    /// defining porosity (additive inverse to one of indicator operator())
    AnalyticalConst2D<T, T> zero(0.);
    AnalyticalConst2D<T, T> one(1.);
    AnalyticalIdentity2D<T, T> tmp(one - indicator);
    IndicatorCircle2D<T> circleI(_pos, indicator.getRadius() + .2 * _converter.getLatticeL());
    /// and dynamics field (translation and rotation)
    _sLattice.addExternalField(_superGeometry, circleI,
                               DESCRIPTOR<T>::ExternalField::velNumerator,
                               DESCRIPTOR<T>::ExternalField::sizeOfVelNum,
                               velocity, indicator);
    _sLattice.addExternalField(_superGeometry, circleI,
                               DESCRIPTOR<T>::ExternalField::velDenominator,
                               DESCRIPTOR<T>::ExternalField::sizeOfVelDenom,
                               indicator);
    /// defining porous field
    _sLattice.multiplyExternalField(_superGeometry, circleI,
                                    DESCRIPTOR<T>::ExternalField::porosityIsAt,
                                    1, tmp);
  }
};

void setBoundaryValues(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                       LBconverter<T> const& converter, int iT,
                       SuperGeometry2D<T>& superGeometry,
                       ParticleDynamics& particle
		      )
{

  OstreamManager clout(std::cout,"setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst2D<T,T> one(1.);
    sLattice.defineExternalField(superGeometry, 1, DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, one);
  }

//  SmoothIndicatorCircle2D<T,T> circle(particle.getPos(), radiusCylinder, 1., 1.2*converter.getLatticeL());
  SmoothIndicatorTriangle2D<T,T> circle(particle.getPos(), radiusCylinder, 1., 1.2*converter.getLatticeL(), angle);


  particle.resetParticleField(circle);
//  particle.computeBoundaryForce(circle);
  //particle.updateParticleDynamics("verlet");
  particle.addParticleField(circle);

  int iTwrite = converter.numTimeSteps(.1);
  if (iT%iTwrite==0) {
		SuperLatticePhysDragIndicator2D<T,DESCRIPTOR> drag(sLattice, superGeometry, circle, converter);

    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag(_drag,input);
    clout << "drag=" << _drag[0] << "; lift=" << _drag[1] << endl;
  }

  // No of time steps for smooth start-up
  int iTmaxStart = converter.numTimeSteps(maxPhysT*0.4);
  int iTupdate = 5;

  if (iT%iTupdate==0 && iT<= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T(iT)};
    T frac[1] = {};
    StartScale(frac,iTvec);
    T maxVelocity = converter.getLatticeU()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU(superGeometry, 3, maxVelocity, distance2Wall);

    sLattice.defineU(superGeometry, 3, poiseuilleU);
  }
}

/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                LBconverter<T> const& converter, int iT,
                SuperGeometry2D<T>& superGeometry, Timer<T>& timer)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTKwriter2D<T> vtkWriter("cylinder2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity2D<T, DESCRIPTOR> externalPor(sLattice, converter);
  SuperLatticePhysExternalParticleVelocity2D<T, DESCRIPTOR> externalVel(sLattice, converter);
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( externalPor );
  vtkWriter.addFunctor( externalVel );

  const int vtkIter  = converter.numTimeSteps(.3);
  const int statIter = converter.numTimeSteps(.1);

  if (iT == 0) {
    /// Writes the converter log file
    writeLogFile(converter, "cylinder2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  if (iT%vtkIter==0) {
    vtkWriter.write(iT);
  }

  /// Writes the vtk files
  if (iT%vtkIter == 0) {
    vtkWriter.write(iT);

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
    BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction(normVel);
    BlockGifWriter<T> gifWriter;
    //gifWriter.write(planeReduction, 0, 0.7, iT, "vel"); //static scale
    gifWriter.write(planeReduction, iT, "vel"); // scaled
  }

  /// Writes output on the console
  if (iT%statIter == 0) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    sLattice.getStatistics().print(iT,converter.physTime(iT));

    /// Drag, lift, pressure drop
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
    AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolatePressure(pressure, true);
//    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag(sLattice, superGeometry, 5, converter);

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerCylinderX - radiusCylinder;
    point1[1] = centerCylinderY;

    point2[0] = centerCylinderX + radiusCylinder;
    point2[1] = centerCylinderY;

    T p1, p2;
    intpolatePressure(&p1,point1);
    intpolatePressure(&p2,point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

//    int input[3] = {};
//    T _drag[drag.getTargetDim()];
//    drag(_drag,input);
//    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << endl;
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  //clout << argv[0] << endl; -> cylinder2d
  //clout << argv[1] << endl; -> angle

  angle = atof(argv[1]);
  //radius = atof(argv[2]);
  clout << "current angle: " << angle << endl;

  LBconverter<T> converter(
    (int) 2,                               // dim
    (T)   L,                               // latticeL_
    (T)   0.02/M,                          // latticeU_
    (T)   0.2*2.*radiusCylinder/Re,        // charNu_
    (T)   2*radiusCylinder,                // charL_ = 1
    (T)   0.2                              // charU_ = 1
  );
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  Vector<T,2> extend(lengthX,lengthY);
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, L, noOfCuboids);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);

  BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>());
  PorousParticleBGKdynamics<T, DESCRIPTOR> designDynamics(converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>());

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition(sLattice);
  // createInterpBoundaryCondition2D<T,DESCRIPTOR>(sBoundaryCondition);
  createLocalBoundaryCondition2D<T,DESCRIPTOR>(sBoundaryCondition);

  sOffLatticeBoundaryCondition2D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
  createBouzidiBoundaryCondition2D<T, DESCRIPTOR> (sOffBoundaryCondition);

  prepareLattice(sLattice, converter, bulkDynamics, designDynamics, sBoundaryCondition, sOffBoundaryCondition, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer(converter.numTimeSteps(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  ParticleDynamics particle(sLattice, converter, superGeometry);
  particle.setPos(centerCylinderX, centerCylinderY);

  for (int iT = 0; iT < converter.numTimeSteps(maxPhysT); ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry, particle);

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer);

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}
