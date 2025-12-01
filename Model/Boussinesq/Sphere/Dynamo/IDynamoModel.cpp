/**
 * @file IDynamoModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoModel.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Induction.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Momentum.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/Transport.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/gitHash.hpp"
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolNSpectrumWriter.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

VectorFormulation::Id IDynamoModel::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IDynamoModel::version() const
{
   return "BoussinesqSphereDynamo:" + std::string(gitHash);
}

void IDynamoModel::addEquations(SharedSimulation spSim)
{
   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Transport>(
      this->spBackend());

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Momentum>(
      this->spBackend());

   // Add induction equation
   spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Induction>(
      this->spBackend());

   #ifdef QUICC_USE_MLIR_GRAPH
   // Add Graph

   // At this point we only support the explicit model
   // If the jitter is not set up, we default to the old implementation
   if (spSim->ss().has(SpatialScheme::Feature::SpectralOrdering123))
   {
      return;
   }

   std::string graphStr = R"mlir(
// type aliases
!real = tensor<?x?x?xf64>
!complex = tensor<?x?x?xcomplex<f64>>

func.func private @bwdScalar(%S: !complex) -> !real {
    // backward scalar path
    %S1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "P"}
    %S1T = quiccir.transpose %S1 permutation = [1, 2, 0] : !complex -> !complex
    %S2 = quiccir.al.prj %S1T : !complex -> !complex attributes{kind = "P"}
    %S2T = quiccir.transpose %S2 permutation = [1, 2, 0] : !complex -> !complex
    %S3 = quiccir.fr.prj %S2T : !complex -> !real attributes{kind = "P"}
    return %S3 : !real
}

func.func private @bwdGradScalar(%S: !complex) -> (!real, !real, !real) {
    // grad R
    %SdR1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "D1"}
    %SdR1T = quiccir.transpose %SdR1 permutation = [1, 2, 0] : !complex -> !complex
    %SdR2 = quiccir.al.prj %SdR1T : !complex -> !complex attributes{kind = "P"}
    %SdR2T = quiccir.transpose %SdR2 permutation = [1, 2, 0] : !complex -> !complex
    %SdR = quiccir.fr.prj %SdR2T : !complex -> !real attributes{kind = "P"}
    // grad Theta
    %SdTh1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %SdTh1T = quiccir.transpose %SdTh1 permutation = [1, 2, 0] : !complex -> !complex
    %SdTh2 = quiccir.al.prj %SdTh1T : !complex -> !complex attributes{kind = "D1"}
    %SdTh2T = quiccir.transpose %SdTh2 permutation = [1, 2, 0] : !complex -> !complex
    %SdTh = quiccir.fr.prj %SdTh2T : !complex -> !real attributes{kind = "P"}
    // grad Phi
    %SdPh1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %SdPh1T = quiccir.transpose %SdPh1 permutation = [1, 2, 0] : !complex -> !complex
    %SdPh2 = quiccir.al.prj %SdPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %SdPh2T = quiccir.transpose %SdPh2 permutation = [1, 2, 0] : !complex -> !complex
    %SdPh = quiccir.fr.prj %SdPh2T : !complex -> !real attributes{kind = "P"}
    return %SdR, %SdTh, %SdPh : !real, !real, !real
}

func.func private @bwdVector(%Tor: !complex, %Pol: !complex) -> (!real, !real, !real) {
    // R
    %PolR1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %PolR1T = quiccir.transpose %PolR1 permutation = [1, 2, 0] : !complex -> !complex
    %PolR2 = quiccir.al.prj %PolR1T : !complex -> !complex attributes{kind = "Ll"}
    %PolR2T = quiccir.transpose %PolR2 permutation = [1, 2, 0] : !complex -> !complex
    %R = quiccir.fr.prj %PolR2T : !complex -> !real attributes{kind = "P"}
    // Theta
    %TorTh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "P"}
    %TorTh1T = quiccir.transpose %TorTh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh2 = quiccir.al.prj %TorTh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %TorTh2T = quiccir.transpose %TorTh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh3 = quiccir.fr.prj %TorTh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolTh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %PolTh1T = quiccir.transpose %PolTh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh2 = quiccir.al.prj %PolTh1T : !complex -> !complex attributes{kind = "D1"}
    %PolTh2T = quiccir.transpose %PolTh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh3 = quiccir.fr.prj %PolTh2T : !complex -> !real attributes{kind = "P"}
    //
    %Theta = quiccir.add %TorTh3, %PolTh3 : !real, !real -> !real
    // Phi
    %TorPh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "P"}
    %TorPh1T = quiccir.transpose %TorPh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh2 = quiccir.al.prj %TorPh1T : !complex -> !complex attributes{kind = "D1"}
    %TorPh2T = quiccir.transpose %TorPh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh3 = quiccir.fr.prj %TorPh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolPh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %PolPh1T = quiccir.transpose %PolPh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh2 = quiccir.al.prj %PolPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %PolPh2T = quiccir.transpose %PolPh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh3 = quiccir.fr.prj %PolPh2T : !complex -> !real attributes{kind = "P"}
    //
    %Phi = quiccir.sub %PolPh3, %TorPh3 : !real, !real -> !real
    return %R, %Theta, %Phi : !real, !real, !real
}

func.func private @bwdCurl(%Tor: !complex, %Pol: !complex) -> (!real, !real, !real) {
    // R
    %TorR1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %TorR1T = quiccir.transpose %TorR1 permutation = [1, 2, 0] : !complex -> !complex
    %TorR2 = quiccir.al.prj %TorR1T : !complex -> !complex attributes{kind = "Ll"}
    %TorR2T = quiccir.transpose %TorR2 permutation = [1, 2, 0] : !complex -> !complex
    %R = quiccir.fr.prj %TorR2T : !complex -> !real attributes{kind = "P"}
    // Theta
    %TorTh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %TorTh1T = quiccir.transpose %TorTh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh2 = quiccir.al.prj %TorTh1T : !complex -> !complex attributes{kind = "D1"}
    %TorTh2T = quiccir.transpose %TorTh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh3 = quiccir.fr.prj %TorTh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolTh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "SphLapl"}
    %PolTh1T = quiccir.transpose %PolTh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh2 = quiccir.al.prj %PolTh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %PolTh2T = quiccir.transpose %PolTh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh3 = quiccir.fr.prj %PolTh2T : !complex -> !real attributes{kind = "P"}
    //
    %Theta = quiccir.sub %TorTh3, %PolTh3 : !real, !real -> !real
    // Phi
    %TorPh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %TorPh1T = quiccir.transpose %TorPh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh2 = quiccir.al.prj %TorPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %TorPh2T = quiccir.transpose %TorPh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh3 = quiccir.fr.prj %TorPh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolPh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "SphLapl"}
    %PolPh1T = quiccir.transpose %PolPh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh2 = quiccir.al.prj %PolPh1T : !complex -> !complex attributes{kind = "D1"}
    %PolPh2T = quiccir.transpose %PolPh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh3 = quiccir.fr.prj %PolPh2T : !complex -> !real attributes{kind = "P"}
    //
    %Phi = quiccir.add %PolPh3, %TorPh3 : !real, !real -> !real
    return %R, %Theta, %Phi : !real, !real, !real
}

func.func private @fwdScalar(%S: !real) -> !complex {
    // forward scalar Nl path
    %S1 = quiccir.fr.int %S : !real -> !complex attributes{kind = "P"}
    %S1T = quiccir.transpose %S1 permutation = [2, 0, 1] : !complex -> !complex
    %S2 = quiccir.al.int %S1T : !complex -> !complex attributes{kind = "P"}
    %S2T = quiccir.transpose %S2 permutation = [2, 0, 1] : !complex -> !complex
    %S3 = quiccir.jw.int %S2T : !complex -> !complex attributes{kind = "I2"}

    return %S3 : !complex
}

func.func private @fwdVel(%R: !real, %Theta: !real, %Phi: !real) -> (!complex, !complex){
    // Tor
    %ThetaTor1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaTor1T = quiccir.transpose %ThetaTor1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor2 = quiccir.al.int %ThetaTor1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %ThetaTor2T = quiccir.transpose %ThetaTor2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor3 = quiccir.jw.int %ThetaTor2T : !complex -> !complex attributes{kind = "I2_Zero"}
    //
    %PhiTor1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiTor1T = quiccir.transpose %PhiTor1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor2 = quiccir.al.int %PhiTor1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %PhiTor2T = quiccir.transpose %PhiTor2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor3 = quiccir.jw.int %PhiTor2T : !complex -> !complex attributes{kind = "I2_Zero"}
    //
    %Tor = quiccir.sub %ThetaTor3, %PhiTor3 : !complex, !complex -> !complex

    // Pol
    %RPol1 = quiccir.fr.int %R : !real -> !complex attributes{kind = "P"}
    %RPol1T = quiccir.transpose %RPol1 permutation = [2, 0, 1] : !complex -> !complex
    %RPol2 = quiccir.al.int %RPol1T : !complex -> !complex attributes{kind = "P"}
    %RPol2T = quiccir.transpose %RPol2 permutation = [2, 0, 1] : !complex -> !complex
    %RPol3 = quiccir.jw.int %RPol2T : !complex -> !complex attributes{kind = "I4DivR1_Zero"}
    //
    %ThetaPol1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaPol1T = quiccir.transpose %ThetaPol1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol2 = quiccir.al.int %ThetaPol1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %ThetaPol2T = quiccir.transpose %ThetaPol2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol3 = quiccir.jw.int %ThetaPol2T : !complex -> !complex attributes{kind = "I4DivR1D1R1_Zero"}
    //
    %PhiPol1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiPol1T = quiccir.transpose %PhiPol1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol2 = quiccir.al.int %PhiPol1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %PhiPol2T = quiccir.transpose %PhiPol2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol3 = quiccir.jw.int %PhiPol2T : !complex -> !complex attributes{kind = "I4DivR1D1R1_Zero"}
    //
    %tmp = quiccir.add %ThetaPol3, %PhiPol3 : !complex, !complex -> !complex
    %Pol = quiccir.sub %tmp, %RPol3 : !complex, !complex -> !complex
    return %Tor, %Pol : !complex, !complex
}

func.func private @fwdMag(%R: !real, %Theta: !real, %Phi: !real) -> (!complex, !complex){
    // Tor
    %RTor1 = quiccir.fr.int %R : !real -> !complex attributes{kind = "P"}
    %RTor1T = quiccir.transpose %RTor1 permutation = [2, 0, 1] : !complex -> !complex
    %RTor2 = quiccir.al.int %RTor1T : !complex -> !complex attributes{kind = "P"}
    %RTor2T = quiccir.transpose %RTor2 permutation = [2, 0, 1] : !complex -> !complex
    %RTor3 = quiccir.jw.int %RTor2T : !complex -> !complex attributes{kind = "I2DivR1_Zero"}
    //
    %ThetaTor1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaTor1T = quiccir.transpose %ThetaTor1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor2 = quiccir.al.int %ThetaTor1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %ThetaTor2T = quiccir.transpose %ThetaTor2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor3 = quiccir.jw.int %ThetaTor2T : !complex -> !complex attributes{kind = "I2DivR1D1R1_Zero"}
    //
    %PhiTor1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiTor1T = quiccir.transpose %PhiTor1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor2 = quiccir.al.int %PhiTor1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %PhiTor2T = quiccir.transpose %PhiTor2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor3 = quiccir.jw.int %PhiTor2T : !complex -> !complex attributes{kind = "I2DivR1D1R1_Zero"}
    //
    %tmpTor = quiccir.sub %RTor3, %ThetaTor3 : !complex, !complex -> !complex
    %Tor = quiccir.sub %tmpTor, %PhiTor3 : !complex, !complex -> !complex

    // Pol
    %ThetaPol1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaPol1T = quiccir.transpose %ThetaPol1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol2 = quiccir.al.int %ThetaPol1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %ThetaPol2T = quiccir.transpose %ThetaPol2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol3 = quiccir.jw.int %ThetaPol2T : !complex -> !complex attributes{kind = "I2_Zero"}
    //
    %PhiPol1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiPol1T = quiccir.transpose %PhiPol1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol2 = quiccir.al.int %PhiPol1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %PhiPol2T = quiccir.transpose %PhiPol2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol3 = quiccir.jw.int %PhiPol2T : !complex -> !complex attributes{kind = "I2_Zero"}
    //
    %Pol = quiccir.sub %ThetaPol3, %PhiPol3 : !complex, !complex -> !complex
    return %Tor, %Pol : !complex, !complex
}

func.func private @nlScalar(%UR: !real, %UTheta: !real, %UPhi: !real,
    %TdR: !real, %TdTheta: !real, %TdPhi: !real) -> !real {
    // U dot grad T
    %DotT = quiccir.dot(%UR, %UTheta, %UPhi), (%TdR, %TdTheta, %TdPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        !real
        attributes{kind = "transport"}
    // U dot R
    %DotR = quiccir.mul.const %UR : !real -> !real attributes{kind = "transport"}
    %TPhysNl = quiccir.sub %DotT, %DotR : !real, !real -> !real
    return %TPhysNl : !real
}

func.func private @nlVel(%UR: !real, %UTheta: !real, %UPhi: !real,
    %CurlUR: !real, %CurlUTheta: !real, %CurlUPhi: !real,
    %BR: !real, %BTheta: !real, %BPhi: !real,
    %CurlBR: !real, %CurlBTheta: !real, %CurlBPhi: !real,
    %T: !real) -> (!real, !real, !real) {
    // Cross Vel
    %CrossU:3 = quiccir.cross(%CurlUR, %CurlUTheta, %CurlUPhi), (%UR, %UTheta, %UPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        (!real, !real, !real)
        attributes{kind = "inertia"}
    // Cross Mag
     %CrossB:3 = quiccir.cross(%BR, %BTheta, %BPhi), (%CurlBR, %CurlBTheta, %CurlBPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        (!real, !real, !real)
        attributes{kind = "lorentz"}
    // Sum by component
    %Cross0 = quiccir.add %CrossU#0, %CrossB#0 : !real, !real -> !real
    %Cross1 = quiccir.add %CrossU#1, %CrossB#1 : !real, !real -> !real
    %Cross2 = quiccir.add %CrossU#2, %CrossB#2 : !real, !real -> !real
    // Add buoyancy
    %Buoy = quiccir.mul.const %T : !real -> !real attributes{kind = "buoyancy"}
    %RTmp = quiccir.sub %Cross0, %Buoy : !real, !real -> !real
    // Add coriolis
    %Cor0 = quiccir.mul.const %UPhi : !real -> !real attributes{kind = "coriolis_sin"}
    %RNl = quiccir.sub %RTmp, %Cor0 : !real, !real -> !real
    %Cor1 = quiccir.mul.const %UPhi : !real -> !real attributes{kind = "coriolis_cos"}
    %ThetaNl = quiccir.sub %Cross1, %Cor1 : !real, !real -> !real
    %Cor2a = quiccir.mul.const %UR : !real -> !real attributes{kind = "coriolis_sin"}
    %PhiTmp = quiccir.add %Cross2, %Cor2a : !real, !real -> !real
    %Cor2b = quiccir.mul.const %UTheta : !real -> !real attributes{kind = "coriolis_cos"}
    %PhiNl = quiccir.add %PhiTmp, %Cor2b : !real, !real -> !real
    return %RNl, %ThetaNl, %PhiNl : !real, !real, !real
}

func.func private @nlMag(%UR: !real, %UTheta: !real, %UPhi: !real,
    %MagR: !real, %MagTheta: !real, %MagPhi: !real) -> (!real, !real, !real) {
    // Cross
    %Cross:3 = quiccir.cross (%MagR, %MagTheta, %MagPhi), (%UR, %UTheta, %UPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        (!real, !real, !real)
        attributes{kind = "induction"}
    return %Cross#0, %Cross#1, %Cross#2  : !real, !real, !real
}

func.func @entry(%T: !complex, %TorVel: !complex, %PolVel: !complex,
   %TorMag: !complex, %PolMag: !complex) -> (!complex, !complex, !complex, !complex, !complex,
      !real, !real, !real, !real, !real, !real) {
    %TPhys = call @bwdScalar(%T) : (!complex) -> !real
    %TGrad:3 = call @bwdGradScalar(%T) : (!complex) -> (!real, !real, !real)
    %Vel:3 = call @bwdVector(%TorVel, %PolVel) : (!complex, !complex) -> (!real, !real, !real)
    %CurlVel:3 = call @bwdCurl(%TorVel, %PolVel) : (!complex, !complex) -> (!real, !real, !real)
    %Mag:3 = call @bwdVector(%TorMag, %PolMag) : (!complex, !complex) -> (!real, !real, !real)
    %CurlMag:3 = call @bwdCurl(%TorMag, %PolMag) : (!complex, !complex) -> (!real, !real, !real)
    %TPhysNl = call @nlScalar(%Vel#0, %Vel#1, %Vel#2, %TGrad#0, %TGrad#1, %TGrad#2) : (!real, !real, !real, !real, !real, !real) -> !real
    %VelNl:3 = call @nlVel(%Vel#0, %Vel#1, %Vel#2, %CurlVel#0, %CurlVel#1, %CurlVel#2, %Mag#0, %Mag#1, %Mag#2, %CurlMag#0, %CurlMag#1, %CurlMag#2, %TPhys) : (!real, !real, !real, !real, !real, !real, !real, !real, !real, !real, !real, !real, !real) -> (!real, !real, !real)
    %MagNl:3 = call @nlMag(%Vel#0, %Vel#1, %Vel#2, %Mag#0, %Mag#1, %Mag#2) : (!real, !real, !real, !real, !real, !real) -> (!real, !real, !real)
    %TNl = call @fwdScalar(%TPhysNl) : (!real) -> !complex
    %TorVelNl, %PolVelNl = call @fwdVel(%VelNl#0, %VelNl#1, %VelNl#2) : (!real, !real, !real) -> (!complex, !complex)
    %TorMagNl, %PolMagNl = call @fwdMag(%MagNl#0, %MagNl#1, %MagNl#2) : (!real, !real, !real) -> (!complex, !complex)
    return %TNl, %TorVelNl, %PolVelNl, %TorMagNl, %PolMagNl, %Vel#0, %Vel#1, %Vel#2, %Mag#0, %Mag#1, %Mag#2: !complex, !complex, !complex, !complex, !complex, !real, !real, !real, !real, !real, !real
}
   )mlir";

   MHDFloat T = 1.0 / spSim->eqParams()->nd(NonDimensional::Ekman::id());
   MHDFloat Ra = spSim->eqParams()->nd(NonDimensional::Rayleigh::id());
   MHDFloat Pr = spSim->eqParams()->nd(NonDimensional::Prandtl::id());
   MHDFloat Pm = spSim->eqParams()->nd(NonDimensional::MagneticPrandtl::id());
   Graph::PhysicalParameters<MHDFloat> physParams;
   physParams.transport = 1.0;
   physParams.inertia = 1.0;
   physParams.buoyancy = Pm * Pm * Ra * T / Pr;
   physParams.coriolis = T * Pm;
   physParams.lorentz = T * Pm;
   physParams.induction = 1.0;
   spSim->addGraph(graphStr, physParams);
   #endif
}

std::map<std::string, std::map<std::string, int>>
IDynamoModel::configTags() const
{
   std::map<std::string, int> onOff;
   onOff.emplace("enable", 1);

   std::map<std::string, int> options;
   options.emplace("enable", 0);
   options.emplace("numbered", 0);
   options.emplace("only_every", 1);

   std::map<std::string, std::map<std::string, int>> tags;
   // temperature
   tags.emplace("temperature_energy", onOff);
   tags.emplace("temperature_l_spectrum", options);
   tags.emplace("temperature_m_spectrum", options);
   tags.emplace("temperature_n_spectrum", options);
   // kinetic
   tags.emplace("kinetic_energy", onOff);
   tags.emplace("kinetic_l_spectrum", options);
   tags.emplace("kinetic_m_spectrum", options);
   tags.emplace("kinetic_n_spectrum", options);
   // magnetic
   tags.emplace("magnetic_energy", onOff);
   tags.emplace("magnetic_l_spectrum", options);
   tags.emplace("magnetic_m_spectrum", options);
   tags.emplace("magnetic_n_spectrum", options);
   // diagnostic
   tags.emplace("angular_momentum", onOff);
   tags.emplace("nusselt", onOff);

   return tags;
}

void IDynamoModel::addAsciiOutputFiles(SharedSimulation spSim)
{
   // Create Nusselt writer
   this->enableAsciiFile<Io::Variable::SphereNusseltWriter>("nusselt", "",
      PhysicalNames::Temperature::id(), spSim);

   // Create temperature energy writer
   this->enableAsciiFile<Io::Variable::SphereScalarEnergyWriter>(
      "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature L energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarLSpectrumWriter>(
      "temperature_l_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature M energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarMSpectrumWriter>(
      "temperature_m_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature N power spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarNSpectrumWriter>(
      "temperature_n_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create kinetic energy writer
   this->enableAsciiFile<Io::Variable::SphereTorPolEnergyWriter>(
      "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolLSpectrumWriter>(
      "kinetic_l_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolMSpectrumWriter>(
      "kinetic_m_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic N power spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolNSpectrumWriter>(
      "kinetic_n_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create magnetic energy writer
   this->enableAsciiFile<Io::Variable::SphereTorPolEnergyWriter>(
      "magnetic_energy", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create magnetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolLSpectrumWriter>(
      "magnetic_l_spectrum", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create magnetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolMSpectrumWriter>(
      "magnetic_m_spectrum", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create kinetic N power spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolNSpectrumWriter>(
      "magnetic_n_spectrum", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create angular momentum writer
   this->enableAsciiFile<Io::Variable::SphereAngularMomentumWriter>(
      "angular_momentum", "", PhysicalNames::Velocity::id(), spSim);
}

} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
