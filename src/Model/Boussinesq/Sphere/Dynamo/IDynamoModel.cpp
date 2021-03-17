/**
 * @file IDynamoModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a sphere (Toroidal/Poloidal formulation)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/IDynamoModel.hpp"

// Project includes
//
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Transport.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Momentum.hpp"
#include "QuICC/Model/Boussinesq/Sphere/Dynamo/Induction.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/MagPrandtl.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/CflTorsional.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/SphereExactStateIds.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

   VectorFormulation::Id IDynamoModel::SchemeFormulation()
   {
      return VectorFormulation::TORPOL;
   }

   void IDynamoModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Transport>(this->spBackend());

      // Add Navier-Stokes equation
      spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Momentum>(this->spBackend());

      // Add induction equation
      spSim->addEquation<Equations::Boussinesq::Sphere::Dynamo::Induction>(this->spBackend());
   }

   void IDynamoModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(true)
      {
         // Shared pointer to equation
         Equations::SharedSphereExactScalarState spScalar;
         Equations::SharedSphereExactVectorState spVector;

         Equations::SHMapType tSH;
         std::pair<Equations::SHMapType::iterator,bool> ptSH;

         // Add temperature initial state generator
         spScalar = spGen->addEquation<Equations::SphereExactScalarState>();
         spScalar->setIdentity(PhysicalNames::Temperature::id());
         switch(1)
         {
            case 0:
               spScalar->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(3,3), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0,2.0)));
               spScalar->setHarmonicOptions(tSH);
               break;

            case 1:
               spScalar->setStateType(Equations::SphereExactStateIds::BENCHTEMPC1);
               break;
         }

         // Add velocity initial state generator
         spVector = spGen->addEquation<Equations::SphereExactVectorState>();
         spVector->setIdentity(PhysicalNames::Velocity::id());
         switch(3)
         {
            case 0:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               break;

            case 1:
               // Poloidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 2:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 3:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHVELC2);
               break;
         }

         // Add magnetic initial state generator
         spVector = spGen->addEquation<Equations::SphereExactVectorState>();
         spVector->setIdentity(PhysicalNames::Magnetic::id());
         switch(3)
         {
            case 0:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               break;

            case 1:
               // Poloidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 2:
               // Toroidal
               spVector->setSpectralType(Equations::SphereExactStateIds::HARMONIC);
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setHarmonicOptions(FieldComponents::Spectral::POL, tSH);
               break;

            case 3:
               spVector->setStateType(Equations::SphereExactStateIds::BENCHMAGC2);
               break;
         }

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add velocity random initial state generator
         spVector = spGen->addEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::Velocity::id());
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add magnetic random initial state generator
         spVector = spGen->addEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::Magnetic::id());
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-4, 1e-4, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-4, 1e-4, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::Temperature::id());
         spScalar->setSpectrum(-1e-4, 1e-4, 1e4, 1e4, 1e4);
      }

      // Add output file
      auto spOut = std::make_shared<Io::Variable::StateFileWriter>(spGen->ss().tag(), spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spOut->expect(PhysicalNames::Magnetic::id());
      spGen->addHdf5OutputFile(spOut);
   }

   void IDynamoModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, true);
      spScalar->setIdentity(PhysicalNames::Temperature::id());

      // Add velocity field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::Velocity::id());

      // Add magnetic field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::Magnetic::id());

      // Add output file
      auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(spVis->ss().tag());
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spOut->expect(PhysicalNames::Magnetic::id());
      spVis->addHdf5OutputFile(spOut);
   }

   void IDynamoModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create Nusselt writer
      auto spNusselt = std::make_shared<Io::Variable::SphereNusseltWriter>("", spSim->ss().tag());
      spNusselt->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spNusselt);

      // Create temperature energy writer
      auto spTemp = std::make_shared<Io::Variable::SphereScalarEnergyWriter>("temperature", spSim->ss().tag());
      spTemp->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spTemp);

#if 0
      // Create temperature L energy spectrum writer
      auto spTempL = std::make_shared<Io::Variable::SphereScalarLSpectrumWriter>("temperature", spSim->ss().tag());
      spTempL->expect(PhysicalNames::Temperature::id());
      //spTempL->numberOutput();
      spSim->addAsciiOutputFile(spTempL);

      // Create temperature M energy spectrum writer
      auto spTempM = std::make_shared<Io::Variable::SphereScalarMSpectrumWriter>("temperature", spSim->ss().tag());
      spTempM->expect(PhysicalNames::Temperature::id());
      //spTempM->numberOutput();
      spSim->addAsciiOutputFile(spTempM);
#endif

      // Create kinetic energy writer
      auto spKinetic = std::make_shared<Io::Variable::SphereTorPolEnergyWriter>("kinetic", spSim->ss().tag());
      spKinetic->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spKinetic);

#if 0
      // Create kinetic L energy spectrum writer
      auto spKineticL = std::make_shared<Io::Variable::SphereTorPolLSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticL->expect(PhysicalNames::Velocity::id());
      //spKineticL->numberOutput();
      spSim->addAsciiOutputFile(spKineticL);

      // Create kinetic M energy spectrum writer
      auto spKineticM = std::make_shared<Io::Variable::SphereTorPolMSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticM->expect(PhysicalNames::Velocity::id());
      //spKineticM->numberOutput();
      spSim->addAsciiOutputFile(spKineticM);
#endif

      // Create magnetic energy writer
      auto spMagnetic = std::make_shared<Io::Variable::SphereTorPolEnergyWriter>("magnetic", spSim->ss().tag());
      spMagnetic->expect(PhysicalNames::Magnetic::id());
      spSim->addAsciiOutputFile(spMagnetic);

#if 0
      // Create magnetic L energy spectrum writer
      auto spMagneticL = std::make_shared<Io::Variable::SphereTorPolLSpectrumWriter>("magnetic", spSim->ss().tag());
      spMagneticL->expect(PhysicalNames::Magnetic::id());
      //spMagneticL->numberOutput();
      spSim->addAsciiOutputFile(spMagneticL);

      // Create magnetic M energy spectrum writer
      auto spMagneticM = std::make_shared<Io::Variable::SphereTorPolMSpectrumWriter>("magnetic", spSim->ss().tag());
      spMagneticM->expect(PhysicalNames::Magnetic::id());
      //spMagneticM->numberOutput();
      spSim->addAsciiOutputFile(spMagneticM);
#endif
   }

}
}
}
}
}
