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
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkVelC2.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkMagC2.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/SpectralKernels/IoValidation.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"

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
      // Shared pointer to equation
      Equations::SharedSphereExactScalarState spScalar;
      Equations::SharedSphereExactVectorState spVector;

      Spectral::Kernel::Complex3DMapType tSH;
      std::pair<Spectral::Kernel::Complex3DMapType::iterator,bool> ptSH;

      // Add temperature initial state generator
      spScalar = spGen->addEquation<Equations::SphereExactScalarState>(this->spBackend());
      spScalar->setIdentity(PhysicalNames::Temperature::id());
      switch(3)
      {
         case 0:
            {
               spScalar->setPhysicalNoise(1e-15);
            }
            break;

         case 1:
            {
               spScalar->setPhysicalConstant(1.0);
            }
            break;

         case 2:
            {
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(3,3), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0,2.0)));
               spScalar->setSpectralModes(tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Physical::Kernel::Sphere::BenchmarkTempC1>();
               spKernel->init(0.0, 1e-5);
               spScalar->setPhysicalKernel(spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::IoValidation>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               spKernel->init(100);
               spScalar->setSrcKernel(spKernel);
            }
            break;

         case 5:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spScalar->setSrcKernel(spKernel);
            }
            break;
      }

      // Add velocity initial state generator
      spVector = spGen->addEquation<Equations::SphereExactVectorState>(this->spBackend());
      spVector->setIdentity(PhysicalNames::Velocity::id());
      switch(3)
      {
         // Toroidal only
         case 0:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Poloidal only
         case 1:
            {
               // Toroidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Toroidal & Poloidal
         case 2:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Physical::Kernel::Sphere::BenchmarkVelC2>();
               spKernel->init();
               spVector->setPhysicalKernel(spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::IoValidation>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               spKernel->init(100);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;

         case 5:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;
      }

      // Add magnetic initial state generator
      spVector = spGen->addEquation<Equations::SphereExactVectorState>(this->spBackend());
      spVector->setIdentity(PhysicalNames::Magnetic::id());
      switch(3)
      {
         // Toroidal only
         case 0:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Poloidal only
         case 1:
            {
               // Toroidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Toroidal & Poloidal
         case 2:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Physical::Kernel::Sphere::BenchmarkMagC2>();
               spKernel->init();
               spVector->setPhysicalKernel(spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::IoValidation>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               spKernel->init(100);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;

         case 5:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;
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
      spScalar = spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
      spScalar->setFields(true, true);
      spScalar->setIdentity(PhysicalNames::Temperature::id());

      // Add velocity field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::Velocity::id());

      // Add magnetic field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
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

#if 1
      // Create temperature L energy spectrum writer
      auto spTempL = std::make_shared<Io::Variable::SphereScalarLSpectrumWriter>("temperature", spSim->ss().tag());
      spTempL->expect(PhysicalNames::Temperature::id());
      spTempL->numberOutput();
      spTempL->onlyEvery(100);
      spSim->addAsciiOutputFile(spTempL);

      // Create temperature M energy spectrum writer
      auto spTempM = std::make_shared<Io::Variable::SphereScalarMSpectrumWriter>("temperature", spSim->ss().tag());
      spTempM->expect(PhysicalNames::Temperature::id());
      spTempM->numberOutput();
      spTempM->onlyEvery(100);
      spSim->addAsciiOutputFile(spTempM);
#endif

#if 1
      // Create temperature N power spectrum writer
      auto spTempN = std::make_shared<Io::Variable::SphereScalarNSpectrumWriter>("temperature", spSim->ss().tag());
      spTempN->expect(PhysicalNames::Temperature::id());
      spTempN->numberOutput();
      spTempN->onlyEvery(100);
      spSim->addAsciiOutputFile(spTempN);
#endif

      // Create kinetic energy writer
      auto spKinetic = std::make_shared<Io::Variable::SphereTorPolEnergyWriter>("kinetic", spSim->ss().tag());
      spKinetic->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spKinetic);

#if 1
      // Create kinetic L energy spectrum writer
      auto spKineticL = std::make_shared<Io::Variable::SphereTorPolLSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticL->expect(PhysicalNames::Velocity::id());
      spKineticL->numberOutput();
      spKineticL->onlyEvery(100);
      spSim->addAsciiOutputFile(spKineticL);

      // Create kinetic M energy spectrum writer
      auto spKineticM = std::make_shared<Io::Variable::SphereTorPolMSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticM->expect(PhysicalNames::Velocity::id());
      spKineticM->numberOutput();
      spKineticM->onlyEvery(100);
      spSim->addAsciiOutputFile(spKineticM);
#endif

#if 1
      // Create kinetic N power spectrum writer
      auto spKineticN = std::make_shared<Io::Variable::SphereTorPolNSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticN->expect(PhysicalNames::Velocity::id());
      spKineticN->numberOutput();
      spKineticN->onlyEvery(100);
      spSim->addAsciiOutputFile(spKineticN);
#endif

      // Create magnetic energy writer
      auto spMagnetic = std::make_shared<Io::Variable::SphereTorPolEnergyWriter>("magnetic", spSim->ss().tag());
      spMagnetic->expect(PhysicalNames::Magnetic::id());
      spSim->addAsciiOutputFile(spMagnetic);

#if 1
      // Create magnetic L energy spectrum writer
      auto spMagneticL = std::make_shared<Io::Variable::SphereTorPolLSpectrumWriter>("magnetic", spSim->ss().tag());
      spMagneticL->expect(PhysicalNames::Magnetic::id());
      spMagneticL->numberOutput();
      spMagneticL->onlyEvery(100);
      spSim->addAsciiOutputFile(spMagneticL);

      // Create magnetic M energy spectrum writer
      auto spMagneticM = std::make_shared<Io::Variable::SphereTorPolMSpectrumWriter>("magnetic", spSim->ss().tag());
      spMagneticM->expect(PhysicalNames::Magnetic::id());
      spMagneticM->numberOutput();
      spMagneticM->onlyEvery(100);
      spSim->addAsciiOutputFile(spMagneticM);
#endif

#if 1
      // Create kinetic N power spectrum writer
      auto spMagneticN = std::make_shared<Io::Variable::SphereTorPolNSpectrumWriter>("magnetic", spSim->ss().tag());
      spMagneticN->expect(PhysicalNames::Magnetic::id());
      spMagneticN->numberOutput();
      spMagneticN->onlyEvery(100);
      spSim->addAsciiOutputFile(spMagneticN);
#endif

#if 1
      // Create angular momentum writer
      auto spAngMom = std::make_shared<Io::Variable::SphereAngularMomentumWriter>("", spSim->ss().tag());
      spAngMom->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spAngMom);
#endif
   }

}
}
}
}
}
