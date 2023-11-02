/**
 * @file IDynamoBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/IDynamoBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/CflTorsional.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Worland/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

std::vector<std::string> IDynamoBackend::fieldNames() const
{
   std::vector<std::string> names = {
      PhysicalNames::Velocity().tag(),
      PhysicalNames::Temperature().tag(),
      PhysicalNames::Magnetic().tag(),
   };

   return names;
}

std::vector<std::string> IDynamoBackend::paramNames() const
{
   std::vector<std::string> names = {NonDimensional::MagneticPrandtl().tag(),
      NonDimensional::Ekman().tag(), NonDimensional::Prandtl().tag(),
      NonDimensional::Rayleigh().tag()};

   return names;
}

std::vector<bool> IDynamoBackend::isPeriodicBox() const
{
   std::vector<bool> periodic = {false, false, false};

   return periodic;
}

std::map<std::string, MHDFloat> IDynamoBackend::automaticParameters(
   const std::map<std::string, MHDFloat>& cfg) const
{
   auto E = cfg.find(NonDimensional::Ekman().tag())->second;
   auto Pm = cfg.find(NonDimensional::MagneticPrandtl().tag())->second;

   std::map<std::string, MHDFloat> params = {
      {NonDimensional::CflInertial().tag(), 0.1 * E / Pm},
      {NonDimensional::CflTorsional().tag(), 0.1 * std::sqrt(E)},
      {NonDimensional::CflAlfvenScale().tag(), Pm / E},
      {NonDimensional::CflAlfvenDamping().tag(), (1.0 + Pm) / 2.0}};

   return params;
}

int IDynamoBackend::nBc(const SpectralFieldId& fId) const
{
   int nBc = 0;

   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::TOR) ||
       fId == std::make_pair(PhysicalNames::Temperature::id(),
                 FieldComponents::Spectral::SCALAR) ||
       fId == std::make_pair(PhysicalNames::Magnetic::id(),
                 FieldComponents::Spectral::TOR) ||
       fId == std::make_pair(PhysicalNames::Magnetic::id(),
                 FieldComponents::Spectral::POL))
   {
      nBc = 1;
   }
   else if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                      FieldComponents::Spectral::POL))
   {
      nBc = 2;
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

void IDynamoBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l, const Resolution& res,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
   auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

   auto bcId = bcs.find(rowId.first)->second;

   SparseSM::Worland::Boundary::Operator bcOp(nN, nN, a, b, l);

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR) &&
       rowId == colId)
   {
      if (l > 0)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::Value>();
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::R1D1DivR1>();
         }
         else
         {
            throw std::logic_error("Boundary conditions for Velocity Toroidal "
                                   "component not implemented");
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL) &&
            rowId == colId)
   {
      if (l > 0)
      {
         if (this->useSplitEquation())
         {
            if (isSplitOperator)
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
            }
            else if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::D1>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::D2>();
            }
            else
            {
               throw std::logic_error("Boundary conditions for Velocity "
                                      "Poloidal component not implemented");
            }
         }
         else
         {
            if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
               bcOp.addRow<SparseSM::Worland::Boundary::D1>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
               bcOp.addRow<SparseSM::Worland::Boundary::D2>();
            }
            else
            {
               throw std::logic_error("Boundary conditions for Velocity "
                                      "Poloidal component not implemented");
            }
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::TOR) &&
            rowId == colId)
   {
      if (l > 0)
      {
         if (bcId == Bc::Name::Insulating::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::Value>();
         }
         else
         {
            throw std::logic_error("Boundary conditions for Magnetic Toroidal "
                                   "component not implemented");
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::POL) &&
            rowId == colId)
   {
      if (l > 0)
      {
         if (bcId == Bc::Name::Insulating::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::InsulatingSphere>();
         }
         else
         {
            throw std::logic_error("Boundary conditions for Magnetic Poloidal "
                                   "component not implemented");
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR) &&
            rowId == colId)
   {
      if (bcId == Bc::Name::FixedTemperature::id())
      {
         bcOp.addRow<SparseSM::Worland::Boundary::Value>();
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         bcOp.addRow<SparseSM::Worland::Boundary::D1>();
      }
      else
      {
         throw std::logic_error(
            "Boundary conditions for Temperature not implemented (" +
            std::to_string(bcId) + ")");
      }
   }

   mat.real() += bcOp.mat();
}

void IDynamoBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const Resolution& res, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
   auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

   auto bcId = bcs.find(fieldId.first)->second;

   int s = this->nBc(fieldId);
   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         SparseSM::Worland::Stencil::R1D1DivR1 bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerkin boundary conditions for Velocity "
                                "Toroidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         SparseSM::Worland::Stencil::ValueD1 bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         SparseSM::Worland::Stencil::ValueD2 bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerin boundary conditions for Velocity "
                                "Poloidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::TOR))
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerkin boundary conditions for Magnetic "
                                "Toroidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::POL))
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         SparseSM::Worland::Stencil::InsulatingSphere bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerin boundary conditions for Magnetic "
                                "Poloidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      if (bcId == Bc::Name::FixedTemperature::id())
      {
         SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         SparseSM::Worland::Stencil::D1 bc(nN, nN - s, a, b, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error(
            "Galerkin boundary conditions for Temperature not implemented");
      }
   }

   if (makeSquare)
   {
      SparseSM::Worland::Id qId(nN - s, nN, a, b, l);
      mat = qId.mat() * mat;
   }
}

void IDynamoBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, const Resolution& res, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, lr)(0);

   auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
   auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

   auto S = mat;
   this->stencil(S, colId, lc, res, false, bcs, nds);

   auto s = this->nBc(rowId);
   SparseSM::Worland::Id qId(nNr - s, nNr, a, b, lr, 0, s);
   mat = qId.mat() * (mat * S);
}

} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
