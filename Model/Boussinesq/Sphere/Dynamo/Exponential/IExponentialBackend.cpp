/**
 * @file IExponentialBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/IExponentialBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/QuasiInverseOnly.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/JacobianMagnetic.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
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

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

namespace {
   const auto mag_tor = std::make_pair(PhysicalNames::Magnetic::id(),
                   FieldComponents::Spectral::TOR);
   const auto jmag_tor = std::make_pair(PhysicalNames::JacobianMagnetic::id(),
                   FieldComponents::Spectral::TOR);
   const auto mag_pol = std::make_pair(PhysicalNames::Magnetic::id(),
                   FieldComponents::Spectral::POL);
   const auto jmag_pol = std::make_pair(PhysicalNames::JacobianMagnetic::id(),
                   FieldComponents::Spectral::POL);
   const auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   const auto jvel_tor = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::TOR);
   const auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   const auto jvel_pol = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::POL);
   const auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);
   const auto jtemp = std::make_pair(PhysicalNames::JacobianTemperature::id(),
                   FieldComponents::Spectral::SCALAR);
}

std::vector<std::string> IExponentialBackend::fieldNames() const
{
   std::vector<std::string> names = {
      PhysicalNames::Velocity().tag(),
      PhysicalNames::Temperature().tag(),
      PhysicalNames::Magnetic().tag(),
      PhysicalNames::JacobianVelocity().tag(),
      PhysicalNames::JacobianTemperature().tag(),
      PhysicalNames::JacobianMagnetic().tag()};

   return names;
}

int IExponentialBackend::nBc(const SpectralFieldId& fId) const
{
   int nBc = 0;

   if (fId == vel_tor ||
       fId == jvel_tor ||
       fId == temp ||
       fId == jtemp ||
       fId == mag_tor ||
       fId == jmag_tor ||
       fId == mag_pol ||
       fId == jmag_pol)
   {
      nBc = 1;
   }
   else if (fId == vel_pol ||
            fId == jvel_pol)
   {
      nBc = 2;
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

void IExponentialBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l,
   std::shared_ptr<details::BlockOptions> opts, const int nN,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto bcId = bcs.find(rowId.first)->second;

   SparseSM::Worland::Boundary::Operator bcOp(nN, nN, a, b, l);

   if ((rowId == vel_tor || rowId == jvel_tor) &&
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
   else if ((rowId == vel_pol || rowId == jvel_pol) &&
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
   else if ((rowId == mag_tor || rowId == jmag_tor) &&
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
   else if ((rowId == mag_pol || rowId == jmag_pol) &&
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
   else if ((rowId == temp || rowId == jtemp) &&
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
   else
   {
      throw std::logic_error("Unknown field for boundary conditions");
   }

   mat.real() += bcOp.mat();
}

void IExponentialBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const int nN, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto bcId = bcs.find(fieldId.first)->second;

   int s = this->nBc(fieldId);
   if(bcId == Bc::Name::QuasiInverseOnly::id())
   {
      SparseSM::Worland::Id qid(nN, nN - s, a, b, l);
      mat = qid.mat();
   }
   else
   {
      if (fieldId == vel_tor || fieldId == jvel_tor)
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
      else if (fieldId == vel_pol || fieldId == jvel_pol)
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
      else if (fieldId == mag_tor || fieldId == jmag_tor)
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
      else if (fieldId == mag_pol || fieldId == jmag_pol)
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
      else if (fieldId == temp || fieldId == jtemp)
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
   }

   if (makeSquare)
   {
      SparseSM::Worland::Id qId(nN - s, nN, a, b, l);
      mat = qId.mat() * mat;
   }
}

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
