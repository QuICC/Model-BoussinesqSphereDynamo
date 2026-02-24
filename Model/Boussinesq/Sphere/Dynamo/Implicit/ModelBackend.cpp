/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Implicit/ModelBackend.hpp"
#include "QuICC/Bc/Name/QuasiInverseOnly.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/QuasiInverse.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I2Qm.hpp"
#include "QuICC/SparseSM/Worland/I2Qp.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/I4Qm.hpp"
#include "QuICC/SparseSM/Worland/I4Qp.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Implicit {

namespace implDetails {

/**
 * @brief Specific options for current model
 */
struct BlockOptionsImpl : public details::BlockOptions
{
   /**
    * @brief default ctor
    */
   BlockOptionsImpl() = default;

   /**
    * @brief default dtor
    */
   virtual ~BlockOptionsImpl() = default;

   /// Jones-Worland alpha
   Internal::MHDFloat a;
   /// Jones-Worland beta
   Internal::MHDFloat b;
   /// Harmonic order m
   int m;
   /// Use truncated quasi-inverse?
   bool truncateQI;
   /// Boundary condition
   std::size_t bcId;
   /// Split operator for influence matrix?
   bool isSplitOperator;
};
} // namespace implDetails

ModelBackend::ModelBackend() :
    IDynamoBackend(),
    mcTruncateQI(true)
{}

void ModelBackend::enableSplitEquation(const bool flag)
{
   if (flag)
   {
      throw std::logic_error(
         "Split equation for implicit model is not implemented");
   }
   else
   {
      IDynamoBackend::enableSplitEquation(flag);
   }
}

bool ModelBackend::isComplex(const SpectralFieldId& fId) const
{
   return true;
}

ModelBackend::SpectralFieldIds ModelBackend::implicitFields(
   const SpectralFieldId& fId) const
{
   SpectralFieldId velTor = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::TOR);
   SpectralFieldId velPol = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::POL);
   SpectralFieldId temp = std::make_pair(PhysicalNames::Temperature::id(),
      FieldComponents::Spectral::SCALAR);

   SpectralFieldIds fields;
   if (fId == velTor || fId == velPol || fId == temp)
   {
      fields = {velTor, velPol, temp};
   }
   else
   {
      fields = {fId};
   }

   // Make sure fields are sorted
   std::sort(fields.begin(), fields.end());
   return fields;
}

void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId,
   const Resolution& res) const
{
   // Operators are real
   info.isComplex = this->isComplex(fId);

   // Operators are real
   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::POL))
   {
      info.isSplitEquation = this->useSplitEquation();
   }
   else
   {
      info.isSplitEquation = false;
   }

   // Implicit coupled fields
   info.im = this->implicitFields(fId);

   // Explicit linear terms
   info.exL.clear();

   // Explicit nonlinear terms
   info.exNL.clear();

   // Explicit nextstep terms
   info.exNS.clear();

   // Index mode
   info.indexMode =
      static_cast<int>(Equations::CouplingIndexType::SLOWEST_SINGLE_RHS);
}

details::BlockDefinition ModelBackend::implicitBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator) const
{
   details::BlockDefinition blkDef;
   blkDef.rowId = rowId;
   blkDef.colId = colId;
   blkDef.isComplex = this->isComplex(rowId);
   blkDef.isGalerkin = this->useGalerkin();

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      blkDef.descr.push_back({});
      auto& d = blkDef.descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::worland_default_t::ALPHA;
      opts->b = Polynomial::Worland::worland_default_t::DBETA;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplitOperator;
      d.opts = opts;

      return d;
   };

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR))
   {
      if (rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               bMat = Pm * i2lapl.mat();
            }

            return bMat;
         };

         // Imaginary part of operator
         auto imagOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));

               SparseSM::Worland::I2 i2(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               bMat = o.m * T * Pm * invlapl * i2.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = imagOp;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::POL))
      {
         // Real part of first lower diagonal
         auto realOpLower = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));

               SparseSM::Worland::I2Qm corQm(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               auto norm = coriolis(l, o.m);
               bMat =
                  -static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQm.mat();
            }

            return bMat;
         };

         // Create first lower diagonal operator
         auto& dLow = getDescription();
         dLow.nRowShift = 1;
         dLow.nColShift = 0;
         dLow.realOp = realOpLower;
         dLow.imagOp = nullptr;

         // Real part of first upper diagonal
         auto realOpUpper = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = MHD_MP(1.0) / (dl * (dl + MHD_MP(1.0)));
               SparseSM::Worland::I2Qp corQp(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               auto norm = -coriolis(l + 1, o.m);
               bMat =
                  -static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQp.mat();
            }

            return bMat;
         };

         // Create first upper diagonal operator
         auto& dUp = getDescription();
         dUp.nRowShift = 0;
         dUp.nColShift = 1;
         dUp.realOp = realOpUpper;
         dUp.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL))
   {
      if (rowId == colId)
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               SparseSM::Worland::I4Lapl2 i4lapl2(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               bMat = Pm * i4lapl2.mat();
            }

            return bMat;
         };

         // Imaginary part of block
         auto imagOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               const auto dl = static_cast<Internal::MHDFloat>(l);
               const auto invlapl = MHD_MP(1.0) / (dl * (dl + MHD_MP(1.0)));
               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();
               SparseSM::Worland::I4Lapl coriolis(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);

               // Correct Laplacian for 4th order system according to:
               // McFadden,Murray,Boisvert,
               // Elimination of Spurious Eigenvalues in the
               // Chebyshev Tau Spectral Method,
               // JCP 91, 228-239 (1990)
               // We simply drop the last column
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
               bMat = static_cast<MHDFloat>(o.m * T * Pm * invlapl) *
                  coriolis.mat() * qid.mat();
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = imagOp;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::TOR))
      {
         // Create real part of block
         auto realOpLower = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));
               SparseSM::Worland::I4Qm corQm(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               auto norm = coriolis(l, o.m);
               bMat =
                  static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQm.mat();
            }

            return bMat;
         };

         // Create first lower diagonal operator
         auto& dLow = getDescription();
         dLow.nRowShift = 1;
         dLow.nColShift = 0;
         dLow.realOp = realOpLower;
         dLow.imagOp = nullptr;

         // Create real part of block
         auto realOpUpper = [](const int nNr, const int nNc, const int l,
                               std::shared_ptr<details::BlockOptions> opts,
                               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                     opts);

               auto coriolis = [](const int l, const int m)
               {
                  return (l - MHD_MP(1.0)) * (l + MHD_MP(1.0)) *
                         Internal::Math::sqrt(
                            ((l - m) * (l + m)) /
                            ((MHD_MP(2.0) * l - MHD_MP(1.0)) *
                               (MHD_MP(2.0) * l + MHD_MP(1.0))));
               };

               const auto T =
                  1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
               const auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())
                                  ->second->value();

               const auto dl = static_cast<MHDFloat>(l);
               const auto invlapl = 1.0 / (dl * (dl + 1.0));
               SparseSM::Worland::I4Qp corQp(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               auto norm = -coriolis(l + 1, o.m);
               bMat =
                  static_cast<MHDFloat>(norm * T * Pm * invlapl) * corQp.mat();
            }

            return bMat;
         };

         // Create first upper diagonal operator
         auto& dUp = getDescription();
         dUp.nRowShift = 0;
         dUp.nColShift = 1;
         dUp.realOp = realOpUpper;
         dUp.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::TOR))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::POL))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto Pr =
               nds.find(NonDimensional::Prandtl::id())->second->value();
            const auto Pm =
               nds.find(NonDimensional::MagneticPrandtl::id())->second->value();

            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            SparseMatrix bMat = (Pm / Pr) * i2lapl.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else
   {
      throw std::logic_error("Equations are not setup properly");
   }

   return blkDef;
}

details::BlockDefinition ModelBackend::timeBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   details::BlockDefinition blkDef;
   blkDef.rowId = rowId;
   blkDef.colId = colId;
   blkDef.isComplex = this->isComplex(rowId);
   blkDef.isGalerkin = this->useGalerkin();

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      blkDef.descr.push_back({});
      auto& d = blkDef.descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::worland_default_t::ALPHA;
      opts->b = Polynomial::Worland::worland_default_t::DBETA;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I4Lapl spasm(nNr, nNc, o.a, o.b, l,
               2 * o.truncateQI);

            // Correct Laplacian for 4th order system according to:
            // McFadden,Murray,Boisvert,
            // Elimination of Spurious Eigenvalues in the
            // Chebyshev Tau Spectral Method,
            // JCP 91, 228-239 (1990)
            // We simply drop the last column
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
            bMat = spasm.mat() * qid.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l, 1 * o.truncateQI);
         SparseMatrix bMat = spasm.mat();

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return blkDef;
}

details::BlockDefinition ModelBackend::qiBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   details::BlockDefinition blkDef;
   blkDef.rowId = rowId;
   blkDef.colId = colId;
   blkDef.isComplex = this->isComplex(rowId);
   blkDef.isGalerkin = this->useGalerkin();

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      blkDef.descr.push_back({});
      auto& d = blkDef.descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::worland_default_t::ALPHA;
      opts->b = Polynomial::Worland::worland_default_t::DBETA;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat(nNr, nNc);
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat(nNr, nNc);
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I4 spasm(nNr, nNc, o.a, o.b, l,
               2 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (l > 0)
         {
            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
            bMat = qid.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l, 1 * o.truncateQI);
         SparseMatrix bMat = spasm.mat();

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return blkDef;
}

void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix,
   const std::size_t opId,
   const Equations::CouplingInformation::FieldId_range imRange,
   const int matIdx, const std::size_t bcType, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 1);
   int m = eigs.at(0);
   auto maxL = res.counter().dim(Dimensions::Simulation::SIM2D,
                  Dimensions::Space::SPECTRAL, m) -
               1;

   auto getNns = [&](const SpectralFieldId& rowId, const SpectralFieldId& colId, const int j0, const int maxJ)
   {
      // Store 1D sizes
      std::vector<int> ns;
      for (int j = j0; j <= maxJ; j++)
      {
         auto nN = this->baseNn(j, res);
         ns.emplace_back(nN);
      }

      return ns;
   };

   // Time operator
   if (opId == ModelOperator::Time::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         auto colId = rowId;
         const auto& fields = this->implicitFields(rowId);
         auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         auto nNs = getNns(rowId, colId, m, maxL);
         buildBlock(rModelMatrix, descr, fields, matIdx, bcType,
            m, maxL, nNs, bcs, nds, false, -1);
      }
   }
   // Linear operator
   else if (opId == ModelOperator::ImplicitLinear::id())
   {
      bool isSplit = false;

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               implicitBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            auto nNs = getNns(rowId, colId, m, maxL);
            buildBlock(rModelMatrix, descr, fields, matIdx,
               bcType, m, maxL, nNs, bcs, nds, isSplit, -1);
         }
      }
   }
   // Boundary operator
   else if (opId == ModelOperator::Boundary::id())
   {
      bool isSplit = false;

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               boundaryBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            auto nNs = getNns(rowId, colId, m, maxL);
            buildBlock(rModelMatrix, descr, fields, matIdx,
               bcType, m, maxL, nNs, bcs, nds, isSplit, -1);
         }
      }
   }
   // Quasi-Inverse operator
   else if (opId == ModelOperator::QuasiInverse::id())
   {
      BcMap qiBcs;
      for(auto& [k,v]: bcs)
      {
         qiBcs.emplace(k, Bc::Name::QuasiInverseOnly::id());
      }

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         auto colId = rowId;
         const auto& fields = this->implicitFields(rowId);
         auto descr = qiBlockBuilder(rowId, colId, res, eigs, qiBcs, nds);
         auto nNs = getNns(rowId, colId, m, maxL);
         buildBlock(rModelMatrix, descr, fields, matIdx, bcType,
            m, maxL, nNs, qiBcs, nds, false, -1);
      }
   }
   else
   {
      throw std::logic_error("Requested operator type is not implemented");
   }
}

details::BlockDefinition ModelBackend::boundaryBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplit) const
{
   details::BlockDefinition blkDef;
   blkDef.rowId = rowId;
   blkDef.colId = colId;
   blkDef.isComplex = this->isComplex(rowId);
   blkDef.isGalerkin = this->useGalerkin();

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      blkDef.descr.push_back({});
      auto& d = blkDef.descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::worland_default_t::ALPHA;
      opts->b = Polynomial::Worland::worland_default_t::DBETA;
      opts->m = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplit;
      d.opts = opts;

      return d;
   };

   if (rowId == colId)
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int l,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         SparseMatrix bMat(nNr, nNc);

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else
   {
      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = nullptr;
      d.imagOp = nullptr;
   }

   return blkDef;
}

void ModelBackend::galerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(this->useGalerkin());
   assert(eigs.size() == 1);
   int m = eigs.at(0);
   auto maxL = res.counter().dim(Dimensions::Simulation::SIM2D,
                  Dimensions::Space::SPECTRAL, m) -
               1;

   std::vector<int> nNs;
   std::vector<int> ls;
   for(int l = m; l <= maxL; l++)
   {
      ls.emplace_back(l);
      auto nN = this->baseNn(l, res);
      nNs.emplace_back(nN);
   }

   // Compute system size
   const auto& fields = this->implicitFields(fieldId);
   const auto sysRows =
      systemInfo(fieldId, fieldId, fields, nNs, makeSquare, false)
         .blockRows;
   const auto sysCols =
      systemInfo(fieldId, fieldId, fields, nNs, true, false)
         .blockCols;

   auto nL = res.counter().dim(Dimensions::Simulation::SIM2D,
      Dimensions::Space::SPECTRAL, m);

   if (mat.size() == 0)
   {
      mat.resize(sysRows, sysCols);
   }
   assert(mat.rows() == sysRows);
   assert(mat.cols() == sysCols);

   assert(eigs.size() == 1);

   int rowShift = 0;
   int colShift = 0;
   for (int i = 0; i < ls.size(); i++)
   {
      SparseMatrix S;
      this->stencil(S, fieldId, ls.at(i), nNs.at(i), makeSquare, bcs, nds);
      this->addBlock(mat, S, rowShift, colShift);

      rowShift += S.rows();
      colShift += S.cols();
   }
}

void ModelBackend::explicitBlock(DecoupledZSparse& mat,
   const SpectralFieldId& fId, const std::size_t opId,
   const SpectralFieldId fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   // Explicit linear operator
   if (opId == ModelOperator::ExplicitLinear::id())
   {
      // Nothing to be done
      throw std::logic_error("There are no explicit linear operators");
   }
   // Explicit nonlinear operator
   else if (opId == ModelOperator::ExplicitNonlinear::id())
   {
      throw std::logic_error("There are no explicit nonlinear operators");
   }
   // Explicit nextstep operator
   else if (opId == ModelOperator::ExplicitNextstep::id())
   {
      throw std::logic_error("There are no explicit nextstep operators");
   }
}

} // namespace Implicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
