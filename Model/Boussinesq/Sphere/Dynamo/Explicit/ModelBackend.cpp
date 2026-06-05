/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Explicit/ModelBackend.hpp"
#include "QuICC/Bc/Name/QuasiInverseOnly.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitQuasiInverse.hpp"
#include "QuICC/ModelOperator/QuasiInverse.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Explicit {

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
   /// Harmonic degree l
   int l;
   /// Use truncated quasi-inverse?
   bool truncateQI;
   /// Boundary condition
   std::size_t bcId;
   /// Split operator for influence matrix?
   bool isSplitOperator;
   /// Use split equation for influence matrix?
   bool useSplitEquation;
};
} // namespace implDetails

namespace {
   const auto mag_tor = std::make_pair(PhysicalNames::Magnetic::id(),
                   FieldComponents::Spectral::TOR);
   const auto mag_pol = std::make_pair(PhysicalNames::Magnetic::id(),
                   FieldComponents::Spectral::POL);
   const auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   const auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   const auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);
}

ModelBackend::ModelBackend() :
    IDynamoBackend(),
    mcTruncateQI(true)
{}

bool ModelBackend::isComplex(const SpectralFieldId& fId) const
{
   return false;
}

ModelBackend::SpectralFieldIds ModelBackend::implicitFields(
   const SpectralFieldId& fId) const
{
   SpectralFieldIds fields = {fId};

   return fields;
}

void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId,
   const Resolution& res) const
{
   // Operators are real
   info.isComplex = this->isComplex(fId);

   // Splitting 4th poloidal equation into two systems
   if (fId == vel_pol)
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
      static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
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
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplitOperator;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (rowId == vel_tor &&
       rowId == colId)
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
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            const auto Pm =
               nds.find(NonDimensional::MagneticPrandtl::id())->second->value();
            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = Pm * i2lapl.mat();
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
   else if (rowId == vel_pol &&
            rowId == colId)
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
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto Pm =
               nds.find(NonDimensional::MagneticPrandtl::id())->second->value();

            if (o.useSplitEquation)
            {
               if (o.isSplitOperator)
               {
                  SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
                     1 * o.truncateQI);
                  bMat = i2lapl.mat();
               }
               else
               {
                  SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
                     1 * o.truncateQI);
                  bMat = Pm * i2lapl.mat();
               }
            }
            else
            {
               SparseSM::Worland::I4Lapl2 i4lapl2(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               bMat = Pm * i4lapl2.mat();
            }
         }

         return bMat;
      };

      // Create diagonal block
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (rowId == mag_tor &&
            rowId == colId)
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
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = i2lapl.mat();
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
   else if (rowId == mag_pol &&
            rowId == colId)
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
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l,
               1 * o.truncateQI);
            bMat = i2lapl.mat();
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
   else if (rowId == temp &&
            rowId == colId)
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
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == vel_tor)
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
   else if (fieldId == vel_pol)
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
            if (o.useSplitEquation)
            {
               SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               bMat = spasm.mat();
            }
            else
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
   else if (fieldId == mag_tor)
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
   else if (fieldId == mag_pol)
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
   else if (fieldId == temp)
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
   const NonDimensional::NdMap& nds, const bool isSplitOperator) const
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
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplitOperator;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == vel_tor)
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
   else if (fieldId == vel_pol)
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
            if (o.useSplitEquation)
            {
               SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l,
                  1 * o.truncateQI);
               if(o.isSplitOperator)
               {
                  SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -2);
                  bMat = spasm.mat() * qid.mat();
               }
               else
               {
                  bMat = spasm.mat();
               }
            }
            else
            {
               SparseSM::Worland::I4 spasm(nNr, nNc, o.a, o.b, l,
                  2 * o.truncateQI);
               bMat = spasm.mat();
            }
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
   else if (fieldId == mag_tor)
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
   else if (fieldId == mag_pol)
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
   else if (fieldId == temp)
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
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplit;
      opts->useSplitEquation = this->useSplitEquation();
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

   return blkDef;
}

details::BlockDefinition
ModelBackend::splitBoundaryValueBlockBuilder(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   details::BlockDefinition blkDef;
   blkDef.rowId = rowId;
   blkDef.colId = colId;
   blkDef.isComplex = true;
   blkDef.isGalerkin = this->useGalerkin();

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      blkDef.descr.push_back({});
      auto& d = blkDef.descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->a = Polynomial::Worland::worland_default_t::ALPHA;
      opts->b = Polynomial::Worland::worland_default_t::DBETA;
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == vel_pol)
   {
      // Boundary value operator
      auto bcValOp = [](const int nNr, const int nNc, const int l,
                        std::shared_ptr<details::BlockOptions> opts,
                        const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat(nNr, 1);

         Eigen::Triplet<MHDFloat> val = {0, 0, 1.0};
         std::vector<Eigen::Triplet<MHDFloat>> triplets = {val};
         bMat.setFromTriplets(triplets.begin(), triplets.end());

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = bcValOp;
      d.imagOp = bcValOp;
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
   int l = eigs.at(0);

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
         auto nNs = getNns(rowId, colId, l, l);
         buildBlock(rModelMatrix, descr, fields, matIdx, bcType,
            l, l, nNs, bcs, nds, false, -1);
      }
   }
   // Linear operator
   else if (opId == ModelOperator::ImplicitLinear::id() ||
            opId == ModelOperator::SplitImplicitLinear::id())
   {
      bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               implicitBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            auto nNs = getNns(rowId, colId, l, l);
            buildBlock(rModelMatrix, descr, fields, matIdx,
               bcType, l, l, nNs, bcs, nds, isSplit, -1);
         }
      }
   }
   // Boundary operator
   else if (opId == ModelOperator::Boundary::id() ||
            opId == ModelOperator::SplitBoundary::id())
   {
      bool isSplit = (opId == ModelOperator::SplitBoundary::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               boundaryBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            auto nNs = getNns(rowId, colId, l, l);
            buildBlock(rModelMatrix, descr, fields, matIdx,
               bcType, l, l, nNs, bcs, nds, isSplit, -1);
         }
      }
   }
   // Split equation boundary value
   else if (opId == ModelOperator::SplitBoundaryValue::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr = splitBoundaryValueBlockBuilder(rowId, colId, res, eigs,
               bcs, nds);
            auto nNs = getNns(rowId, colId, l, l);
            buildBlock(rModelMatrix, descr, fields,
               matIdx, bcType, l, l, nNs, bcs, nds, false, 1);
         }
      }
   }
   // Quasi-Inverse operator
   else if (opId == ModelOperator::QuasiInverse::id() ||
            opId == ModelOperator::SplitQuasiInverse::id())
   {
      bool isSplit = (opId == ModelOperator::SplitQuasiInverse::id());

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
         auto descr = qiBlockBuilder(rowId, colId, res, eigs, qiBcs, nds, isSplit);
         auto nNs = getNns(rowId, colId, l, l);
         buildBlock(rModelMatrix, descr, fields, matIdx, bcType,
            l, l, nNs, qiBcs, nds, false, -1);
      }
   }
   else
   {
      throw std::logic_error("Requested operator type is not implemented");
   }
}

void ModelBackend::galerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 1);
   int l = eigs.at(0);
   auto nN = this->baseNn(l, res);
   this->stencil(mat, fieldId, l, nN, makeSquare, bcs, nds);
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

} // namespace Explicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
