/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "FiniteDiff/Sphere/UniformRadialGrid.hpp"
#include "Model/Boussinesq/Sphere/Dynamo/FdExplicit/ModelBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/Dimensions.hpp"
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
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "FiniteDiff/Sphere/Id.hpp"
#include "FiniteDiff/Sphere/SLapl.hpp"
#include "FiniteDiff/Sphere/Boundary/Value.hpp"
#include "FiniteDiff/Sphere/Boundary/D1.hpp"
#include "FiniteDiff/Sphere/Boundary/D2.hpp"
#include "FiniteDiff/Sphere/Boundary/R1D1DivR1.hpp"
#include "FiniteDiff/Sphere/Boundary/InsulatingSphere.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace FdExplicit {

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

   /// Harmonic degree l
   int l;
   /// Boundary condition
   std::size_t bcId;
   /// Split operator for influence matrix?
   bool isSplitOperator;
   /// Use split equation for influence matrix?
   bool useSplitEquation;
   /// Radial finite difference grid
   Internal::Array igrid;
};
} // namespace implDetails

ModelBackend::ModelBackend() :
    IDynamoBackend()
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
      static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
}

void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId,
   const Resolution& res, const Equations::Tools::ICoupling& coupling,
   const BcMap& bcs) const
{
   // Loop overall matrices/eigs
   for (int idx = 0; idx < info.tauN.size(); ++idx)
   {
      auto eigs = coupling.getIndexes(res, idx);

      int tN, gN, rhs;
      ArrayI shift(3);

      this->blockInfo(tN, gN, shift, rhs, fId, res, eigs.at(0), bcs);

      info.tauN(idx) = tN;
      info.galN(idx) = gN;
      info.galShift.row(idx) = shift;
      info.rhsCols(idx) = rhs;

      // Compute system size
      int sN = 0;
      for (auto f: this->implicitFields(fId))
      {
         this->blockInfo(tN, gN, shift, rhs, f, res, eigs.at(0), bcs);
         sN += gN;
      }

      if (sN == 0)
      {
         sN = info.galN(idx);
      }

      info.sysN(idx) = sN;
   }
}

std::vector<details::BlockDescription> ModelBackend::implicitBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->igrid = this->getGrid(res.sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->isSplitOperator = isSplitOperator;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR) &&
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
            FiniteDiff::Sphere::SLapl opLapl;
            Internal::SparseMatrix lapl; 
            opLapl.compute(lapl, l, o.igrid);
            FiniteDiff::Sphere::Id opId(1,1);
            Internal::SparseMatrix qid;
            opId.compute(qid, l ,o.igrid);
            bMat = Pm * qid * lapl;
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
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL) &&
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
               FiniteDiff::Sphere::SLapl opLapl;
               Internal::SparseMatrix lapl; 
               opLapl.compute(lapl, l, o.igrid);

               FiniteDiff::Sphere::Id opId(1,1);
               Internal::SparseMatrix qid;
               opId.compute(qid, l ,o.igrid);

               if (o.isSplitOperator)
               {
                  bMat = qid * lapl;
               }
               else
               {
                  bMat = Pm * qid * lapl;
               }
            }
            else
            {
               throw std::logic_error("Finite difference bilaplacian is not implemented");
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
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::TOR) &&
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
            FiniteDiff::Sphere::SLapl opLapl;
            Internal::SparseMatrix lapl; 
            opLapl.compute(lapl, l, o.igrid);

            FiniteDiff::Sphere::Id opId(1,1);
            Internal::SparseMatrix qid;
            opId.compute(qid, l ,o.igrid);

            bMat = qid * lapl;
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
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::POL) &&
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
            FiniteDiff::Sphere::SLapl opLapl;
            Internal::SparseMatrix lapl; 
            opLapl.compute(lapl, l, o.igrid);

            FiniteDiff::Sphere::Id opId(1,1);
            Internal::SparseMatrix qid;
            opId.compute(qid, l ,o.igrid);

            bMat = qid * lapl;
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
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR) &&
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

         FiniteDiff::Sphere::SLapl opLapl;
         Internal::SparseMatrix lapl; 
         opLapl.compute(lapl, l, o.igrid);

         FiniteDiff::Sphere::Id opId(1,1);
         Internal::SparseMatrix qid;
         opId.compute(qid, l ,o.igrid);

         SparseMatrix bMat = (Pm / Pr) * qid * lapl;

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

   return descr;
}

std::vector<details::BlockDescription> ModelBackend::timeBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->igrid = this->getGrid(res.sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
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
            FiniteDiff::Sphere::Id opId(1,1);
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
         }
         else
         {
            FiniteDiff::Sphere::Id opId;
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
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
            if (o.useSplitEquation)
            {
               FiniteDiff::Sphere::Id opId(1, 1);
               Internal::SparseMatrix mat; 
               opId.compute(mat, l, o.igrid);
               bMat = mat;
            }
            else
            {
               throw std::logic_error("Finite differences bilaplacian is not implemented!");
            }
         }
         else
         {
            FiniteDiff::Sphere::Id opId;
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
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
            FiniteDiff::Sphere::Id opId(1, 1);
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
         }
         else
         {
            FiniteDiff::Sphere::Id opId;
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
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
            FiniteDiff::Sphere::Id opId(1, 1);
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
         }
         else
         {
            FiniteDiff::Sphere::Id opId;
            Internal::SparseMatrix mat; 
            opId.compute(mat, l, o.igrid);
            bMat = mat;
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

         FiniteDiff::Sphere::Id opId(1,1);
         Internal::SparseMatrix mat; 
         opId.compute(mat, l, o.igrid);
         SparseMatrix bMat = mat;

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return descr;
}

std::vector<details::BlockDescription> ModelBackend::boundaryBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplit) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
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

   return descr;
}

std::vector<details::BlockDescription>
ModelBackend::splitBoundaryValueBlockBuilder(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->l = eigs.at(0);
      opts->bcId = bcs.find(colId.first)->second;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::POL))
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

   return descr;
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

   // Time operator
   if (opId == ModelOperator::Time::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         auto colId = rowId;
         const auto& fields = this->implicitFields(rowId);
         auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx, bcType,
            res, l, l, bcs, nds, false);
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
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, l, l, bcs, nds, isSplit);
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
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, l, l, bcs, nds, isSplit);
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
            buildFixedBlock(rModelMatrix, 1, true, descr, rowId, colId, fields,
               matIdx, bcType, res, l, l, bcs, nds, false);
         }
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
   this->stencil(mat, fieldId, l, res, makeSquare, bcs, nds);
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

const Internal::Array& ModelBackend::getGrid(const int size) const
{
   if(!this->mpGrid)
   {
      this->mpGrid = std::make_shared<Internal::Array>();
      FiniteDiff::Sphere::UniformRadialGrid quad;
      quad.computeGrid(*this->mpGrid, size);
   }

   assert(size == this->mpGrid->size());
   return *this->mpGrid;
}

void ModelBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l,
   std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto bcId = bcs.find(rowId.first)->second;

   SparseMatrix bcOp(nN, nN);
   SparseMatrix bcMat(nN, nN);

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR) &&
       rowId == colId)
   {
      if (l > 0)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            FiniteDiff::Sphere::Boundary::Value valOp;
            valOp.compute(bcMat, 0, l, this->getGrid(nN));
            bcOp += bcMat;

            valOp.compute(bcMat, nN-1, l, this->getGrid(nN));
            bcOp += bcMat;
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            FiniteDiff::Sphere::Boundary::Value valOp;
            valOp.compute(bcMat, 0, l, this->getGrid(nN));
            bcOp += bcMat;

            FiniteDiff::Sphere::Boundary::R1D1DivR1 sfOp;
            sfOp.compute(bcMat, nN-1, l, this->getGrid(nN));
            bcOp += bcMat;
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
               FiniteDiff::Sphere::Boundary::Value valOp;
               valOp.compute(bcMat, 0, l, this->getGrid(nN));
               bcOp += bcMat;

               valOp.compute(bcMat, nN-1, l, this->getGrid(nN));
               bcOp += bcMat;

               SparseMatrix flip(bcOp.rows(), bcOp.cols());
               std::vector<Eigen::Triplet<MHDFloat>> triplets;
               triplets.emplace_back(0, bcOp.cols()-1, 1.0);
               triplets.emplace_back(bcOp.rows()-1, 0, 1.0);
               for(int i = 1; i < bcOp.rows()-1; i++)
               {
                  triplets.emplace_back(i, i, 1.0);
               }
               flip.setFromTriplets(triplets.begin(), triplets.end());
               bcOp = flip*bcOp;

            }
            else if (bcId == Bc::Name::NoSlip::id())
            {
               FiniteDiff::Sphere::Boundary::Value valOp;
               valOp.compute(bcMat, 0, l, this->getGrid(nN));
               bcOp += bcMat;

               FiniteDiff::Sphere::Boundary::D1 d1Op;
               d1Op.compute(bcMat, nN-1, l, this->getGrid(nN));
               bcOp += bcMat;
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               FiniteDiff::Sphere::Boundary::Value valOp;
               valOp.compute(bcMat, 0, l, this->getGrid(nN));
               bcOp += bcMat;

               FiniteDiff::Sphere::Boundary::D2 d2Op;
               d2Op.compute(bcMat, nN-1, l, this->getGrid(nN));
               bcOp += bcMat;
            }
            else
            {
               throw std::logic_error("Boundary conditions for Velocity "
                                      "Poloidal component not implemented");
            }
         }
         else
         {
            throw std::logic_error("4th order system is not implemented for FD scheme");
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
            FiniteDiff::Sphere::Boundary::Value valOp;
            valOp.compute(bcMat, 0, l, this->getGrid(nN));
            bcOp += bcMat;

            valOp.compute(bcMat, nN-1, l, this->getGrid(nN));
            bcOp += bcMat;
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
            FiniteDiff::Sphere::Boundary::Value valOp;
            valOp.compute(bcMat, 0, l, this->getGrid(nN));
            bcOp += bcMat;

            FiniteDiff::Sphere::Boundary::InsulatingSphere magOp;
            magOp.compute(bcMat, nN-1, l, this->getGrid(nN));
            bcOp += bcMat;
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
         FiniteDiff::Sphere::Boundary::Value valOp;
         valOp.compute(bcMat, 0, l, this->getGrid(nN));
         bcOp += bcMat;

         valOp.compute(bcMat, nN-1, l, this->getGrid(nN));
         bcOp += bcMat;
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         FiniteDiff::Sphere::Boundary::Value valOp;
         valOp.compute(bcMat, 0, l, this->getGrid(nN));
         bcOp += bcMat;

         FiniteDiff::Sphere::Boundary::D1 d1Op;
         d1Op.compute(bcMat, nN-1, l, this->getGrid(nN));
         bcOp += bcMat;
      }
      else
      {
         throw std::logic_error(
            "Boundary conditions for Temperature not implemented (" +
            std::to_string(bcId) + ")");
      }
   }

   mat.real() += bcOp;
}

void ModelBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, std::shared_ptr<details::BlockOptions> opts,
   const Resolution& res, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   throw std::logic_error("Galerkin stencil is not implemented for FD scheme");
}

} // namespace FdExplicit
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
