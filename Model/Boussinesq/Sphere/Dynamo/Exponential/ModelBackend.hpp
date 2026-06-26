/**
 * @file ModelBackend.hpp
 * @brief Model backend
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_MODELBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_MODELBACKEND_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Model/Boussinesq/Sphere/Dynamo/Exponential/IExponentialBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

namespace Exponential {

/**
 * @brief Interface for model backend
 */
class ModelBackend : public IExponentialBackend
{
public:
   /**
    * @brief Constructor
    */
   ModelBackend();

   /**
    * @brief Destructor
    */
   virtual ~ModelBackend() = default;

   /**
    * @brief Get equation information
    *
    * @param info Equation information
    * @param fId  Field ID
    * @param res  Resolution object
    */
   virtual void equationInfo(EquationInfo& info, const SpectralFieldId& fId,
      const Resolution& res) const override;

   /**
    * @brief Build model matrix
    *
    * @param rModelMatrix  Input/Output matrix to fill with operators
    * @param opId          Type of model matrix
    * @param imRange       Coupled fields
    * @param matIdx        Matrix index
    * @param bcType        Boundary condition scheme (Tau vs Galerkin)
    * @param res           Resolution object
    * @param eigs          Indexes of other dimensions
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void modelMatrix(DecoupledZSparse& rModelMatrix,
      const std::size_t opId,
      const Equations::CouplingInformation::FieldId_range imRange,
      const int matIdx, const std::size_t bcType, const Resolution& res,
      const std::vector<MHDFloat>& eigs, const BcMap& bcs,
      const NonDimensional::NdMap& nds) const override;

   /**
    * @brief Build galerkin stencil
    *
    * @param mat     Input/Output matrix to fill with stencil
    * @param fId     Field ID
    * @param matIdx  Matrix index
    * @param res     Resolution object
    * @param eigs          Indexes of other dimensions
    * @param makeSquare Truncate stencil to obtain square matrix?
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId,
      const int matIdx, const Resolution& res,
      const std::vector<MHDFloat>& eigs, const bool makeSquare,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

   /**
    * @brief Build explicit block
    *
    * @param rModelMatrix  Input/Output matrix to fill with operators
    * @param fId           Equation field ID
    * @param opId          Type of explicit operator
    * @param fieldId       Coupled field ID
    * @param matIdx        Matrix index
    * @param res           Resolution object
    * @param eigs          Indexes of other dimensions
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId,
      const std::size_t opId, const SpectralFieldId fieldId, const int matIdx,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

protected:
   /**
    * @brief Operators are complex?
    *
    * @param fId  Field ID
    */
   bool isComplex(const SpectralFieldId& fId) const final;

   /**
    * @brief Build model matrix
    *
    * @param fId  Field ID
    */
   SpectralFieldIds implicitFields(const SpectralFieldId& fId) const final;

   /**
    * @brief Get coupled fields
    *
    * @param fId  Field ID
    */
   SpectralFieldIds explicitFields(const SpectralFieldId& fId) const;

   /**
    * @brief Build implicit matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    * @param isSplitOperator  Set operator of split system
    */
   details::BlockDefinition implicitBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const;

   /**
    * @brief Build time matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    */
   details::BlockDefinition timeBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const;

   /**
    * @brief Build explicit linear matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    */
   details::BlockDefinition explicitLinearBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const;

   /**
    * @brief Build quasi-inverse matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    */
   details::BlockDefinition qiBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const;

   /**
    * @brief Build boundary matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    * @param isSplitOperator  Set operator of split system
    */
   details::BlockDefinition boundaryBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const;

   /**
    * @brief Build boundary matrix block description
    *
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param res     Resolution object
    * @param eigs    Slow indexes
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    */
   details::BlockDefinition splitBoundaryValueBlockBuilder(
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const;

private:
   /**
    * @brief Truncate quasi-inverse operators?
    */
   const bool mcTruncateQI;
};

} // namespace Exponential
} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_EXPONENTIAL_MODELBACKEND_HPP
