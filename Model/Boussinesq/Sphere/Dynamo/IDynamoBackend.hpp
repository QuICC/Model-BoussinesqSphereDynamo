/**
 * @file IDynamoBackend.hpp
 * @brief Base model backend for dynamo model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project includes
//
#include "QuICC/Model/ISphericalModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace Dynamo {

/**
 * @brief Base model backend for Dynamo model
 */
class IDynamoBackend : public ISphericalModelBackend
{
public:
   /**
    * @brief Constructor
    */
   IDynamoBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~IDynamoBackend() = default;

   /**
    * @brief Get vector of names for the physical fields
    */
   virtual std::vector<std::string> fieldNames() const override;

   /**
    * @brief Get vector of names for the nondimensional parameters
    */
   virtual std::vector<std::string> paramNames() const override;

   /**
    * @brief Get vector of bools about periodic box
    */
   virtual std::vector<bool> isPeriodicBox() const override;

   /**
    * @brief Get automatically computed parameters based on input parameters
    *
    * @param cfg  Input parameters
    */
   virtual std::map<std::string, MHDFloat> automaticParameters(
      const std::map<std::string, MHDFloat>& cfg) const override;

protected:
   /**
    * @brief Number of boundary conditions
    *
    * @fId  Field ID for which to get number of BC
    */
   int nBc(const SpectralFieldId& fId) const override;

   /**
    * @brief Apply tau line for boundary condition
    *
    * @param mat     Input/Output matrix to apply tau line to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param l       Harmonic degree
    * @param opts    Options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    * @param isSplitOperator  Is second operator of split 4th order system?
    */
   void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int l,
      std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const override;

   /**
    * @brief Boundary condition stencil
    *
    * @param mat        Input/Output matrix to store galerkin stencil
    * @param fID        Field ID
    * @param l          Harmonic degree
    * @param res        Resolution object
    * @param makeSquare Truncate operator to make square
    * @param bcs        Boundary conditions
    * @param nds        Nondimensional parameters
    */
   void stencil(SparseMatrix& mat, const SpectralFieldId& fId, const int l,
      const Resolution& res, const bool makeSquare, const BcMap& bcs,
      const NonDimensional::NdMap& nds) const;

   /**
    * @brief Apply galerkin stencil for boundary condition
    *
    * @param mat     Input/Output matrix to apply stencil to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param lr      Row space harmonic degree
    * @param lc      Column space harmonic degree
    * @param opts    Options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    */
   void applyGalerkinStencil(SparseMatrix& decMat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int lr, const int lc,
      std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

private:
};

} // namespace Dynamo
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_DYNAMO_IDYNAMOBACKEND_HPP
