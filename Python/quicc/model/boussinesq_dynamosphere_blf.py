"""Module provides the functions to generate the Boussinesq rotating thermal dynamo in a sphere with Chebyshev expansion (Toroidal/Poloidal formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.sphere_chebyshev as geo
import quicc.base.base_model as base_model
from quicc.geometry.spherical.sphere_boundary_chebyshev import no_bc


class BoussinesqDynamoSphere(base_model.BaseModel):
    """Class to setup the Boussinesq rotating thermal dynamo in a sphere with Chebyshev expansion (Toroidal/Poloidal formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["magnetic_prandtl", "taylor", "prandtl", "rayleigh"]

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature", "magnetic"]

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

#        # Coupled solve for velocity field and temperature
#        if field_row in [("velocity","tor"), ("velocity","pol"), ("temperature","")]:
#           fields = [("velocity","tor"), ("velocity","pol"), ("temperature","")]
#        # Independent solve for magnetic field
#        else:
#           fields = [field_row]

        # Large coupled solve
        fields =  [("velocity","tor"), ("velocity","pol"), ("temperature",""), ("magnetic","tor"), ("magnetic","pol")]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit linear dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            if field_row == ("temperature",""):
                fields = [("velocity","pol")]
            else:
                fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row == ("temperature",""):
                fields = [("temperature","")]
            else:
                fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, eigs, bcs, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row in [("velocity","tor"), ("magnetic","tor"), ("magnetic","pol"), ("temperature","")]:
                shift_r = 1
            elif field_row == ("velocity","pol"):
                shift_r = 2
            else:
                shift_r = 0

            gal_n = (res[0] - shift_r)

        else:
            gal_n = tau_n
            shift_r = 0

        block_info = (tau_n, gal_n, (shift_r,0,0), 1)
        return block_info

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""
        
        assert(eigs[0].is_integer())

        m = int(eigs[0])

        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return geo.stencil(res[0], res[1], m, bc, make_square)

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_SINGLE_RHS

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-10, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("magnetic","tor"):
                        bc = {0:-10, 'rt':0}
                    elif field_col == ("magnetic","pol"):
                        bc = {0:-13, 'rt':0, 'c':{'l':l}}
                    elif field_col == ("temperature",""):
                        bc = {0:-10, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == field_row:
                            bc = {0:10}
                    elif field_row == ("velocity","pol") and field_col == field_row:
                            bc = {0:20}
                    elif field_row == ("magnetic","tor") and field_col == field_row:
                        bc = {0:10}
                    elif field_row == ("magnetic","pol") and field_col == field_row:
                        bc = {0:13, 'c':{'l':l}}
                    elif field_row == ("temperature","") and field_col == field_row:
                            bc = {0:10}

            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-12, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-21, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:12}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:21}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 1
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 2
                elif field_row == ("magnetic","tor"):
                    bc['rt'] = 1
                elif field_row == ("magnetic","pol"):
                    bc['rt'] = 1
                elif field_row == ("temperature",""):
                    bc['rt'] = 1

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {0:-10, 'rt':1}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-20, 'rt':2}
                    elif field_col == ("magnetic","tor"):
                        bc = {0:-10, 'rt':1}
                    elif field_col == ("magnetic","pol"):
                        bc = {0:-13, 'c':{'l':l}, 'rt':1}
                    elif field_col == ("temperature",""):
                        bc = {0:-10, 'rt':1}

                elif bcId == 1:
                    if field_col == ("velocity","tor"):
                        bc = {0:-12, 'rt':1}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-21, 'rt':2}
        
        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 1
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 2
                elif field_row == ("magnetic","tor"):
                    bc['rt'] = 1
                elif field_row == ("magnetic","pol"):
                    bc['rt'] = 1
                elif field_row == ("temperature",""):
                    bc['rt'] = 1

        else:
            bc = no_bc()

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == ("velocity","pol"):
            mat = geo.i2r2(res[0], res[1], m, bc, -1.0, with_sh_coeff = 'laplh', restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == field_row:
            mat = geo.i2r2(res[0], res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        assert(eigs[0].is_integer())

        Pr = eq_params['prandtl']
        Pm = eq_params['magnetic_prandtl']
        Ra = eq_params['rayleigh']
        T = eq_params['taylor']**0.5

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = geo.i2r2lapl(res[0], res[1], m, bc, Pm, with_sh_coeff = 'laplh', l_zero_fix = 'zero', restriction = restriction)
                bc[0] = min(bc[0], 0)
                mat = mat + geo.i2r2(res[0], res[1], m, bc, 1j*m*T*Pm, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i2r2coriolis(res[0], res[1], m, bc, -T*Pm, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("magnetic","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[1], m, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = geo.i4r4coriolis(res[0], res[1], m, bc, T*Pm, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("velocity","pol"):
                mat = geo.i4r4lapl2(res[0], res[1], m, bc, Pm, with_sh_coeff = 'laplh', l_zero_fix = 'zero', restriction = restriction)
                bc[0] = min(bc[0], 0)
                mat = mat + geo.i4r4lapl(res[0], res[1], m, bc, 1j*m*T*Pm, l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("magnetic","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = geo.i4r4(res[0], res[1], m, bc, -Pm**2*Ra*T/Pr, with_sh_coeff = 'laplh', l_zero_fix = 'zero', restriction = restriction)

        elif field_row == ("magnetic","tor"):
            if field_col == ("velocity","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("velocity","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","tor"):
                mat = geo.i2r2lapl(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("magnetic","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[1], m, bc)

        elif field_row == ("magnetic","pol") and field_col == field_row:
            if field_col == ("velocity","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("velocity","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","pol"):
                mat = geo.i2r2lapl(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'zero', restriction = restriction)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[1], m, bc)

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("velocity","pol"):
                if self.linearize:
                    mat = geo.i2r2(res[0], res[1], m, bc, with_sh_coeff = 'laplh', restriction = restriction)

                else:
                    mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","tor"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("magnetic","pol"):
                mat = geo.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = geo.i2r2lapl(res[0], res[1], m, bc, Pm/Pr, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        assert(eigs[0].is_integer())

        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = geo.i2r2(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("velocity","pol"):
            mat = geo.i4r4lapl(res[0], res[1], m, bc, with_sh_coeff = 'laplh', l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("magnetic","tor"):
            mat = geo.i2r2(res[0], res[1], m, bc, l*(l+1.0), l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("magnetic","pol"):
            mat = geo.i2r2(res[0], res[1], m, bc, l*(l+1.0), l_zero_fix = 'set', restriction = restriction)

        elif field_row == ("temperature",""):
            mat = geo.i2r2(res[0], res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        assert(eigs[0].is_integer())
        m = int(eigs[0])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row in [("velocity","tor"), ("velocity","tor")] and field_row == field_col:
            mat = geo.zblk(res[0], res[1], m, bc, l_zero_fix = 'zero', restriction = restriction)
        else:
            mat = geo.zblk(res[0], res[1], m, bc, restriction = restriction)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat
