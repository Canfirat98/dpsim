// Attribute und Funktionen in einem anderen Branch

#pragma once

/*
#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <bitset>

#include <DPsim.h>


#include <dpsim/Config.h>
#include <dpsim/Solver.h>
#include <dpsim/SolverParameters.h>
#include <dpsim/DataLogger.h>
#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Solver/MNASwitchInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <dpsim-models/SimSignalComp.h>
#include <dpsim-models/SimPowerComp.h>
*/

#include <dpsim/SolverParameters.h>

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
#include <sundials/sundials_types.h>


/* std::size_t is the largest data type. No container can store
 * more than std::size_t elements. Define the number of switches
 * as the log_2 of this value so that we end up with maximally
 * std::size_t matrices. The overhead of statically defining this
 * value should be minimal.
 **/
#define SWITCH_NUM sizeof(std::size_t)*8

namespace DPsim {
	/// Parameter class for Differential Algebraic Equation(DAE) solver systems
	class SolverParametersDAE : public SolverParameters {
	public: 
		SolverParametersDAE() { }
		
		/// Destructor
		virtual ~SolverParametersDAE() {};

		// Solver parameters
		/// Relative tolerance
		realtype mRelativeTolerance = 1e-4;
		///
		bool mUseUserSuppJacMatrix = false;

		// ### nonlinear solver interface optional input parameters ###
		/// specifies maximum no. of nonlinear iterations
		int mIDAMaxNonlinIters = 4;
		/// specifies the maximum number of nonlinear solver convergence failures in one step
		int mIDAMaxConvFails = 10;
		/// Coeff. in the nonlinear convergence test
		realtype mIDANonlinConvCoef = 0.33;

		// ### Solver optional input parameters ###
		/// Maximum order for BDF method (0<mIDAMaxBDFOrder<=5)
		unsigned int mIDAMaxBDFOrder = 5;
		/// Maximum number of steps to be taken by the solver in its attempt to reach the next output time (-1 = unlimited)
		int mIDAMaxNumSteps = -1;
		/// indicates whether or not to use variable integration step size
		bool mVariableStepSize = true;
		/// min integration time
		realtype mMinStepSize = 1e-6;
		/// max integration time
		realtype mMaxStepSize = 1e-3;
		/// Maximum number of error test failures in attempting one step
		unsigned int mIDAMaxErrTestFails = 100;
		/// indicates whether or not to suppress algebraic variables in the local error test
		bool mIDASetSuppressAlg = SUNTRUE;

		// Variables to store ida statistics
        /// number of steps taken by ida
		long int mNumberStepsIDA = 0;
		/// number of calls to the user's res function
        long int mNumberCallsResidualFunctions = 0;
		/// cumulative number of nonlinear iterations performed
		long int mNonLinearIters;
		/// cumulative number of calls made to the linear solver’s setup function (total so far)
		long int mNumLinSolvSetups = 0;
		/// cumulative number of local error test failures that have occurred (total so far)
		long int mNumErrTestFails = 0;
		/// number of failed steps due to a nonlinear solver failure
		long int mNumStepSolveFails = 0;
		/// integration step size taken on the last internal step
		realtype mLastIntegrationStep;
		/// Step size to be attempted on the next step
		realtype mNextIntegrationStep;
		/// Actual initial step size
		realtype mActualInitialStepSize;
		/// Current internal time reached by the solver
		realtype mCurrentInternalTime;
		/// Order used during the last step
		int mOrderLastStep;
		/// Order to be attempted on the next step
		int mOrderNextStep;
		/// vector of solution error weights at the current time
		N_Vector mErrorWeights;
		/// vector of estimated load errors at the current time
		N_Vector mEstLocalErrors;


		/// ### Setters ###
		/// Set relative Tolerance
		void setRelativeTolerance(Real relTol) { mRelativeTolerance = relTol; }
		/// use user supplied jacobian matrix
		void setUseUserSuppliedJacMatrix(bool useUserSuppJacMatrix) { mUseUserSuppJacMatrix = useUserSuppJacMatrix;}

		// ### Setters for nonlinear solver interface optional input functions ###
		/// Set maximum no. of convergence failures
		void setIDAMaxConvFails(int IDAMaxConvFails) { mIDAMaxConvFails = IDAMaxConvFails; }
		/// Set the safety factor in the nonlinear convergence test
		void setIDANonlinConvCoef(Real IDANonlinConvCoef) { mIDANonlinConvCoef = IDANonlinConvCoef; }
		/// Set the maximum number of nonlinear solver iterations in one solve attempt
		void setIDAMaxNonlinIters(int IDAMaxNonlinIters) { mIDAMaxNonlinIters = IDAMaxNonlinIters; }

		// ### Setters for main solver optional input functions ###
		/// Set maximum order for BDF method
		void setIDAMaxBDFOrder(unsigned int IDAMaxBDFOrder) {
			if (IDAMaxBDFOrder>5 || IDAMaxBDFOrder==0)
				IDAMaxBDFOrder = 5;
			mIDAMaxBDFOrder = IDAMaxBDFOrder;
		}
		/// Set maximum no. of internal steps before tout
		void setIDAMaxNumSteps(int IDAMaxNumSteps) { mIDAMaxNumSteps = IDAMaxNumSteps; }
		/// Set flag to specifi whether or not to use variable integration step size
		void setVariableStepSize(bool variableStepSize) { mVariableStepSize = variableStepSize; }
		/// Set min integration time
		void setMinStepSize(Real minStepSize) { mMinStepSize = minStepSize; }
		/// Set max integration time 
		void setMaxStepSize(Real maxStepSize) { mMaxStepSize = maxStepSize; }
		/// Set max integration time
		void setMaxErrTestFails(unsigned int IDAMaxErrTestFails) { mIDAMaxErrTestFails = IDAMaxErrTestFails; }
		/// Set Suppress alg. vars. from error test
		void setIDASetSuppressAlg(bool IDASetSuppressAlg) { mIDASetSuppressAlg = IDASetSuppressAlg; }

    };
}
