/* Copyright 2017-2022 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <list>

#include <dpsim/Solver.h>
#include <dpsim/SolverParametersDAE.h>
#include <dpsim/DataLogger.h>
#include <dpsim-models/Solver/DAEInterface.h>
#include <dpsim/Scheduler.h>
#include <dpsim-models/SimPowerComp.h>
#include <dpsim-models/Logger.h>
#include <dpsim-models/AttributeList.h>


#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
// #include <sunlinsol/sunlinsol_superlumt.h> /* access to SuperLUMT linear solver */
#include <sundials/sundials_types.h>

namespace DPsim {

	/// Solver class which uses Differential Algebraic Equation(DAE) systems
	template <typename VarType>
	class DAESolver : public Solver {
	private:
		/// ### Solver Parameters
		std::shared_ptr<SolverParametersDAE> mSolverParams;

		// General simulation parameters/variables
		///
        CPS::SystemTopology mSystem;
		/// Offsets vector for adding new equations to the residual vector
		std::vector<Int> mOffsets;
		/// Initial sim time
		Real mInitTime;
		/// Current simulation time
		Real mSimTime;
		/// time intervat at which a computed solution is desired
		Real mTimestep;
		/// Number of equations in problem
		Int mNEQ;
		/// Components of the Problem
		CPS::DAEInterface::List mDAEComponents;
		/// Nodes of the Problem
		typename CPS::SimNode<VarType>::List mNodes;

		// IDA simulator variables
		/// Reference to a common simulation context object
		SUNContext mSunctx {nullptr};
		/// Memory block allocated by IDA
		void *mIDAMemoryBlock = NULL;
		/// Vector of problem variables
		N_Vector mStateVector = NULL;
		/// Derivates of the state vector with respect to time
		N_Vector mDerivativeStateVector = NULL;
		/// Vector to indicate differential(0.0) and algebraic(1.0) variable
		N_Vector mStateIDsVector = NULL;
		/// Time IDA reached while solving
		realtype mTimeReachedSolver;
		/// vector of absolute tolerances
		N_Vector mAbsoluteTolerances;
		/// Template Jacobian Matrix
		SUNMatrix mJacobianMatrix = NULL;
		/// Linear solver object
		SUNLinearSolver mLinearSolver = NULL;
		///
        std::vector<CPS::DAEInterface::ResFn> mResidualFunctions;
		///
        std::vector<CPS::DAEInterface::JacobianFn> mJacobianFunctions;

		
		

		/// Residual Function of entire System
		static int residualFunctionWrapper(realtype current_time, N_Vector state, 
									N_Vector dstate_dt, N_Vector resid, void *user_data);
		int residualFunction(realtype current_time, N_Vector state, N_Vector dstate_dt, 
							 N_Vector resid, void *user_data);

		/// Jacobian matrix calculation
		static int jacobianFunctionWrapper(realtype current_time, realtype cj, N_Vector state, 
								N_Vector dstate_dt, N_Vector res_vec, SUNMatrix JJ, void *user_data,
    							N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);
		int calculateJacobianMatrix(realtype current_time, realtype cj, N_Vector state, 
								N_Vector dstate_dt, SUNMatrix JJ);

		/// Initialization of individual components
		void initializeComponents();

		// print and get errors
		/// log weighted errors
		void printEstimatedAndWeightedErrors(std::string& ret);

	public:
		/// Create solve object with given parameters
		DAESolver(const String& name, std::shared_ptr<SolverParametersDAE> solverParams, 
			CPS::SystemTopology system, Real dt, Real t0, 
			CPS::Logger::Level logLevel = CPS::Logger::Level::info);
		/// Deallocate all memory
		~DAESolver();
		
		/// Initialization 
		void initialize();
		/// Solve system for the current time
		Real step(Real time);
		///
		void updateVoltageAndCurrents();
		/// log IDA statistics
		void logStatistics(CPS::Logger::Level minLogLevel = CPS::Logger::Level::debug);
		

	


		// ### Solver tasks ###
		CPS::Task::List getTasks();

		class SolveStep : public CPS::Task {
		public:
			SolveStep(DAESolver& solver) :
				Task(solver.mName + ".SolveStep"), mSolver(solver) {
				mModifiedAttributes.push_back(Scheduler::external);
			}
			void execute(Real time, Int timeStepCount) {
				if (time == mSolver.mInitTime)
					// when time==initTime only log the initial values
					mSolver.updateVoltageAndCurrents();
				else
					
    				mSolver.step(time);
				
			}
		private:
			DAESolver& mSolver;
		};

		///
		class LogTask : public CPS::Task {
		public:
			LogTask(DAESolver& solver) :
				Task(solver.mName + ".Log"), mSolver(solver) {
				mModifiedAttributes.push_back(Scheduler::external);
			}

			void execute(Real time, Int timeStepCount) {
				mSolver.log(time, timeStepCount);
			}

		private:
			DAESolver& mSolver;
		};
	};
}
