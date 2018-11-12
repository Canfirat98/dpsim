/** Simulation
 *
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017-2018, Institute for Automation of Complex Power Systems, EONERC
 *
 * DPsim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <list>

#include <dpsim/Definitions.h>
#include <dpsim/Config.h>
#include <cps/Logger.h>
#include <cps/SystemTopology.h>
#include <cps/Task.h>

namespace DPsim {
	/// Holds switching time and which system should be activated.
	struct SwitchConfiguration {
		Real switchTime;
		UInt systemIndex;
	};

	/// Base class for more specific solvers such as MNA, ODE or IDA.
	class Solver {
	public:
		typedef std::shared_ptr<Solver> Ptr;
		typedef std::vector<Ptr> List;

		virtual ~Solver() { }

		enum class Type { MNA, DAE };

		virtual CPS::Task::List getTasks() = 0;
		/// Log results
		virtual void log(Real time) { };
	};
}
