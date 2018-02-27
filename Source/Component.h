/** Base component
 *
 * @file
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
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
#include <memory>

#include "SharedFactory.h"
#include "SystemModel.h"
#include "Definitions.h"
#include "MathLibrary.h"

namespace DPsim {
	class Terminal;
	class Node;

	/// Base class for all components that might be added to the matrix.
	class Component {
	public:
		struct Attribute {
			enum Type {
				Real,
				Integer,
				String, // value should be *String, not *char!
				Complex
			} mType;
			void* mValue;

			typedef std::map<DPsim::String, Attribute> Map;
		};
	protected:
		/// Component logger
		Logger mLog;
		/// Component logger control for internal variables
		Logger::Level mLogLevel;
		/// Human readable name
		String mName;
		/// Component node 1
		Matrix::Index mNode1 = 0;
		/// Component node 2
		Matrix::Index mNode2 = 0;
		/// Component node 3
		Matrix::Index mNode3 = 0;
		/// Determines the number of Terminals which can be connected to nodes
		Int mNumTerminals = 0;
		/// List of Terminals
		std::vector<std::shared_ptr<Terminal>> mTerminals;
		/// Determines if the component has a virtual node
		Int mNumVirtualNodes = 0;
		/// List of virtual nodes
		std::vector<std::shared_ptr<Node>> mVirtualNodes;
		/// Map of all attributes that should be exported to the Python interface
		Attribute::Map attrMap;
	public:
		typedef std::shared_ptr<Component> Ptr;
		typedef std::vector<Ptr> List;
		/// Unique identifier
		String mUID;			
		///
		Component(String uid, String name, Logger::Level logLevel = Logger::Level::NONE);
		///
		Component(String name, Logger::Level logLevel = Logger::Level::NONE);
		/// Creates a new component with basic features: name, nodes and the log level for this component		
		Component(String name, Matrix::Index node1, Matrix::Index node2, Logger::Level logLevel = Logger::Level::NONE);
		/// Creates a new component with basic features: name, nodes and the log level for this component
		Component(String name, Matrix::Index node1, Matrix::Index node2, Matrix::Index node3, Logger::Level loglevel = Logger::Level::NONE);
		//
		virtual ~Component() { }
		///
		String getName() { return mName; }
		///
		String getType();
		///
		std::map<String, Attribute>& getAttrMap() { return attrMap; }
		/// get value of node1
		Matrix::Index getNode1() { return mNode1; }
		/// get value of node2
		Matrix::Index getNode2() { return mNode2; }
		/// get value of node3
		Matrix::Index getNode3() { return mNode3; }
		/// Returns true if virtual node number is greater than zero.
		Bool hasVirtualNodes() { return mNumVirtualNodes > 0; }
		/// Returns true if virtual node number is greater than zero.
		Int getVirtualNodesNum() { return mNumVirtualNodes; }
		///
		Int getTerminalsNum() { return mNumTerminals; }
		///
		void setVirtualNodeAt(std::shared_ptr<Node> virtualNode, Int nodeNum);
		/// Set Terminals of the component
		virtual void setTerminals(std::vector<std::shared_ptr<Terminal>> terminals);
		///
		virtual void setTerminalAt(std::shared_ptr<Terminal> terminal, Int terminalPosition);
		///
		virtual void setNodes(std::vector<std::shared_ptr<Node>> nodes);
		///
		virtual void initializePowerflow(Real frequency) { }
		// #### MNA section ####
		/// Initializes variables of components
		virtual void initialize(SystemModel& system) { }
		/// Stamps conductance matrix
		virtual void applySystemMatrixStamp(SystemModel& system) = 0;
		/// Stamps current source vector
		virtual void applyRightSideVectorStamp(SystemModel& system) { }
		/// Upgrade values on the current source vector
		virtual void step(SystemModel& system, Real time) { }
		/// Upgrade variable values based on the solution of the step
		virtual void postStep(SystemModel& system) { }		
	};
}
