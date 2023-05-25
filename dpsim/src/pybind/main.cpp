/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <iomanip>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <dpsim/Simulation.h>
#include <dpsim/RealTimeSimulation.h>
#include <dpsim/SolverParameters.h>
#include <dpsim/SolverParametersMNA.h>
#include <dpsim-models/IdentifiedObject.h>
#include <DPsim.h>

#include <dpsim-models/CSVReader.h>

#include <dpsim/pybind/DPComponents.h>
#include <dpsim/pybind/EMTComponents.h>
#include <dpsim/pybind/SPComponents.h>
#include <dpsim/pybind/SignalComponents.h>
#include <dpsim/pybind/BaseComponents.h>
#include <dpsim/pybind/Utils.h>
#include <dpsim/pybind/Attributes.h>

PYBIND11_DECLARE_HOLDER_TYPE(T, CPS::AttributePointer<T>);

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(dpsimpy, m) {
    m.doc() = R"pbdoc(
	DPsim Python bindings
	-----------------------
	The Python bindings provide access to most of the DPsim features implemented in C++.
	It is possible to run powerflow, quasi-static, dynamic phasor and electromagnetic transient simulations
	and to parameterize all components of the network from Python.
    )pbdoc";

	py::enum_<CPS::Logger::Level>(m, "LogLevel")
		.value("trace", CPS::Logger::Level::trace)
		.value("debug", CPS::Logger::Level::debug)
		.value("info", CPS::Logger::Level::info)
		.value("warn", CPS::Logger::Level::warn)
		.value("err", CPS::Logger::Level::err)
		.value("critical", CPS::Logger::Level::critical)
		.value("off", CPS::Logger::Level::off);

	py::class_<CPS::Math>(m, "Math")
		.def_static("single_phase_variable_to_three_phase", &CPS::Math::singlePhaseVariableToThreePhase)
		.def_static("single_phase_parameter_to_three_phase", &CPS::Math::singlePhaseParameterToThreePhase)
		.def_static("single_phase_power_to_three_phase", &CPS::Math::singlePhasePowerToThreePhase);

	m.attr("RMS3PH_TO_PEAK1PH") = RMS3PH_TO_PEAK1PH;
	m.attr("PEAK1PH_TO_RMS3PH") = PEAK1PH_TO_RMS3PH;
	m.attr("P_SNUB_TRANSFORMER") = P_SNUB_TRANSFORMER;
	m.attr("Q_SNUB_TRANSFORMER") = Q_SNUB_TRANSFORMER;

	addAttributes(m);

	py::class_<DPsim::SolverParameters, std::shared_ptr<DPsim::SolverParameters>>(m, "SolverParameters")
		.def(py::init<>())
		.def("set_init_from_nodes_and_terminals", &DPsim::SolverParameters::setInitFromNodesAndTerminals)
		.def("get_init_from_nodes_and_terminals", &DPsim::SolverParameters::getInitFromNodesAndTerminals);

	py::class_<DPsim::SolverParametersMNA, DPsim::SolverParameters, std::shared_ptr<DPsim::SolverParametersMNA>>(m, "SolverParametersMNA", py::multiple_inheritance())
		.def(py::init<>())
		.def("set_domain", &DPsim::SolverParametersMNA::setDomain)
		.def("set_solver_and_component_behaviour", &DPsim::SolverParametersMNA::setSolverAndComponentBehaviour)
		.def("do_split_subnets", &DPsim::SolverParametersMNA::doSplitSubnets)
		.def("do_frequency_parallelization", &DPsim::SolverParametersMNA::doFrequencyParallelization)
		.def("do_system_matrix_recomputation", &DPsim::SolverParametersMNA::doSystemMatrixRecomputation)
		.def("set_direct_linear_solver_implementation", &DPsim::SolverParametersMNA::setDirectLinearSolverImplementation)
		.def("do_steady_state_init", &DPsim::SolverParametersMNA::doSteadyStateInit)
		.def("set_stead_st_ini_time_limit", &DPsim::SolverParametersMNA::setSteadStIniTimeLimit)
		.def("set_stead_st_ini_acc_limit", &DPsim::SolverParametersMNA::setSteadStIniAccLimit)
		.def("get_split_subnets", &DPsim::SolverParametersMNA::getSplitSubnets)
		.def("get_freq_parallel", &DPsim::SolverParametersMNA::getFreqParallel)
		.def("get_system_matrix_recomputation", &DPsim::SolverParametersMNA::getSystemMatrixRecomputation)
		.def("get_steady_state_init", &DPsim::SolverParametersMNA::getSteadyStateInit)
		.def("get_steady_state_init_time_limit", &DPsim::SolverParametersMNA::getSteadyStateInitTimeLimit)
		.def("get_steady_state_init_acc_limit", &DPsim::SolverParametersMNA::getSteadyStateInitAccLimit);

	py::class_<DPsim::DirectLinearSolverConfiguration>(m, "DirectLinearSolverConfiguration")
		.def(py::init<>())
		.def("set_fill_in_reduction_method", &DPsim::DirectLinearSolverConfiguration::setFillInReductionMethod)
		.def("set_scaling_method", &DPsim::DirectLinearSolverConfiguration::setScalingMethod)
		.def("set_partial_refactorization_method", &DPsim::DirectLinearSolverConfiguration::setPartialRefactorizationMethod)
		.def("set_btf", &DPsim::DirectLinearSolverConfiguration::setBTF)
		.def("get_scaling_method", &DPsim::DirectLinearSolverConfiguration::getScalingMethod)
		.def("get_fill_in_reduction_method", &DPsim::DirectLinearSolverConfiguration::getFillInReductionMethod)
		.def("get_partial_refactorization_method", &DPsim::DirectLinearSolverConfiguration::getPartialRefactorizationMethod)
		.def("get_btf", &DPsim::DirectLinearSolverConfiguration::getBTF);

    py::class_<DPsim::Simulation>(m, "Simulation")
	    .def(py::init<std::string, CPS::Logger::Level>(), "name"_a, "loglevel"_a = CPS::Logger::Level::off)
		.def("name", &DPsim::Simulation::name)
		.def("set_time_step", &DPsim::Simulation::setTimeStep)
		.def("set_final_time", &DPsim::Simulation::setFinalTime)
		.def("set_simulation_parameters", &DPsim::Simulation::setSimulationParameters, "time_step"_a, "final_time"_a)
		.def("add_logger", &DPsim::Simulation::addLogger)
		.def("set_system", &DPsim::Simulation::setSystem)
		.def("run", &DPsim::Simulation::run)
		.def("set_solver", &DPsim::Simulation::setSolverType)
		.def("set_domain", &DPsim::Simulation::setDomain)
		.def("set_solver_parameters", &DPsim::Simulation::setSolverParameters, "domain"_a, "type"_a, "solverParameters"_a)
		.def("start", &DPsim::Simulation::start)
		.def("next", &DPsim::Simulation::next)
		.def("get_idobj_attr", &DPsim::Simulation::getIdObjAttribute, "comp"_a, "attr"_a)
		.def("add_interface", &DPsim::Simulation::addInterface, "interface"_a)
		.def("log_idobj_attribute", &DPsim::Simulation::logIdObjAttribute, "comp"_a, "attr"_a)
		.def("log_attribute", &DPsim::Simulation::logAttribute, "name"_a, "attr"_a)
		.def("set_tearing_components", &DPsim::Simulation::setTearingComponents)
		.def("add_event", &DPsim::Simulation::addEvent)
		.def("set_direct_linear_solver_configuration", &DPsim::Simulation::setDirectLinearSolverConfiguration)
		.def("log_lu_times", &DPsim::Simulation::logLUTimes)
		.def("set_relative_rolerance", &DPsim::Simulation::setRelativeTolerance, "rel_tol"_a)
		.def("set_ida_max_conv_fails", &DPsim::Simulation::setIDAMaxConvFails,   "ida_max_conv_fails"_a)
		.def("set_ida_non_lin_conv_coef", &DPsim::Simulation::setIDANonlinConvCoef, "ida_non_lin_conv_coef"_a)
		.def("set_ida_max_non_lin_iters", &DPsim::Simulation::setIDAMaxNonlinIters, "ida_max_non_lin_iters"_a)
		.def("set_ida_max_bdf_order", &DPsim::Simulation::setIDAMaxBDFOrder, "ida_max_bdf_order"_a)
		.def("set_ida_max_num_steps", &DPsim::Simulation::setIDAMaxNumSteps, "ida_max_num_steps"_a)
		.def("set_variable_step_size", &DPsim::Simulation::setVariableStepSize, "variable_step_size"_a)
		.def("set_min_step_size", &DPsim::Simulation::setMinStepSize, "min_step_size"_a)
		.def("set_max_step_size", &DPsim::Simulation::setMaxStepSize, "max_step_size"_a)
		.def("set_max_err_test_fails", &DPsim::Simulation::setMaxErrTestFails, "max_err_test_fails"_a)
		.def("set_ida_suppress_alg", &DPsim::Simulation::setIDASetSuppressAlg, "ida_suppress_alg"_a)
		.def("set_min_step_size", &DPsim::Simulation::setMinStepSize, "max_step_size"_a)
		.def("set_max_step_size", &DPsim::Simulation::setMaxStepSize, "max_step_size"_a);

	py::class_<DPsim::RealTimeSimulation, DPsim::Simulation>(m, "RealTimeSimulation")
		.def(py::init<std::string, CPS::Logger::Level>(), "name"_a, "loglevel"_a = CPS::Logger::Level::info)
		.def("name", &DPsim::RealTimeSimulation::name)
		.def("set_time_step", &DPsim::RealTimeSimulation::setTimeStep)
		.def("set_final_time", &DPsim::RealTimeSimulation::setFinalTime)
		.def("add_logger", &DPsim::RealTimeSimulation::addLogger)
		.def("set_system", &DPsim::RealTimeSimulation::setSystem)
		.def("run", static_cast<void (DPsim::RealTimeSimulation::*)(CPS::Int startIn)>(&DPsim::RealTimeSimulation::run))
		.def("set_solver", &DPsim::RealTimeSimulation::setSolverType)
		.def("set_domain", &DPsim::RealTimeSimulation::setDomain);

	py::class_<CPS::SystemTopology, std::shared_ptr<CPS::SystemTopology>>(m, "SystemTopology")
        .def(py::init<CPS::Real, CPS::TopologicalNode::List, CPS::IdentifiedObject::List>())
		.def(py::init<CPS::Real, CPS::Matrix, CPS::TopologicalNode::List, CPS::IdentifiedObject::List>())
		.def(py::init<CPS::Real>())
		.def("add", &DPsim::SystemTopology::addComponent)
		.def("add", &DPsim::SystemTopology::addComponents)
		.def("node", py::overload_cast<std::string_view>(&DPsim::SystemTopology::node<CPS::TopologicalNode>))
		.def("node", py::overload_cast<CPS::UInt>(&DPsim::SystemTopology::node<CPS::TopologicalNode>))
		.def("connect_component", py::overload_cast<CPS::SimPowerComp<CPS::Real>::Ptr, CPS::SimNode<CPS::Real>::List>(&DPsim::SystemTopology::connectComponentToNodes<CPS::Real>))
		.def("connect_component", py::overload_cast<CPS::SimPowerComp<CPS::Complex>::Ptr, CPS::SimNode<CPS::Complex>::List>(&DPsim::SystemTopology::connectComponentToNodes<CPS::Complex>))
		.def("component", &DPsim::SystemTopology::component<CPS::TopologicalPowerComp>)
		.def("add_tear_component", &DPsim::SystemTopology::addTearComponent)
#ifdef WITH_GRAPHVIZ
		.def("_repr_svg_", &DPsim::SystemTopology::render)
		.def("render_to_file", &DPsim::SystemTopology::renderToFile)
#endif
		.def_readwrite("nodes", &DPsim::SystemTopology::mNodes)
		.def_readwrite("components", &DPsim::SystemTopology::mComponents)
		.def_readonly("tear_components", &DPsim::SystemTopology::mTearComponents)
		.def("list_idobjects", &DPsim::SystemTopology::listIdObjects)
		.def("init_with_powerflow", &DPsim::SystemTopology::initWithPowerflow);

	py::class_<DPsim::Interface, std::shared_ptr<DPsim::Interface>>(m, "Interface");

	py::class_<DPsim::DataLogger, std::shared_ptr<DPsim::DataLogger>>(m, "Logger")
        .def(py::init<std::string>())
		.def_static("set_log_dir", &CPS::Logger::setLogDir)
		.def_static("get_log_dir", &CPS::Logger::logDir)
		.def("log_attribute", py::overload_cast<const CPS::String&, CPS::AttributeBase::Ptr, CPS::UInt, CPS::UInt>(&DPsim::DataLogger::logAttribute), "name"_a, "attr"_a, "max_cols"_a = 0, "max_rows"_a = 0)
		/// Compatibility method. Might be removed later when the python examples have been fully adapted.
		.def("log_attribute", py::overload_cast<const std::vector<CPS::String>&, CPS::AttributeBase::Ptr>(&DPsim::DataLogger::logAttribute), "names"_a, "attr"_a)
		/// Compatibility method. Might be removed later when the python examples have been fully adapted.
		.def("log_attribute", [](DPsim::DataLogger &logger, const CPS::String &name, const CPS::String &attr, const CPS::IdentifiedObject &comp, CPS::UInt rowsMax, CPS::UInt colsMax) {
			logger.logAttribute(name, comp.attribute(attr), rowsMax, colsMax);
		}, "name"_a, "attr"_a, "comp"_a, "rows_max"_a = 0, "cols_max"_a = 0)
		/// Compatibility method. Might be removed later when the python examples have been fully adapted.;
		.def("log_attribute", [](DPsim::DataLogger &logger, const std::vector<CPS::String> &names, const CPS::String &attr, const CPS::IdentifiedObject &comp) {
			logger.logAttribute(names, comp.attribute(attr));
		});

	py::class_<CPS::IdentifiedObject, std::shared_ptr<CPS::IdentifiedObject>>(m, "IdentifiedObject")
		.def("name", &CPS::IdentifiedObject::name)
		/// CHECK: It would be nicer if all the attributes of an IdObject were bound as properties so they show up in the documentation and auto-completion.
		/// I don't know if this is possible to do because it depends on if the attribute map is filled before or after the code in this file is run.
		/// Manually adding the attributes would of course be possible but very tedious to do for all existing components / attributes
		.def("attr", &CPS::IdentifiedObject::attribute, "name"_a)
		.def("print_attribute_list", &printAttributes)
		.def("print_attribute", &printAttribute, "attribute_name"_a)
		.def("__str__", &getAttributeList);

	//Enums

	py::enum_<DPsim::Solver::Behaviour>(m, "SolverBehaviour")
		.value("Initialization", DPsim::Solver::Behaviour::Initialization)
		.value("Simulation", DPsim::Solver::Behaviour::Simulation);

	py::enum_<CPS::Domain>(m, "Domain")
		.value("SP", CPS::Domain::SP)
		.value("DP", CPS::Domain::DP)
		.value("EMT", CPS::Domain::EMT);

	py::enum_<CPS::PhaseType>(m, "PhaseType")
		.value("A", CPS::PhaseType::A)
		.value("B", CPS::PhaseType::B)
		.value("C", CPS::PhaseType::C)
		.value("ABC", CPS::PhaseType::ABC)
		.value("Single", CPS::PhaseType::Single);

	py::enum_<CPS::PowerflowBusType>(m, "PowerflowBusType")
		.value("PV", CPS::PowerflowBusType::PV)
		.value("PQ", CPS::PowerflowBusType::PQ)
		.value("VD", CPS::PowerflowBusType::VD)
		.value("None", CPS::PowerflowBusType::None);

	py::enum_<CPS::GeneratorType>(m, "GeneratorType")
		.value("PVNode", CPS::GeneratorType::PVNode)
		.value("TransientStability", CPS::GeneratorType::TransientStability)
		.value("IdealVoltageSource", CPS::GeneratorType::IdealVoltageSource)
		.value("SG3OrderVBR", CPS::GeneratorType::SG3OrderVBR)
		.value("SG4OrderVBR", CPS::GeneratorType::SG4OrderVBR)
		.value("SG6aOrderVBR", CPS::GeneratorType::SG6aOrderVBR)
		.value("SG6bOrderVBR", CPS::GeneratorType::SG6bOrderVBR)
		.value("FullOrderVBR", CPS::GeneratorType::FullOrderVBR)
		.value("FullOrder", CPS::GeneratorType::FullOrder)
		.value("NONE", CPS::GeneratorType::None);

	py::enum_<DPsim::Solver::Type>(m, "Solver")
		.value("MNA", DPsim::Solver::Type::MNA)
		.value("DAE", DPsim::Solver::Type::DAE)
		.value("NRP", DPsim::Solver::Type::NRP);

	/*
	py::enum_<DPsim::DirectLinearSolverImpl>(m, "DirectLinearSolverImpl")
		.value("Undef", DPsim::DirectLinearSolverImpl::Undef)
		.value("DenseLU", DPsim::DirectLinearSolverImpl::DenseLU)
		.value("SparseLU", DPsim::DirectLinearSolverImpl::SparseLU)
		.value("KLU", DPsim::DirectLinearSolverImpl::KLU)
		.value("CUDADense", DPsim::DirectLinearSolverImpl::CUDADense)
		.value("CUDASparse", DPsim::DirectLinearSolverImpl::CUDASparse)
		.value("CUDAMagma", DPsim::DirectLinearSolverImpl::CUDAMagma);
	*/

	py::enum_<DPsim::SCALING_METHOD>(m, "scaling_method")
		.value("no_scaling", DPsim::SCALING_METHOD::NO_SCALING)
		.value("sum_scaling", DPsim::SCALING_METHOD::SUM_SCALING)
		.value("max_scaling", DPsim::SCALING_METHOD::MAX_SCALING);

	py::enum_<DPsim::FILL_IN_REDUCTION_METHOD>(m, "fill_in_reduction_method")
		.value("amd", DPsim::FILL_IN_REDUCTION_METHOD::AMD)
		.value("amd_nv", DPsim::FILL_IN_REDUCTION_METHOD::AMD_NV)
		.value("amd_ra", DPsim::FILL_IN_REDUCTION_METHOD::AMD_RA)
		.value("colamd", DPsim::FILL_IN_REDUCTION_METHOD::COLAMD);

	py::enum_<DPsim::PARTIAL_REFACTORIZATION_METHOD>(m, "partial_refactorization_method")
		.value("no_partial_refactorization", DPsim::PARTIAL_REFACTORIZATION_METHOD::NO_PARTIAL_REFACTORIZATION)
		.value("factorization_path", DPsim::PARTIAL_REFACTORIZATION_METHOD::FACTORIZATION_PATH)
		.value("refactorization_restart", DPsim::PARTIAL_REFACTORIZATION_METHOD::REFACTORIZATION_RESTART);

	py::enum_<DPsim::USE_BTF>(m, "use_btf")
		.value("no_btf", DPsim::USE_BTF::NO_BTF)
		.value("do_btf", DPsim::USE_BTF::DO_BTF);

	py::enum_<CPS::CSVReader::Mode>(m, "CSVReaderMode")
		.value("AUTO", CPS::CSVReader::Mode::AUTO)
		.value("MANUAL", CPS::CSVReader::Mode::MANUAL);

	py::enum_<CPS::CSVReader::DataFormat>(m, "CSVReaderFormat")
		.value("HHMMSS", CPS::CSVReader::DataFormat::HHMMSS)
		.value("SECONDS", CPS::CSVReader::DataFormat::SECONDS)
		.value("HOURS", CPS::CSVReader::DataFormat::HOURS)
		.value("MINUTES", CPS::CSVReader::DataFormat::MINUTES);

#ifdef WITH_CIM
	py::class_<CPS::CIM::Reader>(m, "CIMReader")
		.def(py::init<std::string, CPS::Logger::Level, CPS::Logger::Level>(), "name"_a, "loglevel"_a = CPS::Logger::Level::info, "comploglevel"_a = CPS::Logger::Level::off)
		.def("loadCIM", (CPS::SystemTopology (CPS::CIM::Reader::*)(CPS::Real, const std::list<CPS::String> &, CPS::Domain, CPS::PhaseType, CPS::GeneratorType)) &CPS::CIM::Reader::loadCIM);
#endif

	py::class_<CPS::CSVReader>(m, "CSVReader")
		.def(py::init<std::string, const std::string &, std::map<std::string, std::string> &, CPS::Logger::Level>())
		.def("assignLoadProfile", &CPS::CSVReader::assignLoadProfile);

	//Base Classes

	py::class_<CPS::TopologicalPowerComp, std::shared_ptr<CPS::TopologicalPowerComp>, CPS::IdentifiedObject>(m, "TopologicalPowerComp");
	py::class_<CPS::SimPowerComp<CPS::Complex>, std::shared_ptr<CPS::SimPowerComp<CPS::Complex>>, CPS::TopologicalPowerComp>(m, "SimPowerCompComplex")
		.def("connect", &CPS::SimPowerComp<CPS::Complex>::connect)
		.def("set_intf_current", &CPS::SimPowerComp<CPS::Complex>::setIntfCurrent)
		.def("set_intf_voltage", &CPS::SimPowerComp<CPS::Complex>::setIntfVoltage)
		.def("get_terminal", &CPS::SimPowerComp<CPS::Complex>::terminal, "index"_a);
	py::class_<CPS::SimPowerComp<CPS::Real>, std::shared_ptr<CPS::SimPowerComp<CPS::Real>>, CPS::TopologicalPowerComp>(m, "SimPowerCompReal")
		.def("connect", &CPS::SimPowerComp<CPS::Real>::connect)
		.def("set_intf_current", &CPS::SimPowerComp<CPS::Real>::setIntfCurrent)
		.def("set_intf_voltage", &CPS::SimPowerComp<CPS::Real>::setIntfVoltage)
		.def("get_terminal", &CPS::SimPowerComp<CPS::Real>::terminal, "index"_a);
	py::class_<CPS::TopologicalNode, std::shared_ptr<CPS::TopologicalNode>, CPS::IdentifiedObject>(m, "TopologicalNode")
		.def("initial_single_voltage", &CPS::TopologicalNode::initialSingleVoltage, "phase_type"_a = CPS::PhaseType::Single);

	py::class_<CPS::TopologicalTerminal, std::shared_ptr<CPS::TopologicalTerminal>, CPS::IdentifiedObject>(m, "TopologicalTerminal")
		.def("set_power", py::overload_cast<CPS::Complex>(&CPS::TopologicalTerminal::setPower))
		.def("set_power", py::overload_cast<CPS::MatrixComp>(&CPS::TopologicalTerminal::setPower));

	py::class_<CPS::SimTerminal<CPS::Complex>, std::shared_ptr<CPS::SimTerminal<CPS::Complex>>, CPS::TopologicalTerminal>(m, "SimTerminalComplex");
	py::class_<CPS::SimTerminal<CPS::Real>, std::shared_ptr<CPS::SimTerminal<CPS::Real>>, CPS::TopologicalTerminal>(m, "SimTerminalReal");


	//Events
	py::module mEvent = m.def_submodule("event", "events");
	py::class_<DPsim::Event, std::shared_ptr<DPsim::Event>>(mEvent, "Event");
	py::class_<DPsim::SwitchEvent, std::shared_ptr<DPsim::SwitchEvent>, DPsim::Event>(mEvent, "SwitchEvent", py::multiple_inheritance())
		.def(py::init<CPS::Real,const std::shared_ptr<CPS::Base::Ph1::Switch>,CPS::Bool>());
	py::class_<DPsim::SwitchEvent3Ph, std::shared_ptr<DPsim::SwitchEvent3Ph>, DPsim::Event>(mEvent, "SwitchEvent3Ph", py::multiple_inheritance())
		.def(py::init<CPS::Real,const std::shared_ptr<CPS::Base::Ph3::Switch>,CPS::Bool>());

	//Components
	py::module mBase = m.def_submodule("base", "base models");
	addBaseComponents(mBase);

	py::module mDP = m.def_submodule("dp", "dynamic phasor models");
	addDPComponents(mDP);

	py::module mEMT = m.def_submodule("emt", "electromagnetic-transient models");
	addEMTComponents(mEMT);

	py::module mSP = m.def_submodule("sp", "static phasor models");
	mSP.attr("SimNode") = mDP.attr("SimNode");
	addSPComponents(mSP);

	py::module mSignal = m.def_submodule("signal", "signal models");
	addSignalComponents(mSignal);


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
