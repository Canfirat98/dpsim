/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

//-----------Power system-----------//
//Voltage level as Base Voltage
Real Vnom = 230e3;

//-----------Generator-----------//
Real nomPower = 500e6;
Real nomPhPhVoltRMS = 22e3;
Real nomFreq = 60;
Real nomOmega= nomFreq* 2*PI;
Real H = 5;
Real Xpd=0.31;
Real Rs = 0.003*0;
Real D = 1*50;
// Initialization parameters
Real initMechPower= 300e6;
Real initActivePower = 300e6;
Real setPointVoltage=nomPhPhVoltRMS + 0.05*nomPhPhVoltRMS;

//-----------Transformer-----------//
Real t_ratio=Vnom/nomPhPhVoltRMS;

//PiLine parameters calculated from CIGRE Benchmark system
Real lineResistance = 6.7;
Real lineInductance = 47./nomOmega;
// Real lineCapacitance = 3.42e-4/nomOmega;
Real lineCapacitance = 0;
Real lineConductance =0;

// Parameters for powerflow initialization
// Slack voltage: 1pu
Real Vslack = Vnom;

void EMT_1ph_SynGenTrStab_SteadyState(String simName, Real timeStep, Real finalTime, bool startFaultEvent, bool endFaultEvent, Real startTimeFault, Real endTimeFault, Real cmdInertia) {
	//  // ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components
	auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);

	//Synchronous generator ideal model
	auto genPF = SP::Ph1::SynchronGenerator::make("Generator", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	genPF->setParameters(nomPower, nomPhPhVoltRMS, initActivePower, setPointVoltage*t_ratio, PowerflowBusType::PV);
	genPF->setBaseVoltage(Vnom);
	genPF->modifyPowerFlowBusType(PowerflowBusType::PV);

	//Grid bus as Slack
	auto extnetPF = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::debug);
	extnetPF->setParameters(Vslack);
	extnetPF->setBaseVoltage(Vnom);
	extnetPF->modifyPowerFlowBusType(PowerflowBusType::VD);

	//Line
	auto linePF = SP::Ph1::PiLine::make("PiLine", Logger::Level::debug);
	linePF->setParameters(lineResistance, lineInductance, lineCapacitance, lineConductance);
	linePF->setBaseVoltage(Vnom);

	// Topology
	genPF->connect({ n1PF });
	linePF->connect({ n1PF, n2PF });
	extnetPF->connect({ n2PF });
	auto systemPF = SystemTopology(60,
			SystemNodeList{n1PF, n2PF},
			SystemComponentList{genPF, linePF, extnetPF});

	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->logAttribute("v1", n1PF->attribute("v"));
	loggerPF->logAttribute("v2", n2PF->attribute("v"));

	// Simulation
	Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
	simPF.doInitFromNodesAndTerminals(false);
	simPF.addLogger(loggerPF);
	simPF.run();

	// ----- Dynamic simulation ------
	String simNameEMT = simName + "_EMT";
	Logger::setLogDir("logs/"+simNameEMT);

	// Nodes
	auto n1EMT = SimNode<Real>::make("n1", PhaseType::ABC);
	auto n2EMT = SimNode<Real>::make("n2", PhaseType::ABC);

	// Components
	auto genEMT = EMT::Ph3::SynchronGeneratorTrStab::make("SynGen", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	genEMT->setStandardParametersPU(nomPower, nomPhPhVoltRMS, nomFreq, Xpd*std::pow(t_ratio,2), cmdInertia*H, Rs, D );

	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower= genPF->getApparentPower();
	genEMT->setInitialValues(initApparentPower, initMechPower);

	//Grid bus as Slack
	auto extnetEMT = EMT::Ph3::NetworkInjection::make("Slack", Logger::Level::debug);

	// Line
	auto lineEMT = EMT::Ph3::PiLine::make("PiLine", Logger::Level::debug);
	lineEMT->setParameters(Math::singlePhaseParameterToThreePhase(lineResistance),
	                      Math::singlePhaseParameterToThreePhase(lineInductance),
					      Math::singlePhaseParameterToThreePhase(lineCapacitance),
						  Math::singlePhaseParameterToThreePhase(lineConductance));

	// Topology
	genEMT->connect({ n1EMT });
	lineEMT->connect({ n1EMT, n2EMT });
	extnetEMT->connect({ n2EMT });
	auto systemEMT = SystemTopology(60,
			SystemNodeList{n1EMT, n2EMT},
			SystemComponentList{genEMT, lineEMT, extnetEMT});

	// Initialization of dynamic topology
	systemEMT.initWithPowerflow(systemPF);


	// Logging
	auto loggerEMT = DataLogger::make(simNameEMT);
	loggerEMT->logAttribute("v1", n1EMT->attribute("v"));
	loggerEMT->logAttribute("v2", n2EMT->attribute("v"));
	//gen
	loggerEMT->logAttribute("Ep", genEMT->attribute("Ep"));
	loggerEMT->logAttribute("v_gen", genEMT->attribute("v_intf"));
	loggerEMT->logAttribute("i_gen", genEMT->attribute("i_intf"));
	loggerEMT->logAttribute("wr_gen", genEMT->attribute("w_r"));
	loggerEMT->logAttribute("delta_r_gen", genEMT->attribute("delta_r"));
	loggerEMT->logAttribute("P_elec", genEMT->attribute("P_elec"));
	loggerEMT->logAttribute("P_mech", genEMT->attribute("P_mech"));
	//line
	loggerEMT->logAttribute("v_line", lineEMT->attribute("v_intf"));
	loggerEMT->logAttribute("i_line", lineEMT->attribute("i_intf"));
	//slack
	loggerEMT->logAttribute("v_slack", extnetEMT->attribute("v_intf"));
	loggerEMT->logAttribute("i_slack", extnetEMT->attribute("i_intf"));



	Simulation simEMT(simNameEMT, Logger::Level::debug);
	simEMT.setSystem(systemEMT);
	simEMT.setTimeStep(timeStep);
	simEMT.setFinalTime(finalTime);
	simEMT.setDomain(Domain::EMT);
	simEMT.addLogger(loggerEMT);

	simEMT.run();
}

int main(int argc, char* argv[]) {


	//Simultion parameters
	String simName="EMT_SynGenTrStab_SMIB_SteadyState";
	Real finalTime = 0.1;
	Real timeStep = 0.000001;
	Bool startFaultEvent=false;
	Bool endFaultEvent=false;
	Real startTimeFault=10;
	Real endTimeFault=10.1;
	Real cmdInertia= 1.0;

	EMT_1ph_SynGenTrStab_SteadyState(simName, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, cmdInertia);
}
