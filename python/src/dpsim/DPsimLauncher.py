import warnings
import sys
sys.path.insert(0,'/home/mmo/git/Can/dpsim/build')
import dpsimpy
import json
import cimpy
from cimpy import CIM2DPsim
from cimpy.CIM2DPsim import Domain
from cimpy.CIM2DPsim import SGModels

TIME_STEP_PF = 0.1
FINAL_TIME_PF = 0.1

class DPsimLauncher:
    def __init__(self, file_path = '/home/mmo-cya/dpsim/python/src/dpsim/Config_network.json'):
        # data contained in the json file
        self.data = {}
        # cimpy topology
        self.cimpy_topology = None
        # log level used in simulations
        self.loglevel = dpsimpy.LogLevel.info
        # domain used in simulations
        self.domain = Domain.DP
        #
        self.dpsimpy_components = self.dpsimpy_components = dpsimpy.dp.ph1
        self.dpsimpy_components.SimNode = getattr(dpsimpy.dp, "SimNode") 
        # simulation start time
        self.start_time = 0
        # simulation end time
        self.end_time = 1
        # simulation step time
        self.time_step = 1e-3
        # simulation name
        self.sim_name = "sim_name"
        # syn gen model used in dynamic simulations
        self.syn_gen_model = SGModels.VBR_4Order
        #
        self.init_from_nodes_and_terminals = True
        
        # power flow simulations
        self.init_from_power_flow = True
        self.sim_name_pf = ""
        self.system_pf = None
        self.sim_pf = None
        
        # dynamic simulation
        self.system = None
        self.sim = None
        
        # read json file
        self.read_json(file_path)
                
    def read_json(self, file_path):
        """
        Diese Funktion liest eine Konfigurationsdatei im JSON-Format ein, in der der Dateipfad zu dem CIM-Netzwerk und deren Parameter zur Simulation
        mittels DPsim enthalten sind. Das Netzwerk wird eingelesen und ins DPsim-Format übersetzt. Der ausgewählte Fehler wird dem Netzwerk hinzugefügt und anschließend wird 
        die Simulation ausgeführt. Ausgegeben werden das neue DPsim Netzwerk und die Simulationsergebnisse mit dem gewünschten Solver.
        """
        
        # Read json file
        try:
            with open(file_path, 'r') as json_file:
                self.data = json.load(json_file)
        except FileNotFoundError:
            raise Exception(f"Die Datei '{file_path}' can not be found!")
        except json.JSONDecodeError:
            raise Exception(f"The file '{file_path}' is not a valid JSON File!")

        # Convert input topology to cimpy
        self.cimpy_topology = cimpy.cim_import(self.data['CIMFiles'], 'cgmes_v2_4_15')
    
        # get loglevel
        if "LogLevel" in self.data["SimulationParameters"].keys():
            if self.data["SimulationParameters"]["LogLevel"] == "debug":
                self.loglevel = dpsimpy.LogLevel.debug
            elif self.data["SimulationParameters"]["LogLevel"] == "info":
                self.loglevel = dpsimpy.LogLevel.info
            elif self.data["SimulationParameters"]["LogLevel"] == "off":
                self.loglevel = dpsimpy.LogLevel.off            
    
        # get domain
        self.dpsimpy_components = None
        if "Domain" in self.data["SolverParameters"].keys():
            if self.data["SolverParameters"]["Domain"] == "RMS":
                self.domain = Domain.SP
                self.dpsimpy_components = dpsimpy.sp.ph1
                self.dpsimpy_components.SimNode = getattr(dpsimpy.sp, "SimNode")
            elif self.data["SolverParameters"]["Domain"] == "DP":
                self.domain = Domain.DP
                self.dpsimpy_components = dpsimpy.dp.ph1
                self.dpsimpy_components.SimNode = getattr(dpsimpy.dp, "SimNode") 
            elif self.data["SolverParameters"]["Domain"] == "EMT":
                self.domain = Domain.EMT
                self.dpsimpy_components = dpsimpy.emt.ph3
                self.dpsimpy_components.SimNode = getattr(dpsimpy.emt, "SimNode") 
            else:
                raise Exception('ERROR: domain {} is not supported in dpsimpy'.format(self.data["SolverParameters"]["Domain"]))

        # get start time
        if "StartTime" in self.data["SimulationParameters"].keys():
            self.start_time = self.data["SimulationParameters"]["StartTime"]

        # get end time
        if "EndTime" in self.data["SimulationParameters"].keys():
            self.end_time = self.data["SimulationParameters"]["EndTime"]

        # get time step
        if "TimeStep" in self.data["SimulationParameters"].keys():
            self.time_step = self.data["SimulationParameters"]["TimeStep"]  

        # get simulation name
        if "SimulationName" in self.data["SimulationParameters"].keys():
            self.sim_name = self.data["SimulationParameters"]["SimulationName"] 

        # check if pf simulation for initializartion has to be performed
        if "InitFromPF" in (self.data['SimulationParameters']):
            self.init_from_power_flow = self.data['SimulationParameters']["InitFromPF"]
            
        # check if pf simulation for initializartion has to be performed
        if "InitFromNodesAndTerminals" in (self.data['SimulationParameters']):
            self.init_from_nodes_and_terminals = self.data['SimulationParameters']["InitFromNodesAndTerminals"]

        # get syn gen model to use in dynamic simulation
        # TODO: ADD MORE SYNGEN MODELS
        if "SGModel" in self.data["SimulationParameters"].keys():
            if self.data["SimulationParameters"] == "3Order":
                self.syn_gen_model = SGModels.VBR_3Order
            elif self.data["SimulationParameters"] == "4Order":
                self.syn_gen_model = SGModels.VBR_4Order
            elif self.data["SimulationParameters"] == "5Order":
                self.syn_gen_model = SGModels.VBR_5Order
            elif self.data["SimulationParameters"] == "6Order":
                self.syn_gen_model = SGModels.VBR_6Order
            else:
                # TODO: WARN MESSAGE
                pass
    
    def start_pf_simulation(self):
        if self.init_from_power_flow:
            self.run_pf_simulation()
        else:
            #TODO
            #load initial pf results from cim files
            pass
        
    def start_dynamic_simulation(self):
        self.name_dyn = self.sim_name + "_Dyn"
        dpsimpy.Logger.set_log_dir("logs/" + self.name_dyn)

         # create dpsim topology for dynamic simulation
        self.system = CIM2DPsim.CIM2DPsim(self.cimpy_topology, domain=self.domain, log_level=self.loglevel, frequency = 60, gen_model=self.syn_gen_model)
        
        ### Simulation
        self.sim = dpsimpy.Simulation(self.name_dyn, self.loglevel)
        self.sim.set_system(self.system)
        if self.domain == Domain.SP:
            self.sim.set_domain(dpsimpy.Domain.SP)
            self.system.init_with_powerflow(systemPF=self.system_pf, domain=dpsimpy.Domain.SP)
        elif self.domain == Domain.DP:
            self.sim.set_domain(dpsimpy.Domain.DP)
            self.system.init_with_powerflow(systemPF=self.system_pf, domain=dpsimpy.Domain.DP)
        elif self.domain == Domain.EMT:
            self.sim.set_domain(dpsimpy.Domain.EMT)
            self.system.init_with_powerflow(systemPF=self.system_pf, domain=dpsimpy.Domain.EMT)
            
        # add events
        self.read_events()
        
        #
        self.sim.do_init_from_nodes_and_terminals(self.init_from_nodes_and_terminals)
        # TODO: ADD SOLVERIMPLEMENTATION to json file
        self.sim.set_direct_solver_implementation(dpsimpy.DirectLinearSolverImpl.SparseLU)
        # TODO: DPSim: AUTOMATIC DETECTION OF system_matrix_recomputation?
        self.sim.do_system_matrix_recomputation(True)
        #
        #self.sim.set_start_time(self.start_time)
        self.sim.set_time_step(self.time_step)
        self.sim.set_final_time(self.end_time)
        
        # create logger add variables to be logged
        logger = dpsimpy.Logger(self.name_dyn)
        logger = self.read_variables_to_log(self.system, logger)
        self.sim.add_logger(logger)
        
        # run simulation
        self.sim.run()
        
    def run_pf_simulation(self):
        # 
        self.sim_name_pf = self.sim_name + "_PF"
        dpsimpy.Logger.set_log_dir("logs/" + self.sim_name_pf)

        # create dpsim topology for pf simulation for initialization
        self.system_pf = CIM2DPsim.CIM2DPsim(self.cimpy_topology, domain=Domain.PF, log_level=self.loglevel, frequency = 60)
    
        #set reference node 
        reference_comp=None
        for node, comp_list in self.system_pf.components_at_node.items():
            if (node.name==self.data['SimulationParameters']['ReferenceNode']):
                for comp in comp_list:
                    if (isinstance(comp, dpsimpy.sp.ph1.SynchronGenerator) or isinstance(comp, dpsimpy.sp.ph1.NetworkInjektion)):
                        reference_comp=comp
                        break
                if (reference_comp is None):  
                    raise Exception('No SynchronGenerator or ExternalNetworkInjection is connected to node: {}!'.format(node.name))    
        if (reference_comp is None):  
            raise Exception('No node named: {} was found!'.format(node.name))
    
        reference_comp.modify_power_flow_bus_type(dpsimpy.PowerflowBusType.VD)
    
        # create logger add variables to be logged
        logger = dpsimpy.Logger(self.sim_name_pf)
        logger = self.read_variables_to_log(self.system_pf, logger)
        
        # start power flow simulation
        self.sim_pf = dpsimpy.Simulation(self.sim_name_pf, self.loglevel)
        self.sim_pf.set_system(self.system_pf)
        self.sim_pf.set_domain(dpsimpy.Domain.SP)
        self.sim_pf.set_solver(dpsimpy.Solver.NRP)
        #self.sim_pf.set_solver_component_behaviour(dpsimpy.SolverBehaviour.Initialization)
        self.sim_pf.do_init_from_nodes_and_terminals(False)
        self.sim_pf.add_logger(logger)
        self.sim_pf.set_time_step(TIME_STEP_PF)
        self.sim_pf.set_final_time(FINAL_TIME_PF)
        self.sim_pf.run()        
    
    def read_variables_to_log(self, system, logger):
        # TODO: Add more variables and components
        for component_type, var_dict in self.data["LoggerVariables"].items():
            if component_type=="Node":
                for variable in var_dict.keys():
                    if variable=="Voltages":
                        if (var_dict["Voltages"][0]=="all"):
                            for node in system.nodes:
                                logger.log_attribute(node.name+'.V', 'v', node)
                        else:
                            for node_name in var_dict["Voltages"]:
                                # search for node in sim_pf
                                node_found = False
                                for node in system.nodes:
                                    if (node.name == node_name):
                                        node_found = True
                                        break
                                if (node_found==False):
                                    # TODO: change exception by warning
                                    raise Exception('Node {} was not found!'.format(node_name))

                                logger.log_attribute(node_name+'.V', 'v', node)
                    elif variable=="Power":
                        if (var_dict["Power"][0]=="all"):
                            for node in system.nodes:
                                logger.log_attribute(node.name+'.S', 's', node)
                        else:
                            for node_name in var_dict["Power"]:
                                # search for node in sim_pf
                                node_found = False
                                for node in system.nodes:
                                    if (node.name == node_name):
                                        node_found = True
                                        break
                                if (node_found==False):
                                    # TODO: change exception by warning
                                    raise Exception('Node {} was not found!'.format(node_name))
                                logger.log_attribute(node_name+'.S', 's', node)
                        
        return logger
        
    def read_events(self):
        # Füge der Topologie den gewünschten Event hinzu
        if "Events" in self.data:
            if self.data["Events"]['EventType'] == "ShortCircuit":
                event_params = self.data["Events"]['EventParameters']
                
                node_name = ""
                if "NodeName" in event_params.keys():
                    node_name = event_params['NodeName']
                else:
                    raise Exception('Paramenter "NodeName" is not in the json file!')
                
                node = None
                for node_ in self.system.nodes:
                    if event_params['NodeName'] == node_name:
                        node = node_
                        break
                if (node is None):
                    # TODO: change exception by warning
                    raise Exception('Node {} was not found!'.format(node_name))   
                    
                open_resistance = 1e12
                if "FaultOpenResistance" in event_params.keys():
                    open_resistance = event_params['FaultOpenResistance']
                else:
                    warnings.warn('Paramenter "FaultOpenResistance" is not in the json file!\n FaultOpenResistance set to 1e12')
                    
                closed_resistance = 1e-3  
                if "FaultClosedResistance" in event_params.keys():
                    closed_resistance = event_params['FaultClosedResistance']
                else:
                    warnings.warn('Paramenter "FaultClosedResistance" is not in the json file!\n FaultClosedResistance set to 1e-3')
                    
                # Füge switch mit Erdung hinzu
                switch = self.dpsimpy_components.Switch('Fault_' + node_name, self.loglevel)
                switch.set_parameters(open_resistance, closed_resistance)
                switch.open()
                self.system.add(switch)
                self.system.connect_component(switch, [self.dpsimpy_components.SimNode.gnd, node])

                # Event hinzufügen
                # TODO: READ TIMES FROM JSON FILE
                sw_event_1 = dpsimpy.event.SwitchEvent(1.0, switch, True)
                self.sim.add_event(sw_event_1)
                sw_event_2 = dpsimpy.event.SwitchEvent(1.1, switch, False)
                self.sim.add_event(sw_event_2)
                
            elif self.data["Events"]['EventType'] == "LoadStep":
                event_params = self.data["Events"]['EventParameters']
                
                node_name = ""
                if "NodeName" in event_params.keys():
                    node_name = event_params['NodeName']
                else:
                    raise Exception('Paramenter "NodeName" is not in the json file!')
                
                node = None
                for node_ in self.system.nodes:
                    if event_params['NodeName'] == node_name:
                        node = node_
                        break
                if (node is None):
                    # TODO: change exception by warning
                    raise Exception('Node {} was not found!'.format(node_name))   
                
                open_resistance = 1e12
                if "FaultOpenResistance" in event_params.keys():
                    open_resistance = event_params['FaultOpenResistance']
                else:
                    warnings.warn('Paramenter "FaultOpenResistance" is not in the json file!\n FaultOpenResistance set to 1e12')
                    
                closed_resistance = 1e-3  
                if "FaultClosedResistance" in event_params.keys():
                    closed_resistance = event_params['FaultClosedResistance']
                else:
                    warnings.warn('Paramenter "FaultClosedResistance" is not in the json file!\n FaultClosedResistance set to 1e-3')
                    
                # Füge switch mit Erdung hinzu
                switch = self.dpsimpy_components.Switch('Fault_' + node_name, self.loglevel)
                switch.set_parameters(open_resistance, closed_resistance)
                switch.open()
                switch.connect([self.dpsimpy_components.SimNode.gnd, node])
                self.system.add(switch)

                # Event hinzufügen
                sw_event_1 = dpsimpy.event.SwitchEvent(1.0, switch, True)
                self.sim.add_event(sw_event_1)
                sw_event_2 = dpsimpy.event.SwitchEvent(1.1, switch, False)
                self.sim.add_event(sw_event_2)
                
            else:
                raise Exception("Bitte gebe einen gültigen Fehlerfall an. Zum Beispiel <ThreePhaseFault>.")

        """
        elif event['EventType'] == "LineDisconnection":
            event_param = event['EventParameters']

            for comp in system.components:
                if event_param['LineName'] == comp.name:
                    # comp ist die zu entfernende Freileitung
                    # Entferne Freileitung aus Komponentenliste    
                    system.components = [x for x in system.components if x != comp]

                    # Entferne Freileitung aus den Knoten-Listen
                    for key, value in system.components_at_node.items():
                        if comp in value:
                            obj_dict = system.components_at_node
                            obj_dict[key] = [c for c in obj_dict[key] if c != comp]
                            system.components_at_node = obj_dict

                    ### DPsim SP simulation
                    sim_parameters = data['SimulationParameters']
                    name = sim_parameters['SimulationName']
                    dpsimpy.Logger.set_log_dir("logs/" + name)

                    ### Simulation
                    sim = dpsimpy.Simulation(name, dpsimpy.LogLevel.debug)
                    sim.set_system(system)
                    sim.set_domain(dpsimpy.Domain.SP)

                    sim.do_init_from_nodes_and_terminals(True)
                    sim.set_direct_solver_implementation(dpsimpy.DirectLinearSolverImpl.SparseLU)
                    sim.do_system_matrix_recomputation(True)

                    # Initialisiere Spannungen an Knoten vom dynamischen System
                    self.system.init_with_powerflow(self.system_pf, domain=self.domain)

                    logger = dpsimpy.Logger(name)
                    sim.add_logger(logger)
                    for node_i in system.nodes:
                        logger.log_attribute(node_i.name + '.V', 'v', node_i)

                    sim.set_time_step(sim_parameters['TimeStep'])
                    sim.set_final_time(sim_parameters['EndTime'])

                    sim.run()

                    return system

            print(f"Die gewünschte Freileitung '{event_param['LineName']}' ist nicht in der Topologie enthalten.")
            return system
        """
   


        