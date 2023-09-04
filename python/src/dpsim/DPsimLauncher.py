import sys
sys.path.insert(0,'/home/mmo-cya/dpsim/build')
import dpsimpy
import json
import cimpy
from cimpy import CIM2DPsim
from cimpy.CIM2DPsim import Domain
from cimpy.CIM2DPsim import SGModels


def DPsimLauncher(path = '/home/mmo-cya/dpsim/python/src/dpsim/Config_network.json' ):
    """
    Diese Funktion liest eine Konfigurationsdatei im JSON-Format ein, in der der Dateipfad zu dem CIM-Netzwerk und deren Parameter zur Simulation
    mittels DPsim enthalten sind. Das Netzwerk wird eingelesen und ins DPsim-Format übersetzt. Der ausgewählte Fehler wird dem Netzwerk hinzugefügt und anschließend wird 
    die Simulation ausgeführt. Ausgegeben werden das neue DPsim Netzwerk und die Simulationsergebnisse mit dem gewünschten Solver.
    """
    # Passe den Dateipfad entsprechend an
    file_path = path

    # Einlesen der Konfigurations Datei im JSON Format
    try:
        with open(file_path, 'r') as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        print(f"Die Datei '{file_path}' wurde nicht gefunden.")
    except json.JSONDecodeError:
        print(f"Die Datei '{file_path}' ist keine gültige JSON-Datei.")

    # Importieren der Netz-Topologie im cimpy-Format
    imported_result = cimpy.cim_import(data['CIMFiles'], 'cgmes_v2_4_15')

    # Übersetzen der cimpy-Topologie in die DPsim-Topologie für die PowerFlow Simulation
    system_PF = CIM2DPsim.CIM2DPsim(imported_result, Domain.PF, log_level=dpsimpy.LogLevel.info)

    # Führe PowerFlow Simulation aus um die Anfangswerte für die dynamische Simulation zu bestimmen
    name_pf = data['SimulationParameters']['SimulationName'] + "_PF"
    dpsimpy.Logger.set_log_dir("logs/" + name_pf)

    #set slack 
    # system_PF.component("Gen_0001").modify_power_flow_bus_type(dpsimpy.PowerflowBusType.VD)
    sim_pf = dpsimpy.Simulation(name_pf, dpsimpy.LogLevel.debug)
    sim_pf.set_system(system_PF)
    sim_pf.set_domain(dpsimpy.Domain.SP)
    sim_pf.set_solver(dpsimpy.Solver.NRP)
    #sim_pf.set_solver_component_behaviour(dpsimpy.SolverBehaviour.Initialization)
    sim_pf.do_init_from_nodes_and_terminals(False)

    logger = dpsimpy.Logger(name_pf)
    sim_pf.add_logger(logger)
    for node in system_PF.nodes:
        logger.log_attribute(node.name+'.V', 'v', node)
        logger.log_attribute(node.name+'.S', 's', node)
        
    sim_pf.set_time_step(0.1)
    sim_pf.set_final_time(0.5)
    sim_pf.run()


   
    # Übersetzen der cimpy-Topologie in die DPsim-Topologie für die dynamische Simulation
    system = CIM2DPsim.CIM2DPsim(imported_result, Domain.SP, log_level=dpsimpy.LogLevel.debug, gen_model= SGModels.VBR_6Order)

    # Füge der Topologie den gewünschten Event hinzu
    event = data['Events'] 
    if event['EventType'] == "ThreePhaseFault":
        event_param = event['EventParameters']

        for node in system.nodes:
            if event_param['NodeName'] == node.name:
                # Knoten ist in der Topologie enthalten und Switch kann an diesen Knoten hinzugefügt werden
                # Füge switch mit Erdung hinzu
                gnd = dpsimpy.sp.SimNode.gnd
                # Switch
                switch = dpsimpy.sp.ph1.Switch('Fault', dpsimpy.LogLevel.debug)
                switch.set_parameters(event_param['FaultOpenResistance'], event_param['FaultClosedResistance'])
                switch.open()
                switch.connect([gnd, node])
                system.add(switch)
                
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
                system.init_with_powerflow(system_PF)
                
                logger = dpsimpy.Logger(name)
                sim.add_logger(logger)
                for node_i in system.nodes:
                    logger.log_attribute(node_i.name+'.V', 'v', node_i)

                sim.set_time_step(sim_parameters['TimeStep'])
                sim.set_final_time(sim_parameters['EndTime'])

                # Event hinzufügen
                sw_event_1 = dpsimpy.event.SwitchEvent(1.0, switch, True)
                sim.add_event(sw_event_1)

                sw_event_2 = dpsimpy.event.SwitchEvent(1.1, switch, False)
                sim.add_event(sw_event_2)
              
                sim.run()
            
                return system
        
        print(f"Der gewünschte Knoten '{event_param['NodeName']}' ist nicht in der Topologie enthalten.")
        return system
    
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
                system.init_with_powerflow(system_PF)
                
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
    
    print("Bitte gebe einen gültigen Fehlerfall an. Zum Beispiel <ThreePhaseFault>.")
    return system

   

