import sys
sys.path.insert(0,'/home/mmo-cya/dpsim/build')
import dpsimpy
import json
import cimpy
from cimpy import CIMpyToDPsim
from cimpy.CIMpyToDPsim import Domain

def DPsimLauncher(path = '/home/mmo-cya/dpsim/python/src/dpsim/Config_network.json' ):
    """
    Diese Funktion liest eine Konfiguartions Datei im JSON Format ein, in der der Dateipfad zu dem CIM-Netzwerk und deren Parameter zur Simulation
    mittels DPsim enthalten sind. Das Netzwerk wird eingelesen und im DPsim Format übersetzt. Der ausgewählte Fehler wird dem Netzwerk hinzugefügt und anschließend wird 
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

    # Übersetzen der cimpy-Topologie in die DPsim-Topologie
    system , K = CIMpyToDPsim.CIMpyToDPsim(imported_result, Domain.PF)
   
    # Füge der Topologie den gewünschten Event hinzu
    event = data['Events'] 
    if event['EventType'] == "ThreePhaseFault":
        event_param = event['EventParameters']

        for node in system.nodes:
            if event_param['NodeName'] == node.name:
                # Knoten ist in der Topologie enthalten und Switch kann an diesen Knoten hinzugefügt werden
                # Füge switch mit Erdung hinzu
                print("HERE")
                gnd = dpsimpy.sp.SimNode.gnd
                # Switch
                switch = dpsimpy.sp.ph1.Switch('Fault', dpsimpy.LogLevel.debug)
                switch.set_parameters(event_param['FaultOpenResistance'], event_param['FaultClosedResistance'])
                switch.open()

                #switch.connect([gnd, node])
                #system.add(switch)
                
                ### DPsim SP simulation
                name = "IEEE14"
                dpsimpy.Logger.set_log_dir("logs/" + name)

                    ### Simulation
                sim_parameters = data['SimulationParameters']
                sim = dpsimpy.Simulation(name, dpsimpy.LogLevel.debug)
                sim.set_system(system)
                sim.set_domain(dpsimpy.Domain.SP)
                sim.set_solver(dpsimpy.Solver.NRP)
                sim.set_solver_component_behaviour(dpsimpy.SolverBehaviour.Initialization)
                sim.set_time_step(sim_parameters['TimeStep'])
                sim.set_final_time(sim_parameters['EndTime'])
                
                sw_event_1 = dpsimpy.event.SwitchEvent(1.0, switch, True)
                sim.add_event(sw_event_1)

                sw_event_2 = dpsimpy.event.SwitchEvent(1.1, switch, False)
                sim.add_event(sw_event_2)
              
                sim.run()
            
                return system
        
        print(f"Der gewünschte Knoten '{event_param['NodeName']}' ist nicht in der Topologie enthalten.")
        return system
    
    print("Bitte gebe einen gültigen Fehlerfall an. Zum Beispiel <ThreePhaseFault>.")
    return system

   

