import json
import PySimpleGUI as sg


# check the type and correct the unit (ns->s,) to make sure correct values are stored in the dictionaries sp and vp
def check_value(k, v):
    if k in ['SCHTransp', 'ThermMod', 'Noise', 'Parasitics', '2D', 'storeN']:
        try:
            v = bool(v)
        except:
            sg.popup('boolean expected')
    elif k in ['odeSolver', 'modFormat', 'vcselDescr']:
        try:
            v = str(v)
        except:
            sg.popup('str expected')
    elif k == 'ni':
        try:
            v = int(v)
        except:
            sg.popup('int expected')
    else:
        try:
            v = float(v)
            if k in ['tmax', 'dt', 'dtFD', 'tb']:   # converts ns to s
                v = round(v * 1e-9, 13)
            elif k in ['Ion', 'Ioff', 'Iss']:       # converts mA to A
                v = round(v * 1e-3, 5)
        except:
            sg.popup('float expected')

    return k, v
 

# read file "file_name.json", convert to dictionary dict of dicts params, and split into individual dictionaries sp and vp
def load_params(file_name):
    with open(file_name, 'r') as openfile:
        d = json.load(openfile)    # dictionary

    return d


# aggregate dictionaries sp and vp into dict of dicts params, convert to json, and save to "file_name.json"
def save_params(sp, vp, file_name):
    params = {'simParams': sp, 'vcselParams': vp}    # dict of dict
    with open(file_name, 'w') as outfile:
        json.dump(params, outfile, indent=4)


# update dictionaries sp and vp on event, check type, adjust units
def update_dict(values, d):
    for k, v in d.items():             # update entire dictionary
        k, d[k] = check_value(k, values[k])

    return d


# disable elements not used for certain parameters
def disable_elements(window, values):
    window['Noise'].Update(disabled = False)  

    # simulation parameters
    if values['odeSolver'] == 'Finite Diff.':
        window['dtFD'].Update(disabled = False)
    else:
        window['dtFD'].Update(disabled = True)
    # modulation pattern
    if values['modFormat'] == 'random bits':
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = False)
        window['eye'].Update(disabled = False)
        window['MTF'].Update(disabled = True)

    elif values['modFormat'] == 'small signal':
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = False)
        window['tb'].Update(disabled = True)
        window['eye'].Update(disabled = True)
        window['MTF'].Update(disabled = False)

    else:
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)
        window['eye'].Update(disabled = True)
        window['MTF'].Update(disabled = True)

    # simulated effects
    if values['Noise'] == True and values['modFormat'] == 'step':
        window['RIN'].Update(disabled = False)
    else:
        window['RIN'].Update(disabled = True)



# update gui (incl. disabling elements) after initializing or loading params file
def update_gui(window, values, sp, vp):
    for k, v in sp.items():             # update entire simParams dictionary
        if k in ['tmax', 'dt', 'dtFD', 'tb']:
            values[k] = sp[k]*1e9 # convert s to ns
        elif k in ['Ion', 'Ioff', 'Iss']:
            values[k] = sp[k]*1e3 # convert A to mA
        else:
            values[k] = sp[k]
        window[k].update(values[k])
    for k, v in vp.items():             # update entire vcselParams dictionary
        values[k] = vp[k]
        window[k].update(values[k])

    disable_elements(window, values)

    return values


def main():

    # 0. load last parameters into sp and vp dictionaries -------------------------------------------------------

    params = load_params('last_params.json')  # loads params file and populates dictionaries sp and vp
    sp, vp = params['simParams'], params['vcselParams']
    del params

    ttips = load_params('ttips.json')  # loads params file and populates dictionaries sp and vp


    # 1. Tabs layout -------------------------------------------------------

    tab1_layout = [
        [sg.Frame('Simulated effects',[
            [sg.Checkbox('carrier transport into the quantum wells', key='SCHTransp', default=sp['SCHTransp'], size=(45,1), tooltip=ttips['SCHTransp'], enable_events=True)],
            [sg.Checkbox('thermal effects', key='ThermMod', default=sp['ThermMod'], disabled=True, tooltip=ttips['ThermMod'], enable_events=True)],
            [sg.Checkbox('noise', key='Noise', default=sp['Noise'], disabled=True, tooltip=ttips['Noise'], enable_events=True)],
            [sg.Checkbox('electrical parasitics', key='Parasitics', default=sp['Parasitics'], disabled=True, tooltip=ttips['Parasitics'], enable_events=True)],
            [sg.Checkbox('2D', key='2D', default=sp['2D'], disabled=True, tooltip=ttips['2D'], enable_events=True)],
            ])],
        [sg.Frame('Simulation parameters',[
            [sg.InputText(sp['ni'], key='ni', size=(7,1), tooltip=ttips['ni'], enable_events=True), sg.Text('radial resolution (typically 7-20 terms)', size=(40,1))],
            [sg.Checkbox('store carrier terms Ni', key='storeN', default=sp['storeN'], disabled=True, tooltip=ttips['storeN'], enable_events=True)],
            [sg.Combo(values=('RK45', 'RK23', 'DOP853', 'Radau', 'LSODA', 'BDF', 'Finite Diff.'), key='odeSolver', default_value=sp['odeSolver'], size=(12,1), tooltip=ttips['odeSolver'], enable_events=True), sg.Text('ODE solver')],
            [sg.InputText(round(sp['tmax']*1e9,0), key='tmax', size=(7,1), tooltip=ttips['tmax'], enable_events=True), sg.Text('simulated time (ns)')],
            [sg.InputText(round(sp['dt']*1e9,3), key='dt', size=(7,1), tooltip=ttips['dt'], enable_events=True), sg.Text('time resolution for storing and plotting (ns)')],
            [sg.InputText(round(sp['dtFD']*1e9,3), key='dtFD', size=(7,1), tooltip=ttips['dtFD'], enable_events=True, readonly=True, disabled_readonly_background_color='grey'), sg.Text('time step for FD solution (ns)')],
            ])],
        [sg.Frame('Current modulation pattern',[
            [sg.Combo(values=('step', 'pulse', 'random bits', 'small signal'), key='modFormat', default_value=sp['modFormat'], size=(12,1), tooltip=ttips['modFormat'], enable_events=True), sg.Text('modulation pattern')],
            [sg.InputText(round(sp['Ion']*1e3,1), key='Ion', size=(7,1), tooltip=ttips['Ion'], enable_events=True), sg.Text('ON current (mA)', size=(40,1))],
            [sg.InputText(round(sp['Ioff']*1e3,1), key='Ioff', size=(7,1), tooltip=ttips['Ioff'], enable_events=True, disabled_readonly_background_color='grey'), sg.Text('OFF current (mA)')],
            [sg.InputText(round(sp['Iss']*1e3,2), key='Iss', size=(7,1), tooltip=ttips['Iss'], enable_events=True, readonly=True, disabled_readonly_background_color='grey'), sg.Text('small signal current step (mA)')],
            [sg.InputText(round(sp['tb']*1e9,1), key='tb', size=(7,1), tooltip=ttips['tb'], enable_events=True, readonly=True, disabled_readonly_background_color='grey'), sg.Text('bit duration (ns)')],
            ])],
        [sg.Frame('Visualization',[
            [sg.Checkbox('2D mode profiles', key='2Dmodes', default=sp['2Dmodes'], size=(45,1), tooltip=ttips['2Dmodes'], enable_events=True)],
            [sg.Checkbox('Steady-state (LI) characteristic Popt(I)', key='LIplot', default=sp['LIplot'], size=(45,1), tooltip=ttips['LIplot'], enable_events=True)],
            [sg.Checkbox('Dynamic characteristic Popt(t)', key='Ptplot', default=sp['Ptplot'], tooltip=ttips['Ptplot'], enable_events=True)],
            [sg.Checkbox('2D optical & carrier profiles within the cavity', key='2Dprofiles', disabled=True, default=sp['2Dprofiles'], size=(45,1), tooltip=ttips['2Dprofiles'], enable_events=True)],
            [sg.Checkbox('Relative Intensity Noise (RIN) spectrum', key='RIN', default=sp['RIN'], disabled=True, tooltip=ttips['RIN'], enable_events=True)],
            [sg.Checkbox('Modulation Transfer Function (MTF)', key='MTF', default=sp['MTF'], disabled=True, tooltip=ttips['MTF'], enable_events=True)],
            [sg.Checkbox('Eye diagram', key='eye', default=sp['eye'], disabled=True, tooltip=ttips['eye'], enable_events=True)],
            ])],
        ]

    tab2_layout = [
        [sg.Frame('Generic description',[
            [sg.Multiline(vp['vcselDescr'], key='vcselDescr', size=(52,2), tooltip=ttips['vcselDescr'], enable_events=True)],
            ])],
        [sg.Frame('Cavity parameters',[
            [sg.InputText(vp['r_ox'], key='r_ox', size=(7,1), tooltip=ttips['r_ox'], enable_events=True), sg.Text('r_ox: oxide aperture radius (cm)', size=(40,1))],
            [sg.InputText(vp['Leff'], key='Leff', size=(7,1), tooltip=ttips['Leff'], enable_events=True), sg.Text('Leff: effective cavity length (cm)')],
            [sg.InputText(vp['nqw'], key='nqw', size=(7,1), tooltip=ttips['nqw'], enable_events=True), sg.Text('nqw: number of quantum wells')],
            [sg.InputText(vp['dqw'], key='dqw', size=(7,1), tooltip=ttips['dqw'], enable_events=True), sg.Text('dqw: single QW thickness (cm)')],
            [sg.InputText(vp['db'], key='db', size=(7,1), tooltip=ttips['db'], enable_events=True), sg.Text('db: SCH thickness (cm)')],
            ])],
        [sg.Frame('Equivalent waveguide parameters',[
            [sg.InputText(vp['wl0'], key='wl0', size=(7,1), tooltip=ttips['wl0'], enable_events=True), sg.Text('wl0: emission wavelength @300k (nm)', size=(40,1))],
            [sg.InputText(vp['nc'], key='nc', size=(7,1), tooltip=ttips['nc'], enable_events=True), sg.Text('nc: core equivalent refractive index')],
            [sg.InputText(vp['ng'], key='ng', size=(7,1), tooltip=ttips['ng'], enable_events=True), sg.Text('ng: group refractive index')],
            [sg.InputText(vp['delta_n'], key='delta_n', size=(7,1), tooltip=ttips['delta_n'], enable_events=True), sg.Text('delta_n: equivalent fractional refractive index change')],
            ])],
        [sg.Frame('General physical parameters',[
            [sg.InputText(vp['Rt'], key='Rt', size=(7,1), tooltip=ttips['Rt'], enable_events=True), sg.Text('Rt: top mirror reflectivity', size=(40,1))],
            [sg.InputText(vp['Rb'], key='Rb', size=(7,1), tooltip=ttips['Rb'], enable_events=True), sg.Text('Rb: bottom mirror reflectivity')],
            [sg.InputText(vp['alpha_i'], key='alpha_i', size=(7,1), tooltip=ttips['alpha_i'], enable_events=True), sg.Text('alpha_i: internal losses (cm-1)')],
            [sg.InputText(vp['tau_N'], key='tau_N', size=(7,1), tooltip=ttips['tau_N'], enable_events=True), sg.Text('tau_N: carrier-lifetime (s)')],
            [sg.InputText(vp['beta'], key='beta', size=(7,1), tooltip=ttips['beta'], enable_events=True), sg.Text('beta: spontaneous recombination coefficient')],
            ])],
        [sg.Frame('Gain parameters',[
            [sg.InputText(vp['gln'], key='gln', size=(7,1), tooltip=ttips['gln'], enable_events=True), sg.Text('gln: logarithmic gain coefficient (cm-1)', size=(40,1))],
            [sg.InputText(vp['Ntr'], key='Ntr', size=(7,1), tooltip=ttips['Ntr'], enable_events=True), sg.Text('Ntr: transparency carrier density (cm-3)')],
            [sg.InputText(vp['epsilon'], key='epsilon', size=(7,1), tooltip=ttips['epsilon'], enable_events=True), sg.Text('epsilon: gain compression factor')],
            [sg.InputText(vp['Gam_r'], key='Gam_r', size=(7,1), tooltip=ttips['Gam_r'], enable_events=True), sg.Text('Gam_r: relative confinement factor')],
            ])],
        [sg.Frame('Carrier transport and diffusion parameters',[
            [sg.InputText(vp['tau_esc'], key='tau_esc', size=(7,1), tooltip=ttips['tau_esc'], enable_events=True), sg.Text('tau_esc: thermionic emission lifetime (s)', size=(40,1))],
            [sg.InputText(vp['tau_cap'], key='tau_cap', size=(7,1), tooltip=ttips['tau_cap'], enable_events=True), sg.Text('tau_cap: ambipolar diffusion time (s)')],
            [sg.InputText(vp['eta_i'], key='eta_i', size=(7,1), tooltip=ttips['eta_i'], enable_events=True), sg.Text('eta_i: current injection efficiency')],
            [sg.InputText(vp['rs'], key='rs', size=(7,1), tooltip=ttips['rs'], enable_events=True), sg.Text('rs: current spreading coefficient (cm)')],
            [sg.InputText(vp['DN'], key='DN', size=(7,1), tooltip=ttips['DN'], enable_events=True), sg.Text('DN: ambipolar diffusion coeff. (cm2/s)')],
            ])],
        ]

    tab1Col_layout = [[sg.Column(tab1_layout, size = (416,754), scrollable=True, vertical_scroll_only=True)]]
    tab2Col_layout = [[sg.Column(tab2_layout, size = (416,754), scrollable=True, vertical_scroll_only=True)]]

    layout = [[sg.TabGroup([[sg.Tab('Simulation parameters', tab1Col_layout), sg.Tab('VCSEL parameters', tab2Col_layout)]])],    
            [sg.Button('initialize params', size=(12,1)), sg.Input(key='load file', enable_events=True, visible=False), sg.FileBrowse('load params', size=(12,1), file_types=(('JSON files', '*.json'),), target='load file'), sg.Input(key='save file', enable_events=True, visible=False), sg.SaveAs('save params as', size=(12,1), file_types=(('JSON files', '*.json'),), target='save file'), sg.Button('run simulation', size=(12,1))]]    


    # 2. window -------------------------------------------------------

    window = sg.Window('VISTAS', layout)

    # 3. events loop -------------------------------------------------------

    while True:
        
        event, values = window.read()
        print(event)

        if event == sg.WIN_CLOSED:
            save_params(sp, vp, 'last_params.json') # aggregates sp and vp in a dict of dicts and saves to a json file
            break

        elif event == 'initialize params':
            params = load_params('default_params.json')  # loads params file and populates dictionaries sp and vp
            sp, vp = params['simParams'], params['vcselParams']
            del params
            save_params(sp, vp, 'last_params.json')
            update_gui(window, values, sp, vp)

        elif event == 'load file':
            file_name = values['load file']
            if file_name != '':
                params = load_params(file_name)  # loads params file and populates dictionaries sp and vp
                sp, vp = params['simParams'], params['vcselParams']
                del params
                update_gui(window, values, sp, vp)   

        elif event == 'save file':
            file_name = values['save file']
            if file_name == '': file_name = 'last_params.json'
            save_params(sp, vp, file_name) # aggregates sp and vp in a dict of dicts and saves to a json file

        elif event == 'run simulation':
            print('RUN')
            save_params(sp, vp, 'last_params.json') # aggregates sp and vp in a dict of dicts and saves to a json file
            break

        else:   # case 'gui-element change'
            sp = update_dict(values, sp)        # update dictionary
            vp = update_dict(values, vp)        # update dictionary
            disable_elements(window, values)    # update GUI

if __name__ == "__main__":
    main()