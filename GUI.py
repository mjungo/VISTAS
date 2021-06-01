import json
import PySimpleGUI as sg

from VISTAS1D import VISTAS1D


# GUI/params mgmt: check the type and correct the unit (ns->s, mA->A) to make sure correct values are stored in the dictionaries sp and vp
def check_value(k, v):

    if k in ['SCHTransp', 'ThermMod', 'Noise', 'Parasitics', '2D', 'storeN', 'Uxyplot', 'PIplot', 'Ptplot', 'NSxyplot', 'RINplot', 'Hfplot', 'eyeplot']:
        try:
            v = bool(v)
        except:
            sg.popup('Bool input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)

    elif k in ['odeSolver', 'modFormat', 'vcselDescr']:
        try:
            v = str(v)
        except:
            sg.popup('String input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)

    elif k == 'ni':
        try:
            v = int(v)
        except:
            sg.popup('Integer input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)
    else:
        try:
            v = float(v)
            if k in ['tmax', 'dt', 'dtFD', 'tb']:   # converts ns to s
                v = round(v * 1e-9, 13)
            elif k in ['Ion', 'Ioff', 'Iss']:       # converts mA to A
                v = round(v * 1e-3, 5)
        except:
            sg.popup('Float input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)

    return k, v
 

# GUI/params mgmt: read file "file_name.json", convert to dictionary dict of dicts params, and split into individual dictionaries sp and vp
def load_params(file_name):
    with open(file_name, 'r') as openfile:
        d = json.load(openfile)    # dictionary

    return d


# GUI/params mgmt: aggregate dictionaries sp and vp into dict of dicts params, convert to json, and save to "file_name.json"
def save_params(sp, vp, sr, file_name):
    params = {'simParams': sp, 'vcselParams': vp, 'simResults': sr}    # dict of dict
    with open(file_name, 'w') as outfile:
        json.dump(params, outfile, indent=4)


# GUI/params mgmt: update dictionaries sp and vp on event, check type, adjust units
def update_dict(values, d):
    for k, v in d.items():             # update entire dictionary
        k, d[k] = check_value(k, values[k])

    return d


# GUI/params mgmt: disable elements not used for certain combinations of parameters
def disable_elements(window, values):

    # simulation parameters
    if values['odeSolver'] == 'Finite Diff.':
        window['dtFD'].Update(disabled = False)
    else:
        window['dtFD'].Update(disabled = True)

    # modulation pattern
    if values['modFormat'] == 'step':
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)
        window['eyeplot'].Update(disabled = True)
        window['Hfplot'].Update(disabled = True)
    
    elif values['modFormat'] == 'pulse':
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)
        window['eyeplot'].Update(disabled = True)
        window['Hfplot'].Update(disabled = True)

    elif values['modFormat'] == 'random bits':
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = False)
        window['eyeplot'].Update(disabled = False)
        window['Hfplot'].Update(disabled = True)

    elif values['modFormat'] == 'small signal':
        window['tmax'].Update(disabled = True)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = False)
        window['tb'].Update(disabled = True)
        window['eyeplot'].Update(disabled = True)
        window['Hfplot'].Update(disabled = False)

    # simulated effects
    if values['Noise'] == True and values['modFormat'] == 'step':
        window['RINplot'].Update(disabled = False)
    else:
        window['RINplot'].Update(disabled = True)


# GUI/params mgmt: update gui (incl. disabling elements) after initializing or loading params file
def update_gui(window, values, sp, vp):
    for k, v in sp.items():             # update entire simParams dictionary
        if k in ['tmax', 'dt', 'dtFD', 'tb']:
            values[k] = sp[k] * 1e9     # convert s to ns
        elif k in ['Ion', 'Ioff', 'Iss']:
            values[k] = sp[k] * 1e3     # convert A to mA
        else:
            values[k] = sp[k]
        window[k].update(values[k])
    for k, v in vp.items():             # update entire vcselParams dictionary
        values[k] = vp[k]
        window[k].update(values[k])

    disable_elements(window, values)

    return values


# GUI/params mgmt: main function, incl. layout and events loop 
def GUI():

    # 0. load last parameters into sp and vp dictionaries -------------------------------------------------------
    params = load_params('last_params.json')  # loads params file and populates dictionaries sp and vp
    sp, vp = params['simParams'], params['vcselParams'] # simulation results not loaded (for that, chose "load params" and select any file but last_params.json or default_params.json)
    sr = {}
    del params

    ttips = load_params('ttips.json')  # loads params file and populates dictionaries sp and vp

    # 1. Tabs layout -------------------------------------------------------
    spTab_layout = [
        [sg.Frame('Simulated effects',[
            [sg.Checkbox('carrier transport into the quantum wells', key='SCHTransp', default=sp['SCHTransp'], size=(45,1), disabled=True, tooltip=ttips['SCHTransp'], enable_events=True)],
            [sg.Checkbox('thermal effects', key='ThermMod', default=sp['ThermMod'], disabled=True, tooltip=ttips['ThermMod'], enable_events=True)],
            [sg.Checkbox('noise', key='Noise', default=sp['Noise'], disabled=True, tooltip=ttips['Noise'], enable_events=True)],
            [sg.Checkbox('electrical parasitics', key='Parasitics', default=sp['Parasitics'], disabled=True, tooltip=ttips['Parasitics'], enable_events=True)],
            [sg.Checkbox('2D', key='2D', default=sp['2D'], disabled=True, tooltip=ttips['2D'], enable_events=True)],
            ])],
        [sg.Frame('Simulation parameters',[
            [sg.InputText(sp['ni'], key='ni', size=(7,1), tooltip=ttips['ni'], enable_events=True), sg.Text('radial resolution (typically 7-20 terms)', size=(40,1))],
            [sg.Checkbox('store carrier terms Ni', key='storeN', default=sp['storeN'], disabled=True, tooltip=ttips['storeN'], enable_events=True)],
            [sg.Combo(values=('RK45', 'RK23', 'DOP853', 'Radau', 'LSODA', 'BDF', 'Finite Diff.'), key='odeSolver', default_value=sp['odeSolver'], size=(12,1), tooltip=ttips['odeSolver'], enable_events=True), sg.Text('ODE solver')],
            [sg.InputText(round(float(0 if sp['tmax'] is None else sp['tmax'])*1e9,0), key='tmax', size=(7,1), tooltip=ttips['tmax'], enable_events=True, readonly=(sp['modFormat']=='small signal'), disabled_readonly_background_color='grey'), sg.Text('simulated time (ns)')],
            [sg.InputText(round(float(0 if sp['dt'] is None else sp['dt'])*1e9,3), key='dt', size=(7,1), tooltip=ttips['dt'], enable_events=True), sg.Text('time resolution for storing and plotting (ns)')],
            [sg.InputText(round(float(0 if sp['dtFD'] is None else sp['dtFD'])*1e9,3), key='dtFD', size=(7,1), tooltip=ttips['dtFD'], enable_events=True, readonly=(sp['odeSolver']!='Finite Diff.'), disabled_readonly_background_color='grey'), sg.Text('time step for FD solution (ns)')],
            ])],
        [sg.Frame('Current modulation pattern',[
            [sg.Combo(values=('step', 'pulse', 'random bits', 'small signal'), key='modFormat', default_value=sp['modFormat'], size=(12,1), tooltip=ttips['modFormat'], enable_events=True), sg.Text('modulation pattern')],
            [sg.InputText(round(float(0 if sp['Ion'] is None else sp['Ion'])*1e3,1), key='Ion', size=(7,1), tooltip=ttips['Ion'], enable_events=True), sg.Text('ON current (mA)', size=(40,1))],
            [sg.InputText(round(float(0 if sp['Ioff'] is None else sp['Ioff'])*1e3,1), key='Ioff', size=(7,1), tooltip=ttips['Ioff'], enable_events=True, readonly=(sp['modFormat']=='small signal' or sp['modFormat']=='step'), disabled_readonly_background_color='grey'), sg.Text('OFF current (mA)')],
            [sg.InputText(round(float(0 if sp['Iss'] is None else sp['Iss'])*1e3,2), key='Iss', size=(7,1), tooltip=ttips['Iss'], enable_events=True, readonly=(sp['modFormat']!='small signal'), disabled_readonly_background_color='grey'), sg.Text('small signal current step (mA)')],
            [sg.InputText(round(float(0 if sp['tb'] is None else sp['tb'])*1e9,1), key='tb', size=(7,1), tooltip=ttips['tb'], enable_events=True, disabled=(sp['modFormat']!='random bits'), disabled_readonly_background_color='grey'), sg.Text('bit duration (ns)')],
            ])],
        [sg.Frame('Results visualization',[
            [sg.Checkbox('2D mode profiles', key='Uxyplot', default=sp['Uxyplot'], size=(45,1), tooltip=ttips['Uxyplot'], enable_events=True)],
            [sg.Checkbox('Steady-state (LI) characteristic Popt(I)', key='PIplot', default=sp['PIplot'], size=(45,1), tooltip=ttips['PIplot'], enable_events=True)],
            [sg.Checkbox('Dynamic characteristic Popt(t)', key='Ptplot', default=sp['Ptplot'], tooltip=ttips['Ptplot'], enable_events=True)],
            [sg.Checkbox('2D optical & carrier profiles within the cavity', key='NSxyplot', default=sp['NSxyplot'], disabled=True, size=(45,1), tooltip=ttips['NSxyplot'], enable_events=True)],
            [sg.Checkbox('Relative Intensity Noise (RIN) spectrum', key='RINplot', default=sp['RINplot'], disabled=(sp['Noise']!=True), tooltip=ttips['RINplot'], enable_events=True)],
            [sg.Checkbox('Frequency response H(f)', key='Hfplot', default=sp['Hfplot'], disabled=(sp['modFormat']!='small signal'), tooltip=ttips['Hfplot'], enable_events=True)],
            [sg.Checkbox('Eye diagram', key='eyeplot', default=sp['eyeplot'], disabled=(sp['modFormat']!='random bits'), tooltip=ttips['eyeplot'], enable_events=True)],
            ])],
        ]

    vpTab_layout = [
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

    spTabCol_layout = [[sg.Column(spTab_layout, size = (416,754), scrollable=True, vertical_scroll_only=True)]]
    vpTabCol_layout = [[sg.Column(vpTab_layout, size = (416,754), scrollable=True, vertical_scroll_only=True)]]

    layout = [[sg.TabGroup([[sg.Tab('Simulation parameters', spTabCol_layout), sg.Tab('VCSEL parameters', vpTabCol_layout)]])],    
            [sg.Button('initialize params', size=(12,1)), sg.Input(key='load file', enable_events=True, visible=False), sg.FileBrowse('load file', size=(12,1), file_types=(('JSON files', '*.json'),), target='load file'), sg.Input(key='save file', enable_events=True, visible=False), sg.SaveAs('save file', size=(12,1), file_types=(('JSON files', '*.json'),), target='save file'), sg.Button('run simulation', size=(12,1))]]

    # 2. window -------------------------------------------------------
    window = sg.Window('VISTAS', layout)

    # 3. events loop -------------------------------------------------------
    while True:
        
        event, values = window.read()
        #print(event)

        if event == sg.WIN_CLOSED:
            save_params(sp, vp, {}, 'last_params.json')         # aggregates sp and vp in a dict of dicts and saves to a json file
            break

        elif event == 'initialize params':
            params = load_params('default_params.json')         # loads params file and populates dictionaries sp and vp
            sp, vp = params['simParams'], params['vcselParams'] # simulation results not loaded
            sr = {}
            del params
            save_params(sp, vp, {}, 'last_params.json')         # simulation results not saved to "last_params.json"
            update_gui(window, values, sp, vp)

        elif event == 'load file':
            file_name = values['load file']
            if file_name != '':
                params = load_params(file_name)             # loads params file and populates dictionaries sp and vp
                sp, vp, sr = params['simParams'], params['vcselParams'], params['simResults']
                # for post-processing, simulations results should be converted from list to numpy array and saved in the respective variables (S, N, ur, etc.)
                # for k, v in sr.items(): # extracts the dictionary items to variables named after the key
                #     exec("%s = %d" % (k, v))
                #     exec("%s = ['%s']" % (k, v))
                del params
                update_gui(window, values, sp, vp)   

        elif event == 'save file':
            file_name = values['save file']
            if file_name == '': file_name = 'last_params.json'
            save_params(sp, vp, sr, file_name)  # aggregates sp and vp and sr in a dict of dicts and saves to a json file

        elif event == 'run simulation':
            save_params(sp, vp, {}, 'last_params.json') # simulation results not saved to "last_params.json"
            sr = VISTAS1D(sp, vp)

        else:   # case 'gui-element change'
            sp = update_dict(values, sp)        # update dictionary
            vp = update_dict(values, vp)        # update dictionary
            disable_elements(window, values)    # update GUI


def main():
    GUI()


if __name__ == "__main__":
    main()