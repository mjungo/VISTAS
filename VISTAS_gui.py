##################################################################################
# MIT License
#
# Copyright (c) 2021 Marc Jungo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##################################################################################

import json
import PySimpleGUI as sg

from VISTAS_algorithm import VISTAS1D
from VISTAS_visualization import *


# GUI/params mgmt: check the type and correct the unit (ns->s, mA->A) to make sure correct values are stored in the dictionaries sp and vp
def check_value(k, v):

    if k in ['SCHtransp', 'ThermMod', 'Noise', 'Parasitics', '2D', 'Modes2Dplot', 'PIplot', 'Ptplot', 'NwS2Dplot', 'RINplot', 'Hfplot', 'Eyeplot']:
        try:
            v = bool(v)
        except:
            #sg.popup('Bool input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)

    elif k in ['odeSolver', 'modFormat', 'vcselDescr']:
        try:
            v = str(v)
        except:
            #sg.popup('String input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)

    elif k == 'nNw':
        try:
            v = int(v)
        except:
            #sg.popup('Integer input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
            v = None    # in case last_params.json is saved before the value is corrected, it saves it as None (jsaon: null)
    else:
        try:
            v = float(v)
            if k in ['tmax', 'dt', 'dtFD', 'tb']:   # converts ns to s
                v = round(v * 1e-9, 13)
            elif k in ['Ion', 'Ioff', 'Iss']:       # converts mA to A
                v = round(v * 1e-3, 5)
        except:
            #sg.popup('Float input expected', button_type=5, auto_close = True, auto_close_duration = 1.5)
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
    for k, v in d.items():                              # update entire dictionary...
        if not k in ['nH', 'nRIN', 'nSeg', 'fmaxplot']: # ... except for the "hyperparameters" that are in the dictionary and json params files, but not in the GUI 
            k, d[k] = check_value(k, values[k])
    return d


# GUI/params mgmt: disable elements not used for certain combinations of parameters
def manage_params_combos(window, values):

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
    
    elif values['modFormat'] == 'pulse':
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)

    elif values['modFormat'] == 'random bits':
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = False)

    elif values['modFormat'] == 'small signal':
        if values['Hfplot'] == True:
            window['tmax'].Update(disabled = True)
        else:
            window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = False)
        window['tb'].Update(disabled = True)

    elif values['modFormat'] == 'steady state':
        if values['RINplot'] == True:
            window['tmax'].Update(disabled = True)
        else:
            window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)      

    # forces FD solution and disables freq resp calculation for simulations including noise
    if values['Noise'] == True:
        window['RINplot'].Update(disabled = False)
        window['Hfplot'].Update(False)
        values['Hfplot'] = False
        window['Hfplot'].Update(disabled = True)
        window['odeSolver'].Update('Finite Diff.')
        values['odeSolver'] = 'Finite Diff.'
        window['dtFD'].Update(disabled = False)
    else:
        window['RINplot'].Update(False)
        values['RINplot'] = False
        window['RINplot'].Update(disabled = True)
        window['Hfplot'].Update(disabled = False)

    # forces modeFormat to 'small signal' when plotting freq resp is selected
    if values['Hfplot'] == True:
        window['modFormat'].Update('small signal')
        values['modFormat'] = 'small signal'
        window['tmax'].Update(disabled = True)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = False)
        window['tb'].Update(disabled = True)
        window['RINplot'].Update(False)
        values['RINplot'] = False
        window['Eyeplot'].Update(False)
        values['Eyeplot'] = False

    # forces modeFormat to 'steady state' when plotting RIN is selected
    if values['RINplot'] == True:
        window['modFormat'].Update('steady state')
        values['modFormat'] = 'steady state'
        window['tmax'].Update(disabled = True)
        window['Ioff'].Update(disabled = True)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = True)
        window['Hfplot'].Update(False)
        values['Hfplot'] = False
        window['Eyeplot'].Update(False)
        values['Eyeplot'] = False

    # forces modeFormat = 'random bits' when plotting eye is selected
    if values['Eyeplot'] == True:
        window['modFormat'].update('random bits')
        values['modFormat'] = 'random bits'
        window['tmax'].Update(disabled = False)
        window['Ioff'].Update(disabled = False)
        window['Iss'].Update(disabled = True)
        window['tb'].Update(disabled = False)
        window['RINplot'].Update(False)
        values['RINplot'] = False
        window['Hfplot'].Update(False)
        values['Hfplot'] = False

    return values


# GUI/params mgmt: update gui (incl. disabling elements) after initializing or loading params file
def update_gui(window, values, sp, vp):
    for k, v in sp.items():                     # update GUI with values stored in the simParams dictionary...
        if not k in ['nH', 'nRIN', 'nSeg', 'fmaxplot']:   # ... except for the "hyperparameters" that are in the dictionary and json params files, but not in the GUI 
            if k in ['tmax', 'dt', 'dtFD', 'tb']:
                values[k] = sp[k] * 1e9         # convert s to ns
            elif k in ['Ion', 'Ioff', 'Iss']:
                values[k] = sp[k] * 1e3         # convert A to mA
            else:
                values[k] = sp[k]
            window[k].update(values[k])         # update GUI values
    for k, v in vp.items():                     # update entire vcselParams dictionary
        values[k] = vp[k]
        window[k].update(values[k])
    # manage_params_combos(window, values)

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
            [sg.Checkbox('Carrier transport into the quantum wells', key='SCHtransp', default=sp['SCHtransp'], size=(45,1), disabled=False, tooltip=ttips['SCHtransp'], enable_events=True)],
            [sg.Checkbox('Thermal effects', key='ThermMod', default=sp['ThermMod'], disabled=True, tooltip=ttips['ThermMod'], enable_events=True)],
            [sg.Checkbox('Noise', key='Noise', default=sp['Noise'], disabled=False, tooltip=ttips['Noise'], enable_events=True)],
            [sg.Checkbox('Electrical parasitics', key='Parasitics', default=sp['Parasitics'], tooltip=ttips['Parasitics'], enable_events=True)],
            [sg.Checkbox('2D', key='2D', default=sp['2D'], disabled=True, tooltip=ttips['2D'], enable_events=True)],
            ])],
        [sg.Frame('Simulation parameters',[
            [sg.InputText(sp['nNw'], key='nNw', size=(7,1), tooltip=ttips['nNw'], enable_events=True), sg.Text('radial resolution (typically 7-20 terms)', size=(40,1))],
            #[sg.Checkbox('store carrier terms Ni', key='storeN', default=sp['storeN'], disabled=True, tooltip=ttips['storeN'], enable_events=True)],
            [sg.Combo(values=('Finite Diff.', 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'), key='odeSolver', default_value=sp['odeSolver'], size=(12,1), tooltip=ttips['odeSolver'], enable_events=True), sg.Text('ODE solver')],
            [sg.InputText(round(float(0 if sp['tmax'] is None else sp['tmax'])*1e9,0), key='tmax', size=(7,1), tooltip=ttips['tmax'], enable_events=True, readonly=(sp['Hfplot']==True or sp['RINplot']==True), disabled_readonly_background_color='grey'), sg.Text('simulated time (ns)')],
            [sg.InputText(round(float(0 if sp['dt'] is None else sp['dt'])*1e9,3), key='dt', size=(7,1), tooltip=ttips['dt'], enable_events=True), sg.Text('time resolution for storing and plotting (ns)')],
            [sg.InputText(round(float(0 if sp['dtFD'] is None else sp['dtFD'])*1e9,3), key='dtFD', size=(7,1), tooltip=ttips['dtFD'], enable_events=True, readonly=(sp['odeSolver']!='Finite Diff.'), disabled_readonly_background_color='grey'), sg.Text('time step for FD solution (ns)')],
            ])],
        [sg.Frame('Current modulation pattern',[
            [sg.Combo(values=('step', 'pulse', 'random bits', 'small signal', 'steady state'), key='modFormat', default_value=sp['modFormat'], size=(12,1), tooltip=ttips['modFormat'], enable_events=True), sg.Text('modulation pattern')],
            [sg.InputText(round(float(0 if sp['Ion'] is None else sp['Ion'])*1e3,1), key='Ion', size=(7,1), tooltip=ttips['Ion'], enable_events=True), sg.Text('ON current (mA)', size=(40,1))],
            [sg.InputText(round(float(0 if sp['Ioff'] is None else sp['Ioff'])*1e3,1), key='Ioff', size=(7,1), tooltip=ttips['Ioff'], enable_events=True, readonly=(sp['modFormat']=='small signal' or sp['modFormat']=='step'), disabled_readonly_background_color='grey'), sg.Text('OFF current (mA)')],
            [sg.InputText(round(float(0 if sp['Iss'] is None else sp['Iss'])*1e3,2), key='Iss', size=(7,1), tooltip=ttips['Iss'], enable_events=True, readonly=(sp['modFormat']!='small signal'), disabled_readonly_background_color='grey'), sg.Text('small signal current step (mA)')],
            [sg.InputText(round(float(0 if sp['tb'] is None else sp['tb'])*1e9,1), key='tb', size=(7,1), tooltip=ttips['tb'], enable_events=True, disabled=(sp['modFormat']!='random bits'), disabled_readonly_background_color='grey'), sg.Text('bit duration (ns)')],
            ])],
        [sg.Frame('Results visualization',[
            [sg.Checkbox('2D mode profiles', key='Modes2Dplot', default=sp['Modes2Dplot'], size=(45,1), tooltip=ttips['Modes2Dplot'], enable_events=True)],
            [sg.Checkbox('Steady-state (LI) characteristic Popt(I)', key='PIplot', default=sp['PIplot'], size=(45,1), tooltip=ttips['PIplot'], enable_events=True)],
            [sg.Checkbox('Dynamic characteristic Popt(t)', key='Ptplot', default=sp['Ptplot'], tooltip=ttips['Ptplot'], enable_events=True)],
            [sg.Checkbox('2D optical & carrier profiles within the cavity', key='NwS2Dplot', default=sp['NwS2Dplot'], size=(45,1), tooltip=ttips['NwS2Dplot'], enable_events=True)],
            [sg.Checkbox('Relative Intensity Noise (RIN) spectrum', key='RINplot', default=sp['RINplot'], disabled=(sp['Noise']!=True), tooltip=ttips['RINplot'], enable_events=True)],
            [sg.Checkbox('Frequency response H(f)', key='Hfplot', default=sp['Hfplot'], disabled=(sp['Noise']==True), tooltip=ttips['Hfplot'], enable_events=True)],
            [sg.Checkbox('Eye diagram', key='Eyeplot', default=sp['Eyeplot'], tooltip=ttips['Eyeplot'], enable_events=True)],
            ])],
        ]

    vpTab_layout = [
        [sg.Frame('Generic description',[
            [sg.Multiline(vp['vcselDescr'], key='vcselDescr', size=(52,2), tooltip=ttips['vcselDescr'], enable_events=True)],
            ])],
        [sg.Frame('Cavity geometry parameters',[
            [sg.InputText(vp['rox'], key='rox', size=(7,1), tooltip=ttips['rox'], enable_events=True), sg.Text('rox: oxide aperture radius (cm)', size=(40,1))],
            [sg.InputText(vp['Leff'], key='Leff', size=(7,1), tooltip=ttips['Leff'], enable_events=True), sg.Text('Leff: effective cavity length (cm)')],
            [sg.InputText(vp['nw'], key='nw', size=(7,1), tooltip=ttips['nw'], enable_events=True), sg.Text('nw: number of quantum wells')],
            [sg.InputText(vp['dw'], key='dw', size=(7,1), tooltip=ttips['dw'], enable_events=True), sg.Text('dw: single QW thickness (cm)')],
            [sg.InputText(vp['db'], key='db', size=(7,1), tooltip=ttips['db'], enable_events=True), sg.Text('db: SCH thickness (cm)')],
            ])],
        [sg.Frame('Equivalent waveguide parameters',[
            [sg.InputText(vp['wl0'], key='wl0', size=(7,1), tooltip=ttips['wl0'], enable_events=True), sg.Text('wl0: emission wavelength @300k (nm)', size=(40,1))],
            [sg.InputText(vp['nc'], key='nc', size=(7,1), tooltip=ttips['nc'], enable_events=True), sg.Text('nc: core equivalent refractive index')],
            [sg.InputText(vp['ng'], key='ng', size=(7,1), tooltip=ttips['ng'], enable_events=True), sg.Text('ng: group refractive index')],
            [sg.InputText(vp['dn'], key='dn', size=(7,1), tooltip=ttips['dn'], enable_events=True), sg.Text('dn: equivalent fractional refractive index change')],
            ])],
        [sg.Frame('Optical parameters',[
            [sg.InputText(vp['Rt'], key='Rt', size=(7,1), tooltip=ttips['Rt'], enable_events=True), sg.Text('Rt: top mirror reflectivity', size=(40,1))],
            [sg.InputText(vp['Rb'], key='Rb', size=(7,1), tooltip=ttips['Rb'], enable_events=True), sg.Text('Rb: bottom mirror reflectivity')],
            [sg.InputText(vp['alphai'], key='alphai', size=(7,1), tooltip=ttips['alphai'], enable_events=True), sg.Text('alphai: internal losses (cm-1)')],
            [sg.InputText(vp['beta'], key='beta', size=(7,1), tooltip=ttips['beta'], enable_events=True), sg.Text('beta: spontaneous recombination coefficient')],
            ])],
        [sg.Frame('Gain parameters',[
            [sg.InputText(vp['gln'], key='gln', size=(7,1), tooltip=ttips['gln'], enable_events=True), sg.Text('gln: logarithmic gain coefficient (cm-1)', size=(40,1))],
            [sg.InputText(vp['Ntr'], key='Ntr', size=(7,1), tooltip=ttips['Ntr'], enable_events=True), sg.Text('Ntr: transparency carrier density (cm-3)')],
            [sg.InputText(vp['epsilon'], key='epsilon', size=(7,1), tooltip=ttips['epsilon'], enable_events=True), sg.Text('epsilon: gain compression factor')],
            [sg.InputText(vp['GamR'], key='GamR', size=(7,1), tooltip=ttips['GamR'], enable_events=True), sg.Text('GamR: relative confinement factor')],
            ])],
        [sg.Frame('Carrier transport and diffusion parameters',[
            [sg.InputText(vp['tauNb'], key='tauNb', size=(7,1), tooltip=ttips['tauNb'], enable_events=True), sg.Text('tauNb: carrier lifetime in the barriers (s)', size=(40,1))],
            [sg.InputText(vp['tauNw'], key='tauNw', size=(7,1), tooltip=ttips['tauNw'], enable_events=True), sg.Text('tauNw: carrier lifetime in the QWs (s)')],
            [sg.InputText(vp['tauCap'], key='tauCap', size=(7,1), tooltip=ttips['tauCap'], enable_events=True), sg.Text('tauCap: ambipolar diffusion time (s)')],
            [sg.InputText(vp['tauEsc'], key='tauEsc', size=(7,1), tooltip=ttips['tauEsc'], enable_events=True), sg.Text('tauEsc: thermionic emission lifetime (s)')],
            [sg.InputText(vp['etai'], key='etai', size=(7,1), tooltip=ttips['etai'], enable_events=True), sg.Text('etai: current injection efficiency')],
            [sg.InputText(vp['rs'], key='rs', size=(7,1), tooltip=ttips['rs'], enable_events=True), sg.Text('rs: current spreading coefficient (cm)')],
            [sg.InputText(vp['DN'], key='DN', size=(7,1), tooltip=ttips['DN'], enable_events=True), sg.Text('DN: ambipolar diffusion coeff. (cm2/s)')],
            ])],
        [sg.Frame('Electrical parasitic elements',[
            [sg.InputText(vp['Cp'], key='Cp', size=(7,1), tooltip=ttips['Cp'], enable_events=True), sg.Text('CP: pad parasitic capacitance (H)', size=(40,1))],
            [sg.InputText(vp['Rm'], key='Rm', size=(7,1), tooltip=ttips['Rm'], enable_events=True), sg.Text('Rm: bragg reflectors parasitic resistance (ohm)')],
            [sg.InputText(vp['Ca'], key='Ca', size=(7,1), tooltip=ttips['Ca'], enable_events=True), sg.Text('Ca: active region parasitic capacitance (H)')],
            [sg.InputText(vp['Ra'], key='Ra', size=(7,1), tooltip=ttips['Ra'], enable_events=True), sg.Text('Ra: active region parasitic resistance (ohm)')],
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
            if 'S' in locals():
                del rho, nrho, phi, nphi, nNw, J0i, nS, LPlm, lvec, Ur, Icw, NScw, It, f, H, Hp, RIN, S2P, teval, S, Nb, Nw
            if 'sr' in locals():    
                del sp, vp, sr
            else:
                del sp, vp
            params = load_params('default_params.json')         # loads params file and populates dictionaries sp and vp
            sp, vp = params['simParams'], params['vcselParams'] # simulation results not loaded
            del params
            save_params(sp, vp, {}, 'last_params.json')         # simulation results not saved to "last_params.json"
            update_gui(window, values, sp, vp)

        elif event == 'load file':
            if 'S' in locals():
                del rho, nrho, phi, nphi, nNw, J0i, nS, LPlm, lvec, Ur, Icw, NScw, It, f, H, Hp, RIN, S2P, teval, S, Nb, Nw
            if 'sr' in locals():    
                del sp, vp, sr
            else:
                del sp, vp
            file_name = values['load file']
            if file_name != '':
                params = load_params(file_name)             # loads params file and populates dictionaries sp and vp
                sp, vp, sr = params['simParams'], params['vcselParams'], params['simResults']   # for post-processing, simulations results should be converted from list to numpy array and saved in the respective variables (S, Nw, ur, etc.)
                del params
                update_gui(window, values, sp, vp)
                values = manage_params_combos(window, values)   # update GUI based on params combos

        elif event == 'save file':
            file_name = values['save file']
            if 'S' in locals():
                sr = {
                    "phi": phi.tolist(),    # dictionary can't store numpy arrays -> conversion to list
                    "nphi": nphi,
                    "rho": rho.tolist(),
                    "nrho": nrho,
                    "nNw": nNw,
                    "J0i": J0i.tolist(),
                    "nS": nS,
                    "LPlm": LPlm,
                    "lvec": lvec.tolist(),
                    "Ur": Ur.tolist(),
                    "Icw": Icw.tolist(),
                    "NScw": NScw.tolist(),
                    "It": It.tolist(),
                    "teval": teval.tolist(),
                    "f": f.tolist(),
                    "H": H.tolist(),
                    "Hp": Hp.tolist(),
                    "RIN": RIN.tolist(),
                    "S2P": S2P.tolist(),
                    "Nb": Nb.tolist(),
                    "Nw": Nw.tolist(),
                    "S": S.tolist(),
                }
            else:
                sr = {}
            if file_name == '': file_name = 'last_params.json'
            save_params(sp, vp, sr, file_name)  # aggregates sp and vp and sr in a dict of dicts and saves to a json file
            del sr

        elif event == 'run simulation':
            save_params(sp, vp, {}, 'last_params.json') # simulation results not saved to "last_params.json"
            if 'S' in locals():
                del rho, nrho, phi, nphi, nNw, J0i, nS, LPlm, lvec, Ur, Icw, NScw, It, f, H, Hp, RIN, S2P, teval, S, Nb, Nw
            
            rho, nrho, phi, nphi, nNw, J0i, nS, LPlm, lvec, Ur, Icw, NScw, It, f, H, Hp, RIN, S2P, teval, S, Nb, Nw = VISTAS1D(sp, vp)
            
            # extracts tmax from teval in case it is calculated in the main loop (RIN or Hf)
            sp['tmax'] = round(max(teval * 1e9), 0) * 1e-9
            window['tmax'].Update(int(sp['tmax'] * 1e9))
            values['tmax'] = int(sp['tmax'] * 1e9)
            
            # visualization            
            try:     # popup to collect figure number; 1 if invalid input
                nfig = int(sg.popup_get_text(message= 'First figure number', title='', default_text='', size=(2, 1), no_titlebar = False, keep_on_top = True))
            except:
                nfig = 1

            if sp['Modes2Dplot'] == 1:  # 2D mode profiles (cosine azimuthal distribution) Ur(x,y)
                nfig = plotModes2D(Ur, LPlm, lvec, nS, rho*1e-2, nrho, phi, nphi, nfig=nfig)              
        
            if sp['PIplot'] == 1:   # steady-state LI characteristic Popt(I)
                nfig = plotPower(Icw * 1e3, S2P*NScw[NScw.shape[0]-nS:,:], LPlm, xlabel = 'current (mA)', nfig=nfig)  

            if sp['Ptplot'] == 1:   # dynamic response Popt(t)
                nfig = plotPower(teval * 1e9, S2P*S, LPlm, xlabel = 'time (ns)', nfig=nfig)

            if sp['NwS2Dplot'] == 1: # carrier and optical field profiles in the cavity
                nfig = plotNwS2D(Nw[:, -1], J0i, S[:, -1], Ur, rho*1e-2, nrho, phi, nphi, nfig=nfig) # by default the last point of Nw
                
            if sp['modFormat'] == 'random bits' and sp['Eyeplot'] == True:
                nfig = plotEye(teval, sp['dt'], sp['tb'], S2P*S, It, nfig=nfig)

            if sp['modFormat'] == 'small signal' and sp['Hfplot'] == True: # small signal response H(f)
                nfig = plotH(f, H, sp['Parasitics'], Hp, 0, 15, -10, 10, nfig=nfig)

            if sp['Noise'] == True and sp['RINplot'] == True:
                nfig = plotRIN(f, RIN, 0, 15, -160, -100, nfig=nfig)

            plt.show(block = False)

        else:   # case 'gui-element change'
            values = manage_params_combos(window, values)   # update GUI based on params combos
            sp = update_dict(values, sp)                    # update dictionary
            vp = update_dict(values, vp)                    # update dictionary
            #update_gui(window, values, sp, vp)


def main():
    GUI()


if __name__ == "__main__":
    main()