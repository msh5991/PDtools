import os,re
from pymatgen.analysis.phase_diagram import PhaseDiagram, CompoundPhaseDiagram, PDEntry, PDPlotter
from pymatgen.core.composition import Composition

class IgorPDPlotter(PDPlotter):
    def __init__(
        self,
        phasediagram: PhaseDiagram,
        show_unstable: float = 0.2,
        backend: str = "plotly",
        **plotkwargs,
    ):
        super().__init__(
            phasediagram,
            show_unstable,
            backend,
            **plotkwargs,
        )
    
    def igorplot_enghalpy(self, pd, target_pressure, itx):
        """
        plot
            (1) formation enthalpy vs pressure
            (2) enthalpy above convex hull vs pressure
        by Igor Pro

        Args:
            pd: PhaseDiagram object (defined in the pymatgen.analysis.phase_diagram module)
            target_pressure: Target pressure [KBar] in the vc-relax calculation
            itx: Name of output .itx file (exported as Igor text format)
        """

        target_pressure /= 10     # Kbar -> GPa conversion
        complist, Eformlist, Eahlist = self._get_enthalpy(pd)

        if not os.path.isfile(itx):
            with open(itx, 'w') as f:
                f.write('IGOR\n')
                waveline = 'WAVES/D pressure '
                for i in complist:
                    waveline = waveline + i + '_Eform '
                    waveline = waveline + i + '_Eah '
                f.write(waveline)
                f.write('\nBEGIN\nEND\n')
                f.write('X Display ' + complist[0] + '_Eform vs pressure\n')
                for i in complist[1:]:
                    f.write('X AppendToGraph ' + i + '_Eform vs pressure\n')
                f.write('X ModifyGraph gfSize=18,marker=19,msize=3,lsize=1,gFont="Arial",tick=2,mirror=1,btLen=8,zero(left)=3,ZisZ=1,standoff=0\n')
                f.write('X Label bottom "pressure (GPa)"\n')
                f.write('X Label left "Formation Enthalpy (eV)"\n')
                f.write('X Display ' + complist[0] + '_Eah vs pressure\n')
                for i in complist[1:]:
                    f.write('X AppendToGraph ' + i + '_Eah vs pressure\n')
                f.write('X ModifyGraph gfSize=18,marker=19,msize=3,lsize=1,gFont="Arial",tick=2,mirror=1,btLen=8,zero(left)=3,ZisZ=1,standoff=0\n')
                f.write('X Label bottom "pressure (GPa)"\n')
                f.write('X Label left "Enthalpy above convex hull (eV)"\n')
        
        with open(itx, 'r') as f:
            lines = f.readlines()
            insertline = lines.index('END\n')
            txt = str(target_pressure)
            for e in Eformlist:
                txt +=  ' ' + str(e) + ' ' + str(Eahlist[Eformlist.index(e)])
            txt += '\n'
            lines.insert(insertline, txt)
    
        with open(itx, 'w') as f:
            f.writelines(lines)

    def igorplot_2Dpd(self, pd, itx, prefix):
        """
        plot 2-dimensional phase diagram by Igor Pro

        Args:
            pd : PhaseDiagram object (defined in the pymatgen.analysis.phase_diagram module)
            itx : Name of output .itx file (exported as Igor text format)
            prefix : prefix of wavename
        """
        
        if pd.dim != 2:
            raise Exception('The Dimension of PhaseDiagram is not equal to 2.')
        
        with open(itx, 'w') as f:
            complist = []
            f.write('IGOR\n')

            # entries on convex hull
            f.write('WAVES/D '+prefix+'_x '+prefix+'_Eform '+prefix+'_Eah\n')
            f.write('BEGIN\n')
            waves = {}
            for entry in pd.stable_entries:
                x = pd.pd_coords(entry.composition)[0]
                Eform = pd.get_form_energy_per_atom(entry)
                Eah = pd.get_e_above_hull(entry)
                text = str(x) + ' ' + str(Eform) + ' ' + str(Eah) + '\n'
                waves[x] = text
            print("waves",waves)
            waves_sorted = sorted(waves.items(), key=lambda x:x[0])
            print("waves_sorted",waves_sorted)
            for w in waves_sorted:
                f.write(w[1])
            f.write('END\n')

            f.write('X Display ' + prefix +'_Eform vs ' + prefix +'_x\n')
            f.write('X ModifyGraph gfSize=18,marker=19,msize=3,lsize=1,gFont="Arial",mode=4,tick=2,mirror=1,btLen=8,zero(left)=3,standoff=0,ZisZ=1\n')
            f.write('X Label bottom ""\n')
            f.write('X Label left "Formation Enthalpy (eV)"\n')

            # entries above convex hull

            # labels
            waves = {}
            for entry in pd.all_entries:
                x = pd.pd_coords(entry.composition)[0]
                comp = str(entry.original_entry.composition.reduced_formula)
                text = re.sub('([0-9]+)', r'\\B\1\\M', comp)
                waves[x] = '\"'+text+'\"\n'
            waves_sorted = sorted(waves.items(), key=lambda x:x[0])
            f.write('WAVES/T '+prefix+'_label\n')
            f.write('BEGIN\n')
            for w in waves_sorted:
                f.write(w[1])
            f.write('END\n')
            f.write('WAVES/D '+prefix+'_labelpos\n')
            f.write('BEGIN\n')
            for w in waves_sorted:
                f.write(str(w[0])+'\n')
            f.write('END\n')
            f.write('X ModifyGraph tkLblRot(bottom)=90,userticks(bottom)={'+prefix+'_labelpos,'+prefix+'_label}\n')
        return

    def _get_enthalpy(self, pd):
        entries = pd.all_entries
        complist = []
        Eformlist = []
        Eahlist = []
        print('Formation Enthalpy (eV):')
        for i in entries:
            e = pd.get_form_energy_per_atom(i)
            print('    '+str(i.original_entry.composition) + ' : ' + str(e))
            complist.append(str(i.original_entry.composition).replace(' ', ''))
            Eformlist.append(str(e))
        print('Enthalpy above convex hull (eV):')
        for i in entries:
            decomp, e = pd.get_decomp_and_e_above_hull(i)
            print('    '+str(i.original_entry.composition) + ' : ' + str(e))
            Eahlist.append(str(e))
        return complist, Eformlist, Eahlist

def PWout2PD(outfiles, target_pressure, terminal_compositions):
    """
    Make PhaseDiagram object (defined in the pymatgen.analysis.phase_diagram module)
    from the output files of vc-relax calculation by PWSCF

    Args:
        outfiles: List of output files for vc-relax calculation by PWSCF
        target_pressure: Target pressure [KBar] in the vc-relax calculation
        terminal_compositions: Terminal compositions for the phase diagram
    
    Returns:
        PhaseDiagram object
    """

    Ry_to_J = 2.179872
    Ry_to_eV = 13.605693123

    entries = []

    for f_out in outfiles:
        with open(f_out, 'r') as f:
            print('reading '+ f_out + ' ...')
            lines = f.readlines()

            #read compositions
            j1 = [i for i, line in enumerate(lines) if 'ATOMIC_POSITIONS' in line]
            j2 = [i for i, line in enumerate(lines) if 'End final coordinates' in line]
            atm = []
            natm = []
            for i in range(j1[-1]+1, j2[0]):
                a = lines[i].split()[0]
                if a in atm:
                    natm[atm.index(a)] += 1
                else:
                    atm.append(a)
                    natm.append(1)
            comp = ''
            for i in range(len(atm)):
                comp += atm[i]
                comp += str(natm[i])
                comp += ' '

            #read Hscf
            j = [i for i, line in enumerate(lines) if 'Begin final coordinates' in line]
            l = re.search(r'[0-9]+\.[0-9]+ +Ang',lines[j[0] + 1])
            volume = re.sub(r' +Ang', '', l.group())
            j = [i for i, line in enumerate(lines) if '!' in line]
            Escf = re.search(r'[-]*[0-9]+\.[0-9]+',lines[j[-1]])
            Escf = float(Escf.group())
            Hscf = Escf + float(target_pressure) * float(volume) /(Ry_to_J *10000)
            Hscf *= Ry_to_eV
            Escf *= Ry_to_eV

        print('    composition = '+comp)
        print('    Energy      = '+str(Escf)+' eV')
        print('    Enthalpy    = '+str(Hscf)+' eV')

        entries.append(PDEntry(Composition(comp), Hscf))
    
    pd = CompoundPhaseDiagram(entries, terminal_compositions, normalize_terminal_compositions=True)
    return pd