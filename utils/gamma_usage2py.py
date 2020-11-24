from collections import OrderedDict
from pathlib import Path
import regex

def usage2decl(module, program, file):
    state = 0
    params = OrderedDict()

    for line in file.read().splitlines():
        if len(line.strip()) == 0:
            continue

        if state == 0:
            if line.startswith('usage:'):
                # Grab arugments from 'usage' line:
                m = regex.match(r'usage: +([\w\-\\/_]+)( *\<([\w\-_]+)\>)*( *\[([\w\-_]+)\])*', line)
                
                try:
                    nci_path = m.captures(1)[0]
                    required_args = m.captures(3)
                    optional_args = m.captures(5)
                    all_args = required_args + optional_args
                    state = 1
                
                except Exception as e:
                    print('line:',line)
                    print(m)
                    print(e)
                    exit()

        elif state == 1:
            if line.startswith('input parameters:'):
                state = 2
                last_arg = -1

        elif state == 2:
            desc = ''

            if not line.startswith('    '):  # Start parsing new arg
                argname = line.split()[0]

                if argname != 'NOTE:':
                    # Assert we're reading in the correct order
                    #assert(all_args.index(argname) == last_arg+1)
                    # DISABLED for now... their documentation has inconsistent namings...
                    last_arg += 1

                    params[argname] = { 'desc': '' }

                desc = line[len(argname)+2:].lstrip()
                
                is_infile = desc.startswith('(input)')
                is_outfile = desc.startswith('(output)')
            
            else: # continue parsing existing arg
                desc = line.lstrip()

            if argname == 'NOTE:':
                continue

            params[argname]['desc'] += desc + '\n'

            # parse description
            enum_prefix = regex.match(r' *(\d+):.*', desc)

            # detect filepath params
            if is_infile or is_outfile:
                params[argname]['type'] = 'path'
                params[argname]['is_outfile'] = is_outfile

            # detect enum params
            elif enum_prefix:
                params[argname]['type'] = 'enum'

                if 'enum' in params[argname]:
                    params[argname]['enum'].append(int(enum_prefix[1]))
                else:
                    params[argname]['enum'] = [int(enum_prefix[1])]

            # All other variables have no validation info
            else:
                params[argname]['type'] = 'unknown'

    return {
        'module': module,
        'program': program,
        'params': params
    }

if __name__ == '__main__':
    basedir = Path.home() / 'GA/gamma_usage'

    decls = {}

    # Parse program details
    for module_dir in basedir.iterdir():
        module = module_dir.name
        
        decls[module] = {}

        for usage_path in module_dir.iterdir():
            program = usage_path.stem

            if program == 'coord_trans':
                continue

            with usage_path.open('r') as usage_file:
                decl = usage2decl(module, program, usage_file)

            decls[module][program] = decl

            if program == 'SLC_copy':
                print(decl)

    # Generate python stubs
    outdir = Path.home() / 'GA' / 'py_gamma_test_proxy.py'

    def indent(n, lines):
        return [('    '*n + line).rstrip() for line in lines]

    with outdir.open('w') as file:
        writeline = lambda x: file.write(x + '\n')

        writeline('\n'.join([
            'from pathlib import Path',
            'from typing import Sequence, NamedTuple, Dict, Union',
            '',
            'PyGammaCall = NamedTuple["PyGammaCall", [("module", str), ("program", str), ("parameters", Dict[str, object])]]',
            '',
            'class PyGammaTestProxy(object):',
        ]))

        writeline('\n'.join(indent(1, [
            'call_sequence: Sequence[PyGammaCall]',
            'call_count: Dict[str, int]',
            '',
            '',
            'def __init__(self)',
            '    self.reset_proxy()',
            '',
            'def reset_proxy(self)',
            '    self.call_sequence = []',
            '    self.call_count = {}',
            '',
            'def _validate(self, condition, result):',
            '    stat, stdout, stderr = result',
            '',
            '    # TODO: error stats?',
            '    stat = stat if condition else 1',
            '    # TODO: stderr?',
            '',
            '    return stat, stdout, stderr',
        ])))

        for module, programs in decls.items():
            for program, decl in programs.items():

                args = ''

                for argname, param in decl['params'].items():
                    if args != '':
                        args += ', '

                    args += argname.replace('-', '_')

                    if param['type'] == 'path':
                        args += ': str'

                writeline('\n'.join(indent(1, [
                    '',
                    f'def {program}(self, {args}):',
                    '    supplied_args = locals()',
                    '    result = (0, "", "")',
                    '',
                    f'    self.call_sequence.append(("{module}", "{program}", supplied_args))',
                    '',
                    f'    if "{program}" in self.call_count:',
                    f'        self.call_count["{program}"] += 1',
                    '    else:',
                    f'        self.call_count["{program}"]  = 1',
                    '',
                ])))

                # validate args, ensure input files exist, touch output files, etc
                # # TODO: maybe generate stdout where appropriate?
                for argname, param in decl['params'].items():
                    if param['type'] == 'path':
                        # Touch output files
                        if param['is_outfile']:
                            writeline('\n'.join(indent(2, [
                                f'Path({argname}).touch()',
                            ])))

                        # Check input files exist
                        else:
                            writeline('\n'.join(indent(2, [
                                f'result = self._validate(Path({argname}).exists(), result)',
                            ])))

                    elif param['type'] == 'enum':
                        valid_values = param['enum']

                        writeline('\n'.join(indent(2, [
                            f'valid_values = {repr(valid_values)}',
                            f'result = self._validate({argname} in valid_values, result)',
                        ])))

                writeline('\n'.join(indent(2, [
                    '',
                    'return result',
                ])))
