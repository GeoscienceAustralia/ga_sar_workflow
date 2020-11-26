from collections import OrderedDict
from pathlib import Path
import regex

def usage2decl(module, program, file):
    state = 0
    params = OrderedDict()

    for line in file.read().splitlines():
        if line.strip().startswith('**') or len(line.strip()) == 0:
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
                    print('line:', line)
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

                ignore_arg = argname.lower().startswith('note:')
                ignore_arg |= argname.lower().startswith('remark')
                ignore_arg |= argname.lower().startswith('example:')
                ignore_arg |= argname.startswith('...')  # HACK: This is... not ideal, ... (usually?) means variable number of arguments

                if not ignore_arg:
                    # Assert we're reading in the correct order
                    # assert(all_args.index(argname) == last_arg+1)
                    # DISABLED for now... their documentation has inconsistent namings...
                    last_arg += 1

                    # HACK: gamma's help is full of typos, we trust the initial argnames... not what's in the desc
                    try:
                        positional_argname = all_args[last_arg]
                        argname = positional_argname
                    except:
                        pass # print('!!!', last_arg, all_args)

                    params[argname] = {
                        'desc': '',
                        'optional': last_arg >= len(required_args)  # argname in optional_args
                    }

                desc = line[len(argname)+2:].lstrip()

                is_infile = desc.startswith('(input)') or desc.startswith('(input/output)')
                is_outfile = desc.startswith('(output)') or desc.startswith('(input/output)')

            else:  # continue parsing existing arg
                desc = line.lstrip()

            if ignore_arg:
                continue

            params[argname]['desc'] += desc + '\n'

            # parse description
            enum_prefix = regex.match(r' *(\d+):.*', desc)

            # detect filepath params
            if is_infile or is_outfile:
                params[argname]['type'] = 'path'
                params[argname]['is_infile'] = is_infile
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

    # A set of programs our parser doesn't handle properly yet
    blacklist = [
        'coord_trans',  # This isn't a normal gamma command
        'mosaic',
        'PALSAR_antpat',
        'JERS_fix',
        'comb_hsi',
        'soil_moisture',
        'validate',
        'ASAR_XCA'  # This has two programs of the same name?
    ]

    # Parse program details
    for module_dir in basedir.iterdir():
        module = module_dir.name

        decls[module] = {}

        for usage_path in module_dir.iterdir():
            program = usage_path.stem

            # Ignore some unused problematic programs
            if program in blacklist:
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
            'PyGammaCall = NamedTuple("PyGammaCall", [("module", str), ("program", str), ("parameters", Dict[str, object])])',
            '',
            '',
            'class SimpleParFile(object):',
            '    values = {}',
            '',
            '    def __init__(self, path):',
            '        with open(path, \'r\') as file:',
            '            lines = file.read().splitlines()[2:]  # Skip header lines',
            '',
            '            for line in lines:',
            '                value_id = line.split(\':\')[0]',
            '                value_data = line[len(value_id)+2:].strip()',
            '',
            '                self.values[value_id] = value_data',
            '',
            '    def get_value(self, value_id: str, dtype = str, index: int = 0):',
            '        if dtype == str:',
            '            return self.values[value_id]',
            '',
            '        return dtype(self.values[value_id].split()[index])',
            '',
            '',
            'class PyGammaTestProxy(object):',
        ]))

        writeline('\n'.join(indent(1, [
            'ParFile = SimpleParFile',
            '',
            'call_sequence: Sequence[PyGammaCall]',
            'call_count: Dict[str, int]',
            '',
            'def __init__(self):',
            '    self.reset_proxy()',
            '',
            'def reset_proxy(self):',
            '    self.call_sequence = []',
            '    self.call_count = {}',
            '',
            'def _validate(self, condition, result):',
            '    stat, stdout, stderr = result',
            '',
            '    # TODO: error stats?',
            '    stat = stat if condition else -1',
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

                    argname = argname.replace('-', '_').replace('/', '_').replace('\\', '_')
                    argname = 'definition' if argname == 'def' else argname

                    args += argname

                    if param['type'] == 'path':
                        args += ': str'

                    if param['optional']:
                        args += ' = None'

                writeline('\n'.join(indent(1, [
                    '',
                    f'def {program.replace("-", "_")}(self, {args}):',
                    '    supplied_args = locals()',
                    '    result = (0, "", "")',
                    '',
                    f'    self.call_sequence.append(PyGammaCall("{module}", "{program}", supplied_args))',
                    '',
                    f'    if "{program}" in self.call_count:',
                    f'        self.call_count["{program}"] += 1',
                    '    else:',
                    f'        self.call_count["{program}"] = 1',
                    '',
                ])))

                # validate args, ensure input files exist, touch output files, etc
                # # TODO: maybe generate stdout where appropriate?
                for argname, param in decl['params'].items():
                    argname = argname.replace('-', '_').replace('/', '_').replace('\\', '_')
                    argname = 'definition' if argname == 'def' else argname

                    # TODO: assert required args are not None, and allow optionals to be none

                    if param['type'] == 'path':
                        # Touch in/out files if they don't already exist
                        if param['is_outfile'] and param['is_infile']:
                            writeline('\n'.join(indent(2, [
                                f'if {argname} is not None and not Path({argname}).exists():',
                                f'    Path({argname}).touch()',
                            ])))

                        # Touch output files
                        elif param['is_outfile']:
                            writeline('\n'.join(indent(2, [
                                f'if {argname} is not None:',
                                f'    Path({argname}).touch()',
                            ])))

                        # Check input files exist
                        else:
                            writeline('\n'.join(indent(2, [
                                f'if {argname} is not None:',
                                f'    result = self._validate(Path({argname}).exists(), result)',
                            ])))

                    elif param['type'] == 'enum':
                        valid_values = param['enum']

                        writeline('\n'.join(indent(2, [
                            f'valid_values = {repr(valid_values)}',
                            f'result = self._validate({argname} in valid_values, result)',
                        ])))

                writeline('\n'.join(indent(2, [
                    'return result',
                ])))
