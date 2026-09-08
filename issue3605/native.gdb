set pagination off
set confirm off
set breakpoint pending on
set debuginfod enabled off
set print thread-events off
set new-console on
handle SIGINT nostop noprint pass
python
import gdb, json, os, subprocess

injections = 0
errors = []
factors = 0
wanted = int(os.environ['PROBE_FACTOR'])

class NativeFactor(gdb.Breakpoint):
    def stop(self):
        global injections
        self.enabled = False
        if injections == 0:
            gdb.execute('bt 20')
        try:
            subprocess.run([os.environ['PROBE_PYTHON'], os.environ['PROBE_SCRIPT'],
                            'send', '--pid', str(gdb.selected_inferior().pid)],
                           check=True, timeout=10)
            injections += 1
        except Exception as error:
            errors.append(str(error))
            print('Ctrl+C injection failed:', error)
            return True
        return False

native = NativeFactor('__dmumps_fac_front_aux_m_MOD_dmumps_fac_i_ldlt')
native.enabled = False

class FactorDriver(gdb.Breakpoint):
    def stop(self):
        global factors
        factors += 1
        if factors == wanted:
            native.enabled = True
        return False

class Arm(gdb.Breakpoint):
    def stop(self):
        global factors
        factors = 0
        native.enabled = False
        return False

FactorDriver('dmumps_fac_driver_')
Arm('arm_interrupt')
end
run
python
with open(os.environ['PROBE_INJECTIONS'], 'w') as f:
    json.dump(dict(injections=injections, errors=errors), f)
print('Verified native-factorization Ctrl+C injections:', injections)
end
